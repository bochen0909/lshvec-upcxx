
#include "argagg.hpp"
#include "utils.h"
#include "CRandProj.h"
#include "kmer.h"
#include "io.h"
#include "model.h"
#include "serialization.h"
#include "config.h"
#include "pbar.h"
#include "ProjConfig.h"

struct Config : public BaseConfig
{
    float learning_rate;
    uint32_t epoch;
    uint32_t dim;
    uint32_t neg_size;
    uint32_t half_window;
    size_t num_seq;
    std::string hash_file;
    std::string output_prefix;
    bool use_cbow = true;
    bool is_fasta = false;
    bool is_fastq = false;

    void print()
    {
        BaseConfig::print();
        myinfo("config: learning_rate=%f", learning_rate);
        myinfo("config: epoch=%ld", epoch);
        myinfo("config: dimension=%ld", dim);
        myinfo("config: neg_size=%ld", neg_size);
        myinfo("config: half_window=%ld", half_window);
        myinfo("config: hash_file=%s", hash_file.c_str());
        myinfo("config: output_prefix=%s", output_prefix.c_str());
        myinfo("config: use_cbow=%s", use_cbow ? "true" : "false");
    }
};

void show_help(const char *prog, argagg::parser &argparser)
{
    std::cout << prog << " v" << PROJECT_VERSION << "\n";

    std::cout << "Usage: " << prog << " [options] file1, file2 ....\n";
    std::cout << "Allowed options:\n";
    std::cout << argparser;
}

void run(Config &);

int main(int argc, char **argv)
{
    Config config;
    config.program = argv[0];
    config.rank = 0;
    config.mpi_hostname = sparc::get_hostname();
    set_spdlog_pattern(config.mpi_hostname.c_str(), config.rank);

    const char *program = argv[0];
    argagg::parser argparser{{

        {"help", {"-h", "--help"}, "shows this help message", 0},

        {"zip_output", {"-z", "--zip"}, "zip output files", 0},
        {"use_skipgram", {"--use-skipgram"}, "use skipgram (ohterwise cbow)", 0},
        {"use_fasta", {"--fasta"}, "input are fasta files", 0},
        {"use_fastq", {"--fastq"}, "input are fastq files", 0},

        {"output", {"-o", "--output"}, "output model prefix (default 'model')", 1},
        {"hash_file", {"--hash-file"}, "hash file to use", 1},
        {"learning_rate", {
                              "--lr",
                          },
         "initial learning rate (default 0.3)",
         1},

        {"dim", {
                    "--dim",
                },
         "number of dimension (default 300)",
         1},
        {"half_window", {
                            "--half-window",
                        },
         "half window size (default 5)",
         1},
        {"neg_size", {
                         "--neg-size",
                     },
         "negative words size (default 5)",
         1},

        {"epoch", {
                      "--epoch",
                  },
         "number of epochs to train (default 100)",
         1},
        {"n_thread", {
                         "--thread",
                     },
         "thread to use (default 0)",
         1},

    }};

    argagg::parser_results args;
    try
    {
        args = argparser.parse(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        show_help(program, argparser);
        return EXIT_FAILURE;
    }

    if (args["help"])
    {
        show_help(program, argparser);
        return EXIT_SUCCESS;
    }

    for (auto x : {"hash_file"})
        if (!args[x])
        {
            std::cerr << "ERROR: " << x << " is not set.\n";
            show_help(program, argparser);
            return EXIT_FAILURE;
        }
    config.nprocs = args["n_thread"].as<uint32_t>(0);
    if (config.nprocs == 0)
    {
        config.nprocs = sparc::get_number_of_thread();
    }
    config.backend = "smp";
    config.dim = args["dim"].as<uint32_t>(300);
    config.learning_rate = args["learning_rate"].as<float>(0.3);
    config.epoch = args["epoch"].as<uint32_t>(100);
    config.neg_size = args["neg_size"].as<uint32_t>(5);
    config.half_window = args["half_window"].as<uint32_t>(5);

    config.zip_output = args["zip_output"];
    config.is_fasta = args["use_fasta"];
    config.is_fastq = args["use_fastq"];
    config.use_cbow = !args["use_skipgram"];
    config.hash_file = args["hash_file"].as<std::string>();

    if (!sparc::file_exists(config.hash_file.c_str()))
    {
        std::cerr << "hash file does not exists: " << config.hash_file << std::endl;
        return EXIT_FAILURE;
    }

    if (args.pos.empty())
    {
        std::cerr << "no input files are provided" << std::endl;
        return EXIT_FAILURE;
    }
    else
    {
        config.inputpath = args.all_as<std::string>();
        bool b_error = false;
        for (size_t i = 0; i < config.inputpath.size(); i++)
        {
            if (!sparc::file_exists(config.inputpath.at(i).c_str()))
            {
                std::cerr << "Error, input file does not exists:  "
                          << config.inputpath.at(i) << std::endl;
                b_error = true;
            }
        }
        if (b_error)
        {
            return EXIT_FAILURE;
        }
    }

    config.output_prefix = args["output"].as<std::string>("model");
    config.print();

    run(config);

    return 0;
}

template <class BR>
void run_epoch(uint32_t this_epoch, Config &config, BR &reader, SingleNodeModel<float> &model, rpns::CRandProj &hash, float learning_rate)
{
    myinfo("Start epoch %ld, learning_rate=%.6f", this_epoch + 1, learning_rate);

    uint32_t batchsize = config.nprocs * 100;
    uint32_t kmer_size = hash.get_kmer_size();
    bool update_wc = this_epoch == 0;
    float sum_loss = 0;
    size_t num_of_seq = 0;
    char msg[128];
    sprintf(msg, "epoch %u:", this_epoch);
    PUnknownBar ubar(msg);
    PBar bar(msg, config.num_seq);
    while (true)
    {
        std::vector<FastaRecord> v = reader.next(batchsize);
        sparc::shuffle(v);
        num_of_seq += v.size();
        if (v.empty())
        {
            break;
        }

#pragma omp parallel for
        for (size_t i = 0; i < v.size(); i++)
        {
            std::string &seq = v.at(i).seq;

            auto strkmers = generate_kmer_for_fastseq(seq, kmer_size, "N", true);
            std::vector<uint32_t> kmers;
            for (auto &strkmer : strkmers)
            {
                kmers.push_back(hash.hash(strkmer, false));
            }

            float loss = model.update(kmers, learning_rate, update_wc);
            if (this_epoch == 0)
            {
                ubar.tick();
            }
            else
            {
                bar.tick();
            }
#pragma omp critical
            sum_loss += loss;
        }
    }
    if (this_epoch == 0)
    {
        ubar.end();
    }
    else
    {
        bar.end();
    }

    if (this_epoch == 0)
    {
        myinfo("Found that number of sequences is %ld", num_of_seq);
        config.num_seq = num_of_seq;
    }

    myinfo("End epoch %ld, loss=%f", this_epoch + 1, num_of_seq == 0 ? 0 : sum_loss / num_of_seq);
}

void run(Config &config)
{
    omp_set_num_threads(config.nprocs);
    rpns::CRandProj hash;
    hash.load(config.hash_file);

    size_t num_word = 1l << hash.get_hash_size();
    myinfo("kmer_size=%lu, hash_size=%lu, num_word=%lu", hash.get_kmer_size(), hash.get_hash_size(), num_word);
    SingleNodeModel<float> model(num_word, config.dim, config.neg_size, config.use_cbow, config.half_window);
    //model.randomize_init();
    model.uniform_init();
    for (uint32_t i = 0; i < config.epoch; i++)
    {
        float learning_rate = config.learning_rate * (config.epoch - i) / config.epoch;

        if (config.is_fasta)
        {
            BatchReader<FastaTextReaderBase, FastaRecord> reader(config.inputpath);
            run_epoch(i, config, reader, model, hash, learning_rate);
        }
        else if (config.is_fastq)
        {
            BatchReader<FastqTextReaderBase, FastaRecord> reader(config.inputpath);
            run_epoch(i, config, reader, model, hash, learning_rate);
        }
        else
        {
            BatchReader<SeqTextReaderBase, FastaRecord> reader(config.inputpath);
            run_epoch(i, config, reader, model, hash, learning_rate);
        }
    }

    save_vector_bin(config.epoch, config.output_prefix, model, config.zip_output);
    save_wordcounts(config.epoch, config.output_prefix, model, config.zip_output);
    myinfo("Finished!");
}

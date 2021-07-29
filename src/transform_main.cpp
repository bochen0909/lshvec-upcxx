
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
    std::string hash_file;
    std::string output_prefix;
    std::string model_path;
    bool is_fasta = false;
    bool is_fastq = false;

    void print()
    {
        BaseConfig::print();
        myinfo("config: hash_file=%s", hash_file.c_str());
        myinfo("config: output_file=%s", output_prefix.c_str());
        myinfo("config: model_path=%s", model_path.c_str());
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

        {"use_fasta", {"--fasta"}, "input are fasta files", 0},
        {"use_fastq", {"--fastq"}, "input are fastq files", 0},

        {"output", {"-o", "--output"}, "output file (default out.txt)", 1},
        {"vec_path", {"--vec"}, "vector file (binary) path", 1},
        {"hash_file", {"--hash-file"}, "hash file to use", 1},

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

    for (auto x : {"hash_file", "vec_path"})
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

    config.zip_output = args["zip_output"];
    config.is_fasta = args["use_fasta"];
    config.is_fastq = args["use_fastq"];
    config.hash_file = args["hash_file"].as<std::string>();
    config.model_path = args["vec_path"].as<std::string>();

    if (!sparc::file_exists(config.hash_file.c_str()))
    {
        std::cerr << "hash file does not exists: " << config.hash_file << std::endl;
        return EXIT_FAILURE;
    }
    if (!sparc::file_exists(config.model_path.c_str()))
    {
        std::cerr << "vector bin file does not exists: " << config.model_path << std::endl;
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

    config.output_prefix = args["output"].as<std::string>("output.txt");
    config.print();

    run(config);

    return 0;
}

void transform(const std::vector<uint32_t> &kmers, const std::vector<Vector<float>> &vecbin, Vector<float> &vec)
{
    vec.zero();
    if (kmers.empty())
    {
        return;
    }
    for (uint32_t kmer : kmers)
    {
        vec.add(vecbin.at(kmer));
    }
    vec.mul(1.0f / kmers.size());
}

template <typename BR, typename OS>
void run(Config &config, BR &reader, rpns::CRandProj &hash, OS &os)
{

    std::string modelpath = config.model_path;
    myinfo("reading vectors from %s\n", modelpath.c_str());
    std::vector<Vector<float>> wi;
    read_vec_bin(modelpath, wi);
    uint32_t dim = (uint32_t)wi.at(0).get_dim();
    myinfo("finish reading %lu vectors[dim=%u] from %s\n", wi.size(), dim, modelpath.c_str());

    uint32_t batchsize = config.nprocs * 100;
    uint32_t kmer_size = hash.get_kmer_size();
    PUnknownBar ubar("transform:");
    size_t count = 0;
    while (true)
    {
        std::vector<FastaRecord> v = reader.next(batchsize);
        count += v.size();
        if (v.empty())
        {
            break;
        }

        std::vector<Vector<float>> embed(v.size());

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
            Vector<float> vec(dim);
            transform(kmers, wi, vec);
            ubar.tick();
            embed[i] = vec;
        }
        for (auto &x : embed)
        {
            x.write_me(os);
        }
    }
    ubar.end();

    myinfo("Finished transforming %ld sequences", count);
}

template <typename BR>
void run(Config &config, BR &reader, rpns::CRandProj &hash)
{

    std::string filepath = config.output_prefix;
    if (config.zip_output)
    {
        filepath += ".gz";
    }

    if (config.zip_output)
    {
        ogzstream output(filepath.c_str());
        run(config, reader, hash, output);
        output.flush();
    }
    else
    {
        std::ofstream output(filepath);
        run(config, reader, hash, output);
        output.flush();
    }
}

void run(Config &config)
{
    omp_set_num_threads(config.nprocs);
    rpns::CRandProj hash;
    hash.load(config.hash_file);

    if (true)
    {

        if (config.is_fasta)
        {
            BatchReader<FastaTextReaderBase, FastaRecord> reader(config.inputpath);
            run(config, reader, hash);
        }
        else if (config.is_fastq)
        {
            BatchReader<FastqTextReaderBase, FastaRecord> reader(config.inputpath);
            run(config, reader, hash);
        }
        else
        {
            BatchReader<SeqTextReaderBase, FastaRecord> reader(config.inputpath);
            run(config, reader, hash);
        }
    }
}

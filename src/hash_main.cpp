
#include "argagg.hpp"
#include "utils.h"
#include "CRandProj.h"
#include "kmer.h"
#include "io.h"
#include "ProjConfig.h"

auto build_rp(bool is_fasta, bool is_fastq, const std::string &fasta_file_path, size_t kmer_size,
              size_t hash_size, size_t n_thread)
{
    spdlog::info("Build hash from {}", fasta_file_path);
    std::vector<FastaRecord> records;
    if (is_fasta)
    {
        records = read_fasta(fasta_file_path);
    }
    else if (is_fastq)
    {
        records = read_fastq(fasta_file_path);
    }
    else
    {
        records = read_seq_text(fasta_file_path);
    }
    size_t n_required = kmer_size * hash_size * 2;
    std::vector<std::string> kmers;
    while (kmers.size() < n_required)
    {
        size_t line_no = (size_t)(sparc::myrand::uniform<double>() * records.size());
        auto &seq = records.at(line_no).seq;
        if (seq.size() > kmer_size)
        {
            std::string kmer = random_generate_kmer(seq, kmer_size, true);
            if (kmer.find('N') == std::string::npos)
            {
                kmers.push_back(kmer);
            }
        }
    }

    rpns::RandProjBuilder builder(kmer_size, hash_size, n_thread);
    builder.create_hash(kmers);
    return builder.make_rp();
}

void show_help(const char *prog, argagg::parser &argparser)
{
    std::cout << prog << " v" << PROJECT_VERSION << "\n";

    std::cout << "Usage: " << prog << " [options] <input file>\n";
    std::cout << "Allowed options:\n";
    std::cout << argparser;
}

int main(int argc, char **argv)
{

    const char *program = argv[0];
    argagg::parser argparser{{

        {"help", {"-h", "--help"}, "shows this help message", 0},

        {"input", {"-i", "--input"}, "input fasta file", 1},
        {"output", {"-o", "--output"}, "output file path (default 'rp.bin')", 1},
        {"use_fasta", {"--fasta"}, "input are fasta files", 0},
        {"use_fastq", {"--fastq"}, "input are fastq files", 0},
        {"kmer_size", {
                          "-k",
                          "--kmer-size",
                      },
         "size of kmer",
         1},

        {"hash_size", {
                          "--hash-size",
                      },
         "number of bits of hash",
         1},

        {"n_thread", {
                         "--thread",
                     },
         "thread to use (default 0)",
         0},

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

    for (auto x : {"input", "kmer_size", "hash_size"})
        if (!args[x])
        {
            std::cerr << "ERROR: " << x << " is not set.\n";
            show_help(program, argparser);
            return EXIT_FAILURE;
        }
    uint32_t n_thread = args["n_thread"].as<uint32_t>(0);
    uint32_t kmer_size = args["kmer_size"].as<uint32_t>();
    uint32_t hash_size = args["hash_size"].as<uint32_t>();
    std::string input_file = args["input"].as<std::string>();
    std::string output_file = args["output"].as<std::string>("rp.bin");
    bool is_fasta = args["use_fasta"];
    bool is_fastq = args["use_fastq"];
    if (!sparc::file_exists(input_file.c_str()))
    {
        std::cerr << "input file does not exists: " << input_file << std::endl;
        return EXIT_FAILURE;
    }

    if (!is_fasta && (sparc::endswith(input_file, ".fa") || sparc::endswith(input_file, ".fa.gz") || sparc::endswith(input_file, ".fasta.gz") || sparc::endswith(input_file, ".fasta.gz")))
    {
        is_fasta = true;
    }

    if (!is_fastq && (sparc::endswith(input_file, ".fq") || sparc::endswith(input_file, ".fq.gz") || sparc::endswith(input_file, ".fastq.gz") || sparc::endswith(input_file, ".fastq.gz")))
    {
        is_fastq = true;
    }

    if (is_fasta && is_fastq)
    {
        std::cerr << "can not be both fasta and fastq" << std::endl;
        return EXIT_FAILURE;
    }

    if (!args.pos.empty())
    {
        std::cerr << "no positional args are allowed." << std::endl;
        return EXIT_FAILURE;
    }

    auto rp = build_rp(is_fasta, is_fastq, input_file, kmer_size, hash_size, n_thread);
    rp->save(output_file);
    return 0;
}

#include "argagg.hpp"
#include "utils.h"
#include "CRandProj.h"
#include "kmer.h"
#include "io.h"

auto build_rp(const std::string &fasta_file_path, size_t kmer_size,
              size_t hash_size, size_t n_thread)
{
    spdlog::info("Build hash from {}", fasta_file_path);
    auto records = read_fasta(fasta_file_path);
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

void show_help(argagg::parser &argparser)
{
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
        std::cerr << argparser;
        return EXIT_FAILURE;
    }

    if (args["help"])
    {
        show_help(argparser);

        return EXIT_SUCCESS;
    }

    for (auto x : {"input", "kmer_size", "hash_size"})
        if (!args[x])
        {
            std::cerr << "ERROR: " << x << " is not set.\n";
            show_help(argparser);
            return EXIT_FAILURE;
        }
    uint32_t n_thread = args["n_thread"].as<uint32_t>(0);
    uint32_t kmer_size = args["kmer_size"].as<uint32_t>();
    uint32_t hash_size = args["hash_size"].as<uint32_t>();
    std::string input_file = args["input"].as<std::string>();
    std::string output_file = args["output"].as<std::string>("rp.bin");

    if (!sparc::file_exists(input_file.c_str()))
    {
        std::cerr << "input file does not exists: " << input_file << std::endl;
        return EXIT_FAILURE;
    }

    if (!args.pos.empty())
    {
        std::cerr << "no positional args are allowed." << std::endl;
        return EXIT_FAILURE;
    }

    auto rp = build_rp(input_file, kmer_size, hash_size, n_thread);
    rp->save(output_file);
    return 0;
}
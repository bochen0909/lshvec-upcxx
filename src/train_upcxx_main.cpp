/*
 * make_kgraph.cpp
 *
 *  Created on: Jul 23, 2021
 *      Author: Bo Chen
 */

#include <string>
#include <exception>
#include <vector>
#include <algorithm> // std::random_shuffle
#include <random>
#include <ctime> // std::time
#include <cstdlib>
#include <string>
#include <iostream>
#include <unistd.h>
#include <upcxx/upcxx.hpp>

#include "argagg.hpp"
#include "CRandProj.h"
#include "utils.h"
#include "kmer.h"
#include "log.h"
#include "config.h"
#include "upchelper.h"

using namespace std;
using namespace sparc;

#define KMER_SEND_BATCH_SIZE (1000 * 1000)

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
	std::cout << "Usage: " << prog << " [options] file1, file2 ....\n";
	std::cout << "Allowed options:\n";
	std::cout << argparser;
}

void check_arg(argagg::parser_results &args, char *name)
{
	if (!args[name])
	{
		cerr << name << " is missing" << endl;
		exit(-1);
	}
}

int run(const std::vector<std::string> &input, Config &config);

int main(int argc, char **argv)
{
	upcxx::init();

	Config config;
	config.program = argv[0];
	config.rank = upcxx::rank_me();
	config.nprocs = upcxx::rank_n();
	config.mpi_hostname = sparc::get_hostname();
	set_spdlog_pattern(config.mpi_hostname.c_str(), config.rank);

	if (config.rank == 0)
	{
		myinfo("Welcome to LSHVec!");
	}
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
		std::cerr << "no input files are provided" << endl;
		return EXIT_FAILURE;
	}
	else
	{
		config.inputpath = args.all_as<string>();
		bool b_error = false;
		for (size_t i = 0; i < config.inputpath.size(); i++)
		{
			if (!dir_exists(config.inputpath.at(i).c_str()))
			{
				cerr << "Error, input dir does not exists:  "
					 << config.inputpath.at(i) << endl;
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
	std::string outputpath = config.output_prefix;

	upcxx::barrier();
	if (upcxx::rank_me() == 0)
	{
		if (dir_exists(outputpath.c_str()))
		{
			cerr << "Error, output dir exists:  " << outputpath << endl;
			return EXIT_FAILURE;
		}
	}

	if (upcxx::rank_me() == 0)
	{
		if (make_dir(outputpath.c_str()) < 0)
		{
			cerr << "Error, mkdir dir failed for " << outputpath << endl;
			return EXIT_FAILURE;
		}
	}

	std::vector<std::string> myinput = get_my_files(config.inputpath);
	myinfo("#of my input: %ld", myinput.size());
	run(myinput, config);

	return 0;
}

rpns::CRandProj g_hash;
int run(const std::vector<std::string> &input, Config &config)
{

	if (config.rank == 0)
	{
		config.print();
	}

	if (true)
	{
		g_hash.load(config.hash_file);
		if (!g_hash.is_defined())
		{
			std::cerr << "load rp hash failed\n";
			return EXIT_FAILURE;
		}
		else
		{
			if (config.rank == 0)
			{
				myinfo("kmer_size=%ld hash_size=%ld", g_hash.get_kmer_size(), g_hash.get_hash_size());
			}
		}
	}



	for (uint32_t i = 1; i < config.epoch + 1; i++)
	{
	}

	upcxx::barrier();

	//dump_vector(config.get_my_output(false), config.rank);

	upcxx::barrier();

	return 0;
}

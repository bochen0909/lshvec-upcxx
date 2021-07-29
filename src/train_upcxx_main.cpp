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
#include "model_upcxx.h"
#include "pbar.h"
#include "ProjConfig.h"

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

		{"output", {"-o", "--output"}, "output folder (default 'output')", 1},
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

	config.outputpath = args["output"].as<std::string>("output");

	std::string outputpath = config.outputpath;

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

template <class BR>
void run_epoch(uint32_t this_epoch, Config &config, BR &reader, UPCXXModel<float> &model, rpns::CRandProj &hash, float learning_rate)
{
	if (config.rank == 0)
	{
		myinfo("Start epoch %ld, learning_rate=%.6f", this_epoch + 1, learning_rate);
	}

	uint32_t batchsize = 10;
	uint32_t kmer_size = hash.get_kmer_size();
	bool update_wc = this_epoch == 0;
	float sum_loss = 0;
	size_t num_of_seq = 0;
	char msg[128];
	sprintf(msg, "epoch %u:", this_epoch);
	PUnknownBar *ubar = 0;
	PBar *bar = 0;
	if (config.rank == 0)
	{
		ubar = new PUnknownBar(msg);
		bar = new PBar(msg, config.num_seq);
	}
	while (true)
	{
		std::vector<FastaRecord> v = reader.next(batchsize);
		sparc::shuffle(v);
		num_of_seq += v.size();
		if (v.empty())
		{
			break;
		}

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
			if (config.rank == 0)
			{
				if (this_epoch == 0)
				{
					ubar->tick();
				}
				else
				{
					bar->tick();
				}
			}
			sum_loss += loss;
		}
	}
	if (config.rank == 0)
	{
		if (this_epoch == 0)
		{
			ubar->end();
		}
		else
		{
			bar->end();
		}
		delete ubar;
		delete bar;
	}

	upcxx::barrier();

	float mean_loss = num_of_seq == 0 ? 0 : sum_loss / num_of_seq;
	int all_sum_loss = upcxx::reduce_all(mean_loss, upcxx_op_add).wait();

	if (this_epoch == 0)
	{
		config.num_seq = num_of_seq;
	}
	if (config.rank == 0)
	{
		if (this_epoch == 0)
		{
			myinfo("Found that number of sequences is %ld", num_of_seq);
		}
		myinfo("End epoch %ld, loss=%f", this_epoch + 1, mean_loss);
	}
}

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

	uint32_t num_word = (uint32_t)(1l << g_hash.get_hash_size());
	UPCXXModel<float> model(num_word, config.rank, config.nprocs, config.dim, config.neg_size, config.use_cbow, config.half_window);
	model.uniform_init();

	for (uint32_t i = 0; i < config.epoch; i++)
	{
		float learning_rate = config.learning_rate * (config.epoch - i) / config.epoch;

		if (config.is_fasta)
		{
			BatchReader<FastaTextReaderBase, FastaRecord> reader(input);
			run_epoch(i, config, reader, model, g_hash, learning_rate);
		}
		else if (config.is_fastq)
		{
			BatchReader<FastqTextReaderBase, FastaRecord> reader(input);
			run_epoch(i, config, reader, model, g_hash, learning_rate);
		}
		else
		{
			BatchReader<SeqTextReaderBase, FastaRecord> reader(input);
			run_epoch(i, config, reader, model, g_hash, learning_rate);
		}
	}

	upcxx::barrier();

	std::string save_prefix = config.get_my_output() + ".vec.bin";
	myinfo("saving vectors to %s", save_prefix.c_str());
	model.dump_vector_bin(save_prefix, config.zip_output);
	save_prefix = config.get_my_output() + ".wc.txt";
	myinfo("saving word counts to %s", save_prefix.c_str());
	model.dump_word_counts(save_prefix, config.zip_output);
	upcxx::barrier();
	if (config.rank == 0)
	{
		myinfo("Finished!");
	}
	upcxx::finalize();
	return 0;
}

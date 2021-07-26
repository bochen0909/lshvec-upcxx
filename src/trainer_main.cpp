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
#include "utils.h"
#include "kmer.h"
#include "log.h"
#include "config.h"
#include "distr_kgraph_graph.h"
#include "upchelper.h"

using namespace std;
using namespace sparc;

#define KMER_SEND_BATCH_SIZE (1000 * 1000)

DistrGraph *g_graph = 0;

size_t g_n_sent = 0;

struct Config : public BaseConfig
{
	float learning_rate;
	int epoch;
	int dim;
	std::string hash_file;

	void print()
	{
		BaseConfig::print();
		myinfo("config: learning_rate=%f", learning_rate);
		myinfo("config: epoch=%ld", epoch);
		myinfo("config: dimension=%ld", dim);
		myinfo("config: hash_file=%s", hash_file.c_str());
	}
};

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
	argagg::parser argparser{{

		{"help", {"-h", "--help"}, "shows this help message", 0},

		{"zip_output", {"-z", "--zip"}, "zip output files", 0},

		{"output", {"-o", "--output"}, "output folder (default './out')", 1},
		{"hash_file", {"--hash-file"}, "hash file to use", 1},
		{"learning_rate", {
							  "--lr",
						  },
		 "initial learning rate (default 0.3)",
		 0},

		{"dim", {
					"--dim",
				},
		 "number of dimension (default 300)",
		 0},
		{"epoch", {
					  "--epoch",
				  },
		 "number of epochs to train (default 100)",
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
		return EXIT_FAILURE;
	}

	if (args["help"])
	{
		std::cerr << argparser;
		return EXIT_SUCCESS;
	}

	config.dim = args["dim"].as<uint32_t>(300);
	config.learning_rate = args["learning_rate"].as<float>(0.3);
	config.epoch = args["epoch"].as<uint32_t>(100);

	config.zip_output = args["zip_output"];
	config.hash_file = args["hash_file"].as<string>();

	if (!sparc::file_exists(config.hash_file))
	{
		std::cerr << "hash file does not exists: " << config.hash_file << endl;
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

	string outputpath = args["output"].as<string>();
	config.outputpath = outputpath;

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

int process_graph_file(const std::string &filepath)
{

	return 0;
}

CRandProj g_hash;

int run(const std::vector<std::string> &input, Config &config)
{

	if (config.rank == 0)
	{
		config.print();
	}

	if (true)
	{
		g_hash.load(config.hash_file);
		if (!g_hash.is_dfined())
		{
			std::cerr < "load rp hash failed\n";
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



	g_graph = new DistrGraph(config.k_neigbhor);

	for(uint32_t i=1;i<config.epoch+1;i++){

	}
	
	
	upcxx::barrier();


	g_graph->dump_model(config.get_my_output(), config.rank);

	delete g_graph;

	return 0;
}

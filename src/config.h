/*
 * config.h
 *
 *  Created on: Jul 28, 2020
 *      Author:
 */

#ifndef SOURCE_DIRECTORY__SRC_SPARC_BASECONFIG_H_
#define SOURCE_DIRECTORY__SRC_SPARC_BASECONFIG_H_

#include <vector>
#include <string>
#include "log.h"

struct BaseConfig
{

	std::string program;
	std::vector<std::string> inputpath;
	std::string outputpath;
	std::string scratch_dir;
	std::string mpi_hostname;
	std::string mpi_ipaddress;
	int port;
	int rank;
	int nprocs;
	bool zip_output;
	std::vector<int> peers_ports;
	std::vector<int> hash_rank_mapping;
	std::vector<std::string> peers_hosts;

	std::string backend = "upcxx";

	void print()
	{
		myinfo("config: program=%s", program.c_str());
		myinfo("config: backend=%s", backend.c_str());
		myinfo("config: zip_output=%d", zip_output ? 1 : 0);
		for (size_t i = 0; i < inputpath.size(); i++)
		{
			myinfo("config: inputpath[%d]=%s", i, inputpath.at(i).c_str());
		}
		myinfo("config: outputpath=%s", outputpath.c_str());
		myinfo("config: scratch_dir=%s", scratch_dir.c_str());
		myinfo("config: #procs=%d", nprocs);
	}
	std::string get_dbpath()
	{
		return get_dbpath(0);
	}
	std::string get_dbpath(int h)
	{
		char tmp[2048];
		sprintf(tmp, "%s/part_h%d_r%d.db", scratch_dir.c_str(), h, rank);
		return tmp;
	}
	std::string get_my_output()
	{
		char tmp[2048];
		sprintf(tmp, "%s/part_r%06d", outputpath.c_str(), rank);
		return tmp;
	}

	int get_my_port()
	{
		return get_my_port(0);
	}
	int get_my_port(int h)
	{
		return port + h * nprocs + rank;
	}
};

#endif /* SOURCE_DIRECTORY__SRC_SPARC_CONFIG_H_ */

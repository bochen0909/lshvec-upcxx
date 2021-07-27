#include <exception>
#include <stdexcept>
#include <string>
#include <algorithm>

#include "io.h"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace std;

TEST_CASE("read fasta file", "[fasta]")
{
	for (auto &filepath : {"knucleotide.fasta", "knucleotide.fasta.gz"})
	{
		auto v = read_fasta(filepath);
		REQUIRE(v.size() == 3);
		REQUIRE(v[0].seq.size() == 50000);
		REQUIRE(v[1].seq.size() == 75000);
		REQUIRE(v[2].seq.size() == 125000);
	}
}

TEST_CASE("read fasta using stream", "[fasta]")
{
	for (auto &filepath : {"sample.fa", "sample.fa.gz"})
	{
		FastaTextReaderBase reader(filepath);
		std::vector<FastaRecord> v;
		while (true)
		{
			FastaRecord r = reader.next();
			if (r.id.empty())
			{
				break;
			}
			v.push_back(r);
		}
		REQUIRE(v.size() == 2);
		REQUIRE(v[0].seq.size() == 351);
		REQUIRE(v[1].seq.size() == 300);
	}
}

TEST_CASE("read fasta using batch", "[fasta]")
{

	BatchReader<FastaTextReaderBase, FastaRecord> reader({"sample.fa", "sample.fa.gz"});

	std::vector<FastaRecord> V;
	while (true)
	{
		std::vector<FastaRecord> v = reader.next(1);
		if (v.empty())
		{
			break;
		}
		V.insert(V.end(), v.begin(), v.end());
	}
	REQUIRE(V.size() == 4);
	REQUIRE(V[0].seq.size() == 351);
	REQUIRE(V[1].seq.size() == 300);
	REQUIRE(V[2].seq.size() == 351);
	REQUIRE(V[3].seq.size() == 300);
}


TEST_CASE("read fastq using stream", "[fastq]")
{
	for (auto &filepath : {"sample.fq", "sample.fq.gz"})
	{
		FastqTextReaderBase reader(filepath);
		std::vector<FastaRecord> v;
		while (true)
		{
			FastaRecord r = reader.next();
			if (r.id.empty())
			{
				break;
			}
			v.push_back(r);
		}
		REQUIRE(v.size() == 2);
		REQUIRE(v[0].seq.size() == 351);
		REQUIRE(v[1].seq.size() == 300);
	}
}

TEST_CASE("read fastq using batch", "[fastq]")
{

	BatchReader<FastqTextReaderBase, FastaRecord> reader({"sample.fq", "sample.fq.gz"});

	std::vector<FastaRecord> V;
	while (true)
	{
		std::vector<FastaRecord> v = reader.next(1);
		if (v.empty())
		{
			break;
		}
		V.insert(V.end(), v.begin(), v.end());
	}
	REQUIRE(V.size() == 4);
	REQUIRE(V[0].seq.size() == 351);
	REQUIRE(V[1].seq.size() == 300);
	REQUIRE(V[2].seq.size() == 351);
	REQUIRE(V[3].seq.size() == 300);
}


TEST_CASE("read seq text using stream", "[seq text]")
{
	for (auto &filepath : {"sample.txt", "sample.txt.gz"})
	{
		SeqTextReaderBase reader(filepath);
		std::vector<std::string> v;
		while (true)
		{
			std::string r = reader.next();
			if (r.empty())
			{
				break;
			}
			v.push_back(r);
		}
		REQUIRE(v.size() == 2);
		REQUIRE(v[0].size() == 351);
		REQUIRE(v[1].size() == 300);
	}
}

TEST_CASE("read seq text using batch", "[seq text]")
{

	BatchReader<SeqTextReaderBase, std::string> reader({"sample.txt", "sample.txt.gz"});

	std::vector<std::string> V;
	while (true)
	{
		std::vector<std::string> v = reader.next(1);
		if (v.empty())
		{
			break;
		}
		V.insert(V.end(), v.begin(), v.end());
	}
	REQUIRE(V.size() == 4);
	REQUIRE(V[0].size() == 351);
	REQUIRE(V[1].size() == 300);
	REQUIRE(V[2].size() == 351);
	REQUIRE(V[3].size() == 300);
}
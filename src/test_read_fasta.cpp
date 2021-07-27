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

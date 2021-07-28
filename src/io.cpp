
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include "utils.h"
#include "gzstream.h"

#include "io.h"

template <class STREAM>
void read_fasta(STREAM &file, std::vector<FastaRecord> &recrods)
{
	std::string line;
	FastaRecord r;
	while (std::getline(file, line))
	{
		sparc::trim(line);
		if (line[0] == '>')
		{
			if (!r.id.empty())
			{
				recrods.push_back(r);
			}
			r.id = line;
			r.seq.clear();
		}
		else
		{ //seq
			transform_seq(line);
			if (r.id.empty())
			{
				throw std::runtime_error("found fasta data without finding header first");
			}
			r.seq += line;
		}
	}
	if (!r.id.empty())
	{
		recrods.push_back(r);
	}
}

std::vector<FastaRecord> read_fasta(const std::string &filepath)
{
	if (!sparc::file_exists(filepath.c_str()))
	{
		throw std::runtime_error(filepath + " does not  exist");
	}
	std::vector<FastaRecord> records;

	if (sparc::endswith(filepath, ".gz"))
	{
		igzstream file(filepath.c_str());
		read_fasta(file, records);
	}
	else
	{
		std::ifstream file(filepath);
		read_fasta(file, records);
	}
	return records;
}

std::vector<FastaRecord> read_fastq(const std::string &filepath)
{
	if (!sparc::file_exists(filepath.c_str()))
	{
		throw std::runtime_error(filepath + " does not  exist");
	}

	FastqTextReaderBase reader(filepath);
	std::vector<FastaRecord> V;
	while (true)
	{
		FastaRecord r = reader.next();
		if (r.id.empty())
		{
			break;
		}
		V.push_back(r);
	}
	return V;
}

std::vector<FastaRecord> read_seq_text(const std::string &filepath)
{
	if (!sparc::file_exists(filepath.c_str()))
	{
		throw std::runtime_error(filepath + " does not  exist");
	}

	SeqTextReaderBase reader(filepath);
	std::vector<FastaRecord> V;
	while (true)
	{
		FastaRecord r = reader.next();
		if (r.id.empty())
		{
			break;
		}
		V.push_back(r);
	}
	return V;
}

template <>
bool record_is_valid(const FastaRecord &t)
{
	return !t.id.empty();
}

template <>
bool record_is_valid(const std::string &t)
{
	return !t.empty();
}
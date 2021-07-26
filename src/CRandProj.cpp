/*
 * CRandProj.cpp
 *
 *  Created on: Feb 19, 2019
 *      Author: bo
 */

#include <iostream>
#include <map>
#include <algorithm>
#include <fstream>
#include <complex>
#include "spdlog/spdlog.h"
#include "CRandProj.h"
using namespace std;

namespace rpns
{
	int CRandProj::magic_number = 0x201902;

	CRandProj::CRandProj()
	{
		// TODO Auto-generated constructor stub
	}
	CRandProj::CRandProj(int hash_size, int kmer_size)
	{
		spdlog::info("CRandProj: hash_size={}, kmer_size={}", hash_size, kmer_size);
		this->num_spokes = hash_size;
		this->kmer_size = kmer_size;
		this->pow2 = new vector<unsigned long>(num_spokes);
		unsigned long one = 1;
		for (int i = 0; i < num_spokes; i++)
		{
			(*pow2)[i] = (one << i);
		}
	}

	void CRandProj::clear()
	{
		num_spokes = -1;
		kmer_size = -1;
		if (pow2)
			delete pow2;
		if (Wc)
			delete Wc;
		if (Wpt)
			delete Wpt;
		pow2 = NULL;
		Wc = NULL;
		Wpt = NULL;
	}

	CRandProj::~CRandProj()
	{
		this->clear();
	}

	void CRandProj::seq_encoding(const string &seq, vector<complex<double>> &mat)
	{
		static map<char, complex<double>> ltc_q = {{'A', std::complex<double>(-1, 0)}, {'T', std::complex<double>(1, 0)}, {'C', std::complex<double>(0, -1)}, {'G', std::complex<double>(0, 1)}};
		mat.resize(seq.size());
		for (size_t i = 0; i < seq.size(); i++)
		{
			complex<double> c = ltc_q[seq[i]];
			mat[i] = c;
		}
	}

	int CRandProj::set_Wc(const vector<double> &re, const vector<double> &im)
	{
		if ((size_t)num_spokes != re.size() || (size_t)num_spokes != im.size())
		{
			cerr << "Error!, size not match" << endl;
			return -1;
		}
		Wc = new vector<complex<double>>(num_spokes);
		for (int i = 0; i < num_spokes; i++)
		{
			(*Wc)[i] = re[i] + im[i] * std::complex<double>(0, 1);
		}
		return 0;
	}

	int CRandProj::set_Wpt(const vector<double> &re, const vector<double> &im)
	{ // kxb in row format
		if ((size_t)(num_spokes * kmer_size) != re.size() || (size_t)(num_spokes * kmer_size) != im.size())
		{
			cerr << "Error!, size not match" << endl;
			return -1;
		}
		Wpt = new vector<vector<complex<double>>>(kmer_size);
		size_t idx = 0;
		for (int i = 0; i < kmer_size; i++)
		{
			vector<complex<double>> row(num_spokes);
			for (int j = 0; j < num_spokes; j++)
			{
				row[j] = re[idx] + im[idx] * std::complex<double>(0, 1);
				idx++;
			}
			(*Wpt)[i] = row;
		}
		return 0;
	}

	bool CRandProj::is_defined()
	{
		if (num_spokes > 0 && kmer_size > 0 && pow2 && Wc && Wpt)
		{
			if (Wc->size() == (size_t)num_spokes && Wpt->size() == (size_t)kmer_size && Wpt->at(0).size() == (size_t)num_spokes)
			{
				return true;
			}
		}
		return false;
	}

	std::vector<std::complex<double>> &CRandProj::multiply(
		const std::vector<std::complex<double>> &x,
		const vector<vector<complex<double>>> &D,
		std::vector<std::complex<double>> &y)
	{

		//x 1xm, D mxn, y 1xn
		y.resize(D.at(0).size());

		for (size_t i = 0; i < y.size(); i++)
		{
			complex<double> sum = 0;
			for (size_t j = 0; j < x.size(); j++)
			{
				sum += x[j] * D[j][i];
			}
			y[i] = sum;
		}
		return y;
	}

	void CRandProj::conjugate(std::vector<std::complex<double>> &x)
	{
		for (size_t i = 0; i < x.size(); i++)
		{
			x[i] = std::conj(x[i]);
		}
	}

	int CRandProj::get_kmer_size()
	{
		return this->kmer_size;
	}
	int CRandProj::get_hash_size()
	{
		return this->num_spokes;
	}

	unsigned long CRandProj::hash(const string &kmer, bool reverse_compliments)
	{
		if (kmer.size() != (size_t)kmer_size)
		{
			cerr << "kmer length is not right: " << kmer << endl;
			return -1;
		}

		std::vector<std::complex<double>> C;
		this->seq_encoding(kmer, C);

		std::vector<std::complex<double>> L;
		this->multiply(C, *Wpt, L);

		unsigned long B = 0;
		for (int i = 0; i < num_spokes; i++)
		{
			if (L[i].real() > (*Wc)[i].real())
			{
				B += (*pow2)[i];
			}
		}

		if (reverse_compliments)
		{
			unsigned long B1 = B;
			std::reverse(C.begin(), C.end());
			for (size_t i = 0; i < C.size(); i++)
			{
				C[i] *= -1;
			}
			this->multiply(C, *Wpt, L);
			unsigned long B = 0;
			for (int i = 0; i < num_spokes; i++)
			{
				if (L[i].real() > (*Wc)[i].real())
				{
					B += (*pow2)[i];
				}
			}
			return B1 < B ? B1 : B;
		}
		else
		{
			return B;
		}
	}

#define READ_NUMBER(NUM, STREAM) STREAM.read((char *)&(NUM), sizeof((NUM)))
	int CRandProj::load(const std::string &path)
	{
		ifstream file(path.c_str(), ios::in | ios::binary);
		int mn = 0;
		READ_NUMBER(mn, file);
		if (mn != magic_number)
		{
			cerr << path << " is not a hash file" << endl;
			return -1;
		}
		this->clear();

		READ_NUMBER(num_spokes, file);
		READ_NUMBER(kmer_size, file);

		pow2 = new vector<unsigned long>(num_spokes);
		unsigned long n;
		for (int i = 0; i < num_spokes; i++)
		{
			READ_NUMBER(n, file);
			(*pow2)[i] = n;
		}
		complex<double> cd;
		Wc = new vector<complex<double>>(num_spokes);
		for (int i = 0; i < num_spokes; i++)
		{
			READ_NUMBER(cd, file);
			(*Wc)[i] = cd;
		}
		Wpt = new vector<vector<complex<double>>>(kmer_size);
		for (int i = 0; i < kmer_size; i++)
		{
			vector<complex<double>> vec(num_spokes);
			for (int j = 0; j < num_spokes; j++)
			{
				READ_NUMBER(cd, file);
				vec[j] = cd;
			}
			(*Wpt)[i] = vec;
		}
		file.close();

		return 0;
	}

#define WRITE_NUMBER(NUM, STREAM) STREAM.write((char *)&(NUM), sizeof((NUM)))
	int CRandProj::save(const std::string &path)
	{
		ofstream file(path.c_str(), ios::out | ios::binary | ios::trunc);
		if (file.is_open())
		{
			WRITE_NUMBER(magic_number, file);
			WRITE_NUMBER(num_spokes, file);
			WRITE_NUMBER(kmer_size, file);

			for (int i = 0; i < num_spokes; i++)
			{
				WRITE_NUMBER(pow2->at(i), file);
			}
			for (int i = 0; i < num_spokes; i++)
			{
				WRITE_NUMBER(Wc->at(i), file);
			}

			for (int i = 0; i < kmer_size; i++)
			{
				for (int j = 0; j < num_spokes; j++)
				{
					WRITE_NUMBER((*Wpt)[i][j], file);
				}
			}
			file.flush();
			file.close();
		}
		else
			cerr << "Unable to open file: " << path << endl;
		return 0;
	}

	std::ostream &operator<<(std::ostream &os, const wsPC &o)
	{
		os << "w: " << o.w << std::endl;
		os << "s: " << o.s << std::endl;
		os << "C: " << o.C << std::endl;
		os << "P: ";
		for (auto &x : o.P)
		{
			os << x << "\t";
		}
		os
			<< std::endl;
		return os;
	}
}
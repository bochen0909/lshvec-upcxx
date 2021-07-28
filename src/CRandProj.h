/*
 * CRandProj.h
 *
 *  Created on: Feb 19, 2019
 *      Author: bo
 */

#ifndef SRC_CRANDPROJ_H_
#define SRC_CRANDPROJ_H_

#include <string>
#include <vector>
#include <complex>
#include <memory>
#include <unordered_map>
#include <iostream>
#include <omp.h>
#include "Eigen/Dense"
#include "spdlog/spdlog.h"
#include "utils.h"
namespace rpns
{
	class RandProjBuilder;

	class CRandProj
	{
		friend class RandProjBuilder;

	public:
		CRandProj();
		CRandProj(int hash_size, int kmer_size);
		virtual ~CRandProj();

		int save(const std::string &path);
		int load(const std::string &path);

		int set_Wc(const std::vector<double> &re, const std::vector<double> &im);
		int set_Wpt(const std::vector<double> &re, const std::vector<double> &im);

		unsigned long hash(const std::string &, bool reverse_compliments);
		inline void conjugate(std::vector<std::complex<double>> &x);
		int get_kmer_size();
		int get_hash_size();
		static void seq_encoding(const std::string &seq, std::vector<std::complex<double>> &mat);
		bool is_defined();

	protected:
		void clear();
	

		std::vector<std::complex<double>> &multiply(
			const std::vector<std::complex<double>> &x,
			const std::vector<std::vector<std::complex<double>>> &D,
			std::vector<std::complex<double>> &y);

	protected:
		std::vector<unsigned long> *pow2 = NULL;
		int num_spokes = 0;
		int kmer_size = 0;
		std::vector<std::complex<double>> *Wc = NULL;
		std::vector<std::vector<std::complex<double>>> *Wpt = NULL;

		static int magic_number;
	};

	struct wsPC
	{
		int w; //wheel id
		int s; // bit pos
		std::vector<std::complex<double>> P;
		std::complex<double> C;
	};

	std::ostream &operator<<(std::ostream &os, const wsPC &o);

	class RandProjBuilder
	{
	protected:
		int kmer_size;
		int hash_size;
		uint32_t n_thread;
		int n_wheels = 1;
		std::vector<wsPC> Wheels;
		std::vector<std::vector<size_t>> __pseduo_index;

	public:
		RandProjBuilder(uint32_t kmer_size, uint64_t hash_size, uint32_t n_thread) : kmer_size(kmer_size), hash_size(hash_size), n_thread(n_thread)
		{
			if (n_thread < 1)
			{
				this->n_thread = sparc::get_number_of_thread();
			}
			omp_set_num_threads(this->n_thread);

			spdlog::info("RandProjBuilder: hash_size={}, kmer_size={}, n_thread={}", hash_size, kmer_size, this->n_thread);
		}

		void print()
		{
			for (auto &x : Wheels)
			{
				std::cout << x << std::endl;
			}
		}

		void __set_pseduo_index(std::vector<std::vector<size_t>> &pseduo_index)
		{
			this->__pseduo_index = pseduo_index;
		}

		void create_hash(const std::vector<std::string> &kmers)
		{
			std::vector<std::vector<std::complex<double>>> encoded_kmers(kmers.size());
			encoded_kmers.resize(kmers.size());
			for (size_t i = 0; i < kmers.size(); i++)
			{
				CRandProj::seq_encoding(kmers.at(i), encoded_kmers[i]);
			}
			_set_wheels(encoded_kmers);
		}
		std::shared_ptr<CRandProj> make_rp()
		{
			std::shared_ptr<CRandProj> rp = std::make_shared<CRandProj>(this->hash_size, this->kmer_size);

			std::vector<std::complex<double>> *Wc = new std::vector<std::complex<double>>();
			std::vector<std::vector<std::complex<double>>> *Wpt = new std::vector<std::vector<std::complex<double>>>();

			for (size_t i = 0; i < Wheels.size(); i++)
			{
				Wc->push_back(Wheels.at(i).C);
			}
			{

				auto m = Wheels.at(0).P.size();
				auto n = Wheels.size();
				Wpt->resize(m);
				for (auto &x : *Wpt)
				{
					x.resize(n);
				}
				for (size_t i = 0; i < n; i++)
				{
					auto &P = Wheels.at(i).P;
					for (size_t j = 0; j < m; j++)
					{
						(*Wpt)[j][i] = P.at(j);
					}
				}
			}

			rp->Wc = Wc;
			rp->Wpt = Wpt;

			return rp;
		};

	protected:
		auto _affine_hull(std::vector<std::vector<std::complex<double>>> &linear_system)
		{
			// linear_system: d dimensions of n docs in this hyperplane
			size_t n = linear_system.size();
			Eigen::MatrixXcd M(n + 1, n + 1);
			for (size_t i = 0; i < n; i++)
			{
				for (size_t j = 0; j < n; j++)
				{
					M(i, j) = linear_system.at(i).at(j);
				}
			}
			for (size_t i = 0; i < n; i++)
			{
				M(i, n) = -1;
			}
			for (size_t i = 0; i < n + 1; i++)
			{
				M(n, i) = 0;
			}

			Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
			auto V = svd.matrixV().conjugate().transpose();
			if (!__pseduo_index.empty())
			{
				auto &S = svd.singularValues();
				auto &U = svd.matrixU();
				std::cout << "M:\n"
						  << M << "\n"
						  << std::endl;
				std::cout << "U:\n"
						  << U << "\n"
						  << std::endl;
				std::cout << "S:\n"
						  << S << "\n"
						  << std::endl;

				std::cout << "V^*^T:\n"
						  << V << "\n"
						  << std::endl;

				auto &USV = U * S.asDiagonal() * V;
				std::cout << "USV^*^T:\n"
						  << USV << "\n"
						  << std::endl;
			}

			std::vector<std::complex<double>> P;
			for (size_t i = 0; i < n; i++)
			{
				P.push_back(V(n, i));
			}
			std::complex<double> C = V(n, n);

			if (!__pseduo_index.empty())
			{
				std::cout << "C:\n"
						  << C << "\n"
						  << std::endl;
				std::cout << "P:\n";
				for (auto &x : P)
				{
					std::cout << x << "\t";
				}
				std::cout
					<< "\n"
					<< std::endl;
			}
			return std::make_tuple(P, C);
		}

		auto _one_wheel(int w, const std::vector<std::vector<std::complex<double>>> &raw_kmers)
		{
			std::vector<wsPC> S;
			S.resize(hash_size);

#pragma omp parallel for
			for (int s = 0; s < hash_size; s++)
			{
				std::vector<std::vector<std::complex<double>>> kmers(raw_kmers.begin(), raw_kmers.end());
				std::vector<std::vector<std::complex<double>>> L;
				if (this->__pseduo_index.empty())
				{
					sparc::myrand::shuffle(kmers);
					L.insert(L.begin(), kmers.begin(), kmers.begin() + kmer_size);
				}
				else
				{
					for (auto &i : __pseduo_index.at(s))
					{
						L.push_back(kmers.at(i));
					}
				}
				std::vector<std::complex<double>> P;
				std::complex<double> C;
				wsPC wspc;
				wspc.w = w;
				wspc.s = s;
				std::tie(wspc.P, wspc.C) = _affine_hull(L);
				S[s] = wspc;
			}
			return S;
		}

		void _set_wheels(std::vector<std::vector<std::complex<double>>> &kmers)
		{
			for (int i = 0; i < n_wheels; i++)
			{
				auto wheel = _one_wheel(i, kmers);
				Wheels.insert(Wheels.end(), wheel.begin(), wheel.end());
			}
		}

		unsigned long fun_pow2(unsigned long i)
		{
			return i == 0 ? 1 : fun_pow2(i - 1) * 2;
		}
	};
};	   //end of namespace
#endif /* SRC_CRANDPROJ_H_ */

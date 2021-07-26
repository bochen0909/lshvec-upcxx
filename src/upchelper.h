/*
 * UPCHELPER.h
 *
 *  Created on: Jul 23, 2020
 *      Author:
 */

#ifndef SOURCE_DIRECTORY__SRC_UPCHELPER_H_
#define SOURCE_DIRECTORY__SRC_UPCHELPER_H_

#include "upcxx_comp.h"
#include "limits.h"
#include "utils.h"
#include "log.h"
#include "kmer.h"

namespace sparc {
bool check_all_peers_on_same_page(const std::vector<std::string> &files) {

	int rank = upcxx::rank_me();
	{ //check #
		int n = 0;
		if (rank == 0) {
			n = (int) files.size();
		}
		upcxx::future<int> fut = upcxx::broadcast(n,
				0/* broadcast from rank 0*/);

		int bcast_n = fut.wait();

		if (bcast_n != (int) files.size()) {
			myerror("#files is not equal to that of rank 0: r=%d, %ld<>%ld",
					rank, files.size(), n);
			return false;
		}
	}
	size_t N = files.size();
	for (size_t i = 0; i < N; i++) { //check files
		std::string fpath = files.at(i);
		std::string bcast_fpath = upcxx_broadcast_nontrivial<std::string>(std::move(fpath),
				0).wait();
		if (std::find(files.begin(), files.end(), bcast_fpath) == files.end()) {
			myerror("%s cannot find locally for rank %d", bcast_fpath.c_str(),
					rank);
			return false;
		}
	}
	return true;
}
bool check_splitted_files(const std::vector<std::vector<std::string>> &myfiles,
		const std::vector<std::string> &allfiles) {
	{
		int n = 0;
		for (auto &v : myfiles) {
			n += (int) v.size();

		}

		int N = upcxx::reduce_all(n, upcxx_op_add).wait();

		if (upcxx::rank_me() == 0 && allfiles.size() != N) {
			myerror("check_splitted_files, %ld<>%ld", allfiles.size(), N);
			upcxx_fatal_error("check_splitted_files failed");
			return false;
		}
	}

	upcxx::barrier();

	for (size_t i = 0; i < allfiles.size(); i++) { //check files
		auto fpath = allfiles.at(i);
		std::string bcast_fpath = upcxx_broadcast_nontrivial<std::string>(std::move(fpath),
				0/* broadcast from rank 0*/).wait();

		int nfound = 0;
		for (auto &v : myfiles) {
			if (std::find(v.begin(), v.end(), bcast_fpath) != v.end()) {
				nfound++;
			}
		}

		int N = upcxx::reduce_all(nfound, upcxx_op_add).wait();

		if (upcxx::rank_me() == 0 && 1 != N) {
			myerror("check_splitted_files # of [%s] is %d<>1",
					bcast_fpath.c_str(), N);
			return false;
		}

	}
	upcxx::barrier();
	return true;
}
void get_all_files(const std::vector<std::string> &folders,
		std::vector<std::string> &ret) {

	for (size_t i = 0; i < folders.size(); i++) {
		auto v = list_dir(folders.at(i).c_str());
		ret.insert(ret.end(), v.begin(), v.end());
	}
	sort(ret.begin(), ret.end());
	ret.erase(unique(ret.begin(), ret.end()), ret.end());
}

bool get_all_files_check(const std::vector<std::string> &folders,
		std::vector<std::string> &ret) {
	std::vector<std::string> files;
	bool good = false;
	for (int i = 0; i < 3; i++) {
		ret.clear();
		get_all_files(folders, ret);
		good = check_all_peers_on_same_page(ret);
		if (good) {
			break;
		} else {
			myerror("Check files failed for iteration %d", i);
		}
	}

	return good;

}

std::vector<std::vector<std::string>> get_my_files(
		const std::vector<std::string> &folders, int n_bucket) {
	std::vector<std::vector<std::string>> myinput;
	std::vector<std::string> allfiles;
	if (get_all_files_check(folders, allfiles)) {
		if (upcxx::rank_me() == 0) {
			myinfo("#of all files: %ld", allfiles.size());
		}
		if (allfiles.empty()) {
			myerror("no input files found");
			upcxx_fatal_error("no input files found");
		}

		std::vector<std::string> localfiles;

		for (size_t j = 0; j < allfiles.size(); j++) {
			auto &file = allfiles.at(j);
			if (j % upcxx::rank_n() == upcxx::rank_me()) {
				localfiles.push_back(file);
			}
		}
		for (int i = 0; i < n_bucket; i++) {
			std::vector<std::string> thisinput;
			for (size_t j = 0; j < localfiles.size(); j++) {
				auto &file = localfiles.at(j);
				if (j % n_bucket == i) {
					thisinput.push_back(file);
				}
			}
			myinput.push_back(thisinput);
		}

		if (!check_splitted_files(myinput, allfiles)) {
			myerror("check_splitted_files failed");
			upcxx_fatal_error("check_splitted_files failed");
		}

	} else {
		myerror("Check files failed for 3 retries");
		upcxx_fatal_error("Check files failed for 3 retries");

	}
	return myinput;
}

std::vector<std::vector<std::string>> get_my_files(const std::string &folder,
		int n_bucket) {
	std::vector<std::string> v = { folder };
	return get_my_files(v, n_bucket);

}

std::vector<std::string> get_my_files(const std::string &folder) {
	std::vector<std::string> v = { folder };
	return get_my_files(v, 1).at(0);

}

std::vector<std::string> get_my_files(const std::vector<std::string> &folders) {
	return get_my_files(folders, 1).at(0);

}

} //end namespace

#endif /* SOURCE_DIRECTORY__SRC_UPCHELPER_H_ */

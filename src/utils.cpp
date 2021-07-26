/*
 * utils.cpp
 *
 *  Created on: May 16, 2020
 *      Author:
 */

#include "utils.h"

#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstring>
#include <set>
#include <map>
#include <algorithm>
#include <random>
#include <cctype>
#include <locale>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include <netdb.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

namespace sparc {

std::string get_env(const std::string &var) {
	const char *val = std::getenv(var.c_str());
	if (val == nullptr) { // invalid to assign nullptr to std::string
		return "";
	} else {
		return val;
	}
}

std::string get_working_dir() {
	char buff[FILENAME_MAX];
	getcwd(buff, FILENAME_MAX);
	return buff;

}

std::string get_hostname() {
	int name_len;
	char hostname[1024];
	hostname[1023] = '\0';
	gethostname(hostname, 1023);
	return hostname;
}

bool endswith(std::string const &fullString, std::string const &ending) {
	if (fullString.length() >= ending.length()) {
		return (0
				== fullString.compare(fullString.length() - ending.length(),
						ending.length(), ending));
	} else {
		return false;
	}
}

std::string int2str(int i) {
	return std::to_string(i);
}

int str2int(const std::string str) {
	return std::stoi(str);
}

bool path_exists(const std::string &s) {
	struct stat buffer;
	return (stat(s.c_str(), &buffer) == 0);
}

bool file_exists(const char *filename) {
	std::ifstream f(filename);
	if (f.good()) {
		f.close();
		return true;
	} else {
		f.close();
		return false;
	}
}

bool dir_exists(const char *path) {
	struct stat info;

	if (stat(path, &info) != 0)
		return false;
	else if (info.st_mode & S_IFDIR)
		return true;
	else
		return false;
}

int make_dir(const char *path) {
	int status = mkdir(path, S_IRWXU | S_IRUSR | S_IXUSR);
	return status;
}

bool is_dir(const char *path) {
	struct stat s;
	if (stat(path, &s) == 0) {
		if (s.st_mode & S_IFDIR) {
			return true;
		}
	}
	return false;
}
std::vector<std::string> list_dir(const char *path) {
	std::vector<std::string> ret;
	std::string strpath = path;
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir(path)) != NULL) {
		/* print all the files and directories within directory */
		while ((ent = readdir(dir)) != NULL) {
			std::string dirname = ent->d_name;
			if (dirname != "." && dirname != "..") {
				std::string s = strpath + "/" + dirname;
				if (!is_dir(s.c_str())) {
					ret.push_back(s);
				}
			}
		}
		closedir(dir);
	} else {
	}

	return ret;
}

// trim from start (in place)
inline void ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
		return !std::isspace(ch);
	}));
}

// trim from end (in place)
inline void rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
		return !std::isspace(ch);
	}).base(), s.end());
}

// trim from both ends (in place)
void trim(std::string &s) {
	ltrim(s);
	rtrim(s);
}

// trim from start (copying)
std::string ltrim_copy(std::string s) {
	ltrim(s);
	return s;
}

// trim from end (copying)
std::string rtrim_copy(std::string s) {
	rtrim(s);
	return s;
}

// trim from both ends (copying)
std::string trim_copy(std::string s) {
	trim(s);
	return s;
}

std::vector<std::string> split(const std::string &source, const char *delimiter,
		bool keepEmpty) {
	std::vector<std::string> results;
	split(results, source, delimiter, keepEmpty);
	return results;
}

void split(std::vector<std::string> &results, const std::string &source,
		const char *delimiter, bool keepEmpty) {

	size_t prev = 0;
	size_t next = 0;

	while ((next = source.find_first_of(delimiter, prev)) != std::string::npos) {
		if (keepEmpty || (next - prev != 0)) {
			results.push_back(source.substr(prev, next - prev));
		}
		prev = next + 1;
	}

	if (prev < source.size()) {
		results.push_back(source.substr(prev));
	}

}

std::vector<std::string> split_not_thread_safe(std::string str,
		std::string sep) {
	char *cstr = const_cast<char*>(str.c_str());
	char *current;
	std::vector<std::string> arr;
	current = strtok(cstr, sep.c_str());
	while (current != NULL) {
		arr.push_back(current);
		current = strtok(NULL, sep.c_str());
	}
	return arr;
}

void split_not_thread_safe(std::vector<std::string> &arr, std::string str,
		std::string sep) {
	char *cstr = const_cast<char*>(str.c_str());
	char *current;
	current = strtok(cstr, sep.c_str());
	while (current != NULL) {
		arr.push_back(current);
		current = strtok(NULL, sep.c_str());
	}
}

std::string get_ip_adderss() {
	return get_ip_adderss(get_hostname());
}

std::string get_ip_adderss(const std::string &hostname) {
	struct hostent *host_entry;
	char *IPbuffer;

	host_entry = gethostbyname(hostname.c_str());

	IPbuffer = inet_ntoa(*((struct in_addr*) host_entry->h_addr_list[0]));

	return IPbuffer;
}

uint32_t hash_string(const std::string &s) {
	if (1) {
		uint32_t h = 0;
		const char *p = s.c_str();
		while (*p) {
			h += *p;
			p++;
		}
		return h;
	} else {
#define A 54059 /* a prime */
#define B 76963 /* another prime */
#define C 86969 /* yet another prime */
#define FIRSTH 37 /* also prime */
		const char *p = s.c_str();
		unsigned h = FIRSTH;
		while (*p) {
			h = (h * A) ^ (p[0] * B);
			p++;
		}
		return h; // or return h % C;
	}
}

}

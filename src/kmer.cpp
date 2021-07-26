#include <map>
#include <algorithm>    // std::transform
#include <unordered_set>
#include <tuple>
#include <exception>
#include <sstream>
//#include <boost/algorithm/string.hpp>
//#include <boost/lambda/lambda.hpp>
//#include <boost/algorithm/searching/knuth_morris_pratt.hpp>
#include <stdlib.h>
#include <cmath>
#include "kmer.h"
#include "utils.h"
using namespace std;

#define rand01 (((double) rand() / (RAND_MAX)) )

/*
 * Kmer functions
 */

unsigned long kmer_to_number(const string &kmer) {
	if (kmer.size() > 31) {
		throw std::runtime_error("at most 31-mer");
	}
	unsigned long num = 0;
	for (int i = 0; i < int(kmer.size()); i++) {
		if (i > 0) {
			num = num << 2;
		}
		char chr = kmer.at(i);
		if (chr == 'C') {
			num += 0;
		} else if (chr == 'T') {
			num += 1;
		} else if (chr == 'A') {
			num += 2;
		} else if (chr == 'G') {
			num += 3;
		} else {
			throw std::runtime_error("only ATCG can be in kmer");
		}
	}
	return num;
}

string number_to_kmer(unsigned long num, int k) {
	if (k > 31) {
		throw std::runtime_error("at most 31-mer");
	}
	char kmer[k + 1];
	for (int i = 0; i < k; i++) {
		int rem = num % 4;
		char c;
		if (rem == 0) {
			c = 'C';
		} else if (rem == 1) {
			c = 'T';
		} else if (rem == 2) {
			c = 'A';
		} else if (rem == 3) {
			c = 'G';
		}
		kmer[k - i - 1] = c;
		num = num / 4;
	}
	kmer[k] = '\0';
	return string(kmer);
}

static map<char, char> KMER_RC = { { 'A', 'T' }, { 'T', 'A' }, { 'C', 'G' }, {
		'G', 'C' }, { 'N', 'N' } };

inline char rc_transform(char c) {
	return KMER_RC.at(c);
}

inline vector<string> _internal_generate_kmer(const string &seq, int k,
		bool is_canonical) {
	vector<string> kmers;

	for (int i = 0; i < int(seq.size()) - k + 1; i++) {
		if (is_canonical) {
			kmers.push_back(canonical_kmer(seq.substr(i, k)));
		} else {

			kmers.push_back(seq.substr(i, k));
		}
	}

	return kmers;
}

string random_generate_kmer(const string &seq, int k, bool is_canonical) {

	int m = int(seq.size()) - k;
	int i = int(rand01 * m);
	if (rand01 > 0.5) {
		if (is_canonical)
			return canonical_kmer(seq.substr(i, k));
		else
			return seq.substr(i, k);
	} else {
		string seqrc = reverse_complement(seq);
		if (is_canonical)
			return canonical_kmer(seqrc.substr(i, k));
		else
			return seqrc.substr(i, k);

	}
}

vector<string> generate_kmer(const string &seq, int k, const char *err_char,
		bool is_canonical) {
	std::vector<std::string> splits;

	sparc::split(splits, seq, err_char);
	//boost::split(splits, seq, boost::lambda::_1 == err_char);

	std::vector<string> results;

	for (size_t i = 0; i < splits.size(); i++) {

		vector<string> tmp = _internal_generate_kmer(splits.at(i), k,
				is_canonical);
		results.insert(results.end(), tmp.begin(), tmp.end());
		tmp = _internal_generate_kmer(reverse_complement(splits.at(i)), k,
				is_canonical);
		results.insert(results.end(), tmp.begin(), tmp.end());

	}

	std::unordered_set<string> s;
	for (string i : results)
		s.insert(i);
	results.assign(s.begin(), s.end());

	return results;
}

vector<unsigned long> generate_kmer_number(const string &seq, int k,
		const char *err_char, bool is_canonical) {
	vector<string> kmers = generate_kmer(seq, k, err_char, is_canonical);
	vector<unsigned long> numbers(kmers.size());
	for (size_t i = 0; i < kmers.size(); i++) {
		numbers[i] = kmer_to_number(kmers.at(i));
	}
	return numbers;
}

vector<string> generate_kmer_for_fastseq(const string &seq, int k,
		const char *err_char, bool is_canonical) {
	std::vector<std::string> splits;

	sparc::split(splits, seq, err_char);

	//boost::split(splits, seq, boost::lambda::_1 == err_char);

	std::vector<string> results;

	for (size_t i = 0; i < splits.size(); i++) {

		vector<string> tmp = _internal_generate_kmer(splits.at(i), k,
				is_canonical);
		results.insert(results.end(), tmp.begin(), tmp.end());
	}
	return results;
}

string canonical_kmer(const string &seq) {
	string rc = reverse_complement(seq);
	return rc < seq ? rc : seq;
}

string reverse_complement(const string &seq) {
	string rseq;
	rseq.resize(seq.size());
	std::transform(seq.rbegin(), seq.rend(), rseq.begin(), rc_transform);
	return rseq;
}

/*
 * overlap functions
 */

inline std::vector<std::tuple<int, string>> split_evenly(const string &text,
		int n) {
	std::vector<std::tuple<int, string>> v;
	int N = text.length();
	if (N < n) {
		throw std::runtime_error("text is too small");
	}
	int partlen = int(N / n);
	int i = 0;
	while (i < n) {
		int k = i + 1 == n ? N - partlen * (n - 1) : partlen;
		auto a = std::make_tuple(i * partlen, text.substr(i * partlen, k));
		v.push_back(a);
		i++;
	}
	return v;
}

/*
 * combination functions
 */

/* arr[]  ---> Input Array
 data[] ---> Temporary array to store current combination
 start & end ---> Staring and Ending indexes in arr[]
 index  ---> Current index in data[]
 r ---> Size of a combination to be printed
 from https://ideone.com/ywsqBz
 */
void combinationUtil(uint32_t arr[], uint32_t data[], int start, int end,
		int index, std::vector<std::pair<uint32_t, uint32_t>> &results) {
	int r = 2;
	// Current combination is ready to be printed, print it
	if (index == r) {
		if (0) {
			for (int i = 0; i < r; i++)
				printf("%d ", data[i]);
			printf("\n");
		}
		results.push_back(std::make_pair(data[0], data[1]));
		return;
	}

	// replace index with all possible elements. The condition
	// "end-i+1 >= r-index" makes sure that including one element
	// at index will make a combination with remaining elements
	// at remaining positions
	for (int i = start; i <= end && end - i + 1 >= r - index; i++) {
		data[index] = arr[i];
		combinationUtil(arr, data, i + 1, end, index + 1, results);

		// Remove duplicates
		while (arr[i] == arr[i + 1])
			i++;
	}
}

void make_combination2(std::vector<uint32_t> &v,
		std::vector<std::pair<uint32_t, uint32_t> > &results) {

	// A temporary array to store all combination one by one
	uint32_t data[2];
	int n = v.size();
	uint32_t *arr = v.data();

	// Print all combination using temprary array 'data[]'
	combinationUtil(arr, data, 0, n - 1, 0, results);
}

inline void random_choice(uint32_t *arr, int n, uint32_t *data, int m) {
	for (int i = 0; i < m; i++)
		data[i] = arr[rand() % n];
}
//assume reads has no duplicates.
std::vector<std::pair<uint32_t, uint32_t> > generate_edges(
		std::vector<uint32_t> &reads, size_t max_degree) {
	if (reads.size() <= max_degree) {

		std::vector<std::pair<uint32_t, uint32_t> > results;
		make_combination2(reads, results);
		return results;
	} else {
		std::vector<std::pair<uint32_t, uint32_t> > results;
		int m = std::ceil(max_degree / 2.0 * reads.size());
		uint32_t d1[m], d2[m];
		random_choice(reads.data(), reads.size(), d1, m);
		random_choice(reads.data(), reads.size(), d2, m);
		for (int i = 0; i < m; i++) {
			if (d1[i] != d2[i])
				results.push_back(std::make_pair(d1[i], d2[i]));
		}
		return results;
	}
}

//copy from facebook fasttext
// The correct implementation of fnv should be:
// h = h ^ uint32_t(uint8_t(str[i]));
// Unfortunately, earlier version of fasttext used
// h = h ^ uint32_t(str[i]);
// which is undefined behavior (as char can be signed or unsigned).
// Since all fasttext models that were already released were trained
// using signed char, we fixed the hash function to make models
// compatible whatever compiler is used.
uint32_t fnv_hash(const std::string &str) {
	uint32_t h = 2166136261;
	for (size_t i = 0; i < str.size(); i++) {
		h = h ^ uint32_t(int8_t(str[i]));
		h = h * 16777619;
	}
	return h;
}

uint32_t fnv_hash(uint32_t n) {
	return fnv_hash(std::to_string(n));
}

static const char BASE64_CHARS[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		"abcdefghijklmnopqrstuvwxyz"
		"0123456789+/";

static map<char, int> CHARS_BASE64 = { { 'A', 0 }, { 'B', 1 }, { 'C', 2 }, {
		'D', 3 }, { 'E', 4 }, { 'F', 5 }, { 'G', 6 }, { 'H', 7 }, { 'I', 8 }, {
		'J', 9 }, { 'K', 10 }, { 'L', 11 }, { 'M', 12 }, { 'N', 13 },
		{ 'O', 14 }, { 'P', 15 }, { 'Q', 16 }, { 'R', 17 }, { 'S', 18 }, { 'T',
				19 }, { 'U', 20 }, { 'V', 21 }, { 'W', 22 }, { 'X', 23 }, { 'Y',
				24 }, { 'Z', 25 }, { 'a', 26 }, { 'b', 27 }, { 'c', 28 }, { 'd',
				29 }, { 'e', 30 }, { 'f', 31 }, { 'g', 32 }, { 'h', 33 }, { 'i',
				34 }, { 'j', 35 }, { 'k', 36 }, { 'l', 37 }, { 'm', 38 }, { 'n',
				39 }, { 'o', 40 }, { 'p', 41 }, { 'q', 42 }, { 'r', 43 }, { 's',
				44 }, { 't', 45 }, { 'u', 46 }, { 'v', 47 }, { 'w', 48 }, { 'x',
				49 }, { 'y', 50 }, { 'z', 51 }, { '0', 52 }, { '1', 53 }, { '2',
				54 }, { '3', 55 }, { '4', 56 }, { '5', 57 }, { '6', 58 }, { '7',
				59 }, { '8', 60 }, { '9', 61 }, { '+', 62 }, { '/', 63 } };
std::string ulong2base64(unsigned long n) {
	if (n == 0)
		return "A";
	stringstream s;
	int r;
	while (n != 0) {
		r = n % 64;
		s << BASE64_CHARS[r];
		n = n >> 6;
	}
	std::string str = s.str();
	std::reverse(str.begin(), str.end());
	return str;
}

unsigned long base64toulong(const std::string &str) {
	unsigned long n = 0;
	for (size_t i = 0; i < str.size(); i++) {
		n = n << 6;
		char c = str.at(i);
		n += CHARS_BASE64[c];
	}
	return n;
}

static const char base64_table[] =
		"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

int Base64encode_len(int len) {
	return ((len + 2) / 3 * 4) + 1;
}

std::string base64_encode(const unsigned char *src, size_t len) {
	unsigned char *out, *pos;
	const unsigned char *end, *in;

	size_t olen;

	olen = 4 * ((len + 2) / 3); /* 3-byte blocks to 4-byte */

	if (olen < len)
		return std::string(); /* integer overflow */

	std::string outStr;
	outStr.resize(olen);
	out = (unsigned char*) &outStr[0];

	end = src + len;
	in = src;
	pos = out;
	while (end - in >= 3) {
		*pos++ = base64_table[in[0] >> 2];
		*pos++ = base64_table[((in[0] & 0x03) << 4) | (in[1] >> 4)];
		*pos++ = base64_table[((in[1] & 0x0f) << 2) | (in[2] >> 6)];
		*pos++ = base64_table[in[2] & 0x3f];
		in += 3;
	}

	if (end - in) {
		*pos++ = base64_table[in[0] >> 2];
		if (end - in == 1) {
			*pos++ = base64_table[(in[0] & 0x03) << 4];
			*pos++ = '=';
		} else {
			*pos++ = base64_table[((in[0] & 0x03) << 4) | (in[1] >> 4)];
			*pos++ = base64_table[(in[1] & 0x0f) << 2];
		}
		*pos++ = '=';
	}

	return outStr;
}
std::string Base64encode(const std::string &word) {
	return base64_encode((const unsigned char*) (word.c_str()), word.size());
}

static std::map<char, int> dnabase_encoding = { { 'A', 0 }, { 'T', 1 },
		{ 'C', 2 }, { 'G', 3 } };

string kmer_to_base64(const string &kmer) {
	//one byte for 4 bases
	int size = (int) (kmer.size() / 4);
	if (kmer.size() % 4 > 0) {
		size++;
	}
	unsigned char bytes[size + 1];

	for (size_t i = 0; i < kmer.size(); i += 4) {
		int a = 0, b = 0, c = 0, d = 0;
		a = dnabase_encoding[kmer.at(i)];
		if (i + 1 < kmer.size()) {
			b = dnabase_encoding[kmer.at(i + 1)];
		}

		if (i + 2 < kmer.size()) {
			c = dnabase_encoding[kmer.at(i + 2)];
		}
		if (i + 3 < kmer.size()) {
			d = dnabase_encoding[kmer.at(i + 3)];
		}
		int e = (a << 6) + (b << 4) + (c << 2) + d;
		bytes[i / 4] = e;
	}
	bytes[size] = '\0';
	return base64_encode(bytes, size);
}

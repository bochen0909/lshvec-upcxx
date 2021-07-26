#include <vector>
#include <string>
#include <map>

unsigned long kmer_to_number(const std::string &kmer);

std::string number_to_kmer(unsigned long num, int k);

std::vector<std::string> generate_kmer(const std::string &seq, int k,
		const char *err_char, bool is_canonical);

std::vector<unsigned long> generate_kmer_number(const std::string &seq, int k,
		const char *err_char, bool is_canonical);

std::vector<std::string> generate_kmer_for_fastseq(const std::string &seq,
		int k, const char *err_char, bool is_canonical);

std::string random_generate_kmer(const std::string &seq, int k,
		bool is_canonical);

std::string canonical_kmer(const std::string &seq);

std::string reverse_complement(const std::string &seq);

std::vector<std::pair<uint32_t, uint32_t> > generate_edges(
		std::vector<uint32_t> &reads, size_t max_degree);

uint32_t fnv_hash(const std::string &str);

uint32_t fnv_hash(uint32_t n);

std::string ulong2base64(unsigned long n);

unsigned long base64toulong(const std::string &str);

std::string kmer_to_base64(const std::string &kmer);

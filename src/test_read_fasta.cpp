#include <exception>
#include <stdexcept>
#include <string>
#include <algorithm>

#include "acutest.h"
#include "io.h"

using namespace std;

#define TEST_INT_EQUAL(a, b)       \
	TEST_CHECK((a) == (b));        \
	TEST_MSG("Expected: %d", (a)); \
	TEST_MSG("Produced: %d", (b));

#define TEST_DOUBLE_EQUAL(a, b)              \
	TEST_CHECK(std::abs((a) - (b)) < 1e-10); \
	TEST_MSG("Expected: %f", (double)(a));   \
	TEST_MSG("Produced: %f", (double)(b));

#define TEST_STR_EQUAL(a, b)            \
	TEST_CHECK(string(a) == string(b)); \
	TEST_MSG("Expected: %s", (a));      \
	TEST_MSG("Produced: %s", (b));

void test_read(void)
{
	for (auto &filepath : {"knucleotide.fasta", "knucleotide.fasta.gz"})
	{
		auto v = read_fasta(filepath);
		TEST_INT_EQUAL(v.size(), 3);
		TEST_INT_EQUAL(v[0].seq.size(), 50000);
		TEST_INT_EQUAL(v[1].seq.size(), 75000);
		TEST_INT_EQUAL(v[2].seq.size(), 125000);
	}
}

TEST_LIST = {{"test_read", test_read},

			 {NULL, NULL}};

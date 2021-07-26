#include <exception>
#include <stdexcept>
#include <string>
#include <algorithm>

#include "acutest.h"
#include "utils.h"

using namespace std;
using namespace sparc;

#define TEST_INT_EQUAL(a,b) \
		TEST_CHECK( (a) == (b)); \
		TEST_MSG("Expected: %d", (a)); \
		TEST_MSG("Produced: %d", (b));

#define TEST_DOUBLE_EQUAL(a,b) \
		TEST_CHECK( std::abs((a) - (b))<1e-10); \
		TEST_MSG("Expected: %f", (double)(a)); \
		TEST_MSG("Produced: %f", (double)(b));

#define TEST_STR_EQUAL(a,b) \
		TEST_CHECK( string(a) == string(b)); \
		TEST_MSG("Expected: %s", (a)); \
		TEST_MSG("Produced: %s", (b));

void test_trim(void) {
	std::string s = string(" ");
	trim(s);
	TEST_STR_EQUAL("", s)
	s = string("\r");
	trim(s);
	TEST_STR_EQUAL("", s)
	s = string("\n");
	trim(s);
	TEST_STR_EQUAL("", s)
	s = string("\t");
	trim(s);
	TEST_STR_EQUAL("", s)
	s = string("\r\n");
	trim(s);
	TEST_STR_EQUAL("", s)
	s = string(" \t\n");
	trim(s);
	TEST_STR_EQUAL("", s)
	s = string("\t \t      AAA \t  \t \r \n");
	trim(s);
	TEST_STR_EQUAL("AAA", s)
}

void test_split(void) {
	std::vector<std::string> v;
	char *sep = " ";
	{
		v = split("", sep);
		TEST_INT_EQUAL(0, v.size());
		v = split(" ", sep);
		TEST_INT_EQUAL(0, v.size());
		v = split("  a b   c       d  e  ", sep);
		TEST_INT_EQUAL(5, v.size());
	}
	sep = "\t";
	{
		v = split("", sep);
		TEST_INT_EQUAL(0, v.size());
		v = split("\t", sep);
		TEST_INT_EQUAL(0, v.size());
		v = split("   a\t b \t  c  \t     d \t e  ", sep);
		TEST_INT_EQUAL(5, v.size());
	}

	sep = ",";
	{
		v = split("", sep);
		TEST_INT_EQUAL(0, v.size());
		v = split(",", sep);
		TEST_INT_EQUAL(0, v.size());
		v = split("123,12 313,asdsaf,,", sep);
		TEST_INT_EQUAL(3, v.size());
	}

}

TEST_LIST = { {"test_trim", test_trim},

	{	"test_split", test_split},
	{	NULL, NULL}};


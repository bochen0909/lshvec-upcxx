/*
 * io_test.h
 *
 *  Created on: Jun 24, 2021
 *      Author: Bo Chen
 */


#include "utils.h"
#include "CRandProj.h"
#include "acutest.h"

#define TOLERANCE 0.00001f

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

using namespace rpns;
void rp_test ()
{
    std::vector<std::string> kmers{"GCGGA",
                                   "GGCTG",
                                   "TCGCA",
                                   "CGTCA",
                                   "GATCT",
                                   "GCAAA",
                                   "ACTTG",
                                   "ATGAT",
                                   "ACATC",
                                   "TAAAG",
                                   "CCGGT",
                                   "GCAAT",
                                   "GCATG",
                                   "GGGCG",
                                   "ATTTC",
                                   "AGAGC",
                                   "TGAGG",
                                   "TATTG",
                                   "TCTGC",
                                   "TCGAA",
                                   "CAGCC",
                                   "TACCA",
                                   "ATCGT",
                                   "TTCCT",
                                   "CGTTT",
                                   "TAGGC",
                                   "TGGTC",
                                   "GAGTG",
                                   "CCGGT",
                                   "ATCTC",
                                   "ACACG",
                                   "GTTAC",
                                   "CCATC",
                                   "CTAGT",
                                   "CGACA",
                                   "GATAT",
                                   "TGAGC",
                                   "AATAA",
                                   "TAACA",
                                   "AAGTG",
                                   "CCAGG",
                                   "CGGAT",
                                   "GAACG",
                                   "TCTCA",
                                   "TCTCG",
                                   "ATTCT",
                                   "CAGTC",
                                   "CGCGG",
                                   "GCCTG",
                                   "TGGTC",
                                   "CGTTG",
                                   "AATGT",
                                   "CGGTT",
                                   "CTCGT",
                                   "TTGCT",
                                   "GAATT",
                                   "AACAA",
                                   "ACGCC",
                                   "CACTT",
                                   "AGGAT",
                                   "ATTTT",
                                   "GCCTC",
                                   "CCGGC",
                                   "GACAC",
                                   "TTTCT",
                                   "GTCTG",
                                   "AGCCA",
                                   "CAGTA",
                                   "CCGAA",
                                   "ATTCT",
                                   "TACAC",
                                   "TTGGT",
                                   "GACAA",
                                   "GAGTC",
                                   "ATGGA",
                                   "TTGCA",
                                   "CATCG",
                                   "GGTGT",
                                   "TGGCA",
                                   "ACCGG",
                                   "CCGTA",
                                   "TCTGT",
                                   "TTCTG",
                                   "GCTCC",
                                   "CCCTC",
                                   "CCCCG",
                                   "GCCCA",
                                   "GTAGA",
                                   "GTAAG",
                                   "GTGCG",
                                   "CAAGG",
                                   "CTGCG",
                                   "CACCC",
                                   "GTGGA",
                                   "GTCCG",
                                   "GATCT",
                                   "AAGGC",
                                   "GCCAT",
                                   "GGCAG",
                                   "GTGGG",
                                   "AGTCA",
                                   "GTTCT",
                                   "TTAGT",
                                   "GTCGT",
                                   "ATGAT",
                                   "CCAAT",
                                   "TCAGC",
                                   "TCTTA",
                                   "TAAAC",
                                   "GTCCG",
                                   "AACAC",
                                   "CATAA",
                                   "AAGAT",
                                   "TGTTC",
                                   "TCAGT",
                                   "AACCA",
                                   "GCTCT",
                                   "CGGGT",
                                   "TCTCC",
                                   "TGGGT",
                                   "CCCAA",
                                   "GTTGC",
                                   "TGGTC",
                                   "TTCGA",
                                   "CAGGT",
                                   "CTATT",
                                   "GAGGA",
                                   "GTATT",
                                   "TTACT",
                                   "GCAGT",
                                   "GAAGA",
                                   "CTGTA",
                                   "TTCTC",
                                   "AACGA",
                                   "TGATT",
                                   "TAATA",
                                   "CGGGT",
                                   "TTTAC",
                                   "GTAAA",
                                   "TGTTC",
                                   "GGTGT",
                                   "ACCCT",
                                   "GTAGG",
                                   "TACAC",
                                   "TGACA",
                                   "ATTCT",
                                   "CGACC",
                                   "CACAC",
                                   "CAGAG",
                                   "GCAGT",
                                   "AGAGC",
                                   "TTCAG",
                                   "GTTGC",
                                   "GAGAT",
                                   "GGAAA",
                                   "CTGAG",
                                   "ATCGC",
                                   "GAATA",
                                   "CATAA",
                                   "CCGTG",
                                   "AAAAA",
                                   "CCATT",
                                   "GCCAC",
                                   "GTTTG",
                                   "TGCCA",
                                   "GCGCC",
                                   "CTCCA",
                                   "GAGCC",
                                   "GCTTC",
                                   "GCTCT",
                                   "ATCTA",
                                   "GTGGG",
                                   "CTCTC",
                                   "GGGAC",
                                   "GAGCT",
                                   "GACTC",
                                   "TATTA",
                                   "CTGTC",
                                   "TCAGC",
                                   "TGAGG",
                                   "CTCTA",
                                   "ACCAG",
                                   "ATAGT",
                                   "CCCGG",
                                   "GCAGG",
                                   "ATTTG",
                                   "TGCAG",
                                   "TCCAC",
                                   "AACTC",
                                   "CGTTC",
                                   "CTGAA",
                                   "GGGCG",
                                   "CACAT",
                                   "TTCTA",
                                   "CGCGT",
                                   "AGACA",
                                   "GGAGC",
                                   "ACTTC",
                                   "CCTGT",
                                   "GTCGC"};

    std::vector<std::vector<size_t>> pseduo_index{{98, 186, 80, 7, 33}, {95, 162, 127, 156, 10}, {17, 170, 0, 163, 97}};
    int kmer_size = 5;
    int hash_size = 3;
    RandProjBuilder builder(kmer_size, hash_size, 0);
    builder.__set_pseduo_index(pseduo_index);
    builder.create_hash(kmers);
    builder.print();

    auto rp = builder.make_rp();

    std::unordered_map<size_t, int> value_counts;
    for (size_t i = 0; i < kmers.size(); i++)
    {
        auto h = rp->hash(kmers.at(i), true);
        if (i < 5)
        {
            std::cout << kmers.at(i) << " " << h << std::endl;
        }
        value_counts[h]++;
    }
    for (auto &kv : value_counts)
    {
        std::cout << kv.first << "\t" << kv.second << "\n";
    }
    value_counts.clear();
    std::cout << "no rc \n";
    for (size_t i = 0; i < kmers.size(); i++)
    {
        auto h = rp->hash(kmers.at(i), false);
        if (i < 5)
        {
            std::cout << kmers.at(i) << " " << h << std::endl;
        }
        value_counts[h]++;
    }

    for (auto &kv : value_counts)
    {
        std::cout << kv.first << "\t" << kv.second << "\n";
    }
}

void rp_test_save_load ()

{
    std::vector<std::string> kmers{"GCGGA",
                                   "GGCTG",
                                   "TCGCA",
                                   "CGTCA",
                                   "GGTGT",
                                   "TGGCA",
                                   "ACCGG",
                                   "CCGTA",
                                   "TCTGT",
                                   "CTCTC",
                                   "GGGAC",
                                   "GAGCT",
                                   "GACTC",
                                   "TATTA",
                                   "CTGTC",
                                   "TCAGC",
                                   "TGAGG",
                                   "CTCTA",
                                   "ACCAG",
                                   "ATAGT",
                                   "CCCGG",
                                   "GCAGG",
                                   "ATTTG",
                                   "TGCAG",
                                   "TCCAC",
                                   "AACTC",
                                   "CGTTC",
                                   "CTGAA",
                                   "GGGCG",
                                   "CACAT",
                                   "TTCTA",
                                   "CGCGT",
                                   "AGACA",
                                   "GGAGC",
                                   "ACTTC",
                                   "CCTGT",
                                   "GTCGC"};

    int kmer_size = 5;
    int hash_size = 3;
    RandProjBuilder builder(kmer_size, hash_size, 0);
    builder.create_hash(kmers);

    auto rp = builder.make_rp();
    rp->save("rp1.bin");
    std::vector<size_t> h1;
    for (size_t i = 0; i < kmers.size(); i++)
    {
        auto h = rp->hash(kmers.at(i), true);
        h1.push_back(h);
    }

    CRandProj rp2;
    rp2.load("rp1.bin");
    for (size_t i = 0; i < kmers.size(); i++)
    {
        auto h = rp2.hash(kmers.at(i), true);
        TEST_INT_EQUAL(h, h1.at(i));
    }
}


TEST_LIST = {
    {"rp_test", rp_test},
    {"rp_test_save_load", rp_test_save_load},
			 {NULL, NULL}};
#ifndef __SRC_IO_H__
#define __SRC_IO_H__

#include <string> 
#include <fstream> 
#include <vector> 

struct FastaRecord{
    std::string id;
    std::string seq;
};
std::vector<FastaRecord>  read_fasta(const std::string& filepath);


#endif 
//
//  plink.hpp
//  plink
//
//  Created by Shing Wan Choi on 18/08/2016.
//  Copyright © 2016 Shing Wan Choi. All rights reserved.
//

#ifndef plink_hpp
#define plink_hpp

#include <stdio.h>
#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <deque>
#include <limits.h>
#include <vector>
#include "misc.hpp"

class PLINK{
public:
    //Initialize plink object with the bim bed fam file prefix
    PLINK(std::string prefix):m_prefix(prefix){
        m_init = false;
        #if defined(__LP64__)
        m_bit_size = sizeof(uint64_t)*CHAR_BIT;
        #elif defined(_WIN64)
        m_bit_size = sizeof(uint64_t)*CHAR_BIT;
        #else
        m_bit_size = sizeof(uint32_t)*CHAR_BIT;
        #endif
        m_num_bytes=0;
        m_num_sample=0;
        m_num_snp=0;
        m_required_bit = 0;
        m_snp_iter=0;
    };
    ~PLINK();
    void initialize();
    int read_snp(int num_snp);
    void lerase(int num);
    size_t get_num_snp() const{return m_num_snp; };
    size_t get_num_sample() const{return m_num_sample;};
    std::vector<int> get_genotype(int geno) const;
    int get_distance() const{
    		return m_bp_list.back()-m_bp_list.front();
    };
    int get_first_bp() const{ return m_bp_list.front(); };
    int get_last_bp() const{return m_bp_list.back();};
private:
    bool openPlinkBinaryFile(const std::string s, std::ifstream & BIT);
    bool m_init;
    std::string m_prefix;
    std::ifstream m_bed;
    std::ifstream m_bim;
    size_t m_num_sample;
    size_t m_num_snp;
    size_t m_num_bytes;
    size_t m_snp_iter;
    size_t m_bit_size;
    size_t m_required_bit;
    std::deque<std::string> m_chr_list;
    std::deque<std::string> m_snp_list;
    std::deque<size_t> m_cm_list;
    std::deque<size_t> m_bp_list;
    std::deque<std::string> m_ref_allele;
    std::deque<std::string> m_alt_allele;
    std::deque<double> m_maf;
#if defined(__LP64__) || defined(_WIN64)
    // LP64 machine, OS X or Linux
    std::deque<uint64_t*> m_genotype;
    std::deque<uint64_t*> m_missing;
    uint64_t FIVEMASK = ((~0LLU) / 3);
    uint64_t THREE = 3LLU;
    uint64_t THREEMASK = 0x3333333333333333LLU;
    uint64_t OFMASK = 0x0f0f0f0f0f0f0f0fLLU;
#if defined(_WIN64)
	#define __LP64__
#endif
#else
    // 32-bit machine, Windows or Linux or OS X
    std::deque<uint32_t*> m_genotype;
    std::deque<uint32_t*> m_missing;
    const uint32_t FIVEMASK = ((~0LU) / 3);
    const uint32_t THREE = 3LU;
    const uint32_t THREEMASK = 0x33333333;
    	const uint32_t OFMASK = 0x0f0f0f0f;
    	const uint32_t ONEZEROMASK = 0x01010101;
    	const uint32_t AAAAMASK = 0xaaaaaaaa;
#endif
    
};
#endif /* plink_hpp */

//
//  plink.cpp
//  plink
//
//  Created by Shing Wan Choi on 18/08/2016.
//  Copyright © 2016 Shing Wan Choi. All rights reserved.
//

#include "plink.hpp"

PLINK::~PLINK(){
	for(size_t i = 0; i < m_genotype.size(); ++i){
		delete m_genotype[i];
		delete m_missing[i];
	}
}

double PLINK::get_r2(const size_t i, const size_t j){
	double r2 =0.0;
	if(i >= m_genotype.size() || j >=m_genotype.size()) throw std::runtime_error("Out of bound error! In R2 calculation.");
	// Crazy stuff of plink here
	return r2;
}

void PLINK::compute_clump( size_t index, size_t i_start, size_t i_end, std::vector<SNP> &snp_list, const std::deque<size_t> &index_check, const double r2_threshold){
	size_t ref_index = index_check[index];
	for(size_t i = i_start; i < i_end && i < index_check.size(); ++i){
		if(i != index){
			// This is a pair of SNPs
			double r2 = get_r2(i, index);
			if(r2 >= r2_threshold){
				size_t target_index = index_check[i];
				if(snp_list[target_index].get_p_value() < snp_list[ref_index].get_p_value()) snp_list[target_index].add_clump(ref_index);
				else snp_list[ref_index].add_clump(target_index);
			}
		}
	}
}

void PLINK::clump(const size_t index, const std::deque<size_t> &index_check, std::vector<SNP> &snp_list, const double r2_threshold){
	std::vector<std::thread> thread_store;
	if((index_check.size()-1) < m_thread){
		for(size_t i = 0; i < index_check.size(); ++i){
			if(index_check[i]!=index) thread_store.push_back(std::thread(&PLINK::compute_clump, this, index,i, i+1, std::ref(snp_list), std::cref(index_check), r2_threshold));
		}
	}
	else{
		int num_snp_per_thread =(int)(index_check.size()-1) / (int)m_thread;  //round down
		int remain = (int)(index_check.size()-1) % (int)m_thread;
		int cur_start = 0;
		int cur_end = num_snp_per_thread;
		for(size_t i = 0; i < m_thread; ++i){
			thread_store.push_back(std::thread(&PLINK::compute_clump, this, index, cur_start, cur_end+(remain>0), std::ref(snp_list), std::cref(index_check),r2_threshold ));
			cur_start = cur_end+(remain>0);
			cur_end+=num_snp_per_thread+(remain>0);
			if(cur_end>index_check.size()) cur_end =index_check.size();
			remain--;
		}
	}
	for (size_t i = 0; i < thread_store.size(); ++i) thread_store[i].join();
	thread_store.clear();
}

void PLINK::clumping(std::map<std::string, bool> inclusion, std::vector<SNP> &snp_list, const std::map<std::string, size_t> &snp_index, double p_threshold, double r2_threshold, size_t kb_threshold){
	// Go through all SNPs
	std::deque<size_t> bp_check;
	std::deque<size_t> index_check;
	size_t require_bp=0, require_index=0; //this is the index on index_check
	bool requiring=false;
	for(size_t i = 0; i < m_snp_list.size(); ++i){
		std::string rs=m_snp_list[i];
		if(inclusion.find(rs)!=inclusion.end() && snp_index.find(rs)!=snp_index.end()){
			// required SNP
			size_t cur_loc = snp_list[snp_index.at(rs)].get_loc();
			while(requiring && cur_loc-require_bp > kb_threshold){
				// All genotypes are here, now calculate the P-value
				clump(require_index, index_check, snp_list, r2_threshold);
				// now go through the whole remaining index to check if there are any required
				requiring = false;
				for(size_t check = require_index+1; check< index_check.size(); ++check){
					if(snp_list[index_check[check]].get_p_value() < p_threshold){
						requiring = true;
						require_index =check;
						require_bp = snp_list[index_check[check]].get_loc();
					}
				}
				if(requiring){
					// remove SNPs in the front that are too far away
					size_t num_remove=0;
					for(size_t check; check < require_index; ++check){
						if(require_bp-snp_list[index_check[check]].get_loc() > kb_threshold)num_remove++;
						else break;
					}
					if(num_remove!=0) lerase(num_remove);
				}
			}
			read_snp(1, true);
			bp_check.push_back(cur_loc);
			index_check.push_back(snp_index.at(rs));
			if(!requiring && snp_list[index_check.back()].get_p_value() < p_threshold){
				require_bp = cur_loc;
				require_index = index_check.size()-1;
			}
		}
		else m_bed.seekg(m_num_bytes, m_bed.cur); //This should skip 1 SNP
	}
	if(m_genotype.size() != 0){
		//Still got stuff to work on
		// keep clumping
		while(requiring){
			clump(require_index, index_check, snp_list, r2_threshold);
			requiring = false;
			for(size_t check = require_index+1; check< index_check.size(); ++check){
				if(snp_list[index_check[check]].get_p_value() < p_threshold){
					requiring = true;
					require_index =check;
					require_bp = snp_list[index_check[check]].get_loc();
				}
			}
			if(requiring){
				size_t num_remove=0;
				for(size_t check; check < require_index; ++check){
					if(require_bp-snp_list[index_check[check]].get_loc() > kb_threshold)num_remove++;
					else break;
				}
				if(num_remove!=0) lerase(num_remove);
			}
		}
	}
	//completed
}

std::vector<int> PLINK::get_genotype(int geno) const{
	if(geno >= m_genotype.size()){
		std::string error_message = "Asked for "+std::to_string(geno)+" genotype but contain only "+std::to_string(m_genotype.size());
		throw std::runtime_error(error_message);
	}
	std::vector<int> res(m_num_sample);
	for(size_t i = 0; i < m_num_sample; ++i){
		int index =(i*2)/m_bit_size;
		int info = m_genotype[geno][index] >> (m_bit_size-(i+1)*2)& THREEMASK;
		switch(info){
			case 0:
				res[i] = 0;
				break;
			case 1:
				res[i] = -1;
				break;
			case 2:
				res[i] = 1;
				break;
			case 3:
				res[i] = 2;
				break;
			default:
				throw std::runtime_error("Undefined genotype");
		}
	}
	return res;
}

void PLINK::lerase(int num){
    if(num <0){
        std::string error_message = "Number of removed SNPs cannot be less than 1: "+std::to_string(num);
        throw std::runtime_error(error_message);
    }
    if(num >= m_genotype.size()){
        std::string error_message = "Number of removed SNPs exceed number of SNPs available "+std::to_string(num)+" "+std::to_string(m_genotype.size());
        throw std::runtime_error(error_message);
    }
    for(size_t i = 0; i < num; ++i){
        delete m_genotype[i];
        delete m_missing[i];
    }
    m_genotype.erase(m_genotype.begin(), m_genotype.begin()+num);
    m_missing.erase(m_missing.begin(), m_missing.begin()+num);
    m_chr_list.erase(m_chr_list.begin(), m_chr_list.begin()+num);
    //m_snp_list.erase(m_snp_list.begin(), m_snp_list.begin()+num);
    m_cm_list.erase(m_cm_list.begin(), m_cm_list.begin()+num);
    m_bp_list.erase(m_bp_list.begin(), m_bp_list.begin()+num);
    m_ref_allele.erase(m_ref_allele.begin(), m_ref_allele.begin()+num);
    m_alt_allele.erase(m_alt_allele.begin(), m_alt_allele.begin()+num);
    m_maf.erase(m_maf.begin(), m_maf.begin()+num);
}

//The return value should be the number of remaining SNPs
int PLINK::read_snp(int num_snp, bool ld){
    if(!m_init) throw std::runtime_error("Class uninitialize! Must initialize before use!");
    if(num_snp <= 0){
        std::string error_message = "Number of required SNPs cannot be less than 1: "+std::to_string(num_snp);
        throw std::runtime_error(error_message);
    }
    std::string line;
    //First get the information of the SNPs
    size_t cur_iter = 0;
    for(;m_snp_iter < m_num_snp & cur_iter<num_snp; ++m_snp_iter){
    		std::getline(m_bim, line);
    		misc::trim(line);
    	    	if(!line.empty()){
    	    		std::vector<std::string> token = misc::split(line);
    	    		if(token.size() >= 6){
    	    			m_chr_list.push_back(token[0]);
    	    			//m_snp_list.push_back(token[1]);
    	    			int temp = misc::convert<int>(token[2]);
    	    			if(temp < 0){
    	    				std::string error_message = "Negative CM: "+line;
    	    				throw std::runtime_error(error_message);
    	    			}
    	    			m_cm_list.push_back(temp);
    	    			temp = misc::convert<int>(token[3]);
    	    			if(temp < 0){
    	    				std::string error_message = "Negative BP: "+line;
    	    				throw std::runtime_error(error_message);
    	    			}
    	    			m_bp_list.push_back(temp);
    	    			m_ref_allele.push_back(token[4]);
    	    			m_alt_allele.push_back(token[5]);
    	    		}
    	    		else throw std::runtime_error("Malformed bim file");
    	    	}else throw std::runtime_error("Malformed bim file");
        cur_iter++;
        char genotype_list[m_num_bytes];
        m_bed.read(genotype_list, m_num_bytes);
        size_t i_genotype = 0;
        size_t total_allele = 0;
        size_t num_missing = 0;
        long_type *genotype = new long_type[(m_required_bit /(m_bit_size))+1];
        long_type *missing = new long_type[(m_required_bit /(m_bit_size))+1];
        for(size_t byte_runner= 0; byte_runner < m_num_bytes;){
#ifdef __LP64__
            long_type current_genotypes = 0ULL;
#else
            long_type current_genotypes=0UL;
#endif
            for(int byte_set = 0; byte_set < sizeof(long_type)/sizeof(char) && byte_runner < m_num_bytes; ++byte_set){
            		long_type current_byte = static_cast<long_type>(genotype_list[byte_runner]) << ((sizeof(long_type)-1)*CHAR_BIT) >> CHAR_BIT*byte_set;
                current_genotypes |= current_byte;
                byte_runner++;
            }
            long_type five_masked_geno = current_genotypes & FIVEMASK;
            long_type inter = (five_masked_geno & (current_genotypes>>1)) ^ five_masked_geno;
            long_type current_missing = inter | (inter << 1);
            if(!ld) genotype[i_genotype] = current_genotypes;
            else genotype[i_genotype] = (current_genotypes^(five_masked_geno>>1))&(~current_missing);
            missing[i_genotype] = ~current_missing;
            total_allele += __builtin_popcountll(current_genotypes & (~current_missing));
            num_missing +=__builtin_popcountll(current_missing);
            i_genotype++;
        }
        m_genotype.push_back(genotype);
        m_missing.push_back(missing);
        double maf = (double)total_allele/((double)m_required_bit-(double)num_missing);
        maf = (maf > 0.5)? 1.0-maf: maf;
        m_maf.push_back(maf);

    }
    return m_num_snp-m_snp_iter;
}

void PLINK::initialize(){
    std::string fam_name = m_prefix+".fam";
    std::string bim_name = m_prefix+".bim";
    std::string bed_name = m_prefix+".bed";
    // Start processing the fam file
    std::ifstream fam;
    fam.open(fam_name.c_str());
    if(!fam.is_open()){
        std::string error_message = "Cannot open fam file: "+fam_name;
        throw std::runtime_error(error_message);
    }
    std::string line;
    while(std::getline(fam, line))
        if(!misc::trimmed(line).empty()) m_num_sample++;
    fam.close();
    // Check whether if the bed file is correct
    bool snp_major = openPlinkBinaryFile(bed_name, m_bed);
    if(!snp_major) throw std::runtime_error("Currently does not support sample major format");
    // Check whether if the bim file is correct
    m_bim.open(bim_name.c_str());
    if(!m_bim.is_open()){
        std::string error_message = "Cannot open bim file: "+bim_name;
        throw std::runtime_error(error_message);
    }
    
    while(std::getline(m_bim, line)){
        if(!misc::trimmed(line).empty()) m_num_snp++;
        m_snp_list.push_back(misc::split(line)[1]); // This is dangerous as we don't check bim file format
    }
    m_bim.clear();
    m_bim.seekg(0);
    m_num_bytes=ceil((double)m_num_sample/4.0);
    m_required_bit = m_num_sample*2;
    m_snp_iter=0;
    m_init = true;
}

bool PLINK::openPlinkBinaryFile(const std::string s, std::ifstream & BIT){
    BIT.open(s.c_str(), std::ios::in | std::ios::binary);
    if(!BIT.is_open()){
        std::string error_message= "Cannot open the bed file: "+s;
        throw std::runtime_error(error_message);
    }
    // 1) Check for magic number
    // 2) else check for 0.99 SNP/Ind coding
    // 3) else print warning that file is too old
    char ch[3];
    BIT.read(ch,3);
    bool bfile_SNP_major = false;
    bool v1_bfile = true;
    // If v1.00 file format
    // Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
    // check magic number
    if(static_cast<uint32_t>(ch[0])== 108 && static_cast<uint32_t>(ch[1])==27)
        bfile_SNP_major=(static_cast<uint32_t>(ch[2])==1);
    else
        v1_bfile = false;
    // Reset file if < v1
    if ( ! v1_bfile ) {
        std::cerr << "Warning, old BED file <v1.00 : will try to recover..." << std::endl;
        std::cerr << "  but you should --make-bed from PED )" << std::endl;
        BIT.close();
        BIT.clear();
        BIT.open(s.c_str(), std::ios::in | std::ios::binary);
        BIT.read(ch,1);
        uint32_t file_info = static_cast<uint32_t>(ch[0]);
        if(file_info != 1 && file_info!=0){
            std::cerr << std::endl <<   " *** Possible problem: guessing that BED is < v0.99       *** " << std::endl;
            std::cerr <<                " *** High chance of data corruption, spurious results     *** " << std::endl;
            std::cerr <<                " *** Unless you are _sure_ this really is an old BED file *** " << std::endl;
            std::cerr <<                " *** you should recreate PED -> BED                       *** " << std::endl << std::endl;
            bfile_SNP_major = false;
            BIT.close();
            BIT.clear();
            BIT.open(s.c_str(), std::ios::in | std::ios::binary);
        }
        else{
            std::cerr << "Binary PED file is v0.99" << std::endl;
            bfile_SNP_major = (file_info==1);
        }
    }
    return bfile_SNP_major;
}

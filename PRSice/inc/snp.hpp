#ifndef SNP_H
#define SNP_H

#include <string>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>

#include "commander.hpp"
#include "misc.hpp"

class SNP
{
    public:
        SNP();
#if defined(__LP64__) || defined(_WIN64)
        SNP(const std::string rs_id, const std::string chr, const size_t loc, const std::string ref_allele, const std::string alt_allele, const double statistic, const double se, const double p_value, uint64_t *flag, const size_t size_of_flag);
#else
        SNP(const std::string rs_id, const std::string chr, const size_t loc, const std::string ref_allele, const std::string alt_allele, const double statistic, const double se, const double p_value, uint32_t *flag, const size_t size_of_flag);
#endif
        virtual ~SNP();
        std::string get_ref_allele() { return m_ref_allele; }
        std::string get_alt_allele() { return m_alt_allele; }
        std::string get_rs_id() { return m_rs_id; }
        std::string get_chr() { return m_chr; }
        size_t get_loc() { return m_loc; }
        double set_stat() { return m_stat; }
        double get_p_value() { return m_p_value; }
        static std::vector<size_t> sort_by_p(const std::vector<SNP> &input);
        static std::vector<int> get_index(const Commander &c_commander, const std::string &c_input);
        bool check_loc(const std::string &chr, const size_t loc, const std::string &ref_allele, const std::string &alt_allele){
        		if(chr.compare(m_chr)!=0) return false;
          	if(loc!= m_loc) return false;
          	//Check if allele is the same
          	if(ref_allele.compare(m_ref_allele)!=0 && alt_allele.compare(m_ref_allele)!=0) return false; // not possible even after flipping
         	else if(ref_allele.compare(m_ref_allele)!=0 && alt_allele.compare(m_ref_allele)==0){ // flipping can be performed
         		//flipping here
         		if(!m_alt_allele.empty()){
         			std::string temp = m_alt_allele;
         			m_alt_allele = m_ref_allele;
         			m_ref_allele = temp;
         		}
         		// here, we need to flip the test statistic
         		//TODO: work out the how to flip the test statistic
         	}
          	return true;
        };
        void add_clump(const size_t i) { m_clump_target.push_back(i);};
        void add_clump( std::vector<size_t> &i){m_clump_target.insert( m_clump_target.end(), i.begin(), i.end() );};
        bool clumped() const{return m_clumped;};
        void set_clumped() { m_clumped = true;};
        void clump_all(std::vector<SNP> &snp_list){
        		for(size_t i = 0; i < m_clump_target.size(); ++i){
        			snp_list[i].set_clumped();
        			// update flag
        			for(size_t j = 0; j < m_size_of_flag; ++j) m_flags[j] |= snp_list[i].m_flags[j];
        		}
        }
    protected:
    private:
        static size_t index_check(const std::string &c_in);
        static size_t index_check(const std::string &c_in, const std::vector<std::string> &c_header, const std::string &typeOfError);
        std::string m_ref_allele;
        std::string m_alt_allele;
        std::string m_rs_id;
        std::string m_chr;
        size_t m_loc;
        size_t m_size_of_flag;
        double m_stat;
        double m_standard_error;
        double m_p_value;
        bool m_clumped = false;
        std::vector<size_t> m_clump_target; // index of SNPs that are clumped under this SNP
#if defined(__LP64__) || defined(_WIN64)
        uint64_t *m_flags;
#else
        uint32_t *m_flags;
#endif
};

#endif // SNP_H

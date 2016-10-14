#ifndef SNP_H
#define SNP_H

#include <string>
#include <fstream>
#include <stdexcept>
#include <deque>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <boost/ptr_container/ptr_vector.hpp>

#include "commander.hpp"
#include "misc.hpp"

class SNP
{
    public:
#if defined(__LP64__) || defined(_WIN64)
        typedef uint64_t long_type;
#else
        typedef uint32_t long_type;
#endif
        SNP();
        SNP(const std::string rs_id, const std::string chr, const int loc,
        		const std::string ref_allele, const std::string alt_allele,
				const double statistic, const double se, const double p_value,
				long_type *flag, const size_t size_of_flag);
        virtual ~SNP();
        std::string get_ref_allele() const { return m_ref_allele; }
        std::string get_alt_allele() const { return m_alt_allele; }
        std::string get_rs_id() const { return m_rs_id; }
        std::string get_chr() const { return m_chr; }
        int get_loc() const { return m_loc; }
        double get_stat() const { return m_stat; }
        double get_p_value() const { return m_p_value; }
        static std::vector<size_t> sort_by_p(const boost::ptr_vector<SNP> &input);
        static std::vector<int> get_index(const Commander &c_commander,
        		const std::string &c_input);
        bool check_loc(const std::string &chr, const int loc,
        		const std::string &ref_allele, const std::string &alt_allele){
        		if(chr.compare(m_chr)!=0) return false;
          	if(loc!= m_loc) return false;
          	//Check if allele is the same
          	if(ref_allele.compare(m_ref_allele)!=0 && alt_allele.compare(m_ref_allele)!=0 &&
          			ref_allele.compare(complement(m_ref_allele))!=0 && alt_allele.compare(complement(m_ref_allele))!=0 ) return false; // not possible even after flipping
          	if(m_alt_allele.empty()){
          		// can only use the reference allele to do things, more dangerous
          		if((ref_allele.compare(m_ref_allele)!=0 && alt_allele.compare(m_ref_allele)==0) ||
          				(ref_allele.compare(complement(m_ref_allele))!=0 && alt_allele.compare(complement(m_ref_allele))==0)){
          			m_alt_allele = ref_allele;
          			m_ref_allele = alt_allele;
          			m_flipped=true;
          		}
          	}
          	else{
          		// can use both
          		if((ref_allele.compare(m_alt_allele)==0 && alt_allele.compare(m_ref_allele)==0) ||
          			(ref_allele.compare(complement(m_alt_allele))==0 && alt_allele.compare(complement(m_ref_allele))==0)	){
          			// need to flip
          			m_alt_allele = ref_allele;
          			m_ref_allele = alt_allele;
          			m_flipped=true;
             		//TODO: work out the how to flip the test statistic

          		}
          	}
          	return true;
        };
        void add_clump( std::vector<size_t> &i){ m_clump_target.insert( m_clump_target.end(), i.begin(), i.end() ); };
        void add_clump_r2( std::vector<double> &i){ m_clump_r2.insert( m_clump_r2.end(), i.begin(), i.end() ); };
        bool clumped() const { return m_clumped; };
        bool flipped() const { return m_flipped; };
        // indication of whether if the SNP is within the specific region
        bool in(size_t i) const{
        		size_t index = i/m_bit_size;
        		return (m_flags[i/m_bit_size] >> i%m_bit_size) & ONE; // 1 = true, 0 = false
        }
        void set_clumped() { m_clumped = true;};
        void set_flag(long_type *flag){
        		if(m_flags!=nullptr) delete [] m_flags;
        		m_flags = flag;
        };
        void clump_all(boost::ptr_vector<SNP> &snp_list, double r2_threshold){
        		for(size_t i = 0; i < m_clump_target.size(); ++i){
        			if(!snp_list[m_clump_target[i]].clumped()){
        				snp_list[m_clump_target[i]].set_clumped();
        				if(m_clump_r2[i] >= r2_threshold){
        					for(size_t j = 0; j < m_size_of_flag; ++j)  m_flags[j] |= snp_list[m_clump_target[i]].m_flags[j];
        				}
        			}
        		}
        }
        void clump(boost::ptr_vector<SNP> &snp_list){
        		for(size_t i = 0; i < m_clump_target.size(); ++i){
        			if(!snp_list[m_clump_target[i]].clumped()){
        				int sum_total = 0;
        				for(size_t j = 0; j < m_size_of_flag; ++j){
        						// if there is any overlap this should set the snp_list to the new flag
        						snp_list[m_clump_target[i]].m_flags[j]= snp_list[m_clump_target[i]].m_flags[j] ^
        												(m_flags[j] & snp_list[m_clump_target[i]].m_flags[j]);
        						sum_total+=snp_list[m_clump_target[i]].m_flags[j];
        				}
        				if(sum_total==0)  snp_list[m_clump_target[i]].set_clumped();
        			}
        		}
        }
        double score(int geno) const {
			int g = (geno-1 > 0)? (geno-1) : 0;
			if(!m_flipped) g=2-g;
			return (g>0)? (0.5*(double)g)*m_stat: g;
        }
        static bool sort_snp (const SNP& i, const SNP& j){
            if(i.get_chr().compare(j.get_chr()) == 0)
        		if(i.get_loc() == j.get_loc()) return i.get_rs_id().compare(j.get_rs_id()) < 0;
        		else return i.get_loc() < j.get_loc();
        	else return (i.get_chr().compare(j.get_chr()) < 0);
        }
    protected:
    private:

        std::string complement(std::string allele){
			if(allele.compare("A")==0 || allele.compare("a")==0) return "T";
			if(allele.compare("T")==0 || allele.compare("t")==0) return "A";
			if(allele.compare("G")==0 || allele.compare("g")==0) return "C";
			if(allele.compare("C")==0 || allele.compare("c")==0) return "G";
			else return allele; // Cannot flip, so will just return it as is
        }
        static size_t index_check(const std::string &c_in);
        static size_t index_check(const std::string &c_in,
        		const std::vector<std::string> &c_header, const std::string &typeOfError);
        std::string m_ref_allele;
        std::string m_alt_allele;
        std::string m_rs_id;
        std::string m_chr;
        int m_loc;
        size_t m_size_of_flag;
        size_t m_bit_size;
        double m_stat;
        double m_standard_error;
        double m_p_value;
        bool m_clumped = false;
        bool m_flipped=false;
        std::vector<size_t> m_clump_target; // index of SNPs that are clumped under this SNP
        std::vector<double> m_clump_r2; // index of SNPs that are clumped under this SNP
        long_type *m_flags;
        long_type *m_region_clumped; //place holder for now
#if defined(__LP64__) || defined(_WIN64)
        long_type ONE = 1LLU;
#else
        long_type ONE = 1LU;
#endif
};

#endif // SNP_H

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
#include "storage.hpp"
#include "commander.hpp"
#include "misc.hpp"

class SNP
{
public:
#if defined(__LP64__) || defined(_WIN64)
    typedef uint64_t long_type;
	#define ONE  0x1LLU

#else
    typedef uint32_t long_type;
	#define ONE  0x1LU
#endif

    SNP();
    SNP(const std::string rs_id, const std::string chr, const int loc,
    		const std::string ref_allele, const std::string alt_allele,
			const double statistic, const double se, const double p_value);
    virtual ~SNP();


    // for future me, plink encoding is 00 hom alt, 10 het, 11 hom ref and 01 missing
    // in binary. So by converting into integer, we've got 0, 2, 3, 1
    // and we should filter out 1 before anyway
    inline int geno(int geno_in) const
    {
    		int g = (geno_in -1 > 0)? (geno_in -1) : 0;
    		if(!m_flipped) g= 2-g;
    		return g;
    }
    inline double score(int geno) const
    {
        int g = (geno-1 > 0)? (geno-1) : 0;
        if(!m_flipped) g=2-g;
        return (g>0)? (0.5*(double)g)*m_stat: g;
    }
    // location checker
    void set_loc(int loc) { m_loc = loc; };
    bool check_loc(const std::string &chr, const int loc, const std::string &ref_allele, const std::string &alt_allele);

    // clumping
    void add_clump( std::vector<size_t> &i) { m_clump_target.insert( m_clump_target.end(), i.begin(), i.end() ); };
    void add_clump_r2( std::vector<double> &i) { m_clump_r2.insert( m_clump_r2.end(), i.begin(), i.end() ); };
    void set_clumped() { m_clumped = true; };
    void proxy_clump(boost::ptr_vector<SNP> &snp_list, double r2_threshold);
    void clump(boost::ptr_vector<SNP> &snp_list);

    // region
    inline bool in(size_t i) const
    {
    		if(i/m_bit_size >= m_flags.size()) throw std::out_of_range("Out of range for flag");
        return (m_flags[i/m_bit_size] >> i%m_bit_size) & ONE; // 1 = true, 0 = false
    }
    void set_flag(std::vector<long_type> flag) { m_flags = flag; };

    // sorting
    static std::vector<size_t> sort_by_p(const boost::ptr_vector<SNP> &input);
    bool operator < (const SNP& j) const
    {
        if(m_chr.compare(j.get_chr()) == 0)
            return (m_loc == j.get_loc())? (m_rs.compare(j.get_rs_id()) < 0) : (m_loc < j.get_loc());
        else return (m_chr.compare(j.get_chr()) < 0);
    }
    bool operator == (const SNP& j) const { return m_rs.compare(j.get_rs_id()) == 0; }

    // getter
    std::string get_ref_allele() const { return m_ref; };
    std::string get_alt_allele() const { return m_alt; };
    std::string get_rs_id() const { return m_rs; };
    std::string get_chr() const { return m_chr; };
    int get_loc() const { return m_loc; };
    double get_stat() const { return m_stat; };
    double get_p_value() const { return m_p_value; };
    bool clumped() const { return m_clumped; };
    bool flipped() const { return m_flipped; };

    // header check
    static std::vector<int> get_index(const Commander &c_commander, const std::string &c_input);
    static bool valid_snp(std::string allele)
    {
    		if(allele.compare("A")==0 || allele.compare("a")==0) return true;
    		else if(allele.compare("C")==0 || allele.compare("c")==0) return true;
    		else if(allele.compare("T")==0 || allele.compare("t")==0) return true;
    		else if(allele.compare("G")==0 || allele.compare("g")==0) return true;
    		else if(allele.compare("1")==0 || allele.compare("2")==0) return true;
    		else return false;
    }
    static bool ambiguous(std::string ref_allele, std::string alt_allele)
    {
    		return (ref_allele == "A" && alt_allele == "T")
    				|| (ref_allele == "a" && alt_allele == "t")
					|| (ref_allele == "G" && alt_allele == "C")
					|| (ref_allele == "g" && alt_allele == "c");
    }
protected:
private:
    //basic info
    std::string m_ref;
    std::string m_alt;
    std::string m_rs;
    std::string m_chr;
    int m_loc;
    double m_stat;
    double m_se;
    double m_p_value;
    //clump related
    bool m_clumped;
    std::vector<size_t> m_clump_target; // index of SNPs that are clumped under this SNP
    std::vector<double> m_clump_r2; // index of SNPs that are clumped under this SNP
    //region related
    size_t m_bit_size;
    std::vector<long_type> m_flags;
    //others
    bool m_flipped;

    //functions
    std::string complement(std::string allele)
    {
        if(allele.compare("A")==0 || allele.compare("a")==0) return "T";
        if(allele.compare("T")==0 || allele.compare("t")==0) return "A";
        if(allele.compare("G")==0 || allele.compare("g")==0) return "C";
        if(allele.compare("C")==0 || allele.compare("c")==0) return "G";
        else return allele; // Cannot flip, so will just return it as is
    }

    //private function headers
    static size_t index_check(const std::string &c_in);
    static size_t index_check(const std::string &c_in,
                              const std::vector<std::string> &c_header, const std::string &typeOfError);


//    size_t m_size_of_flag;
//    std::vector<long_type> m_region_clumped; //place holder for now

};

#endif // SNP_H

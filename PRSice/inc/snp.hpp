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
#else
    typedef uint32_t long_type;
#endif
    SNP();
    SNP(const std::string rs_id, const std::string chr, const int loc,
        const std::string ref_allele, const std::string alt_allele,
        const double statistic, const double se, const double p_value,
        long_type *flag, const size_t size_of_flag);
    virtual ~SNP();
    static std::vector<size_t> sort_by_p(const boost::ptr_vector<SNP> &input);
    static std::vector<int> get_index(const Commander &c_commander, const std::string &c_input);

    std::string get_ref_allele() const
    {
        return m_ref_allele;
    }
    std::string get_alt_allele() const
    {
        return m_alt_allele;
    }
    std::string get_rs_id() const
    {
        return m_rs_id;
    }
    std::string get_chr() const
    {
        return m_chr;
    }
    int get_loc() const
    {
        return m_loc;
    }
    double get_stat() const
    {
        return m_stat;
    }
    double get_p_value() const
    {
        return m_p_value;
    }
    inline double score(int geno) const
    {
        int g = (geno-1 > 0)? (geno-1) : 0;
        if(!m_flipped) g=2-g;
        return (g>0)? (0.5*(double)g)*m_stat: g;
    }

    void add_clump( std::vector<size_t> &i)
    {
        m_clump_target.insert( m_clump_target.end(), i.begin(), i.end() );
    };
    void add_clump_r2( std::vector<double> &i)
    {
        m_clump_r2.insert( m_clump_r2.end(), i.begin(), i.end() );
    };
    void set_clumped()
    {
        m_clumped = true;
    };
    void set_flag(long_type *flag)
    {
        if(m_flags!=nullptr) delete [] m_flags;
        m_flags = flag;
    };
    void clump_all(boost::ptr_vector<SNP> &snp_list, double r2_threshold);
    void clump(boost::ptr_vector<SNP> &snp_list);

    void set_loc(int loc)
    {
        m_loc = loc;
    };
    bool check_loc(const std::string &chr, const int loc, const std::string &ref_allele,
                   const std::string &alt_allele);
    bool clumped() const
    {
        return m_clumped;
    };
    bool flipped() const
    {
        return m_flipped;
    };
    // indication of whether if the SNP is within the specific region
    inline bool in(size_t i) const
    {
        size_t index = i/m_bit_size;
        return (m_flags[i/m_bit_size] >> i%m_bit_size) & ONE; // 1 = true, 0 = false
    }
    bool operator < (const SNP& j) const
    {
        if(m_chr.compare(j.get_chr()) == 0)
            if(m_loc == j.get_loc()) return m_rs_id.compare(j.get_rs_id()) < 0;
            else return m_loc < j.get_loc();
        else return (m_chr.compare(j.get_chr()) < 0);
    }

    bool operator == (const SNP& j) const
    {
        return m_rs_id.compare(j.get_rs_id()) == 0;
    }


protected:
private:
    std::string complement(std::string allele)
    {
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

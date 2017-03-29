#ifndef SNP_H
#define SNP_H

#include <string>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include "storage.hpp"
#include "commander.hpp"
#include "misc.hpp"

class SNP
{
public:
    SNP();
    SNP(const std::string rs_id, const int chr, const int loc,
    		const std::string ref_allele, const std::string alt_allele);
    SNP(const std::string rs_id, const int chr, const int loc,
    		const std::string ref_allele, const std::string alt_allele,
			const int bound_start, const int bound_end);
    SNP(const std::string rs_id, const int chr, const int loc,
    		const std::string ref_allele, const std::string alt_allele,
			const double statistic, const double se, const double p_value,
			const int category, const double p_threshold);
    virtual ~SNP();

    static bool ambiguous(std::string ref_allele, std::string alt_allele)
    {
    	return (ref_allele == "A" && alt_allele == "T")
    			|| (ref_allele == "a" && alt_allele == "t")
				|| (ref_allele == "G" && alt_allele == "C")
				|| (ref_allele == "g" && alt_allele == "c");
    };

    void not_required() { m_required=false; };
    void required() { m_required=true; };
    void set_statistic(const double statistic, const double se, const double p_value, const int category,
    		const double p_threshold)
    {
    	m_stat = statistic;
    	m_se = se;
    	m_p_value = p_value;
    	m_category = category;
    	m_p_threshold = p_threshold;
    };
    bool is_required() const { return m_required; };
    std::string get_rs() const { return m_rs; };


    bool operator == (const SNP &Ref) const
    {
    	if(m_chr == Ref.m_chr && m_loc==Ref.m_loc && m_rs.compare(Ref.m_rs)==0)
    	{
    		if(m_ref.compare(Ref.m_ref)==0){
    			if(!m_alt.empty() && !Ref.m_alt.empty())
    			{
    				return m_alt.compare(Ref.m_alt)==0;
    			}else return true;
    		}
    		else if(complement(m_ref).compare(Ref.m_ref)==0)
    		{
    			if(!m_alt.empty() && !Ref.m_alt.empty())
    			{
    				return complement(m_alt).compare(Ref.m_alt)==0;
    			}else return true;
    		}
    		else if(!m_alt.empty() && !Ref.m_alt.empty())
    		{
    			if(m_ref.compare(Ref.m_alt)==0 && m_alt.compare(Ref.m_ref)==0)
    			{
    				return true;
    			}
    			if(complement(m_ref).compare(Ref.m_alt)==0 && complement(m_alt).compare(Ref.m_alt)==0)
    			{
    				return true;
    			}
    			return false;
    		}
    		else return false; // cannot flip nor match
    	}
    	else
    	{
    		return false;
    	}
    };

private:
    //basic info
    std::string m_ref;
    std::string m_alt;
    std::string m_rs;
    int m_chr;
    int m_loc;
    int m_category;
    // This indicate where this SNP's bound is at
    // useful for PRSlice and also clumping
    // thinking about it. Even if the location isn't given for
    // PRSet or PRSlice, we can still use the coordinates from
    // the target / reference file
    int m_bound_start;
    int m_bound_end;
    double m_stat;
    double m_se;
    double m_p_value;
    double m_p_threshold;
    //clump related
    bool m_clumped;
    bool m_required;
    inline std::string complement(const std::string &allele) const
    {
    	if(allele.compare("A")==0 || allele.compare("a")==0) return "T";
    	if(allele.compare("T")==0 || allele.compare("t")==0) return "A";
    	if(allele.compare("G")==0 || allele.compare("g")==0) return "C";
    	if(allele.compare("C")==0 || allele.compare("c")==0) return "G";
    	else return allele; // Cannot flip, so will just return it as is
    }

};

#endif // SNP_H

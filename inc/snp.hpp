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

    void set_statistic(const double stat, const double se, const double p_value, const int category,
    		const double p_threshold)
    {
    		statistic.stat = stat;
    		statistic.se = se;
    		statistic.p_value = p_value;
    		threshold.category = category;
    		threshold.p_threshold = p_threshold;
    };
    void set_flipped() { statistic.flipped = true; };
    std::string get_rs() const { return basic.rs; };
    int range_start() const { return m_range_start; };
    int range_end() const { return m_range_end; };
    static std::vector<size_t> sort_by_p(const std::vector<SNP> &input);

    bool operator == (const SNP &Ref) const
    {
    	if(basic.chr == Ref.basic.chr && basic.loc==Ref.basic.loc && basic.rs.compare(Ref.basic.rs)==0)
    	{
    		if(basic.ref.compare(Ref.basic.ref)==0){
    			if(!basic.alt.empty() && !Ref.basic.alt.empty())
    			{
    				return basic.alt.compare(Ref.basic.alt)==0;
    			}else return true;
    		}
    		else if(complement(basic.ref).compare(Ref.basic.ref)==0)
    		{
    			if(!basic.alt.empty() && !Ref.basic.alt.empty())
    			{
    				return complement(basic.alt).compare(Ref.basic.alt)==0;
    			}else return true;
    		}
    		else if(!basic.alt.empty() && !Ref.basic.alt.empty())
    		{
    			if(basic.ref.compare(Ref.basic.alt)==0 && basic.alt.compare(Ref.basic.ref)==0)
    			{
    				return true;
    			}
    			if(complement(basic.ref).compare(Ref.basic.alt)==0 &&
    					complement(basic.alt).compare(Ref.basic.alt)==0)
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

    inline void flipped(){ statistic.flipped = true; };
    inline void fill_info(int chr, int loc, std::string alt)
    {
    		if(basic.chr==-1) basic.chr = chr;
    		if(basic.loc==-1) basic.loc = loc;
    		if(basic.alt.empty()) basic.alt=alt;
    };
    inline bool matching (int chr, int loc, std::string ref, std::string alt, bool &flipped) const{
    		if(chr != -1 && basic.chr != -1 && chr != basic.chr) return false;
    		if(loc != -1 && basic.loc != -1 && loc != basic.loc) return false;
    		if(basic.ref.compare(ref)==0){
    			if(!basic.alt.empty() && !alt.empty())
    			{
    				return basic.alt.compare(alt)==0;
    			}else return true;
    		}
    		else if(complement(basic.ref).compare(ref)==0)
    		{
    			if(!basic.alt.empty() && !alt.empty())
    			{
    				return complement(basic.alt).compare(alt)==0;
    			}else return true;
    		}
    		else if(!basic.alt.empty() && !alt.empty())
    		{
    			if(basic.ref.compare(alt)==0 && basic.alt.compare(ref)==0)
    			{
    				flipped = true;
    				return true;
    			}
    			if(complement(basic.ref).compare(alt)==0 && complement(basic.alt).compare(alt)==0)
    			{
    				flipped = true;
    				return true;
    			}
    			return false;
    		}
    		else return false; // cannot flip nor match
    };

    int chr() const { return basic.chr; };
    int loc() const { return basic.loc; };
    int snp_id() const { return file_info.id; };
    double p_value() const { return statistic.p_value; };
    std::string rs() const { return basic.rs; };
    std::string ref() const { return basic.ref; };
    std::string alt() const { return basic.alt; };
    void set_upper(int upper){ m_range_end = upper; };
    void set_lower(int lower){ m_range_start = lower;};
    void set_flag(std::vector<long_type> flag) { m_flags = flag; };
    void not_required(){ m_required = false; };

    std::string file_name() const { return file_info.file; };
    bool is_required() const { return m_required; };
private:
    //basic info
    struct{
    	std::string ref;
    	std::string alt;
    	std::string rs;
    	int chr;
    	int loc;
    }basic;

    struct{
    	std::string file;
    	int id;
    }file_info;

    struct{
    	double stat;
        double se;
        double p_value;
        bool flipped;
    }statistic;

    struct{
    	int category;
        double p_threshold;
    }threshold;

    // This indicate where this SNP's bound is at
    // useful for PRSlice and also clumping
    // thinking about it. Even if the location isn't given for
    // PRSet or PRSlice, we can still use the coordinates from
    // the target / reference file
    // the bound is [ )
    bool m_required = true;
    int m_range_start;
    int m_range_end;
    //clump related
    bool m_clumped;
    //prset related
    std::vector<long_type> m_flags;

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

#include "snp.hpp"

SNP::SNP()
{
	basic.chr=-1;
	basic.loc=-1; // default is -1 to indicate that it is not provided
    statistic.flipped = false;
    statistic.stat=0.0;
    statistic.se=0.0;
    statistic.p_value=2.0; // set this to 2 such that only SNPs in base file have valid P-value
    threshold.p_threshold=0.0;
    threshold.category=0;
    m_clumped = false;
    m_bit_size = sizeof(long_type)*CHAR_BIT;
}


SNP::SNP(const std::string rs_id, const std::string chr, const int loc,
		const std::string ref_allele, const std::string alt_allele,
		)
	: basic.ref(ref_allele), basic.alt(alt_allele),
	  basic.rs(rs_id), basic.chr(chr), basic.loc(loc)
{
	statistic.se = 0.0;
	statistic.p_value = 0.0;
	statistic.stat = 0.0;
	statistic.flipped = false;
    threshold.p_threshold=0.0;
    threshold.category=0;
    m_bit_size = sizeof(long_type)*CHAR_BIT;
    m_clumped = false;
}


SNP::SNP(const std::string rs_id, const int chr, const int loc,
		const std::string ref_allele, const std::string alt_allele,
		const int range_start, const int range_end)
	: SNP(rs_id, chr, loc, ref_allele, alt_allele), m_range_start(range_start), m_range_end(range_end){}

SNP::SNP(const std::string rs_id, const std::string chr, const int loc,
		const std::string ref_allele, const std::string alt_allele,
		const double statistic, const double se, const double p_value, const int category,
		const double p_threshold)
	: SNP(rs_id, chr, loc, ref_allele, alt_allele), statistic.stat(statistic),
	  statistic.se(se), statistic.p_value(p_value),
	  threshold.p_threshold(p_threshold), threshold.category(category){}

SNP::~SNP(){}

std::vector<size_t> SNP::sort_by_p(const std::vector<SNP> &input)
{
    std::vector<size_t> idx(input.size());
    std::iota(idx.begin(), idx.end(),0);
    std::sort(idx.begin(), idx.end(), [&input](size_t i1, size_t i2)
    {
    	// plink do it with respect to the location instead of statistic
        if(input[i1].statistic.p_value==input[i2].statistic.p_value)
        {
        	if(input[i1].basic.chr==input[i2].basic.chr)
        	{
        		if(input[i1].basic.loc == input[i2].basic.loc)
        		{
        			if(fabs(input[i1].statistic.stat)==fabs(input[i2].statistic.stat))
        			{
        				return input[i1].statistic.se < input[i2].statistic.se;
        			}
        			else return fabs(input[i1].statistic.stat) > fabs(input[2].statistic.stat);
        		}
        		else return input[i1].basic.loc < input[i2].basic.loc;
        	}
			else return input[i1].basic.chr<input[i2].basic.chr;
        }
        else return input[i1].statistic.p_value < input[i2].statistic.p_value;
    });
    return idx;
}


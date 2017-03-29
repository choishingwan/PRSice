#include "snp.hpp"

SNP::SNP()
{
    m_loc=-1; // default is -1 to indicate that it is not provided
    m_stat=0.0;
    m_se=0.0;
    m_p_value=0.0;
    m_p_threshold=0.0;
    m_category=0;
    m_flipped = false;
    m_clumped = false;
    m_required=true;
    m_bit_size = sizeof(long_type)*CHAR_BIT;
}


SNP::SNP(const std::string rs_id, const std::string chr, const int loc,
		const std::string ref_allele, const std::string alt_allele,
		)
	: m_ref(ref_allele), m_alt(alt_allele),
		m_rs(rs_id), m_chr(chr), m_loc(loc)
{
	m_se = 0.0;
	m_p_value = 0.0;
	m_stat = 0.0;
    m_p_threshold=0.0;
    m_category=0;
    m_bit_size = sizeof(long_type)*CHAR_BIT;
    m_clumped = false;
    m_required=true;
}


SNP::SNP(const std::string rs_id, const int chr, const int loc,
		const std::string ref_allele, const std::string alt_allele,
		const int bound_start, const int bound_end)
	: SNP(rs_id, chr, loc, ref_allele, alt_allele), m_bound_start(bound_start), m_bound_end(bound_end){}

SNP::SNP(const std::string rs_id, const std::string chr, const int loc,
		const std::string ref_allele, const std::string alt_allele,
		const double statistic, const double se, const double p_value, const int category,
		const double p_threshold)
	: SNP(rs_id, chr, loc, ref_allele, alt_allele),  m_category(category), m_stat(statistic),
		m_se(se), m_p_value(p_value), m_p_threshold(p_threshold){}

SNP::~SNP(){}

std::vector<size_t> SNP::sort_by_p(const std::vector<SNP> &input)
{
    std::vector<size_t> idx(input.size());
    std::iota(idx.begin(), idx.end(),0);
    std::sort(idx.begin(), idx.end(), [&input](size_t i1, size_t i2)
    {
    	// plink do it with respect to the location instead of statistic
        if(input[i1].m_p_value==input[i2].m_p_value)
        {
        	if(input[i1].m_chr==input[i2].m_chr)
        	{
        		if(input[i1].m_loc == input[i2].m_loc)
        		{
        			if(fabs(input[i1].m_stat)==fabs(input[i2].m_stat))
        			{
        				return input[i1].m_se < input[i2].m_se;
        			}
        			else return fabs(input[i1].m_stat) > fabs(input[2].m_stat);
        		}
        		else return input[i1].m_loc < input[i2].m_loc;
        	}
			else return input[i1].m_chr<input[i2].m_chr;
        }
        else return input[i1].m_p_value < input[i2].m_p_value;
    });
    return idx;
}


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
    m_bit_size = sizeof(long_type)*CHAR_BIT;
    clump_info.clumped=false;
	clump_info.genotype = nullptr;
	clump_info.contain_missing = false;
}


SNP::SNP(const std::string rs_id, const int chr, const int loc,
		const std::string ref_allele, const std::string alt_allele,
		const std::string file_name, const int num_line)
{
	basic.ref=ref_allele;
	basic.alt=alt_allele;
	basic.rs=rs_id;
	basic.chr=chr;
	basic.loc=loc;
	statistic.se = 0.0;
	statistic.p_value = 0.0;
	statistic.stat = 0.0;
	statistic.flipped = false;
    threshold.p_threshold=0.0;
    threshold.category=0;
    m_bit_size = sizeof(long_type)*CHAR_BIT;
    clump_info.clumped=false;
	clump_info.genotype = nullptr;
	clump_info.contain_missing = false;
    file_info.file = file_name;
    file_info.id = num_line;
}




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


void SNP::clump(std::vector<SNP> &snp_list)
{
	for(auto &&target : clump_info.target){
		if(!snp_list[target].clumped())
		{
			int sum_total = 0;
			for(size_t i_flag = 0; i_flag < m_flags.size(); ++i_flag)
			{
				// if there is any overlap this should set the snp_list to the new flag
				snp_list[target].m_flags[i_flag] = snp_list[target].m_flags[i_flag] ^
						(m_flags[i_flag] & snp_list[target].m_flags[i_flag]);
				sum_total+=snp_list[target].m_flags[i_flag];
			}
			if(sum_total==0)  snp_list[target].set_clumped();
		}
	}
	clump_info.clumped=true; // protect from other SNPs tempering its flags
}

void SNP::proxy_clump(std::vector<SNP> &snp_list, double r2_threshold)
{
    for(size_t i_target = 0; i_target < clump_info.target.size(); ++i_target)
    {
        if(!snp_list[clump_info.target[i_target]].clumped())
        {
            snp_list[clump_info.target[i_target]].set_clumped();
            if(clump_info.r2[i_target] >= r2_threshold)
            {
                for(size_t j = 0; j < m_flags.size(); ++j)  m_flags[j] |= snp_list[clump_info.target[i_target]].m_flags[j];
            }
        }
    }
	clump_info.clumped=true; // protect from other SNPs tempering its flags
}

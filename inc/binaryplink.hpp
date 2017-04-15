

#ifndef BINARYPLINK
#define BINARYPLINK

#include "genotype.hpp"

class BinaryPlink: public Genotype{
public:
	BinaryPlink(std::string prefix,  std::string remove_sample, std::string keep_sample,
			bool ignore_fid, int num_auto=22, bool no_x=false, bool no_y=false, bool no_xy=false,
			bool no_mt=false, const size_t thread=1, bool verbose=false);
	 ~BinaryPlink();
private:
	uintptr_t m_bed_offset = 3;
	std::vector<Sample> load_samples(bool ignore_fid);
	std::vector<SNP> load_snps();
	void check_bed();

	void cleanup(){
		fclose(m_bedfile);
		m_bedfile=nullptr;
	};
	inline void read_genotype(uintptr_t* genotype, const uint32_t snp_index, const std::string &file_name)
	{
		if(m_cur_file.empty() || m_cur_file.compare(file_name)!=0)
		{
			if (m_bedfile != nullptr)
			{
				fclose(m_bedfile);
				m_bedfile=nullptr;
			}
			std::string bedname = file_name+".bed";
			m_bedfile = fopen(bedname.c_str(), FOPEN_RB); // assume there is no error
		}
		if (fseeko(m_bedfile, m_bed_offset + (snp_index* ((uint64_t)m_unfiltered_sample_ct4))
				, SEEK_SET))
		{
			throw std::runtime_error("ERROR: Cannot read the bed file!");
		}
		std::memset(m_tmp_genotype, 0x0, m_unfiltered_sample_ctl*2*sizeof(uintptr_t));
		if(load_and_collapse_incl(m_unfiltered_sample_ct, m_founder_ct, m_founder_info, m_final_mask,
				false, m_bedfile, m_tmp_genotype, genotype))
		{
			throw std::runtime_error("ERROR: Cannot read the bed file!");
		}
	}

	void read_score(std::vector< std::vector<Sample_lite> > &current_prs_score, size_t start_index, size_t end_bound);
	FILE* m_bedfile = nullptr;
	std::string m_cur_file;
	uintptr_t m_final_mask;
	uintptr_t *m_tmp_genotype;
};

#endif

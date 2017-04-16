#ifndef BGEN_H
#define BGEN_H
#include "genotype.hpp"

class BGEN: public Genotype
{
    public:
        BGEN(std::string prefix, std::string remove_sample,
                std::string keep_sample, bool ignore_fid, int num_auto = 22,
                bool no_x = false, bool no_y = false, bool no_xy = false,
                bool no_mt = false, const size_t thread = 1, bool verbose =
                        false);
        ~BGEN();
    private:

        std::vector<Sample> load_samples(bool ignore_fid);
        std::vector<SNP> load_snps();
        void cleanup()
        {

        }
        ;
        inline void read_genotype(uintptr_t* genotype, const uint32_t snp_index,
                const std::string &file_name)
        {
        }

        void read_score(
                std::vector<std::vector<Sample_lite> > &current_prs_score,
                size_t start_index, size_t end_bound);
};
#endif

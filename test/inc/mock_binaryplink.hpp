#ifndef MOCK_BINARYPLINK_HPP
#define MOCK_BINARYPLINK_HPP
#include "binaryplink.hpp"
#include "genotype.hpp"
#include <bitset>

class mock_binaryplink : public ::BinaryPlink
{
public:
    mock_binaryplink(GenoFile& geno, Phenotype& pheno, const std::string& delim,
                     Reporter* reporter)
        : BinaryPlink(geno, pheno, delim, reporter)
    {
    }
mock_binaryplink(){}

    std::vector<uintptr_t> sample_for_ld() const { return m_sample_for_ld; }
    std::vector<uintptr_t> calculate_prs() const { return m_calculate_prs; }
    std::vector<Sample_ID> sample_id() const { return m_sample_id; }
    void test_check_bed(const std::string& bed_name, size_t num_marker,
                   uintptr_t& bed_offset){
        check_bed(bed_name, num_marker, bed_offset);
    }

    void set_sample(uintptr_t n_sample) { m_unfiltered_sample_ct = n_sample; }
    void gen_bed_head(const std::string& name, size_t num_sample,
                      size_t num_snp, bool new_version, bool sample_major)
    {
        std::ofstream out(name.c_str(), std::ios::out | std::ios::binary);
        uintptr_t unfiltered_sample_ct4 = (num_sample + 3) / 4;
        std::bitset<8> b;
        char ch[1];
        if(new_version){
            b.reset();
            b.set(2);
            b.set(3);
            b.set(5);
            b.set(6);
            ch[0] = static_cast<char>( b.to_ulong());
            out.write(ch, 1);
            b.reset();
            b.set(0);
            b.set(1);
            b.set(3);
            b.set(4);
            ch[0] =static_cast<char>( b.to_ulong());
            out.write(ch, 1);
        }
        b.reset();
        if (!sample_major) { b.set(0); }
        ch[0] =static_cast<char>( b.to_ulong());
        out.write(ch, 1);
        for(size_t i = 0; i < unfiltered_sample_ct4*num_snp;++i){
            out.write(ch, 1);
        }
        out.close();
    }
};

#endif // MOCK_BINARYPLINK_HPP

#ifndef MOCK_BINARYPLINK_HPP
#define MOCK_BINARYPLINK_HPP
#include "binaryplink.hpp"
#include "genotype.hpp"

class mock_binaryplink: public ::BinaryPlink
{
public:
    mock_binaryplink(GenoFile& geno, Phenotype&pheno, const std::string &delim, Reporter* reporter): BinaryPlink(geno, pheno, delim, reporter){}


    std::vector<uintptr_t> sample_for_ld() const { return m_sample_for_ld; }
    std::vector<uintptr_t> calculate_prs() const { return m_calculate_prs; }
    std::vector<Sample_ID> sample_id() const { return m_sample_id;}
};

#endif // MOCK_BINARYPLINK_HPP

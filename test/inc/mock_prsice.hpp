#ifndef MOCK_PRSICE_HPP
#define MOCK_PRSICE_HPP
#include "catch.hpp"
#include "prsice.hpp"
class mock_prsice : public PRSice
{
public:
    mock_prsice(Reporter* reporter) { set_reporter(reporter); }
    mock_prsice() {};
    mock_prsice(const CalculatePRS& prs_info, const PThresholding& p_info,
                const Phenotype& pheno, const Permutations& perm,
                const std::string& output, Reporter* reporter)
        : PRSice(prs_info, p_info, pheno, perm, output, reporter)
    {
    }
    void set_reporter(Reporter* reporter) { m_reporter = reporter; }

    std::tuple<size_t, bool>
    test_get_pheno_idx(const std::vector<std::string_view>& column,
                       const std::string& pheno)
    {
        return get_pheno_idx(column, pheno);
    }
    void test_parse_pheno_header(std::unique_ptr<std::istream> pheno_file)
    {
        parse_pheno_header(std::move(pheno_file));
    }
    Phenotype get_pheno_info() const { return m_pheno_info; }
};

#endif // MOCK_PRSICE_HPP

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

    static std::tuple<size_t, bool>
    test_get_pheno_idx(const std::vector<std::string_view>& column,
                       const Phenotype& pheno_info, const std::string& pheno)
    {
        return PRSice::get_pheno_idx(column, pheno_info, pheno);
    }
    static void
    test_parse_pheno_header(std::unique_ptr<std::istream> pheno_file,
                            Phenotype& pheno_info, Reporter& reporter)
    {
        parse_pheno_header(std::move(pheno_file), pheno_info, reporter);
    }
    Phenotype get_pheno_info() const { return m_pheno_info; }
    std::tuple<size_t, size_t> get_progress()
    {
        return {m_total_process, m_total_competitive_process};
    }
    std::tuple<size_t, size_t> get_current_progress()
    {
        return {m_analysis_done, m_total_competitive_perm_done};
    }
    void set_progress(size_t analysis_done, size_t competitive_done)
    {
        m_analysis_done = analysis_done;
        m_total_competitive_perm_done = competitive_done;
    }
};

#endif // MOCK_PRSICE_HPP

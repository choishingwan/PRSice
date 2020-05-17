#ifndef MOCK_PRSICE_HPP
#define MOCK_PRSICE_HPP
#include "catch.hpp"
#include "prsice.hpp"
class mock_prsice : public PRSice
{
public:
    mock_prsice(Reporter* reporter) { set_reporter(reporter); }
    mock_prsice(const bool binary, Reporter* reporter)
        : PRSice(CalculatePRS(), PThresholding(), Permutations(), "PRSice", 3,
                 3, binary, reporter)
    {
    }
    mock_prsice() {};
    mock_prsice(const CalculatePRS& prs_info, const PThresholding& p_info,
                const Permutations& perm, const std::string& output,
                const size_t max_fid, const size_t max_iid, const bool binary,
                Reporter* reporter)
        : PRSice(prs_info, p_info, perm, output, max_fid, max_iid, binary,
                 reporter)
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
    std::unordered_map<std::string, std::string>
    test_load_pheno_map(const std::string& delim, const size_t idx,
                        const bool ignore_fid,
                        std::unique_ptr<std::istream> pheno_file)
    {
        return load_pheno_map(delim, idx, ignore_fid, std::move(pheno_file));
    }
    void test_parse_pheno(const std::string& pheno,
                          std::vector<double>& pheno_store, int& max_pheno_code)
    {
        parse_pheno(pheno, pheno_store, max_pheno_code);
    }
    bool
    test_quantitative_pheno_is_valid(const std::vector<double>& pheno_store)
    {
        return quantitative_pheno_is_valid(pheno_store);
    }

    std::tuple<bool, size_t, size_t>
    test_binary_pheno_is_valid(const int max_pheno_code,
                               std::vector<double>& pheno_store)
    {
        return binary_pheno_is_valid(max_pheno_code, pheno_store);
    }


    std::tuple<std::vector<double>, size_t, int>
    test_process_phenotype_info(const std::string& delim, const bool ignore_fid,
                                Genotype& target)
    {
        return process_phenotype_info(delim, ignore_fid, target);
    }
    std::unordered_map<std::string, size_t>& sample_with_phenotypes()
    {
        return m_sample_with_phenotypes;
    }

    std::tuple<std::vector<double>, size_t, size_t, int>
    test_process_phenotype_file(const std::string& file_name,
                                const std::string& delim,
                                const std::size_t pheno_idx,
                                const bool ignore_fid, Genotype& target)
    {
        return process_phenotype_file(file_name, delim, pheno_idx, ignore_fid,
                                      target);
    }
    void test_print_pheno_log(const std::string& name, const size_t sample_ct,
                              const size_t num_not_found,
                              const size_t invalid_pheno,
                              const int max_pheno_code, const bool ignore_fid,
                              std::vector<double>& pheno_store)
    {
        print_pheno_log(name, sample_ct, num_not_found, invalid_pheno,
                        max_pheno_code, ignore_fid, pheno_store);
    }
    void test_gen_pheno_vec(const std::string& pheno_file,
                            const std::string& pheno_name,
                            const std::string& delim, const size_t pheno_idx,
                            const bool ignore_fid, Genotype& target)
    {
        gen_pheno_vec(pheno_file, pheno_name, delim, pheno_idx, ignore_fid,
                      target);
    }
    static bool test_is_valid_covariate(const std::set<size_t>& factor_idx,
                                        const std::vector<size_t>& cov_idx,
                                        std::vector<std::string>& cov_line,
                                        std::vector<size_t>& missing_count)
    {
        PRSice prsice;
        return prsice.is_valid_covariate(factor_idx, cov_idx, cov_line,
                                         missing_count);
    }
    Eigen::VectorXd& phenotype_matrix() { return m_phenotype; }
    void test_update_phenotype_matrix(const std::vector<bool>& valid_samples,
                                      const size_t num_valid, Genotype& target)
    {
        update_phenotype_matrix(valid_samples, num_valid, target);
    }
    std::vector<std::unordered_map<std::string, size_t>>
    test_cov_check_and_factor_level_count(
        const std::set<size_t>& factor_idx,
        const std::vector<std::string>& cov_names,
        const std::vector<size_t>& cov_idx, const std::string& delim,
        const bool ignore_fid, std::unique_ptr<std::istream>& cov_file,
        Genotype& target)
    {
        return cov_check_and_factor_level_count(factor_idx, cov_names, cov_idx,
                                                delim, ignore_fid, cov_file,
                                                target);
    }
};

#endif // MOCK_PRSICE_HPP

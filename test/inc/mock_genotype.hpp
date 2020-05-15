#ifndef MOCK_GENOTYPE_H
#define MOCK_GENOTYPE_H

#include "genotype.hpp"
#include "reporter.hpp"
#include <memory>

class mockGenotype : public Genotype
{
public:
    std::string test_initialize(const GenoFile& geno, const Phenotype& pheno,
                                const std::string& delim,
                                const std::string& type, Reporter* reporter)
    {
        return initialize(geno, pheno, delim, type, reporter);
    }
    std::vector<size_t> sorted_p_index() const { return m_sort_by_p_index; }
    std::vector<std::string>
    test_load_genotype_prefix(std::unique_ptr<std::istringstream> in)
    {
        return load_genotype_prefix(std::move(in));
    }
    bool test_ambiguous(const std::string& a, const std::string& b)
    {
        return ambiguous(a, b);
    }
    size_t test_get_rs_column(const std::string& input)
    {
        if (m_reporter == nullptr)
        {
            Reporter reporter("log", 60, true);
            m_reporter = &reporter;
        }
        return get_rs_column(input);
    }
    std::vector<std::string>
    test_load_snp_list(std::unique_ptr<std::istream> input)
    {
        auto res = load_snp_list(std::move(input));
        std::vector<std::string> result;
        result.insert(result.end(), res.begin(), res.end());
        return result;
    }
    std::vector<std::string> test_load_ref(std::unique_ptr<std::istream> input,
                                           const std::string& delim,
                                           bool ignore_fid)
    {
        m_delim = delim;
        auto tmp = load_ref(std::move(input), ignore_fid);
        std::vector<std::string> result;
        result.insert(result.end(), tmp.begin(), tmp.end());
        return result;
    }
    bool test_parse_chr(const std::vector<std::string_view>& token,
                        const BaseFile& base_file,
                        std::vector<size_t>& filter_count, size_t& chr)
    {
        return parse_chr(token, base_file, filter_count, chr);
    }
    void test_init_chr(int num_auto = 22, bool no_x = false, bool no_y = false,
                       bool no_xy = false, bool no_mt = false)
    {
        init_chr(num_auto, no_x, no_y, no_xy, no_mt);
    }
    std::tuple<std::vector<size_t>, std::unordered_set<std::string>>
    test_transverse_base_file(
        const BaseFile& base_file, const QCFiltering& base_qc,
        const PThresholding& threshold_info,
        const std::vector<IITree<size_t, size_t>>& exclusion_regions,
        const std::streampos file_length, const bool gz_input,
        std::unique_ptr<std::istream> input)
    {
        return transverse_base_file(base_file, base_qc, threshold_info,
                                    exclusion_regions, file_length, gz_input,
                                    std::move(input));
    }
    void test_parse_allele(const std::vector<std::string_view>& token,
                           const BaseFile& base_file, size_t index,
                           std::string& allele)
    {
        parse_allele(token, base_file, index, allele);
    }
    bool test_parse_loc(const std::vector<std::string_view>& token,
                        const BaseFile& base_file, size_t& loc)
    {
        return parse_loc(token, base_file, loc);
    }
    bool test_parse_pvalue(const std::string_view& p_value_str,
                           const double max_threshold,
                           std::vector<size_t>& filter_count, double& pvalue)
    {
        return parse_pvalue(p_value_str, max_threshold, filter_count, pvalue);
    }
    bool test_parse_stat(const std::string_view& stat_str, const bool odd_ratio,
                         std::vector<size_t>& filter_count, double& stat)
    {
        return (parse_stat(stat_str, odd_ratio, filter_count, stat));
    }
    bool test_base_filter_by_value(const std::vector<std::string_view>& token,
                                   const BaseFile& base_file,
                                   const double& threshold,
                                   std::vector<size_t>& filter_count,
                                   size_t type, size_t index)
    {
        return (base_filter_by_value(token, base_file, threshold, filter_count,
                                     type, index));
    }
    std::vector<std::pair<size_t, size_t>> test_get_chrom_boundary()
    {
        return get_chrom_boundary();
    }
    void test_post_sample_read_init() { post_sample_read_init(); }
    bool test_parse_rs_id(const std::vector<std::string_view>& token,
                          const BaseFile& base_file,
                          std::unordered_set<std::string>& processed_rs,
                          std::unordered_set<std::string>& dup_index,
                          std::vector<size_t>& filter_count, std::string& rs_id)
    {
        return parse_rs_id(token, base_file, processed_rs, dup_index,
                           filter_count, rs_id);
    }
    unsigned long long
    test_cal_bar_category(const double& pvalue,
                          const std::vector<double>& barlevels, double& pthres)
    {
        return cal_bar_category(pvalue, barlevels, pthres);
    }
    unsigned long long
    test_calculate_category(const PThresholding& thresholding,
                            const double& pvalue, double& pthres)
    {
        return calculate_category(thresholding, pvalue, pthres);
    }
    void test_init_sample_vectors() { init_sample_vectors(); }
    void test_gen_sample(const size_t fid_idx, const size_t iid_idx,
                         const size_t sex_idx, const size_t dad_idx,
                         const size_t mum_idx, const size_t cur_idx,
                         const std::unordered_set<std::string>& founder_info,
                         const std::string& pheno,
                         std::vector<std::string>& token,
                         std::vector<Sample_ID>& sample_storage,
                         std::unordered_set<std::string>& sample_in_file,
                         std::vector<std::string>& duplicated_sample_id)
    {
        gen_sample(fid_idx, iid_idx, sex_idx, dad_idx, mum_idx, cur_idx,
                   founder_info, pheno, token, sample_storage, sample_in_file,
                   duplicated_sample_id);
    }
    bool test_check_ambig(const std::string& a1, const std::string& a2,
                          const std::string& ref, bool& flipping)
    {
        return check_ambig(a1, a2, ref, flipping);
    }
    bool test_check_chr(const std::string& chr_str, std::string& prev_chr,
                        size_t& chr_num, bool& chr_error, bool& sex_error)
    {
        return (check_chr(chr_str, prev_chr, chr_num, chr_error, sex_error));
    }
    bool test_check_rs(const std::string& snp_id, const std::string& chr_id,
                       std::string& rs_id,
                       std::unordered_set<std::string>& processed_snps,
                       std::unordered_set<std::string>& duplicated_snps,
                       Genotype* genotype)
    {
        return check_rs(snp_id, chr_id, rs_id, processed_snps, duplicated_snps,
                        genotype);
    }
    bool test_process_snp(
        const std::vector<IITree<size_t, size_t>>& exclusion_regions,
        const std::string& mismatch_snp_record_name,
        const std::string& mismatch_source, const std::string& snpid, SNP& snp,
        std::unordered_set<std::string>& processed_snps,
        std::unordered_set<std::string>& duplicated_snps,
        std::vector<bool>& retain_snp, Genotype* genotype)
    {
        return process_snp(exclusion_regions, mismatch_snp_record_name,
                           mismatch_source, snpid, snp, processed_snps,
                           duplicated_snps, retain_snp, genotype);
    }

    bool test_not_in_xregion(
        const std::vector<IITree<size_t, size_t>>& exclusion_regions,
        const SNP& base, const SNP& target)
    {
        return not_in_xregion(exclusion_regions, base, target);
    }
    void add_select_sample(const std::string& in)
    {
        m_sample_selection_list.insert(in);
    }
    void change_sample_selection(bool remove) { m_remove_sample = remove; }
    void add_select_snp(const std::string& in, bool exclude)
    {
        m_snp_selection_list.insert(in);
        m_exclude_snp = exclude;
    }
    void set_keep_nonfounder(bool keep_nonfounder)
    {
        m_keep_nonfounder = keep_nonfounder;
    }
    void load_snp(const std::string& rs)
    {
        m_existed_snps_index[rs] = m_existed_snps.size();
        m_existed_snps.emplace_back(SNP(rs, 1, 1, "A", "C", 0, 0, 1, 1));
    }
    void load_snp(SNP snp)
    {
        m_existed_snps_index[snp.rs()] = m_existed_snps.size();
        m_existed_snps.emplace_back(snp);
    }
    std::vector<SNP>& modify_existed_snps() { return m_existed_snps; }
    uint32_t num_auto() const { return m_autosome_ct; }
    std::vector<int32_t> xymt_codes() const { return m_xymt_codes; }
    std::vector<uintptr_t> haploid_mask() const { return m_haploid_mask; }
    std::string delim() const { return m_delim; }
    std::vector<std::string> genotype_file_names() const
    {
        return m_genotype_file_names;
    }
    void add_file_name(const std::string& in)
    {
        m_genotype_file_names.push_back(in);
    }
    std::vector<SNP> existed_snps() const { return m_existed_snps; }
    std::unordered_map<std::string, size_t> existed_snps_idx() const
    {
        return m_existed_snps_index;
    }
    std::string keep_file() const { return m_keep_file; }
    std::string remove_file() const { return m_remove_file; }
    std::string sample_file() const { return m_sample_file; }
    bool ignore_fid() const { return m_ignore_fid; }
    // helper
    void set_reporter(Reporter* reporter) { m_reporter = reporter; }
    void set_sample(uintptr_t n_sample) { m_unfiltered_sample_ct = n_sample; }
    void set_delim(const std::string& delim) { m_delim = delim; }
    void set_ignore_fid(const bool ignore_fid) { m_ignore_fid = ignore_fid; }
    size_t num_male() const { return m_num_male; }
    size_t num_female() const { return m_num_female; }
    size_t num_ambig_sex() const { return m_num_ambig_sex; }
    size_t num_ambig() const { return m_num_ambig; }
    size_t num_xrange() const { return m_num_xrange; }
    size_t num_nonfounder() const { return m_num_non_founder; }
    size_t base_missed() const { return m_base_missed; }
    uintptr_t num_founder() const { return m_founder_ct; }
    uintptr_t num_sample() const { return m_sample_ct; }
    size_t max_window() const { return m_max_window_size; }
    std::vector<uintptr_t> sample_for_ld() const { return m_sample_for_ld; }
    std::vector<uintptr_t> calculate_prs() const { return m_calculate_prs; }
    void test_update_index_tot(const uintptr_t founder_ctl2,
                               const uintptr_t founder_ctv2,
                               const uintptr_t founder_count,
                               std::vector<uintptr_t>& index_data,
                               std::vector<uintptr_t>& index_tots,
                               std::vector<uintptr_t>& founder_include2,
                               uintptr_t* index_genotype)
    {
        update_index_tot(founder_ctl2, founder_ctv2, founder_count, index_data,
                         index_tots, founder_include2, index_genotype);
    }

    double test_get_r2(const uintptr_t founder_ctl2,
                       const uintptr_t founder_ctv2, uintptr_t* window_data_ptr,
                       std::vector<uintptr_t>& index_data,
                       std::vector<uintptr_t>& index_tots)
    {
        return get_r2(founder_ctl2, founder_ctv2, window_data_ptr, index_data,
                      index_tots);
    }
    void set_sample_vector(const size_t n_sample)
    {
        for (size_t i = 0; i < n_sample; ++i)
        { SET_BIT(i, m_calculate_prs.data()); }
        m_sample_ct = n_sample;
    }

    void set_founder_vector(const size_t n_sample)
    {
        for (size_t i = 0; i < n_sample; ++i)
        { SET_BIT(i, m_sample_for_ld.data()); }
        m_founder_ct = n_sample;
    }
    void set_founder_vector(const std::vector<bool>& founder)
    {
        for (size_t i = 0; i < founder.size(); ++i)
        {
            if (founder[i])
            {
                SET_BIT(i, m_sample_for_ld.data());
                ++m_founder_ct;
            }
        }
    }
    void set_very_small_thresholds() { m_very_small_thresholds = true; }
};

#endif // MOCK_GENOTYPE_H

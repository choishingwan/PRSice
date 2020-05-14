// This file is part of PRSice-2, copyright (C) 2016-2019
// Shing Wan Choi, Paul F. O’Reilly
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#ifndef GENOTYPE_H
#define GENOTYPE_H

#include "IITree.h"
#include "commander.hpp"
#include "genotype_pool.hpp"
#include "misc.hpp"
#include "plink_common.hpp"
#include "reporter.hpp"
#include "snp.hpp"
#include "storage.hpp"
#include "thread_queue.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <deque>
#include <fstream>
#include <functional>
#include <memory>
#include <memoryread.hpp>
#include <mutex>
#include <random>
#include <set>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#ifdef __APPLE__
#include <sys/sysctl.h> // sysctl()
#endif


#define MULTIPLEX_LD 1920
#define MULTIPLEX_2LD (MULTIPLEX_LD * 2)
class Genotype
{
public:
    /*!
     * \brief default constructor of Genotype
     */
    Genotype() {}
    /*!
     * \brief Genotype constructor
     * \param thread = number of thread
     * \param ignore_fid = if FID should be ignored
     * \param keep_nonfounder = if we should keep non-founders
     * \param keep_ambig = if we should keep ambiguous SNPs
     * \param is_ref = if this is a reference genome
     */
    Genotype(const uint32_t thread, const bool ignore_fid,
             const bool keep_nonfounder, const bool keep_ambig,
             const bool is_ref = false)
        : m_thread(thread)
        , m_ignore_fid(ignore_fid)
        , m_is_ref(is_ref)
        , m_keep_nonfounder(keep_nonfounder)
        , m_keep_ambig(keep_ambig)
    {
    }
    virtual ~Genotype();

    // after load samples, samples that we don't want to include
    // will all have their FID, IID and phenotype masked by empty string
    /*!
     * \brief Function to load samples into the genotype class object. Will call
     *        gen_sample_vector. Will only store the sample info if this is the
     *        target.
     * \param keep_file is the file containing samples to be included
     * \param remove_file is the file containing samples to be excluded
     * \param verbose boolean indicate whether we should output the count
     * \param reporter is the logger
     */
    void load_samples(bool verbose = true);
    // We need the exclusion_region parameter because when we read in the base
    // we do allow users to provide a base without the CHR and LOC, which forbid
    // us to do the regional filtering. However, as exclusion and extractions
    // are based on the SNP ID, we can be confident that they were already done
    // in read_base
    void load_snps(const std::string& out,
                   const std::vector<IITree<size_t, size_t>>& exclusion_regions,
                   bool verbose, Genotype* target = nullptr);

    void calc_freqs_and_intermediate(const QCFiltering& filter_info,
                                     const std::string& prefix, bool verbose,
                                     Genotype* target = nullptr);
    /*!
     * \brief Return the number of SNPs, use for unit test
     * \return reuturn the number of SNPs included
     */
    size_t num_snps() const { return m_existed_snps.size(); }
    /*!
     * \brief Function to re-propagate the m_existed_snps_index after reading
     * the reference panel
     */
    void update_snp_index()
    {
        m_existed_snps_index.clear();
        for (size_t i_snp = 0; i_snp < m_existed_snps.size(); ++i_snp)
        { m_existed_snps_index[m_existed_snps[i_snp].rs()] = i_snp; }
    }
    /*!
     * \brief Return the number of sample we wish to perform PRS on
     * \return the number of sample
     */
    size_t num_sample() const { return m_sample_id.size(); }

    /*!
     * \brief Function to prepare clumping. Should sort all the SNPs by their
     * chromosome number and then by their p-value
     * \return true if there are SNPs to sort
     */
    bool sort_by_p()
    {
        if (m_existed_snps.size() == 0) return false;
        m_sort_by_p_index = SNP::sort_by_p_chr(m_existed_snps);
        return true;
    }

    static bool chr_prefix(std::string_view str)
    {
        if (str.length() <= 3) { return false; }
        return ((str[0] & 0xdf) == 'C' && (str[1] & 0xdf) == 'H'
                && (str[2] & 0xdf) == 'R');
    }
    static int32_t get_chrom_code(std::string_view str)
    {
        if (chr_prefix(str)) { str.remove_prefix(3); }
        if (str.length() == 0) return -1;
        if ((str[0] & 0xdf) == 'X')
        {
            if (str.length() == 2 && (str[1] & 0xdf) == 'Y')
                return CHROM_XY;
            else if (str.length() == 1)
                return CHROM_X;
        }
        else if (str.length() == 1 && (str[0] & 0xdf) == 'Y')
        {
            return CHROM_Y;
        }
        else if ((str[0] & 0xdf) == 'M')
        {
            if (str.length() == 1
                || (str.length() == 2 && (str[1] & 0xdf) == 'T'))
                return CHROM_MT;
        }
        else if (str[0] == '0')
        {
            if (str.length() == 2)
            {
                switch (str[1] & 0xdf)
                {
                case 'X': return CHROM_X;
                case 'Y': return CHROM_Y;
                case 'M': return CHROM_MT;
                }
            }
        }
        try
        {
            return misc::Convertor::convert<int32_t>(std::string(str));
        }
        catch (const std::runtime_error&)
        {
            return -1;
        }
    }


    std::string
    print_duplicated_snps(const std::unordered_set<std::string>& snp_name,
                          const std::string& out_prefix);
    bool base_filter_by_value(const std::vector<std::string_view>& token,
                              const BaseFile& base_file,
                              const double& threshold,
                              std::vector<size_t>& filter_count, size_t type,
                              size_t index)
    {
        if (!base_file.has_column[index]) return true;
        if (filter_count.size() != +FILTER_COUNT::MAX)
        { filter_count.resize(+FILTER_COUNT::MAX, 0); }
        double value = 1;
        try
        {
            value = misc::Convertor::convert<double>(
                std::string(token[base_file.column_index[index]]).c_str());
        }
        catch (...)
        {
            ++filter_count[type];
            return false;
        }
        if (value < threshold)
        {
            ++filter_count[type];
            return false;
        }
        return true;
    }
    bool parse_chr(const std::vector<std::string_view>& token,
                   const BaseFile& base_file, std::vector<size_t>& filter_count,
                   size_t& chr)
    {
        if (filter_count.size() != +FILTER_COUNT::MAX)
        { filter_count.resize(+FILTER_COUNT::MAX, 0); }
        chr = ~size_t(0);
        if (!base_file.has_column[+BASE_INDEX::CHR]) return true;
        int32_t chr_code =
            get_chrom_code(token[base_file.column_index[+BASE_INDEX::CHR]]);
        if (chr_code < 0)
        {
            ++filter_count[+FILTER_COUNT::CHR];
            return false;
        }
        if (static_cast<size_t>(chr_code) > m_autosome_ct
            || is_set(m_haploid_mask.data(), static_cast<uint32_t>(chr_code)))
        {
            ++filter_count[+FILTER_COUNT::HAPLOID];
            return false;
        }
        chr = static_cast<size_t>(chr_code);
        return true;
    }


    bool parse_loc(const std::vector<std::string_view>& token,
                   const BaseFile& base_file, size_t& loc)
    {
        loc = ~size_t(0);
        if (!base_file.has_column[+BASE_INDEX::BP]) return true;
        try
        {
            loc = misc::Convertor::convert<size_t>(
                std::string(token[base_file.column_index[+BASE_INDEX::BP]])
                    .c_str());
        }
        catch (...)
        {
            return false;
        }
        return true;
    }
    void init_sample_vectors()
    {
        if (m_unfiltered_sample_ct == 0)
        {
            throw std::runtime_error("Error: No sample found in genotype file");
        }
        const uintptr_t unfiltered_sample_ctl =
            BITCT_TO_WORDCT(m_unfiltered_sample_ct);
        m_sample_for_ld.resize(unfiltered_sample_ctl, 0);
        m_calculate_prs.resize(unfiltered_sample_ctl, 0);
        m_num_male = 0;
        m_num_female = 0;
        m_num_ambig_sex = 0;
        m_num_non_founder = 0;
        m_sample_ct = 0;
        m_founder_ct = 0;
        m_vector_initialized = true;
    }
    void post_sample_read_init()
    {
        if (m_unfiltered_sample_ct == 0)
        {
            throw std::runtime_error("Error: No sample found in genotype file");
        }
        const uintptr_t unfiltered_sample_ctl =
            BITCT_TO_WORDCT(m_unfiltered_sample_ct);
        const uintptr_t unfiltered_sample_ctv2 = 2 * unfiltered_sample_ctl;
        m_tmp_genotype.resize(unfiltered_sample_ctv2, 0);
        m_prs_info.resize(m_sample_ct, PRS());
        m_in_regression.resize(m_calculate_prs.size(), 0);
        m_sample_include2.resize(unfiltered_sample_ctv2, 0);
        m_founder_include2.resize(unfiltered_sample_ctv2, 0);
        // fill it with the required mask (copy from PLINK2)
        init_quaterarr_from_bitarr(m_calculate_prs.data(),
                                   m_unfiltered_sample_ct,
                                   m_sample_include2.data());
        init_quaterarr_from_bitarr(m_sample_for_ld.data(),
                                   m_unfiltered_sample_ct,
                                   m_founder_include2.data());
    }
    void clumping(const Clumping& clump_info, Genotype& reference,
                  size_t threads);
    std::vector<std::pair<size_t, size_t>> get_chrom_boundary();
    template <typename T>
    void
    threaded_clumping(const std::vector<std::pair<size_t, size_t>> snp_range,
                      const Clumping& clump_info, T& progress_observer,
                      std::vector<std::atomic<bool>>& remained_snps,
                      std::atomic<size_t>& num_core, Genotype& reference);
    void efficient_clumping(const Clumping& clump_info, Genotype& reference);

    /*!
     * \brief Before each run of PRSice, we need to reset the in regression
     * flag to false and propagate it later on to indicate if the sample is
     * used in the regression model
     */

    void reset_in_regression_flag()
    {
        std::fill(m_in_regression.begin(), m_in_regression.end(), 0);
    }
    void reset_std_flag()
    {
        std::fill(m_exclude_from_std.begin(), m_exclude_from_std.end(), 0);
    }
    void exclude_from_std(const size_t idx)
    {
        SET_BIT(idx, m_exclude_from_std.data());
    }
    /*!
     * \brief Function to prepare the object for PRSice. Will sort the
     * m_existed_snp vector according to their p-value.
     * \return True if there are SNPs to process
     */
    bool prepare_prsice(const PThresholding& m_p_info);
    /*!
     * \brief This function will return the sample ID
     * \param i is the index of the sample
     * \return the sample ID
     */
    std::string sample_id(const size_t i, const std::string& delim) const
    {
        if (i > m_sample_id.size())
            throw std::out_of_range("Sample name vector out of range");
        if (m_ignore_fid)
            return m_sample_id[i].IID;
        else
            return m_sample_id[i].FID + delim + m_sample_id[i].IID;
    }

    bool in_regression(size_t i) const
    {
        return m_sample_id.at(i).in_regression;
    }

    // as we don't check i is within range, the following 3 functions
    // have a potential of array out of bound
    /*!
     * \brief Check if the i th sample has been included in regression
     * \param i is the index of the sample
     * \return true if it is included
     */
    bool sample_in_regression(size_t i) const
    {
        return IS_SET(m_in_regression.data(), i);
    }
    /*!
     * \brief Inform PRSice that the i th sample has been included in the
     * regression
     *
     * \param i is the sample index
     */
    void set_in_regression(size_t i) { SET_BIT(i, m_in_regression.data()); }
    /*!
     * \brief Return the phenotype stored in the fam file of the i th sample
     * \param i is the index of the  sample
     * \return the phenotype of the sample
     */
    std::string pheno(size_t i) const { return m_sample_id[i].pheno; }
    /*!
     * \brief This function return if the i th sample has NA as phenotype
     * \param i is the sample ID
     * \return true if the phenotype is NA
     */
    bool pheno_is_na(size_t i) const { return m_sample_id[i].pheno == "NA"; }
    /*!
     * \brief Return the fid of the i th sample
     * \param i is the sample index
     * \return FID of the i th sample
     */
    std::string fid(size_t i) const { return m_sample_id.at(i).FID; }
    /*!
     * \brief Return the iid fo the i th sample
     * \param i is the sample index
     * \return IID of the i th sample
     */
    std::string iid(size_t i) const { return m_sample_id.at(i).IID; }
    /*!
     * \brief This function will calculate the required PRS for the i th
     * sample \param score_type is the type of score user want to calculate
     * \param i is the sample index
     * \return the PRS
     */
    inline double calculate_score(const std::vector<PRS>& prs_list,
                                  size_t i) const
    {
        if (i >= prs_list.size())
            throw std::out_of_range("Sample name vector out of range");
        const size_t num_snp = prs_list[i].num_snp;
        double prs = prs_list[i].prs;
        double avg = prs;
        if (num_snp == 0) { avg = 0.0; }
        else
        {
            avg = prs / static_cast<double>(num_snp);
        }
        switch (m_prs_calculation.scoring_method)
        {
        case SCORING::SUM: return prs_list[i].prs;
        case SCORING::STANDARDIZE:
        case SCORING::CONTROL_STD: return (avg - m_mean_score) / m_score_sd;
        default:
            // default is avg
            return avg;
        }
    }
    inline double calculate_score(size_t i) const
    {
        return calculate_score(m_prs_info, i);
    }
    /*!
     * \brief Function for calculating the PRS from the null set
     * \param set_size is the size of the set
     * \param prev_size is the amount of SNPs we have already processed
     * \param background_list is the vector containing the permuted
     * background index \param first_run is a boolean representing if we
     * need to reset the PRS to
     * 0
     * \param require_standardize is a boolean representing if we need to
     * calculate the mean and SD
     */
    void get_null_score(std::vector<PRS>& prs_list, const size_t& set_size,
                        const size_t& prev_size,
                        std::vector<size_t>& background_list,
                        const bool first_run);
    void get_null_score(const size_t& set_size, const size_t& prev_size,
                        std::vector<size_t>& background_list,
                        const bool first_run)
    {
        get_null_score(m_prs_info, set_size, prev_size, background_list,
                       first_run);
    }
    /*!
     * \brief return the largest chromosome allowed
     * \return  the largest chromosome
     */
    uint32_t max_chr() const { return m_max_code; }
    /*!
     * \brief Indicate we will be using a reference file. Use for bgen
     * intermediate output generation
     */
    void expect_reference() { m_expect_reference = true; }
    std::tuple<std::vector<size_t>, std::unordered_set<std::string>>
    read_base(const BaseFile& base_file, const QCFiltering& base_qc,
              const PThresholding& threshold_info,
              const std::vector<IITree<size_t, size_t>>& exclusion_regions);
    std::tuple<std::vector<size_t>, std::unordered_set<std::string>>
    transverse_base_file(
        const BaseFile& base_file, const QCFiltering& base_qc,
        const PThresholding& threshold_info,
        const std::vector<IITree<size_t, size_t>>& exclusion_regions,
        const std::streampos file_length, const bool gz_input,
        std::unique_ptr<std::istream> input);
    void print_base_stat(const std::vector<size_t>& filter_count,
                         const std::unordered_set<std::string>& dup_index,
                         const std::string& out, const double info_score);
    void build_clump_windows(const unsigned long long& clump_distance);
    intptr_t cal_avail_memory(const uintptr_t founder_ctv2);
    void
    build_membership_matrix(std::vector<std::vector<size_t>>& region_membership,
                            const size_t num_sets, const std::string& out,
                            const std::vector<std::string>& region_name,
                            const bool print_snps);
    size_t num_threshold(size_t idx) const
    {
        return m_set_thresholds[idx].size();
    }
    std::vector<std::set<double>> get_set_thresholds() const
    {
        return m_set_thresholds;
    }
    bool get_score(std::vector<size_t>::const_iterator& start_index,
                   const std::vector<size_t>::const_iterator& end_index,
                   double& cur_threshold, uint32_t& num_snp_included,
                   const bool first_run);
    static bool within_region(const std::vector<IITree<size_t, size_t>>& cr,
                              const size_t chr, const size_t loc)
    {
        std::vector<size_t> output;
        if (chr >= cr.size()) return false;
        size_t start = (loc == 0) ? 1 : loc;
        cr[chr].overlap(start, loc + 1, output);
        return !output.empty();
    }

    bool perform_freqs_and_inter(const QCFiltering& filter_info,
                                 const std::string& prefix, Genotype* target);
    static void
    construct_flag(const std::vector<IITree<size_t, size_t>>& gene_sets,
                   const size_t required_size,
                   const bool genome_wide_background, SNP* snp)
    {
        auto& flag = snp->get_flag();
        if (flag.size() != required_size) { flag.resize(required_size); }
        std::fill(flag.begin(), flag.end(), 0);
        SET_BIT(0, flag.data());
        if (genome_wide_background) { SET_BIT(1, flag.data()); }
        if (gene_sets.empty()) return;
        if (snp->chr() >= gene_sets.size()) return;
        std::vector<size_t> out;
        if (!gene_sets[snp->chr()].has_overlap(snp->loc(), out)) return;
        for (auto&& idx : out) { SET_BIT(idx, flag.data()); }
    }

    void add_flags(const std::vector<IITree<size_t, size_t>>& cr,
                   const size_t num_sets, const bool genome_wide_background);

    // Refactoring
    Genotype& keep_nonfounder(bool keep)
    {
        m_keep_nonfounder = keep;
        return *this;
    }
    Genotype& keep_ambig(bool keep)
    {
        m_keep_ambig = keep;
        return *this;
    }
    Genotype& reference()
    {
        m_is_ref = true;
        return *this;
    }
    Genotype& intermediate(bool use)
    {
        m_intermediate = use;
        return *this;
    }
    Genotype& set_prs_instruction(const CalculatePRS& prs)
    {
        m_has_prs_instruction = true;
        m_prs_calculation = prs;
        return *this;
    }
    void snp_extraction(const std::string& extract_snps,
                        const std::string& exclude_snps);
    void clump_progress_observer(Thread_Queue<size_t>& progress_observer,
                                 size_t total_snp, size_t num_thread,
                                 bool verbose);
    class dummy_reporter
    {
        bool m_completed = false;
        bool m_verbose = true;
        size_t m_total_snp = 0;
        size_t m_processed = 0;
        double m_prev_progress = 0;

    public:
        dummy_reporter(size_t total, bool verbose = true)
            : m_verbose(verbose), m_total_snp(total)
        {
        }
        void emplace(size_t&& item)
        {
            m_processed += item;
            auto progress = (static_cast<double>(m_processed)
                             / static_cast<double>(m_total_snp))
                            * 100;
            if (m_verbose && progress - m_prev_progress > 0.01)
            {
                fprintf(stderr, "\rClumping Progress: %03.2f%%", progress);
                m_prev_progress = progress;
            }
        }
        void completed() { m_completed = true; }
    };

    Genotype& set_weight()
    {
        assert(m_has_prs_instruction);
        switch (m_prs_calculation.genetic_model)
        {
        case MODEL::HETEROZYGOUS:
            m_homcom_weight = 0;
            m_het_weight = 1;
            m_homrar_weight = 0;
            break;
        case MODEL::DOMINANT:
            m_homcom_weight = 0;
            m_het_weight = 1;
            m_homrar_weight = 1;
            break;
        case MODEL::RECESSIVE:
            m_homcom_weight = 0;
            m_het_weight = 0;
            m_homrar_weight = 1;
            break;
        default:
            m_homcom_weight = 0;
            m_het_weight = 1;
            m_homrar_weight = 2;
            break;
        }
        return *this;
    }

    void set_thresholds(const QCFiltering& qc)
    {
        m_hard_threshold = qc.hard_threshold;
        m_dose_threshold = qc.dose_threshold;
    }

    void load_genotype_to_memory();
    bool genotyped_stored() const { return m_genotype_stored; }
    const std::unordered_map<std::string, size_t>& included_snps_idx() const
    {
        return m_existed_snps_index;
    }
    const std::vector<SNP>& included_snps() const { return m_existed_snps; }

protected:
    // friend with all child class so that they can also access the
    // protected elements
    friend class BinaryPlink;
    friend class BinaryGen;
    // vector storing all the genotype files
    // std::vector<Sample> m_sample_names;
    FileRead m_genotype_file;
    std::vector<SNP> m_existed_snps;
    std::unordered_map<std::string, size_t> m_existed_snps_index;
    std::unordered_set<std::string> m_sample_selection_list;
    std::unordered_set<std::string> m_snp_selection_list;
    std::vector<std::set<double>> m_set_thresholds;
    std::vector<Sample_ID> m_sample_id;
    std::vector<PRS> m_prs_info;
    std::vector<std::string> m_genotype_file_names;
    std::vector<uintptr_t> m_tmp_genotype;
    // std::vector<uintptr_t> m_chrom_mask;
    std::vector<uintptr_t> m_sample_for_ld;
    std::vector<uintptr_t> m_founder_include2;
    std::vector<uintptr_t> m_calculate_prs;
    std::vector<uintptr_t> m_sample_include2;
    std::vector<uintptr_t> m_exclude_from_std;
    std::vector<uintptr_t> m_in_regression;
    std::vector<uintptr_t> m_haploid_mask;
    std::vector<size_t> m_sort_by_p_index;
    // std::vector<uintptr_t> m_sex_male;
    std::vector<int32_t> m_xymt_codes;
    std::ofstream m_mismatch_snp_record;
    // std::vector<int32_t> m_chrom_start;
    // sample file name. Fam for plink
    std::string m_sample_file;
    std::string m_delim;
    std::string m_keep_file;
    std::string m_remove_file;
    double m_mean_score = 0.0;
    double m_score_sd = 0.0;
    double m_hard_threshold = 0.0;
    double m_dose_threshold = 0.0;
    double m_homcom_weight = 0;
    double m_het_weight = 1;
    double m_homrar_weight = 2;
    size_t m_num_thresholds = 0;
    size_t m_thread = 1; // number of final samples
    size_t m_max_window_size = 0;
    size_t m_num_ambig = 0;
    size_t m_num_maf_filter = 0;
    size_t m_num_geno_filter = 0;
    size_t m_num_info_filter = 0;
    size_t m_num_miss_filter = 0;
    size_t m_num_xrange = 0;
    size_t m_base_missed = 0;
    size_t m_num_male = 0;
    size_t m_num_female = 0;
    size_t m_num_ambig_sex = 0;
    size_t m_num_non_founder = 0;
    uintptr_t m_unfiltered_sample_ct = 0; // number of unfiltered samples
    uintptr_t m_unfiltered_marker_ct = 0;
    uintptr_t m_sample_ct = 0;
    uintptr_t m_founder_ct = 0;
    uintptr_t m_marker_ct = 0;
    uint32_t m_max_category = 0;
    uint32_t m_autosome_ct = 0;
    uint32_t m_max_code = 0;
    std::random_device::result_type m_seed = 0;
    uint32_t m_num_ref_target_mismatch = 0;
    bool m_genotype_stored = false;
    bool m_use_proxy = false;
    bool m_has_prs_instruction = false;
    bool m_ignore_fid = false;
    bool m_intermediate = false;
    bool m_is_ref = false;
    bool m_keep_nonfounder = false;
    bool m_keep_ambig = false;
    bool m_remove_sample = true;
    bool m_exclude_snp = true;
    bool m_hard_coded = false;
    bool m_expect_reference = false;
    bool m_memory_initialized = false;
    bool m_very_small_thresholds = false;
    bool m_vector_initialized = false;
    Reporter* m_reporter = nullptr;
    CalculatePRS m_prs_calculation;

    std::string initialize(const GenoFile& geno, const Phenotype& pheno,
                           const std::string& delim, const std::string& type,
                           Reporter* reporter)
    {
        m_sample_file = "";
        m_ignore_fid = pheno.ignore_fid;
        m_keep_file = geno.keep;
        m_remove_file = geno.remove;
        m_delim = delim;
        m_reporter = reporter;
        init_chr(geno.num_autosome);
        const bool use_list = !(geno.file_list.empty());
        const std::string file_name =
            use_list ? geno.file_list : geno.file_name;
        std::vector<std::string> token = misc::split(file_name, ",");
        if (token.empty())
        {
            throw std::runtime_error("Error: Invalid user input: " + file_name);
        }
        const bool external_sample = (token.size() == 2);
        if (token.size() > 2)
        {
            throw std::runtime_error("Error: Undefine user input: "
                                     + file_name);
        }
        // check for empty input
        for (auto&& t : token)
        {
            if (t.empty())
            {
                throw std::runtime_error("Error: Invalid user input: "
                                         + file_name);
            }
        }
        std::string message = "Initializing Genotype ";
        if (use_list)
        {
            auto input = misc::load_stream(token[0]);
            // we move input to function, can no longer use it
            m_genotype_file_names = load_genotype_prefix(std::move(input));
            message.append("info from file: " + token[0] + " (" + type + ")\n");
        }
        else
        {
            m_genotype_file_names = set_genotype_files(token[0]);
            message.append("file: " + token[0] + " (" + type + ")\n");
        }
        if (external_sample)
        {
            m_sample_file = token[1];
            message.append("With external fam file: " + m_sample_file + "\n");
        }
        return message;
    }
    bool has_parent(const std::unordered_set<std::string>& founder_info,
                    const std::vector<std::string>& token,
                    const std::string& fid, const size_t idx);
    void gen_sample(const size_t fid_idx, const size_t iid_idx,
                    const size_t sex_idx, const size_t dad_idx,
                    const size_t mum_idx, const size_t cur_idx,
                    const std::unordered_set<std::string>& founder_info,
                    const std::string& pheno, std::vector<std::string>& token,
                    std::vector<Sample_ID>& sample_storage,
                    std::unordered_set<std::string>& sample_in_file,
                    std::vector<std::string>& duplicated_sample_id);
    void recalculate_categories(const PThresholding& p_info);
    void print_mismatch(const std::string& out, const std::string& type,
                        const SNP& target, const SNP& new_snp);
    /*
        bool parse_chr_id(const std::vector<std::string_view>& token,
                          const BaseFile& base_file,
                          std::unordered_set<std::string>& processed_idx,
                          std::unordered_set<std::string>& dup_rs,
                          std::vector<size_t>& filter_count, std::string&
       chr_id)
        {
            if (!base_file.has_column[+BASE_INDEX::CHR]
                && !base_file.has_column[+BASE_INDEX::BP])
            {
                throw std::runtime_error(
                    "Error: Need chr and bp column to construct chr id");
            }
            chr_id = SNP::chr_id(token[+BASE_INDEX::CHR],
       token[+BASE_INDEX::BP]); return snp_dup_selection_check(chr_id,
       processed_idx, dup_rs, filter_count);
        }
        */
    bool snp_dup_selection_check(const std::string& id,
                                 std::unordered_set<std::string>& processed_idx,
                                 std::unordered_set<std::string>& dup_rs,
                                 std::vector<size_t>& filter_count)
    {
        if (filter_count.size() != +FILTER_COUNT::MAX)
        { filter_count.resize(+FILTER_COUNT::MAX, 0); }
        if (processed_idx.find(id) != processed_idx.end())
        {
            ++filter_count[+FILTER_COUNT::DUP_SNP];
            dup_rs.insert(id);
            return false;
        }
        auto&& selection = m_snp_selection_list.find(id);
        if ((!m_exclude_snp && selection == m_snp_selection_list.end())
            || (m_exclude_snp && selection != m_snp_selection_list.end()))
        {
            ++filter_count[+FILTER_COUNT::SELECT];
            return false;
        }
        processed_idx.insert(id);
        return true;
    }
    bool parse_rs_id(const std::vector<std::string_view>& token,
                     const BaseFile& base_file,
                     std::unordered_set<std::string>& processed_idx,
                     std::unordered_set<std::string>& dup_rs,
                     std::vector<size_t>& filter_count, std::string& rs_id)
    {
        if (!base_file.has_column[+BASE_INDEX::RS])
        { throw std::runtime_error("Error: RS ID column not provided!"); }
        rs_id = token[base_file.column_index[+BASE_INDEX::RS]];
        return (snp_dup_selection_check(rs_id, processed_idx, dup_rs,
                                        filter_count));
    }

    void parse_allele(const std::vector<std::string_view>& token,
                      const BaseFile& base_file, size_t index,
                      std::string& allele)
    {
        allele = (base_file.has_column[index])
                     ? token[base_file.column_index[index]]
                     : "";
        misc::to_upper(allele);
    }

    bool parse_pvalue(const std::string_view& p_value_str,
                      const double max_threshold,
                      std::vector<size_t>& filter_count, double& pvalue)
    {
        if (filter_count.size() != +FILTER_COUNT::MAX)
        { filter_count.resize(+FILTER_COUNT::MAX, 0); }
        try
        {
            pvalue = misc::Convertor::convert<double>(std::string(p_value_str));
        }
        catch (...)
        {
            ++filter_count[+FILTER_COUNT::NOT_CONVERT];
            return false;
        }
        if (pvalue < 0.0 || pvalue > 1.0)
        {
            throw std::runtime_error("Error: Invalid p-value: "
                                     + std::string(p_value_str));
        }
        if (pvalue > max_threshold)
        {
            ++filter_count[+FILTER_COUNT::P_EXCLUDED];
            return false;
        }
        return true;
    }
    bool parse_stat(const std::string_view& stat_str, const bool odd_ratio,
                    std::vector<size_t>& filter_count, double& stat)
    {
        if (filter_count.size() != +FILTER_COUNT::MAX)
        { filter_count.resize(+FILTER_COUNT::MAX, 0); }
        try
        {
            stat = misc::Convertor::convert<double>(std::string(stat_str));
            if (odd_ratio && misc::logically_equal(stat, 0.0))
            {
                ++filter_count[+FILTER_COUNT::NOT_CONVERT];
                return false;
            }
            else if (odd_ratio && stat < 0.0)
            {
                ++filter_count[+FILTER_COUNT::NEGATIVE];
                return false;
            }
            else if (odd_ratio)
                stat = log(stat);
            return true;
        }
        catch (...)
        {
            ++filter_count[+FILTER_COUNT::NOT_CONVERT];
            return false;
        }
    }
    /*!
     * \brief Calculate the threshold bin based on the p-value and bound
     * info \param pvalue the input p-value \param bound_start is the start
     * of p-value threshold \param bound_inter is the step size of p-value
     * threshold \param bound_end is the end of p-value threshold \param
     * pthres return the name of p-value threshold this SNP belongs to
     * \param no_full indicate if we want the p=1 threshold
     * \return the category where this SNP belongs to
     */
    unsigned long long calculate_category(const PThresholding& thresholding,
                                          const double& pvalue, double& pthres)
    {
        // NOTE: Threshold is x < p <= end and minimum category is 0
        unsigned long long category = 0;
        double division;
        if (std::fpclassify(thresholding.inter) != FP_NORMAL)
        { throw std::runtime_error("Error: interval is too small."); }
        else if (pvalue > thresholding.upper && !thresholding.no_full)
        {
            division = (thresholding.upper + 0.1 - thresholding.lower)
                       / thresholding.inter;
            if (division > std::numeric_limits<unsigned long long>::max())
            {
                throw std::runtime_error(
                    "Error: Number of threshold required are too large");
            }
            category = static_cast<unsigned long long>(std::ceil(division));
            pthres = 1.0;
        }
        else
        {
            if (pvalue > thresholding.lower
                || misc::logically_equal(pvalue, thresholding.lower))
            {
                division = (pvalue - thresholding.lower) / thresholding.inter;
                if (division > std::numeric_limits<unsigned long long>::max())
                {
                    throw std::runtime_error("Error: Number of threshold "
                                             "required are too large");
                }
                category = static_cast<unsigned long long>(std::ceil(division));
            }
            pthres = category * thresholding.inter + thresholding.lower;
        }
        return category;
    }

    unsigned long long cal_bar_category(const double& pvalue,
                                        const std::vector<double>& barlevels,
                                        double& pthres)
    {
        for (unsigned long long i = 0; i < barlevels.size(); ++i)
        {
            if (pvalue < barlevels[i]
                || misc::logically_equal(pvalue, barlevels[i]))
            {
                pthres = barlevels[i];
                return i;
            }
        }
        pthres = 1.0;
        return static_cast<unsigned long long>(barlevels.size());
    }

    /*!
     * \brief Replace # in the name of the genotype file and generate list
     * of file for subsequent analysis \param prefix contains the name of
     * the genotype file \return a vector of string containing names of the
     * genotype files
     */
    std::vector<std::string> set_genotype_files(const std::string& prefix);
    /*!
     * \brief Read in the genotype list file and add the genotype file names
     * to the vector \param file_name is the name of the list file \return a
     * vector of string containing names of all genotype files
     */
    std::vector<std::string>
    load_genotype_prefix(std::unique_ptr<std::istream> in);
    /*!
     * \brief Initialize vector related to chromosome information
     * e.g.haplotype. Currently not really useful except for setting the
     * max_code which is later used to transform chromosome strings to
     * chromosome code \param num_auto is the number of autosome, we fix it
     * to 22 for human \param no_x indicate if chrX is missing for this
     * organism \param no_y indicate if chrY is missing for this organism
     * \param no_xy indicate if chrXY is missing for this organism
     * \param no_mt indicate if chrMT is missing for this organism
     */
    void init_chr(int num_auto = 22, bool no_x = false, bool no_y = false,
                  bool no_xy = false, bool no_mt = false);


    /*!
     * \brief Function to read in the sample. Any subclass must implement
     * this function. They \b must initialize the \b m_sample_info \b
     * m_founder_info \b m_founder_ct \b m_sample_ct \b m_prs_info \b
     * m_in_regression and \b m_tmp_genotype (optional) \return vector
     * containing the sample information
     */
    virtual std::vector<Sample_ID> gen_sample_vector()
    {
        return std::vector<Sample_ID>(0);
    }

    virtual void gen_snp_vector(
        const std::vector<IITree<size_t, size_t>>& /*exclusion_regions*/,
        const std::string& /*out_prefix*/, Genotype* /*target*/)
    {
    }
    virtual bool calc_freq_gen_inter(const QCFiltering& /*QC info*/,
                                     const std::string& /*prefix*/,
                                     Genotype* /*target=nullptr*/)
    {
        return false;
    }

    void update_index_tot(const uintptr_t founder_ctl2,
                          const uintptr_t founder_ctv2,
                          const uintptr_t founder_count,
                          std::vector<uintptr_t>& index_data,
                          std::vector<uintptr_t>& index_tots,
                          std::vector<uintptr_t>& founder_include2,
                          uintptr_t* index_genotype)
    {
        // reset the index_data information
        if (index_genotype == nullptr)
        { throw std::runtime_error("Error: Genotype is null!"); }
        assert(index_genotype != nullptr);
        std::fill(index_data.begin(), index_data.end(), 0);
        vec_datamask(founder_count, 0, index_genotype, founder_include2.data(),
                     index_data.data());
        index_tots[0] = popcount2_longs(index_data.data(), founder_ctl2);
        vec_datamask(founder_count, 2, index_genotype, founder_include2.data(),
                     &(index_data[founder_ctv2]));
        index_tots[1] =
            popcount2_longs(&(index_data[founder_ctv2]), founder_ctl2);
        vec_datamask(founder_count, 3, index_genotype, founder_include2.data(),
                     &(index_data[2 * founder_ctv2]));
        index_tots[2] =
            popcount2_longs(&(index_data[2 * founder_ctv2]), founder_ctl2);
    }

    double get_r2(const uintptr_t founder_ctl2, const uintptr_t founder_ctv2,
                  uintptr_t* window_data_ptr,
                  std::vector<uintptr_t>& index_data,
                  std::vector<uintptr_t>& index_tots)
    {
        assert(window_data_ptr != nullptr);
        // is_x is used in PLINK to indicate if the genotype is from the X
        // chromsome, as PRSice ignore any sex chromosome, we can set it as
        // a constant false
        const bool is_x = false;
        uint32_t counts[18];
        double freq11;
        double freq11_expected;
        double freq1x;
        double freq2x;
        double freqx1;
        double freqx2;
        double dxx;
        // calculate the counts
        // these counts are then used for calculation of R2. However, I
        // don't fully understand the algorithm here (copy from PLINK2)
        genovec_3freq(window_data_ptr, index_data.data(), founder_ctl2,
                      &(counts[0]), &(counts[1]), &(counts[2]));
        counts[0] = index_tots[0] - counts[0] - counts[1] - counts[2];
        genovec_3freq(window_data_ptr, &(index_data[founder_ctv2]),
                      founder_ctl2, &(counts[3]), &(counts[4]), &(counts[5]));
        counts[3] = index_tots[1] - counts[3] - counts[4] - counts[5];
        genovec_3freq(window_data_ptr, &(index_data[2 * founder_ctv2]),
                      founder_ctl2, &(counts[6]), &(counts[7]), &(counts[8]));
        counts[6] = index_tots[2] - counts[6] - counts[7] - counts[8];
        if (!em_phase_hethet_nobase(counts, is_x, is_x, &freq1x, &freq2x,
                                    &freqx1, &freqx2, &freq11))
        {
            // if the calculation is sucessful, we can then calculate the R2
            freq11_expected = freqx1 * freq1x;
            dxx = freq11 - freq11_expected;
            // message from PLINK:
            // if r^2 threshold is 0, let everything else through but
            // exclude the apparent zeroes.  Zeroes *are* included if
            // r2_thresh is negative,
            // though (only nans are rejected then).
            if (fabs(dxx) < SMALL_EPSILON
                || fabs(freq11_expected * freq2x * freqx2) < SMALL_EPSILON)
            { return 0.0; }
            else
            {
                return dxx * dxx / (freq11_expected * freq2x * freqx2);
            }
        }
        return -1;
    }

    void update_sample_prs(PRS& sample_prs, const uint32_t geno,
                           const std::vector<double>& scores,
                           const std::vector<size_t>& counts)
    {
        sample_prs.num_snp += counts[geno];
        sample_prs.prs += scores[geno];
    }
    void initialize_sample_prs(PRS& sample_prs, const uint32_t geno,
                               const std::vector<double>& scores,
                               const std::vector<size_t>& counts)
    {
        sample_prs.num_snp = counts[geno];
        sample_prs.prs = scores[geno];
    }
    template <class T>
    void process_sample_prs(std::vector<uintptr_t>& genotype,
                            std::vector<PRS>& prs_list,
                            const std::vector<double>& scores,
                            const std::vector<size_t>& counts, T load_prs)
    {
        uintptr_t* lbptr = genotype.data();
        uintptr_t ulii;
        uint32_t uii;
        uint32_t ujj;
        uint32_t ukk;
        uii = 0;
        ulii = 0;
        do
        {
            // ulii contain the numeric representation of the current
            // genotype
            ulii = ~(*lbptr++);
            if (uii + BITCT2 > m_unfiltered_sample_ct)
            {
                // this is PLINK, not sure exactly what this is about
                ulii &= (ONELU << ((m_unfiltered_sample_ct & (BITCT2 - 1)) * 2))
                        - ONELU;
            }
            // ujj sample index of the current genotype block
            ujj = 0;
            while (ujj < BITCT)
            {
                // go through the whole genotype block
                // ukk is the current genotype
                ukk = (ulii >> ujj) & 3;
                // and the sample index can be calculated as uii+(ujj/2)
                if (uii + (ujj / 2) >= m_sample_ct) { break; }
                auto& sample_prs = prs_list[uii + (ujj / 2)];
                // now we will get all genotypes (0, 1, 2, 3)
                (this->*load_prs)(sample_prs, ukk, scores, counts);
                // ulii &= ~((3 * ONELU) << ujj);
                // as each sample is represented by two byte, we will add 2
                // to the index
                ujj += 2;
            }
            // uii is the number of samples we have finished so far
            uii += BITCT2;
        } while (uii < m_sample_ct);
    }

    void read_prs(std::vector<uintptr_t>& genotype, std::vector<PRS>& prs_list,
                  const size_t ploidy, const double stat,
                  const double adj_score, const double miss_score,
                  const size_t miss_count, const double homcom_weight,
                  const double het_weight, const double homrar_weight,
                  const bool not_first)
    {
        std::vector<double> scores = {homcom_weight * stat - adj_score,
                                      het_weight * stat - adj_score, miss_score,
                                      homrar_weight * stat - adj_score};
        std::vector<size_t> counts = {ploidy, ploidy, miss_count, ploidy};
        if (not_first)
        {
            process_sample_prs(genotype, prs_list, scores, counts,
                               &Genotype::update_sample_prs);
        }
        else
        {
            process_sample_prs(genotype, prs_list, scores, counts,
                               &Genotype::initialize_sample_prs);
        }
    }

    /*!
     * \brief Function to read in the genotype in PLINK binary format. Any
     * subclass must implement this function to assist the processing of
     * their specific file type. The first argument is the genotype vector
     * use to store the PLINK binary whereas the second parameter is the
     * streampos, allowing us to seekg to the specific location of the file
     * as indicate in the thrid parameter
     */
    virtual inline void read_genotype(uintptr_t* /*genotype*/,
                                      const SNP& /*snp*/)
    {
    }
    virtual inline void read_genotype(FileRead& /*genotype_file*/,
                                      uintptr_t* /*tmp_store*/,
                                      uintptr_t* /*genotype*/,
                                      const SNP& /*snp*/)
    {
    }
    virtual void
    read_score(std::vector<PRS>& /*prs_list*/,
               const std::vector<size_t>::const_iterator& /*start*/,
               const std::vector<size_t>::const_iterator& /*end*/,
               bool /*reset_zero*/, bool ultra = false)
    {
    }
    void read_score(const std::vector<size_t>::const_iterator& start,
                    const std::vector<size_t>::const_iterator& end,
                    bool reset_zero, bool ultra = false)
    {
        read_score(m_prs_info, start, end, reset_zero, ultra);
    }
    void standardize_prs();
    // for loading the sample inclusion / exclusion set
    /*!
     * \brief Function to load in the sample extraction exclusion list
     * \param input the file name
     * \param ignore_fid whether we should ignore the FID (use 2 column or
     * 1) \return an unordered_set use for checking if the sample is in the
     * file
     */
    std::unordered_set<std::string>
    load_ref(std::unique_ptr<std::istream> input, bool ignore_fid);
    bool
    not_in_xregion(const std::vector<IITree<size_t, size_t>>& exclusion_regions,
                   const SNP& base, const SNP& target);
    bool check_rs(const std::string& snpid, const std::string& chrid,
                  std::string& rsid,
                  std::unordered_set<std::string>& processed_snps,
                  std::unordered_set<std::string>& duplicated_snps,
                  Genotype* genotype);
    bool check_ambig(const std::string& a1, const std::string& a2,
                     const std::string& ref, bool& flipping);

    bool check_chr(const std::string& chr_str, std::string& prev_chr,
                   size_t& chr_num, bool& chr_error, bool& sex_error);
    bool
    process_snp(const std::vector<IITree<size_t, size_t>>& exclusion_regions,
                const std::string& mismatch_snp_record_name,
                const std::string& mismatch_source, const std::string& snpid,
                SNP& snp, std::unordered_set<std::string>& processed_snps,
                std::unordered_set<std::string>& duplicated_snps,
                std::vector<bool>& retain_snp, Genotype* genotype);
    void shrink_snp_vector(const std::vector<bool>& retain)
    {
        m_existed_snps.erase(
            std::remove_if(m_existed_snps.begin(), m_existed_snps.end(),
                           [&retain, this](const SNP& s) {
                               return !retain[(&s - &*begin(m_existed_snps))];
                           }),
            m_existed_snps.end());
        m_existed_snps.shrink_to_fit();
    }
    void shrink_snp_vector(const std::vector<std::atomic<bool>>& retain)
    {
        m_existed_snps.erase(
            std::remove_if(m_existed_snps.begin(), m_existed_snps.end(),
                           [&retain, this](const SNP& s) {
                               return !retain[(&s - &*begin(m_existed_snps))];
                           }),
            m_existed_snps.end());
        m_existed_snps.shrink_to_fit();
    }
    bool filter_snp(const uint32_t ref_ct, const uint32_t het_ct,
                    const uint32_t alt_ct, const uint32_t ref_founder_ct,
                    const uint32_t het_founder_ct,
                    const uint32_t alt_founder_ct, const double geno,
                    const double maf, uint32_t& missing_founder_ct)
    {
        uint32_t total_alleles = ref_ct + het_ct + alt_ct;
        double cur_geno =
            1.0 - total_alleles / (static_cast<double>(m_sample_ct));
        if (total_alleles == m_sample_ct)
        {
            ++m_num_miss_filter;
            return true;
        }
        if (geno < cur_geno)
        {
            ++m_num_geno_filter;
            return true;
        }

        uint32_t total_founder_alleles =
            ref_founder_ct + het_founder_ct + alt_founder_ct;
        missing_founder_ct =
            static_cast<uint32_t>(m_founder_ct) - total_founder_alleles;
        if (missing_founder_ct == m_founder_ct)
        {
            ++m_num_miss_filter;
            return true;
        }
        double cur_maf =
            (static_cast<double>(2.0 * alt_founder_ct + het_founder_ct))
            / (static_cast<double>(2.0 * total_founder_alleles));
        cur_maf = (cur_maf > 0.5) ? 1 - cur_maf : cur_maf;
        if (alt_founder_ct == total_founder_alleles
            || ref_founder_ct == total_founder_alleles || cur_maf < maf)
        {
            ++m_num_maf_filter;
            return true;
        }
        return false;
    }
    /*!
     * \brief Function to load in SNP extraction exclusion list
     * \param input the file name of the SNP list
     * \param reporter the logger
     * \returnan unordered_set use for checking if the SNP is in the file
     */
    std::unordered_set<std::string>
    load_snp_list(std::unique_ptr<std::istream> input);
    size_t get_rs_column(const std::string& input);
    /** Misc information **/
    /*!
     * \brief Function to check if the two alleles are ambiguous
     * \param ref_allele the reference allele
     * \param alt_allele the alternative allele
     * \return true if it is ambiguous
     */
    inline bool ambiguous(const std::string& ref_allele,
                          const std::string& alt_allele) const
    {
        // the allele should all be in upper case but whatever
        // true if equal
        if (ref_allele.empty() || alt_allele.empty()) return false;
        if (ref_allele == alt_allele) return true;
        return (ref_allele == "A" && alt_allele == "T")
               || (alt_allele == "A" && ref_allele == "T")
               || (ref_allele == "G" && alt_allele == "C")
               || (alt_allele == "G" && ref_allele == "C");
    }


    // no touchy area (PLINK Code)
    uint32_t em_phase_hethet(double known11, double known12, double known21,
                             double known22, uint32_t center_ct,
                             double* freq1x_ptr, double* freq2x_ptr,
                             double* freqx1_ptr, double* freqx2_ptr,
                             double* freq11_ptr, uint32_t* onside_sol_ct_ptr);
    uint32_t em_phase_hethet_nobase(uint32_t* counts, uint32_t is_x1,
                                    uint32_t is_x2, double* freq1x_ptr,
                                    double* freq2x_ptr, double* freqx1_ptr,
                                    double* freqx2_ptr, double* freq11_ptr);
    double calc_lnlike(double known11, double known12, double known21,
                       double known22, double center_ct_d, double freq11,
                       double freq12, double freq21, double freq22,
                       double half_hethet_share, double freq11_incr);
    std::mutex m_clump_mutex;
};

#endif /* SRC_GENOTYPE_HPP_ */

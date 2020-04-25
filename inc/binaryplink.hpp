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


#ifndef BINARYPLINK
#define BINARYPLINK

#include "commander.hpp"
#include "genotype.hpp"
#include "misc.hpp"
#include "reporter.hpp"
#include <functional>
class BinaryPlink : public Genotype
{
public:
    /*!
     * \brief Constructor of BinaryPlink
     * \param file_list is the file contain all genotype file prefix
     * \param file is the prefix of genotype file
     * \param thread is the allowed number of thread
     * \param ignore_fid indicate if we should assume FID is absent from all
     * input
     * \param keep_nonfounder indicate if user wish to include non-founders in
     * their model
     * \param keep_ambig indicate if user wish to include ambiguous SNPs
     * \param is_ref indicate if this object should be reference (T) or target
     * (F)
     * \param reporter is the logger
     */
    BinaryPlink() {}
    BinaryPlink(const GenoFile& geno, const Phenotype& pheno,
                const std::string& delim, Reporter* reporter);
    ~BinaryPlink();

protected:
    std::vector<uintptr_t> m_sample_mask;
    std::streampos m_prev_loc = 0;
    std::vector<Sample_ID> gen_sample_vector();
    void
    gen_snp_vector(const std::vector<IITree<size_t, size_t>>& exclusion_regions,
                   const std::string& out_prefix, Genotype* target = nullptr);
    bool calc_freq_gen_inter(const QCFiltering& filter_info, const std::string&,
                             Genotype* target = nullptr);
    void check_bed(const std::string& bed_name, size_t num_marker,
                   uintptr_t& bed_offset);
    size_t transverse_bed_for_snp(
        const std::vector<IITree<size_t, size_t>>& exclusion_regions,
        const std::string mismatch_snp_record_name, const size_t idx,
        const uintptr_t unfiltered_sample_ct4, const uintptr_t bed_offset,
        std::unique_ptr<std::istream> bim,
        std::unordered_set<std::string>& duplicated_snps,
        std::unordered_set<std::string>& processed_snps,
        std::vector<bool>& retain_snp, bool& chr_error, bool& sex_error,
        Genotype* genotype);
    std::unordered_set<std::string>
    get_founder_info(std::unique_ptr<std::istream>& famfile);
    inline void read_genotype(uintptr_t* __restrict genotype,
                              const std::streampos byte_pos,
                              const size_t& file_idx)
    {
        // first, generate the mask to mask out the last few byte that we don't
        // want (if our sample number isn't a multiple of 16, it is possible
        // that there'll be trailling bytes that we don't want
        const uintptr_t final_mask =
            get_final_mask(static_cast<uint32_t>(m_founder_ct));
        const uintptr_t unfiltered_sample_ct4 =
            (m_unfiltered_sample_ct + 3) / 4;

        // now we start reading / parsing the binary from the file
        assert(unfiltered_sample_ct);

        m_genotype_file.read(m_genotype_file_names[file_idx] + ".bed", byte_pos,
                             unfiltered_sample_ct4,
                             reinterpret_cast<char*>(m_tmp_genotype.data()));
        if (m_unfiltered_sample_ct == m_founder_ct)
        {
            copy_quaterarr_nonempty_subset(
                m_tmp_genotype.data(), m_sample_for_ld.data(),
                static_cast<uint32_t>(m_unfiltered_sample_ct),
                static_cast<uint32_t>(m_founder_ct), genotype);
        }
        else
        {
            // not sure why the genotype = m_tmp_genotype.data() worked before
            // but not now.
            for (size_t i = 0; i < m_tmp_genotype.size(); ++i)
            { *genotype++ = m_tmp_genotype[i]; }
            // genotype = m_tmp_genotype.data();
            genotype[(m_unfiltered_sample_ct - 1) / BITCT2] &= final_mask;
        }
    }

    virtual void
    read_score(std::vector<PRS>& prs_list,
               const std::vector<size_t>::const_iterator& start_idx,
               const std::vector<size_t>::const_iterator& end_idx,
               bool reset_zero, bool ultra = false);

    // modified version of the
    // single_marker_freqs_and_hwe function from PLINK (plink_filter.c)
    // we remove the HWE calculation (we don't want to include that yet,
    // as that'd require us to implement a version for bgen)
    inline void single_marker_freqs_and_hwe(
        uintptr_t unfiltered_sample_ctl2, uintptr_t* lptr,
        uintptr_t* sample_include2, uintptr_t* founder_include2,
        uintptr_t sample_ct, uint32_t* ll_ctp, uint32_t* lh_ctp,
        uint32_t* hh_ctp, uintptr_t sample_f_ct, uint32_t* ll_ctfp,
        uint32_t* lh_ctfp, uint32_t* hh_ctfp)
    {
        uint32_t tot_a = 0;
        uint32_t tot_b = 0;
        uint32_t tot_c = 0;
        uint32_t tot_a_f = 0;
        uint32_t tot_b_f = 0;
        uint32_t tot_c_f = 0;
        uintptr_t* lptr_end = &(lptr[unfiltered_sample_ctl2]);
        uintptr_t loader;
        uintptr_t loader2;
        uintptr_t loader3;
#ifdef __LP64__
        uintptr_t cur_decr = 120;
        uintptr_t* lptr_12x_end;
        unfiltered_sample_ctl2 -= unfiltered_sample_ctl2 % 12;
        while (unfiltered_sample_ctl2 >= 120)
        {
        single_marker_freqs_and_hwe_loop:
            lptr_12x_end = &(lptr[cur_decr]);
            count_3freq_1920b((__m128i*) lptr, (__m128i*) lptr_12x_end,
                              (__m128i*) sample_include2, &tot_a, &tot_b,
                              &tot_c);
            count_3freq_1920b((__m128i*) lptr, (__m128i*) lptr_12x_end,
                              (__m128i*) founder_include2, &tot_a_f, &tot_b_f,
                              &tot_c_f);
            lptr = lptr_12x_end;
            sample_include2 = &(sample_include2[cur_decr]);
            founder_include2 = &(founder_include2[cur_decr]);
            unfiltered_sample_ctl2 -= cur_decr;
        }
        if (unfiltered_sample_ctl2)
        {
            cur_decr = unfiltered_sample_ctl2;
            goto single_marker_freqs_and_hwe_loop;
        }
#else
        uintptr_t* lptr_twelve_end =
            &(lptr[unfiltered_sample_ctl2 - unfiltered_sample_ctl2 % 12]);
        while (lptr < lptr_twelve_end)
        {
            count_3freq_48b(lptr, sample_include2, &tot_a, &tot_b, &tot_c);
            count_3freq_48b(lptr, founder_include2, &tot_a_f, &tot_b_f,
                            &tot_c_f);
            lptr = &(lptr[12]);
            sample_include2 = &(sample_include2[12]);
            founder_include2 = &(founder_include2[12]);
        }
#endif
        while (lptr < lptr_end)
        {
            loader = *lptr++;
            loader2 = *sample_include2++;
            loader3 = (loader >> 1) & loader2;
            loader2 &= loader;
            // N.B. because of the construction of sample_include2, only
            // even-numbered bits can be present here.  So popcount2_long is
            // safe.
            tot_a += popcount2_long(loader2);
            tot_b += popcount2_long(loader3);
            tot_c += popcount2_long(loader & loader3);
            loader2 = *founder_include2++;
            loader3 = (loader >> 1) & loader2;
            loader2 &= loader;
            tot_a_f += popcount2_long(loader2);
            tot_b_f += popcount2_long(loader3);
            tot_c_f += popcount2_long(loader & loader3);
        }
        *hh_ctp = tot_c;
        *lh_ctp = tot_b - tot_c;
        *ll_ctp = sample_ct - tot_a - *lh_ctp;
        *hh_ctfp = tot_c_f;
        *lh_ctfp = tot_b_f - tot_c_f;
        *ll_ctfp = sample_f_ct - tot_a_f - *lh_ctfp;
    }
};

#endif

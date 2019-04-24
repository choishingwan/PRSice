
// This file is part of PRSice2.0, copyright (C) 2016-2017
// Shing Wan Choi, Jack Euesden, Cathryn M. Lewis, Paul F. O’Reilly
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
    BinaryPlink(const std::string& file_list, const std::string& file,
                const uint32_t thread, const bool ignore_fid,
                const bool keep_nonfounder, const bool keep_ambig,
                const bool is_ref, Reporter& reporter);
    BinaryPlink() {}
    ~BinaryPlink();

protected:
    std::vector<uintptr_t> m_sample_mask;
    std::string m_cur_file;
    std::ifstream m_bed_file;
    std::streampos m_prev_loc = 0;

    /*!
     * \brief Generate the sample vector
     * \return Vector containing the sample information
     */
    std::vector<Sample_ID> gen_sample_vector();
    /*!
     * \brief Function to generate the SNP vector
     * \param out_prefix is the output prefix
     * \param maf_threshold is the maf threshold
     * \param maf_filter is the boolean indicate if we want to perform maf
     * filtering
     * \param geno_threshold is the geno threshold
     * \param geno_filter is the boolean indicate if we want to perform geno
     * filtering
     * \param exclusion is the exclusion region
     * \param target  contain the target genotype information (if is reference)
     * \return a vector of SNP
     */
    std::vector<SNP>
    gen_snp_vector(const std::string& out_prefix, const double& maf_threshold,
                   const bool maf_filter, const double& geno_threshold,
                   const bool geno_filter, const double& /*hard_threshold*/,
                   const bool /*hard_coded*/, const double& /*info_threshold*/,
                   const bool /*info_filter*/, cgranges_t* exclusion_regions,
                   Genotype* target = nullptr);
    void gen_snp_vector(const std::string& out_prefix,
                        cgranges_t* exclusion_regions,
                        Genotype* target = nullptr);
    void calc_freq_gen_inter(const double& maf_threshold,
                             const double& geno_threshold, const double&,
                             const double& /*info_threshold*/,
                             const bool maf_filter, const bool geno_filter,
                             const bool, const bool /*info_filter*/,
                             Genotype* target = nullptr);
    /*!
     * \brief This function is use to check the bed version. Most importantly,
     *        this should give the correct bed_offset for file reading
     * \param bed_name is the name of the file
     * \param num_marker is the number of markers, use for checking the size of
     *        the bed file
     * \param bed_offset is the offset for bed file reading
     */
    void check_bed(const std::string& bed_name, size_t num_marker,
                   uintptr_t& bed_offset);
    /*!
     * \brief Function for read in the genotype file directly from the binary
     * plink file
     *
     * \param genotype is the vector use to store the binary information from
     * the file
     *
     * \param byte_pos is the stream pos of the SNP in the file, this allow us
     * to jump directly to the target SNP without iterate through the file
     *
     * \param file_name is the name of the file
     */
    inline void read_genotype(uintptr_t* genotype,
                              const std::streampos byte_pos,
                              const std::string& file_name)
    {
        // first, generate the mask to mask out the last few byte that we don't
        // want (if our sample number isn't a multiple of 16, it is possible
        // that there'll be trailling bytes that we don't want
        uintptr_t final_mask =
            get_final_mask(static_cast<uint32_t>(m_founder_ct));

        // we don't need to check if m_cur_file is empty because empty equals
        // only to empty, which shouldn't happen
        if (m_cur_file != file_name) {
            if (m_bed_file.is_open()) {
                m_bed_file.close();
            }
            std::string bedname = file_name + ".bed";
            // open the bed file in binary mode
            m_bed_file.open(bedname.c_str(), std::ios::binary);
            if (!m_bed_file.is_open()) {
                throw std::runtime_error(std::string(
                    "Error: Cannot open bed file: " + file_name + ".bed"));
            }
            // we reset the prev_loc value to 0 as this is a new file
            m_prev_loc = 0;
            // update teh current file name
            m_cur_file = file_name;
        }
        // if our current position does not equal to the required read position,
        // we will try to seek from the beginning of file to the target location
        if ((m_prev_loc != byte_pos)
            && !m_bed_file.seekg(byte_pos, std::ios_base::beg))
        {
            // if the location is not equal and seek fail, we have problem
            // reading the bed file
            throw std::runtime_error("Error: Cannot seek within the bed file!");
        }
        // now we start reading / parsing the binary from the file
        if (load_and_collapse_incl(
                static_cast<uint32_t>(m_unfiltered_sample_ct),
                static_cast<uint32_t>(m_founder_ct), m_founder_info.data(),
                final_mask, false, m_bed_file, m_tmp_genotype.data(), genotype))
        {
            throw std::runtime_error("Error: Cannot read the bed file!");
        }
        // directly read in the current location to avoid possible calculation
        // error
        m_prev_loc = m_bed_file.tellg();
    }
    /*!
     * \brief read_score is the master function for performing the score reading
     * \param index contain the index of SNPs that we should read from
     * \param reset_zero is a boolean indicate if we want to reset the score to
     * 0
     */
    void read_score(const std::vector<size_t>& index, bool reset_zero);
    /*!
     * \brief read_score is the master function for performing the score reading
     * \param start_index is the index of SNP that we should start reading from
     * \param end_bound is the index of the first SNP for us to ignore
     * \param region_index is the index of the region of interest
     * \param reset_zero is a boolean indicate if we want to reset the score to
     * 0
     */
    void read_score(const size_t start_index, const size_t end_bound,
                    const size_t region_index, bool reset_zero);
    /*!
     * \brief This is a slightly modified version of load_and_collapse_incl copy
     * from PLINK2, main difference is the use of ifstream
     *
     * \param unfiltered_sample_ct is the number of unfiltered sample
     * \param sample_ct is the number of sample we want
     * \param sample_include is the binary vector indicate if the specific
     * sample is required
     *
     * \param final_mask is the mask to mask out un used regions
     * \param do_reverse not sure what it does, but should always be false here
     * \param bedfile is the ifstream pointing to the plink binary file
     * \param rawbuf is the vector for storing the raw genotype read, which will
     * then be transformed into the required genotype data format
     *
     * \param mainbuf is the vector for storing the final genotype data
     * \return 0 if the read is sucessful, and anyother number if it is not
     */
    uint32_t load_and_collapse_incl(uint32_t unfiltered_sample_ct,
                                    uint32_t sample_ct,
                                    const uintptr_t* __restrict sample_include,
                                    uintptr_t final_mask, uint32_t do_reverse,
                                    std::ifstream& bedfile,
                                    uintptr_t* __restrict rawbuf,
                                    uintptr_t* __restrict mainbuf)
    {
        assert(unfiltered_sample_ct);
        uint32_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
        // if we don't perform selection, we can directly perform the read on
        // the mainbuf
        if (unfiltered_sample_ct == sample_ct) {
            rawbuf = mainbuf;
        }
        // we try to read in the data and store it in rawbug
        if (!bedfile.read((char*) rawbuf, unfiltered_sample_ct4)) {
            return RET_READ_FAIL;
        }
        if (unfiltered_sample_ct != sample_ct) {
            // if we need to perform selection, we will remove all unwanted
            // sample and push the data forward
            copy_quaterarr_nonempty_subset(rawbuf, sample_include,
                                           unfiltered_sample_ct, sample_ct,
                                           mainbuf);
        }
        else
        {
            // if we dno't need filtering, then we simply mask out the unwanted
            // region (to avoid the leftover, if any)
            mainbuf[(unfiltered_sample_ct - 1) / BITCT2] &= final_mask;
        }
        if (do_reverse) {
            // this will never be callsed in PRSice
            reverse_loadbuf(sample_ct, (unsigned char*) mainbuf);
        }
        return 0;
    }

    inline uint32_t load_raw(uintptr_t unfiltered_sample_ct4,
                             std::ifstream& bedfile, uintptr_t* rawbuf)
    {
        // only use this if all accesses to the data involve
        // 1. some sort of mask, or
        // 2. explicit iteration from 0..(unfiltered_sample_ct-1).
        // otherwise improper trailing bits might cause a segfault, when we
        // should be ignoring them or just issuing a warning.
        if (!bedfile.read((char*) rawbuf, unfiltered_sample_ct4)) {
            return RET_READ_FAIL;
        }
    }

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
        while (unfiltered_sample_ctl2 >= 120) {
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
        if (unfiltered_sample_ctl2) {
            cur_decr = unfiltered_sample_ctl2;
            goto single_marker_freqs_and_hwe_loop;
        }
#else
        uintptr_t* lptr_twelve_end =
            &(lptr[unfiltered_sample_ctl2 - unfiltered_sample_ctl2 % 12]);
        while (lptr < lptr_twelve_end) {
            count_3freq_48b(lptr, sample_include2, &tot_a, &tot_b, &tot_c);
            count_3freq_48b(lptr, founder_include2, &tot_a_f, &tot_b_f,
                            &tot_c_f);
            lptr = &(lptr[12]);
            sample_include2 = &(sample_include2[12]);
            founder_include2 = &(founder_include2[12]);
        }
#endif
        while (lptr < lptr_end) {
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

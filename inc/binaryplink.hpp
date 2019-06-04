
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
                const size_t thread, const bool ignore_fid,
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
    std::vector<Sample_ID> gen_sample_vector(const std::string& delim);
    void gen_snp_vector(const std::vector<IITree<int, int>>& exclusion_regions,
                        const std::string& out_prefix,
                        Genotype* target = nullptr);
    void calc_freq_gen_inter(const double& maf_threshold,
                             const double& geno_threshold, const double&,
                             const bool maf_filter, const bool geno_filter,
                             const bool, const bool,
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
        const uintptr_t final_mask =
            get_final_mask(static_cast<uint32_t>(m_founder_ct));
        const uintptr_t unfiltered_sample_ct4 =
            (m_unfiltered_sample_ct + 3) / 4;
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
            std::string error_message =
                "Error: Failed to read the bed file: " + m_cur_file;
            throw std::runtime_error(error_message);
        }
        // directly read in the current location to avoid possible calculation
        // error
        m_prev_loc =
            static_cast<std::streampos>(unfiltered_sample_ct4) + byte_pos;
    }
    virtual void
    read_score(const std::vector<size_t>::const_iterator& start_idx,
               const std::vector<size_t>::const_iterator& end_idx,
               bool reset_zero, const bool use_ref_maf);
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
        return 0;
    }
};

#endif

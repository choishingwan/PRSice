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
#include <mio.hpp>
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
                const bool is_ref, Reporter* reporter);
    BinaryPlink() {}
    ~BinaryPlink();
    void init_mmap()
    {
        m_genotype_files.resize(m_genotype_file_names.size());
        std::error_code error;
        for (size_t i = 0; i < m_genotype_file_names.size(); ++i)
        {
            m_genotype_files[i].map(m_genotype_file_names[i] + ".bed", error);
            if (error) { throw std::runtime_error(error.message()); }
        }
    }

protected:
    std::vector<uintptr_t> m_sample_mask;
    std::streampos m_prev_loc = 0;
    std::vector<Sample_ID> gen_sample_vector(const std::string& delim);
    void
    gen_snp_vector(const std::vector<IITree<size_t, size_t>>& exclusion_regions,
                   const std::string& out_prefix, Genotype* target = nullptr);
    void calc_freq_gen_inter(const double& maf_threshold,
                             const double& geno_threshold, const double&,
                             const bool maf_filter, const bool geno_filter,
                             const bool, const bool,
                             Genotype* target = nullptr);
    void check_bed(const std::string& bed_name, size_t num_marker,
                   uintptr_t& bed_offset);
    inline void read_genotype(uintptr_t* __restrict genotype,
                              const unsigned long long byte_pos,
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
        // if we don't perform selection, we can directly perform the read on
        // the mainbuf
        auto&& cur_map = m_genotype_files[file_idx];
        const unsigned long long max_file_size = cur_map.mapped_length();

        // read in the genotype information to the genotype vector
        if (byte_pos + unfiltered_sample_ct4 > max_file_size)
        {
            std::string error_message =
                "Erorr: Reading out of bound: " + misc::to_string(byte_pos)
                + " " + misc::to_string(unfiltered_sample_ct4) + " "
                + misc::to_string(max_file_size);
            throw std::runtime_error(error_message);
        }
        char* geno;
        if (m_unfiltered_sample_ct == m_founder_ct)
        { geno = reinterpret_cast<char*>(genotype); }
        else
        {
            geno = reinterpret_cast<char*>(m_tmp_genotype.data());
        }
        for (unsigned long long i = 0; i < unfiltered_sample_ct4; ++i)
        {
            *geno = cur_map[byte_pos + i];
            ++geno;
        }
        if (m_unfiltered_sample_ct != m_founder_ct)
        {
            copy_quaterarr_nonempty_subset(
                m_tmp_genotype.data(), m_founder_info.data(),
                static_cast<uint32_t>(m_unfiltered_sample_ct),
                static_cast<uint32_t>(m_founder_ct), genotype);
        }
        else
        {
            genotype[(m_unfiltered_sample_ct - 1) / BITCT2] &= final_mask;
        }
    }

    virtual void
    read_score(const std::vector<size_t>::const_iterator& start_idx,
               const std::vector<size_t>::const_iterator& end_idx,
               bool reset_zero, const bool use_ref_maf);
};

#endif

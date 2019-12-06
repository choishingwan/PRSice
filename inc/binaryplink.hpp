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
                             Genotype* target = nullptr,
                             bool force_cal = false);
    void check_bed(const std::string& bed_name, size_t num_marker,
                   uintptr_t& bed_offset);
    inline void read_genotype(uintptr_t* __restrict genotype,
                              const long long byte_pos, const size_t& file_idx)
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
        if (m_unfiltered_sample_ct != m_founder_ct)
        {
            copy_quaterarr_nonempty_subset(
                m_tmp_genotype.data(), m_founder_info.data(),
                static_cast<uint32_t>(m_unfiltered_sample_ct),
                static_cast<uint32_t>(m_founder_ct), genotype);
        }
        else
        {
            genotype = m_tmp_genotype.data();
            genotype[(m_unfiltered_sample_ct - 1) / BITCT2] &= final_mask;
        }
    }

    virtual void
    read_score(const std::vector<size_t>::const_iterator& start_idx,
               const std::vector<size_t>::const_iterator& end_idx,
               bool reset_zero, bool ultra = false);
};

#endif

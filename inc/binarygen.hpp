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

#ifndef BinaryGEN_H
#define BinaryGEN_H

#include "bgen_lib.hpp"
#include "binarygen_setters.hpp"
#include "genotype.hpp"
#include "reporter.hpp"
#include <stdexcept>
#include <zlib.h>

/**
 * Potential problem:
 * Where multi-allelic variants exist in these data, they have been
 * split into a series of bi-allelic variants. This implies that
 * several variants may share the same genomic position but with
 * different alternative alleles.
 */
class BinaryGen : public Genotype
{
public:
    BinaryGen() {}
    BinaryGen(const GenoFile& geno, const Phenotype& pheno,
              const std::string& delim, Reporter* reporter);
    ~BinaryGen();

    //
    /*!
     * \brief check if the sample file is of the sample format specified by bgen
     *        or just a simple text file
     * \return
     */
    static bool check_is_sample_format(std::unique_ptr<std::istream>& input);

protected:
    typedef std::vector<std::vector<double>> Data;
    std::vector<genfile::bgen::Context> m_context_map;
    std::vector<genfile::byte_t> m_buffer1, m_buffer2;
    bool m_target_plink = false;
    bool m_ref_plink = false;
    bool m_has_external_sample = false;

    /*!
     * \brief Generate the sample vector
     * \return Vector containing the sample information
     */
    std::vector<Sample_ID> gen_sample_vector() override;
    void handle_pheno_header(std::unique_ptr<std::istream>& sample);
    void
    gen_snp_vector(const std::vector<IITree<size_t, size_t>>& exclusion_regions,
                   const std::string& out_prefix,
                   Genotype* target = nullptr) override;
    bool calc_freq_gen_inter(const QCFiltering& filter_info,
                             const std::string& prefix,
                             Genotype* genotype = nullptr) override;

    genfile::bgen::Context get_context(const size_t& idx);
    size_t get_sex_col(const std::string& header,
                       const std::string& format_line);
    /*!
     * \brief Check if the sample information and ordering of the bgen file
     *        matched the sample / phenotype file
     * \param bgen_name is the name of the bgen file
     * \param context is the context object
     * \return true if the sample is consistent
     */
    void check_sample_consistent(const genfile::bgen::Context& context,
                                 std::istream& stream);

    size_t transverse_bgen_for_snp(
        const std::vector<IITree<size_t, size_t>>& exclusion_regions,
        const std::string mismatch_snp_record_name, const size_t file_idx,
        std::unique_ptr<std::istream> bgen_file,
        std::unordered_set<std::string>& duplicated_snps,
        std::unordered_set<std::string>& processed_snps,
        std::vector<bool>& retain_snp, bool& chr_error, bool& sex_error,
        Genotype* genotype);
    inline void read_genotype(const std::unique_ptr<SNP>& snp,
                              const uintptr_t /*selected_size*/,
                              FileRead& genotype_file,
                              uintptr_t* __restrict /*tmp_genotype*/,
                              uintptr_t* __restrict genotype,
                              uintptr_t* __restrict subset_mask,
                              bool is_ref = false) override
    {
        auto [file_idx, byte_pos] = snp->get_file_info(is_ref);
        const uintptr_t unfiltered_sample_ct4 =
            (m_unfiltered_sample_ct + 3) / 4;
        if ((m_ref_plink && is_ref) || (!is_ref && m_target_plink))
        {
            genotype_file.read(m_genotype_file_names[file_idx], byte_pos,
                               unfiltered_sample_ct4,
                               reinterpret_cast<char*>(genotype));
        }
        else if (!load_and_collapse_incl(byte_pos, file_idx, genotype_file,
                                         genotype, subset_mask))
        {
            throw std::runtime_error("Error: Cannot read the bgen file!");
        }
    }
    bool load_and_collapse_incl(const std::streampos byte_pos,
                                const size_t& file_idx, FileRead& genotype_file,
                                uintptr_t* __restrict mainbuf,
                                uintptr_t* __restrict subset_mask)
    {
        assert(m_unfiltered_sample_ct);
        try
        {
            PLINK_generator setter(subset_mask, mainbuf, m_hard_threshold,
                                   m_dose_threshold);
            genfile::bgen::read_and_parse_genotype_data_block<PLINK_generator>(
                genotype_file, m_genotype_file_names[file_idx] + ".bgen",
                m_context_map[file_idx], setter, &m_buffer1, &m_buffer2,
                byte_pos);
        }
        catch (...)
        {
            return false;
        }
        return true;
    }

    void count_and_read_genotype(const std::unique_ptr<SNP>&) override;
    void read_score(std::vector<PRS>& prs_list,
                    const std::vector<size_t>::const_iterator& start_idx,
                    const std::vector<size_t>::const_iterator& end_idx,
                    bool reset_zero) override;
    void hard_code_score(std::vector<PRS>& prs_list,
                         const std::vector<size_t>::const_iterator& start_idx,
                         const std::vector<size_t>::const_iterator& end_idx,
                         bool reset_zero);
    void dosage_score(std::vector<PRS>& prs_list,
                      const std::vector<size_t>::const_iterator& start_idx,
                      const std::vector<size_t>::const_iterator& end_idx,
                      bool reset_zero);

    /*
     * Different structures use for reading in the bgen info
     */
};

#endif

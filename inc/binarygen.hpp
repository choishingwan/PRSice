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
    std::vector<Sample_ID> gen_sample_vector();
    void handle_pheno_header(std::unique_ptr<std::istream>& sample);
    void
    gen_snp_vector(const std::vector<IITree<size_t, size_t>>& exclusion_regions,
                   const std::string& out_prefix, Genotype* target = nullptr);
    bool calc_freq_gen_inter(const QCFiltering& filter_info,
                             const std::string& prefix,
                             Genotype* genotype = nullptr);

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
    /*!
     * \brief This function will read in the bgen probability data and transform
     * it into PLINK bianry
     *
     * \param genotype is the vector containing the plink binary
     * \param byte_pos is the streampos of the bgen file, for quick seek
     * \param file_name is the file name of the bgen file
     */
    inline void read_genotype(uintptr_t* genotype,
                              const std::streampos byte_pos,
                              const size_t& file_idx)
    {
        const uintptr_t unfiltered_sample_ct4 =
            (m_unfiltered_sample_ct + 3) / 4;
        if (m_ref_plink)
        {
            // when m_ref_plink is set, it suggest we are using the
            // intermediate, which is a binary plink format. Therefore we can
            // directly read in the binary data to genotype, this should already
            // be well formated when we write it into the file
            // intermediate file, so don't need .bgen
            std::string file_name = m_genotype_file_names[file_idx];
            m_genotype_file.read(file_name, byte_pos, unfiltered_sample_ct4,
                                 reinterpret_cast<char*>(genotype));
        }
        // if not, we will try to parse the binary GEN format into a plink
        // format (using the default intermediate flag = false)
        else if (load_and_collapse_incl(byte_pos, file_idx, genotype))
        {
            throw std::runtime_error("Error: Cannot read the bgen file!");
        }
    }

    /*!
     * \brief We modify the load_and_collapse_incl function from PLINK to read
     * in and reformat the bgen data into the binary PLINK data
     *
     * \param byte_pos is the location of the file
     * \param file_name is the file name
     * \param unfiltered_sample_ct is the number of unfiltered sample
     * \param sample_ct is the number of sample we want
     * \param sample_include is the bianry vector containing the boolean
     * indicate if we want the target sample
     *
     * \param final_mask is a mask to remove trailing bytes
     * \param do_reverse is implemented in PLINK but we will not use
     * \param rawbuf is the raw vector to store the information
     * \param mainbuf is the result vector containing the genotype information
     * \param intermediate is boolean indicate if we want to generate an
     * intermediate file
     * \return 0 if sucessful
     */
    uint32_t load_and_collapse_incl(const std::streampos byte_pos,
                                    const size_t& file_idx,
                                    uintptr_t* __restrict mainbuf)
    {
        // check sample size != 0
        assert(m_unfiltered_sample_ct);
        // we first check if we will read in a new file
        auto&& context = m_context_map[file_idx];
        // we initailize the PLINK generator with the m_sample_include
        // vector, the mainbuf (result storage) and also the hard threshold.
        // We don't need to bother about founder or founder info here as all
        // samples in bgen are considered as founder (at least this is the
        // case at the moment)
        // TODO, might want to make this struct a member so that we don't
        // need to re-initialize it? Though it only contain pointers and
        // doesn't have any big structures
        PLINK_generator setter(&m_calculate_prs, mainbuf, m_hard_threshold,
                               m_dose_threshold);
        // we can now use the bgen library to parse the BGEN input and
        // transform it into PLINK format (NOTE: The
        // read_and_parse_genotype_data_block function has been modified
        // such that it will always call .sample_completed() when finish
        // reading each sample. This allow for a more elegant implementation
        // on our side
        std::string file_name = m_genotype_file_names[file_idx] + ".bgen";
        // WARNING: Problem here
        genfile::bgen::read_and_parse_genotype_data_block<PLINK_generator>(
            m_genotype_file, file_name, context, setter, &m_buffer1, &m_buffer2,
            byte_pos);
        // output from load_raw should have already copied all samples
        // to the front without the need of subseting
        // mainbuf should contains the information
        // update the m_prev_loc counter. However, for bgen, it is most likely
        // that we will need to seek for every SNP as there are padded data
        // between each SNP's genotype
        return 0;
    }

    void read_score(std::vector<PRS>& prs_list,
                    const std::vector<size_t>::const_iterator& start_idx,
                    const std::vector<size_t>::const_iterator& end_idx,
                    bool reset_zero, bool ultra = false);
    void hard_code_score(std::vector<PRS>& prs_list,
                         const std::vector<size_t>::const_iterator& start_idx,
                         const std::vector<size_t>::const_iterator& end_idx,
                         bool reset_zero, bool ultra = false);
    void dosage_score(std::vector<PRS>& prs_list,
                      const std::vector<size_t>::const_iterator& start_idx,
                      const std::vector<size_t>::const_iterator& end_idx,
                      bool reset_zero);

    /*
     * Different structures use for reading in the bgen info
     */
};

#endif

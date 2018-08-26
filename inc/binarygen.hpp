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
#ifndef BinaryGEN_H
#define BinaryGEN_H

#include "bgen_lib.hpp"
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
    /*!
     * \brief Constructor of binarygen object
     * \param commander contains the user input
     * \param reporter is the logger
     * \param is_ref indicate if this is the reference panel
     */
    BinaryGen(const Commander& commander, Reporter& reporter,
              const bool is_ref = false);
    ~BinaryGen();

private:
    typedef std::vector<std::vector<double>> Data;
    std::unordered_map<std::string, genfile::bgen::Context> m_context_map;
    std::vector<genfile::byte_t> m_buffer1, m_buffer2;
    std::ifstream m_bgen_file;
    std::string m_cur_file;
    std::string m_intermediate_file;
    std::streampos m_prev_loc = 0;
    bool m_intermediate = false;
    bool m_target_plink = false;
    bool m_ref_plink = false;

    /*!
     * \brief Generate the sample vector
     * \return Vector containing the sample information
     */
    std::vector<Sample_ID> gen_sample_vector();
    //
    /*!
     * \brief check if the sample file is of the sample format specified by bgen
     *        or just a simple text file
     * \return
     */
    bool check_is_sample_format();
    /*!
     * \brief Generate the SNP vector
     * \param commander contain the user input
     * \param exclusion contain the exclusion region information
     * \param target contain the target genotype (only for m_is_ref)
     * \return a vector containing the SNP vector or a empty factor if m_is_ref
     * = T
     */
    std::vector<SNP> gen_snp_vector(const Commander& commander,
                                    Region& exclusion,
                                    Genotype* target = nullptr);
    /*!
     * \brief Read in the context information for the bgen. This will propergate
     * the m_context_map
     * \param prefix is the input bgen file prefix
     */
    void get_context(std::string& prefix);
    /*!
     * \brief Check if the sample information and ordering of the bgen file
     *        matched the sample / phenotype file
     * \param bgen_name is the name of the bgen file
     * \param context is the context object
     * \return true if the sample is consistent
     */
    bool check_sample_consistent(const std::string& bgen_name,
                                 const genfile::bgen::Context& context);

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
                              const std::string& file_name)
    {
        const uintptr_t unfiltered_sample_ct4 =
            (m_unfiltered_sample_ct + 3) / 4;
        if (m_ref_plink) {
            // when m_ref_plink is set, it suggest we are using the
            // intermediate, which is a binary plink format. Therefore we can
            // directly read from the file
            if (file_name != m_cur_file || !m_bgen_file.is_open()) {
                if (m_bgen_file.is_open()) m_bgen_file.close();
                // read in the file in binary format
                m_bgen_file.open(file_name.c_str(), std::ifstream::binary);
                if (!m_bgen_file.is_open()) {
                    std::string error_message =
                        "Error: Cannot open bgen file: " + file_name;
                    throw std::runtime_error(error_message);
                }
                // update the binary file name
                m_cur_file = file_name;
                // set prev_location to 0
                m_prev_loc = 0;
            }
            if (byte_pos != m_prev_loc
                && !m_bgen_file.seekg(byte_pos, std::ios_base::beg))
            {
                // if the location is not equal and seek fail, we have problem
                // reading the bed file
                throw std::runtime_error(
                    "Error: Cannot seek within the bgen file!");
            }
            // directly read in the binary data to genotype, this should already
            // be well formated when we write it into the fiel
            m_bgen_file.read((char*) genotype, unfiltered_sample_ct4);
            // update the location to previous location
            m_prev_loc = byte_pos;
        }
        // if not, we will try to parse the binary GEN format into a plink
        // format (using the default intermediate flag = false)
        else if (load_and_collapse_incl(byte_pos, file_name,
                                        m_founder_info.data(), genotype))
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
                                    const std::string& file_name,
                                    uintptr_t* __restrict mainbuf,
                                    bool intermediate = false)
    {
        // check sample size != 0
        assert(m_unfiltered_sample_ct);
        // we first check if we will read in a new file
        if (file_name != m_cur_file || !m_bgen_file.is_open()) {
            if (m_bgen_file.is_open()) m_bgen_file.close();

            std::string bgen_name = file_name + ".bgen";
            if (intermediate) {
                // if it is the intermeidate file, we don't need to add the
                // suffix
                bgen_name = file_name;
            }
            // we can now open the file in binary mode
            m_bgen_file.open(bgen_name.c_str(), std::ifstream::binary);
            if (!m_bgen_file.is_open()) {
                std::string error_message =
                    "Error: Cannot open bgen file: " + bgen_name;
                throw std::runtime_error(error_message);
            }
            // reset the name of current file
            m_cur_file = bgen_name;
            // and reset the prev_loc counter
            m_prev_loc = 0;
        }
        if (m_prev_loc != byte_pos
            || !m_bgen_file.seekg(byte_pos, std::ios_base::beg))
        {
            // if the location is not equal and seek fail, we have problem
            // reading the bed file
            throw std::runtime_error(
                "Error: Cannot seek within the bgen file!");
        }

        if (!intermediate) {
            // we are not using intermediate file at the moment, so this must be
            // in bgen format obtain the context information of the current bgen
            // file
            auto&& context = m_context_map[file_name];
            // we initailize the PLINK generator with the m_sample_include
            // vector, the mainbuf (result storage) and also the hard threshold.
            // We don't need to bother about founder or founder info here as all
            // samples in bgen are considered as founder (at least this is the
            // case at the moment)
            PLINK_generator setter(&m_sample_include, mainbuf,
                                   m_hard_threshold);
            // we can now use the bgen library to parse the BGEN input and
            // transform it into PLINK format (NOTE: The
            // read_and_parse_genotype_data_block function has been modified
            // such that it will always call .sample_completed() when finish
            // reading each sample. This allow for a more elegant implementation
            // on our side
            genfile::bgen::read_and_parse_genotype_data_block<PLINK_generator>(
                m_bgen_file, context, setter, &m_buffer1, &m_buffer2, false);
            // output from load_raw should have already copied all samples
            // to the front without the need of subseting
        }
        else
        {
            const uintptr_t unfiltered_sample_ct4 =
                (m_unfiltered_sample_ct + 3) / 4;
            // directly read from the intermediate file. No need to worry about
            // the transformation as that should be dealt with when we generate
            // the intermidate file
            m_bgen_file.read((char*) mainbuf, unfiltered_sample_ct4);
        }
        // mainbuf should contains the information
        // update the m_prev_loc counter. However, for bgen, it is most likely
        // that we will need to seek for every SNP as there are padded data
        // between each SNP's genotype
        m_prev_loc = m_bgen_file.tellg();
        return 0;
    }

    void read_score(std::vector<size_t>& index, bool reset_zero);
    void hard_code_score(std::vector<size_t>& index, int32_t homcom_weight,
                         uint32_t het_weight, uint32_t homrar_weight,
                         bool set_zero);
    void dosage_score(std::vector<size_t>& index, uint32_t homcom_weight,
                      uint32_t het_weight, uint32_t homrar_weight,
                      bool set_zero);
    void read_score(size_t start_index, size_t end_bound,
                    const size_t region_index, bool set_zero);
    void hard_code_score(size_t start_index, size_t end_bound,
                         uint32_t homcom_weight, uint32_t het_weight,
                         uint32_t homrar_weight, const size_t region_index,
                         bool set_zero);
    void dosage_score(size_t start_index, size_t end_bound,
                      uint32_t homcom_weight, uint32_t het_weight,
                      uint32_t homrar_weight, const size_t region_index,
                      bool set_zero);


    /*
     * Different structures use for reading in the bgen info
     */
    struct PRS_Interpreter
    {
        ~PRS_Interpreter(){};
        PRS_Interpreter(std::vector<PRS>* sample_prs,
                        std::vector<uintptr_t>* sample_inclusion,
                        MISSING_SCORE missing)
            : m_sample_prs(sample_prs)
            , m_sample_inclusion(sample_inclusion)
            , m_missing_score(missing)
        {
            // m_sample contains only samples extracted
            // m_score.resize(m_sample_prs->size(), 0);

            m_miss_count = (missing != MISSING_SCORE::SET_ZERO);
            m_sample_missing_index.reserve(m_sample_prs->size());
        }
        void set_stat(double stat, uint32_t homcom_weight, uint32_t het_weight,
                      uint32_t homrar_weight, bool flipped, bool not_first)
        {
            m_stat = stat;
            m_flipped = flipped;
            m_homcom_weight = homcom_weight;
            m_het_weight = het_weight;
            m_homrar_weight = homrar_weight;
            m_not_first = not_first;
            if (m_flipped) {
                std::swap(m_homcom_weight, m_homrar_weight);
            }
        }
        void initialise(std::size_t number_of_samples,
                        std::size_t number_of_alleles)
        {
            m_prs_sample_i = 0;
            m_bgen_sample_i = 0;
            m_sample_missing_index.clear();
            // std::fill(m_sample_missing.begin(), m_sample_missing.end(),
            // false);  std::fill(m_score.begin(), m_score.end(), 0.0);
        }
        void set_min_max_ploidy(uint32_t min_ploidy, uint32_t max_ploidy,
                                uint32_t min_entries, uint32_t max_entries)
        {
        }

        bool set_sample(std::size_t i)
        {
            m_entry_i = 0;
            m_is_missing = false;
            m_start_geno = false;
            m_sum = 0.0;
            m_bgen_sample_i = i;
            return IS_SET(m_sample_inclusion->data(), m_bgen_sample_i);
        }

        void set_number_of_entries(std::size_t ploidy,
                                   std::size_t number_of_entries,
                                   genfile::OrderType order_type,
                                   genfile::ValueType value_type)
        {
        }

        // Called once for each genotype (or haplotype) probability per sample.
        void set_value(uint32_t, double value)
        {
            int geno = 2 - m_entry_i;
            auto&& sample_prs = (*m_sample_prs)[m_prs_sample_i];
            // for bgen 1.1, we will still get set_value even for missing data
            // however, all values will be 0
            // Therefore for missing sample, we would've add 1, which is ok if
            // MISSING_SCORE !=SET_ZERO
            switch (geno)
            {
            default:
                sample_prs.num_snp =
                    sample_prs.num_snp * (m_not_first || m_start_geno) + 1;
                sample_prs.prs = sample_prs.prs * (m_not_first || m_start_geno)
                                 + m_homcom_weight * value * m_stat * 0.5;
                break;
            case 1:
                sample_prs.num_snp =
                    sample_prs.num_snp * (m_not_first || m_start_geno) + 1;
                sample_prs.prs = sample_prs.prs * (m_not_first || m_start_geno)
                                 + m_het_weight * value * m_stat * 0.5;
                break;
            case 2:
                sample_prs.num_snp =
                    sample_prs.num_snp * (m_not_first || m_start_geno) + 1;
                sample_prs.prs = sample_prs.prs * (m_not_first || m_start_geno)
                                 + m_homrar_weight * value * m_stat * 0.5;
            }
            m_total_prob += value * geno;
            ++m_entry_i;
            m_sum += value;
            m_start_geno = true;
        }
        // call if sample is missing
        void set_value(uint32_t, genfile::MissingValue value)
        {
            m_is_missing = true;
        }

        void sample_completed()
        {
            if (m_sum == 0.0 || m_is_missing) {
                // this is a missing sample
                m_sample_missing_index.push_back(m_prs_sample_i);
                // remove the problematic count if needed
                (*m_sample_prs)[m_prs_sample_i].num_snp -= m_miss_count;
            }
            // go to next sample
            ++m_prs_sample_i;
        }

        void finalise()
        {
            size_t num_miss = m_sample_missing_index.size();
            size_t num_prs = m_sample_prs->size();
            if (num_prs == num_miss) {
                // all samples are missing
                return;
            }
            double expected_value =
                m_total_prob / (((double) num_prs - (double) num_miss) * 2);
            // worth separating centre score out
            size_t miss_index = m_sample_missing_index.front();
            switch (m_missing_score)
            {
            default:
                // centre score
                for (auto&& index : m_sample_missing_index) {
                    (*m_sample_prs)[index].prs += expected_value;
                }
                break;
            case MISSING_SCORE::CENTER:
                for (size_t i = 0; i < num_prs; ++i) {
                    if (miss_index == i) {
                        ++miss_index;
                        continue;
                    }
                    (*m_sample_prs)[m_prs_sample_i].prs -= expected_value;
                }
                break;
            case MISSING_SCORE::SET_ZERO:
                // do nothing
                break;
            }
        }

    private:
        std::vector<PRS>* m_sample_prs;
        std::vector<uintptr_t>* m_sample_inclusion;
        std::vector<uint32_t> m_sample_missing_index;
        double m_stat = 0.0;
        double m_total_prob = 0.0;
        double m_sum = 0.0;
        uint32_t m_homcom_weight = 0;
        uint32_t m_het_weight = 1;
        uint32_t m_homrar_weight = 2;
        uint32_t m_bgen_sample_i = 0;
        uint32_t m_prs_sample_i = 0;
        uint32_t m_entry_i = 0;
        uint32_t m_miss_count = 0;
        MISSING_SCORE m_missing_score;
        bool m_flipped = false;
        bool m_not_first = false;
        bool m_is_missing = false;
        bool m_start_geno = false;
    };

    /*!
     * \brief The PLINK_generator struct. This will be passed into the bgen
     * library and used for parsing the BGEN data. This will generate a plink
     * binary as an output (stored in the genotype vector passed during
     * construction)
     */
    struct PLINK_generator
    {
        PLINK_generator(std::vector<uintptr_t>* sample, uintptr_t* genotype,
                        double hard_threshold)
            : m_sample(sample)
            , m_genotype(genotype)
            , m_hard_threshold(hard_threshold)
        {
        }
        void initialise(std::size_t number_of_samples,
                        std::size_t number_of_alleles)
        {
            // initialize is called when we meet a new SNP
            // we therefore want to reset the index and shift, use for indicate
            // where we push the byte onto genotype
            m_index = 0;
            m_shift = 0;
            // we also clean the running stat (not RS ID), so that we can go
            // through another round of calculation of mean and sd
            rs.clear();
        }
        /*!
         * \brief set_min_max_ploidy is a function required by BGEN. As we only
         * allow diploid, this function should not concern us (as we will error
         * out otherwise
         *
         * \param min_ploidy
         * \param max_ploidy
         * \param min_entries
         * \param max_entries
         */
        void set_min_max_ploidy(uint32_t min_ploidy, uint32_t max_ploidy,
                                uint32_t min_entries, uint32_t max_entries)
        {
        }
        /*!
         * \brief set_sample is the function to be called when BGEN start to
         * process a new sample
         *
         * \param i is the ID of the sample (index on vector)
         * \return True if we want to process this sample, false if not
         */
        bool set_sample(std::size_t i)
        {
            // we set the sample index to i
            m_sample_i = i;
            // set the genotype to 1 (missing)
            m_geno = 1;
            // and set the entry number to 3 (homozygous effective) bit = 11
            m_entry_i = 3;
            // we also reset the hard_prob to 0
            m_hard_prob = 0.0;
            // and expected value to 0
            m_exp_value = 0.0;
            // then we determine if we want to include sample using by
            // consulting the flag on m_sample
            return IS_SET(m_sample->data(), m_sample_i);
        }
        /*!
         * \brief set_number_of_entries is yet another function required by bgen
         * library that we don't use
         *
         * \param ploidy
         * \param number_of_entries
         * \param order_type
         * \param value_type
         */
        void set_number_of_entries(std::size_t ploidy,
                                   std::size_t number_of_entries,
                                   genfile::OrderType order_type,
                                   genfile::ValueType value_type)
        {
        }

        // Called once for each genotype (or haplotype) probability per sample.
        /*!
         * \brief set_value is called for each genotype or haplotype probability
         * per sample
         * \param value where value is the probability of the genotype
         */
        void set_value(uint32_t, double value)
        {
            // if the current probability is the highest and the value is higher
            // than the required threshold, we will assign it
            if (value > m_hard_prob && value >= m_hard_threshold) {
                /*
                 * Representation of each m_entry to their genotype:
                 *   m_entry   desired binary
                 *   0           00
                 *   1           00
                 *   2           01
                 *   3           11
                 *   should in theory not reach 0
                 */
                m_geno = (m_entry_i == 1) ? 0 : m_entry_i;
                m_hard_prob = value;
            }
            m_exp_value += m_geno * value;
            --m_entry_i;
        }
        /*!
         * \brief set_value will capture the missingness signature and do
         * nothing
         *
         * \param value this is the missing signature used by bgen v1.2+
         */
        void set_value(uint32_t, genfile::MissingValue value) {}
        /*!
         * \brief finalise This function is called when the SNP is processed. Do
         * nothing here
         */
        void finalise() {}
        /*!
         * \brief sample_completed is called when each sample is completed. This
         * will set the binary vector accordingly
         */
        void sample_completed()
        {
            // if m_shift is zero, it is when we meet the index for the first
            // time, therefore we want to initialize it
            if (m_shift == 0) m_genotype[m_index] = 0;
            // we now add the genotype to the vector
            m_genotype[m_index] |= m_geno << m_shift;
            // as the genotype is represented by two bit, we will +=2
            m_shift += 2;
            // if we reach the boundary, we will now add the index and reset the
            // shift
            if (m_shift == BITCT) {
                m_index++;
                m_shift = 0;
            }
            // we can now push in the expected value for this sample. This can
            // then use for the calculation of the info score
            rs.push(m_exp_value);
        }
        /*!
         * \brief info_score is the function use to calculate the info score
         * \return return the  MaCH INFO score
         */
        double info_score() const
        {
            double p = rs.mean() / 2.0;
            double p_all = 2.0 * p * (1.0 - p);
            return (rs.var() / p_all);
        }

    private:
        // is the sample inclusion vector, if bit is set, sample is required
        std::vector<uintptr_t>* m_sample;
        // is the genotype vector
        uintptr_t* m_genotype;
        misc::RunningStat rs;
        double m_hard_threshold = 0.0;
        double m_hard_prob = 0;
        double m_exp_value = 0.0;
        uintptr_t m_geno = 1;
        uint32_t m_shift = 0;
        uint32_t m_index = 0;
        size_t m_sample_i = 0;
        size_t m_entry_i = 3;
    };
};

#endif

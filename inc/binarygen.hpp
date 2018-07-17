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
    BinaryGen(const std::string& prefix, const std::string& sample_file,
              const std::string& multi_input, const size_t thread = 1,
              const bool ignore_fid = false, const bool keep_nonfounder = false,
              const bool keep_ambig = false, const bool is_ref = false, const bool intermediate=false);
    ~BinaryGen();
private:

    std::unordered_map<std::string, genfile::bgen::Context> m_context_map;
    typedef std::vector<std::vector<double>> Data;
    std::string m_cur_file;
    bool m_intermediate = false;

    std::vector<Sample_ID> gen_sample_vector();
    // check if the sample file is of the sample format specified by bgen
    // or just a simple text file
    bool check_is_sample_format();
    std::vector<SNP> gen_snp_vector(const double geno, const double maf,
                                    const double info_score,
                                    const double hard_threshold,
                                    const bool hard_coded, Region& exclusion,
                                    const std::string& out_prefix,
                                    Genotype* target = nullptr);
    void get_context(std::string& prefix);
    bool check_sample_consistent(const std::string& bgen_name,
                                 const genfile::bgen::Context& context);


    inline void read_genotype(uintptr_t* genotype,
                              const std::streampos byte_pos,
                              const std::string& file_name)
    {
        // the bgen library seems over complicated
        // try to use PLINK one. The difficulty is ifstream vs FILE
        // this is for the LD calculation
        uintptr_t final_mask = get_final_mask(m_founder_ct);
        // std::fill(m_tmp_genotype.begin(), m_tmp_genotype.end(), 0);
        // std::memset(m_tmp_genotype, 0x0, m_unfiltered_sample_ctl * 2 *
        // sizeof(uintptr_t));
        if (load_and_collapse_incl(byte_pos, file_name, m_unfiltered_sample_ct,
                                   m_founder_ct, m_founder_info.data(),
                                   final_mask, false, m_tmp_genotype.data(),
                                   genotype))
        {
            throw std::runtime_error("ERROR: Cannot read the bgen file!");
        }
    };

    // borrowed from plink
    uint32_t load_and_collapse_incl(const std::streampos byte_pos,
                                    const std::string& file_name,
                                    uint32_t unfiltered_sample_ct,
                                    uint32_t sample_ct,
                                    const uintptr_t* __restrict sample_include,
                                    uintptr_t final_mask, uint32_t do_reverse,
                                    uintptr_t* __restrict rawbuf,
                                    uintptr_t* __restrict mainbuf)
    {
        assert(unfiltered_sample_ct);
        if (m_cur_file.empty() || file_name.compare(m_cur_file) != 0
                    || !m_bgen_file.is_open())
        {
        	if (m_bgen_file.is_open()) m_bgen_file.close();
        	std::string bgen_name = file_name + ".bgen";
        	m_bgen_file.open(bgen_name.c_str(), std::ifstream::binary);
        	if (!m_bgen_file.is_open()) {
        		std::string error_message =
        				"ERROR: Cannot open bgen file: " + file_name;
        		throw std::runtime_error(error_message);
        	}
        	m_cur_file = file_name;
        }
        auto&& context = m_context_map[file_name];
        // anyway to avoid reinitializing the memory for buffer 1 and 2?
        std::vector<genfile::byte_t> buffer1, buffer2;
        m_bgen_file.seekg(byte_pos, std::ios_base::beg);

        PLINK_generator setter(&m_sample_include, mainbuf, m_hard_threshold);
        genfile::bgen::read_and_parse_genotype_data_block<PLINK_generator>(
        		m_bgen_file, context, setter, &buffer1, &buffer2, false);



        // output from load_raw should have already copied all samples
        // to the front without the need of subseting
        if (do_reverse) {
            reverse_loadbuf(sample_ct, (unsigned char*) mainbuf);
        }

        // mainbuf should contains the information
        return 0;
    }

    void read_score(std::vector<size_t>& index, bool reset_zero);
    void hard_code_score(std::vector<size_t>& index);
    void dosage_score(std::vector<size_t>& index);
    void read_score(size_t start_index, size_t end_bound,
                    const size_t region_index, bool reset_zero);
    void hard_code_score(size_t start_index, size_t end_bound,
                         const size_t region_index);
    void dosage_score(size_t start_index, size_t end_bound,
                      const size_t region_index);


    std::ifstream m_bgen_file;

    /*
     * Different structures use for reading in the bgen info
     */


    struct PRS_Interpreter
    {
        ~PRS_Interpreter(){};
        PRS_Interpreter(std::vector<PRS>* sample_prs, MODEL model,
                        MISSING_SCORE missing,
                        std::vector<uintptr_t>* sample_inclusion, double stat,
                        bool flipped)
            : m_sample_prs(sample_prs)
            , m_sample_inclusion(sample_inclusion)
            , m_model(model)
            , m_missing(missing)
            , m_stat(stat)
            , m_flipped(flipped)
        {
            // m_sample contains only samples extracted
            m_score.resize(m_sample_prs->size(), 0);
        }
        void initialise(std::size_t number_of_samples,
                        std::size_t number_of_alleles)
        {
            m_total_sample_size = number_of_samples;
        }
        void set_min_max_ploidy(uint32_t min_ploidy, uint32_t max_ploidy,
                                uint32_t min_entries, uint32_t max_entries)
        {
        }

        bool set_sample(std::size_t i)
        {
            if (!first && !exclude) {
                if (missing || m_sum == 0) {
                    m_missing_samples.push_back(m_sample_i);
                }
            }
            first = false;
            m_entry_i = 0;
            missing = false;
            m_sum = 0.0;
            m_sample_i = i;
            exclude = !IS_SET(m_sample_inclusion->data(), m_sample_i);
            if (!exclude) {
                // as this is size_t, we can't use -1. Therefore in subsequent
                // use, we should always -1
                m_score_i++;
            }
            // exclude = !m_sample->at(m_sample_i).included;
            m_num_included_samples += !exclude;
            // don't bother reading the score of anyone that is being excluded
            return !exclude;
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

            int geno = (!m_flipped) ? 2 - m_entry_i : m_entry_i;
            if (m_model == MODEL::HETEROZYGOUS && geno == 2)
                geno = 0;
            else if (m_model == MODEL::DOMINANT && geno == 2)
                geno = 1;
            else if (m_model == MODEL::RECESSIVE)
                geno = std::max(geno - 1, 0);
            m_score[m_score_i - 1] += value * geno;
            if (!exclude) m_total_prob += value * geno;
            m_entry_i++;
            m_sum += value;
        }
        // call if sample is missing
        void set_value(uint32_t, genfile::MissingValue value)
        {
            missing = true;
        }

        void finalise()
        {
            // TODO: Double check this function in case there are any problem
            // summarize the last sample's info
            if (!first && !exclude) {
                if (missing || m_sum == 0) {
                    // this is missing
                    m_missing_samples.push_back(m_sample_i);
                }
            }
            // now update the PRS score


            if (m_num_included_samples == m_missing_samples.size()) {
                valid = false;
            }
            size_t i_missing = 0;
            // missing sample should be unique and sorted, right?
            // std::sort(m_missing_samples.begin(), m_missing_samples.end());
            // m_missing_samples.erase(
            //    std::unique(m_missing_samples.begin(),
            //    m_missing_samples.end()), m_missing_samples.end());
            size_t num_miss = m_missing_samples.size();
            double mean =
                m_total_prob
                / (((double) m_num_included_samples - (double) num_miss) * 2);
            size_t idx_sample_include = 0;
            for (size_t i_sample = 0; i_sample < m_sample_prs->size();
                 ++i_sample)
            {
                if (IS_SET(m_sample_inclusion->data(), i_sample)) {
                    auto&& sample = m_sample_prs->at(idx_sample_include++);
                    if (i_missing < num_miss
                        && i_sample == m_missing_samples[i_missing])
                    {
                        if (m_missing == MISSING_SCORE::MEAN_IMPUTE)
                            sample.prs += m_stat * mean;
                        // m_prs_score->at(i_sample + m_vector_pad) +=
                        //    m_stat * mean;
                        if (m_missing != MISSING_SCORE::SET_ZERO)
                            sample.num_snp++;
                        // m_num_snps->at(i_sample + m_vector_pad)++;
                        i_missing++;
                    }
                    else
                    { // not missing sample
                        if (m_missing == MISSING_SCORE::CENTER) {
                            // if centering, we want to keep missing at 0
                            sample.prs -= m_stat * mean;
                            // m_prs_score->at(i_sample + m_vector_pad) -=
                            //    m_stat * mean;
                        }
                        // again, so that it will generate the same result as
                        // genotype file format when we are 100% certain of the
                        // genotypes
                        sample.prs += m_score[i_sample] * m_stat * 0.5;
                        // m_prs_score->at(i_sample + m_vector_pad) +=
                        //    m_score[i_sample] * m_stat * 0.5;
                        sample.num_snp++;
                        // m_num_snps->at(i_sample + m_vector_pad)++;
                    }
                }
                if (idx_sample_include >= m_num_included_samples) break;
            }
        }
        void sample_completed(){

        }
    private:
        size_t m_sample_i = 0;
        size_t m_entry_i = 0;
        size_t m_num_included_samples = 0;
        size_t m_total_sample_size = 0;
        size_t m_score_i = 0;
        double m_sum = 0.0;
        bool missing = false;
        bool first = true;
        bool exclude = false;
        bool valid = true;
        std::vector<PRS>* m_sample_prs;
        std::vector<uintptr_t>* m_sample_inclusion;
        std::vector<double> m_score;
        std::vector<size_t> m_missing_samples;
        MODEL m_model;
        MISSING_SCORE m_missing;
        double m_stat = 0.0;
        double m_total_prob = 0.0;
        bool m_flipped;
    };


    struct PLINK_generator
    {
        PLINK_generator(std::vector<uintptr_t>* sample, uintptr_t* genotype,
                        double hard_threshold)
            : m_sample(sample)
            , m_genotype(genotype)
        	, m_hard_threshold(hard_threshold)
        {}
        void initialise(std::size_t number_of_samples,
                        std::size_t number_of_alleles)
        {
        	m_index = 0;
        	m_shift =0;
        	rs.clear();
        }
        void set_min_max_ploidy(uint32_t min_ploidy, uint32_t max_ploidy,
                                uint32_t min_entries, uint32_t max_entries)
        {}

        bool set_sample(std::size_t i)
        {
        	m_sample_i = i;
        	m_geno = 1;
        	m_entry_i = 0;
        	m_hard_prob = 0.0;
        	m_exp_value = 0.0;
        	m_include_snp = IS_SET(m_sample->data(), m_sample_i);
            return m_include_snp;
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
            if (value > m_hard_prob && value >= m_hard_threshold) {
                // geno = 2 - m_entry_i;
            	m_geno = (m_entry_i == 0) ? 0 : m_entry_i + 1;
                m_hard_prob = value;
            }
            m_entry_i++;
        }
        // call if sample is missing
        void set_value(uint32_t, genfile::MissingValue value) {}

        void finalise()
        {}
        void sample_completed(){
            if (m_shift == 0)
                m_genotype[m_index] = 0; // match behaviour of binaryplink
            m_genotype[m_index] |= m_geno << m_shift;
            m_shift += 2;
            if (m_shift == BITCT) {
                m_index++;
                m_shift = 0;
            }
        }
        double info_score() const{
            double p = rs.mean() / 2.0;
            double p_all = 2.0 * p * (1.0 - p);
            return(rs.var() / p_all);

        }
    private:
        std::vector<uintptr_t>* m_sample;
        uintptr_t* m_genotype;
        misc::RunningStat rs;
        double m_hard_threshold = 0.0;
        double m_hard_prob = 0;
        double m_exp_value = 0.0;
        uintptr_t m_geno = 1;
        uint32_t m_shift = 0;
        uint32_t m_index = 0;
        size_t m_sample_i = 0;
        size_t m_entry_i = 0;
        bool m_include_snp = false;
    };
};

#endif

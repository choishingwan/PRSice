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
              const bool keep_ambig = false, const bool is_ref = false,
              const bool intermediate = false);
    ~BinaryGen();

private:
    typedef std::vector<std::vector<double>> Data;
    std::unordered_map<std::string, genfile::bgen::Context> m_context_map;
    std::vector<genfile::byte_t> m_buffer1, m_buffer2;
    std::ifstream m_bgen_file;
    std::string m_cur_file;
    std::string m_intermediate_file;
    bool m_intermediate = false;
    bool m_target_plink = false;
    bool m_ref_plink = false;
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
        const uintptr_t unfiltered_sample_ct4 =
            (m_unfiltered_sample_ct + 3) / 4;
        // std::fill(m_tmp_genotype.begin(), m_tmp_genotype.end(), 0);
        // std::memset(m_tmp_genotype, 0x0, m_unfiltered_sample_ctl * 2 *
        // sizeof(uintptr_t));
        if (m_ref_plink) {
            // this is a plink file format, so we will directly read the file
            if (m_cur_file.empty() || file_name.compare(m_cur_file) != 0
                || !m_bgen_file.is_open())
            {
                if (m_bgen_file.is_open()) m_bgen_file.close();
                m_bgen_file.open(file_name.c_str(), std::ifstream::binary);
                if (!m_bgen_file.is_open()) {
                    std::string error_message =
                        "Error: Cannot open bgen file: " + file_name;
                    throw std::runtime_error(error_message);
                }
                m_cur_file = file_name;
            }

            m_bgen_file.seekg(byte_pos, std::ios_base::beg);
            m_bgen_file.read((char*) genotype, unfiltered_sample_ct4);
        }
        else if (load_and_collapse_incl(byte_pos, file_name,
                                        m_unfiltered_sample_ct, m_founder_ct,
                                        m_founder_info.data(), final_mask,
                                        false, m_tmp_genotype.data(), genotype))
        {
            throw std::runtime_error("Error: Cannot read the bgen file!");
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
                                    uintptr_t* __restrict mainbuf,
									bool intermediate=false)
    {
        assert(unfiltered_sample_ct);
        if(!intermediate){
        	if (m_cur_file.empty() || file_name.compare(m_cur_file) != 0
        			|| !m_bgen_file.is_open())
        	{
        		if (m_bgen_file.is_open()) m_bgen_file.close();
        		std::string bgen_name = file_name + ".bgen";
        		m_bgen_file.open(bgen_name.c_str(), std::ifstream::binary);
        		if (!m_bgen_file.is_open()) {
        			std::string error_message =
        					"Error: Cannot open bgen file: " + file_name;
        			throw std::runtime_error(error_message);
        		}
        		m_cur_file = file_name;
        	}
        	auto&& context = m_context_map[file_name];
        	m_bgen_file.seekg(byte_pos, std::ios_base::beg);
        	PLINK_generator setter(&m_sample_include, mainbuf, m_hard_threshold);
        	genfile::bgen::read_and_parse_genotype_data_block<PLINK_generator>(
        			m_bgen_file, context, setter, &m_buffer1, &m_buffer2, false);
        	// output from load_raw should have already copied all samples
        	// to the front without the need of subseting
        	if (do_reverse) {
        		reverse_loadbuf(sample_ct, (unsigned char*) mainbuf);
        	}
        }else{
            const uintptr_t unfiltered_sample_ct4 =
                (m_unfiltered_sample_ct + 3) / 4;
        	 if (m_cur_file.empty() || file_name.compare(m_cur_file) != 0
        			 || !m_bgen_file.is_open())
        	 {
        		 if (m_bgen_file.is_open()) m_bgen_file.close();
        		 m_bgen_file.open(file_name.c_str(), std::ifstream::binary);
        		 if (!m_bgen_file.is_open()) {
        			 std::string error_message =
        					 "Error: Cannot open bgen file: " + file_name;
        			 throw std::runtime_error(error_message);
        		 }
        		 m_cur_file = file_name;
        	 }

        	 m_bgen_file.seekg(byte_pos, std::ios_base::beg);
        	 m_bgen_file.read((char*) mainbuf, unfiltered_sample_ct4);
        }
        // mainbuf should contains the information
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
                sample_prs.num_snp = sample_prs.num_snp * (m_not_first || m_start_geno) + 1;
                sample_prs.prs = sample_prs.prs * (m_not_first || m_start_geno)
                                 + m_homcom_weight * value * m_stat * 0.5;
                break;
            case 1:
                sample_prs.num_snp = sample_prs.num_snp * (m_not_first || m_start_geno) + 1;
                sample_prs.prs = sample_prs.prs * (m_not_first || m_start_geno)
                                 + m_het_weight * value * m_stat * 0.5;
                break;
            case 2:
                sample_prs.num_snp = sample_prs.num_snp * (m_not_first || m_start_geno) + 1;
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
        bool m_is_missing  = false;
        bool m_start_geno = false;
    };


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
            m_index = 0;
            m_shift = 0;
            rs.clear();
        }
        void set_min_max_ploidy(uint32_t min_ploidy, uint32_t max_ploidy,
                                uint32_t min_entries, uint32_t max_entries)
        {
        }

        bool set_sample(std::size_t i)
        {
            m_sample_i = i;
            m_geno = 1;
            m_entry_i = 3;
            m_hard_prob = 0.0;
            m_exp_value = 0.0;
            m_include_sample = IS_SET(m_sample->data(), m_sample_i);
            return m_include_sample;
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
                // m_geno = 2 - m_entry_i;
                // 0 2 3 as 1 is missing
                // problem is the genotype should be
                // 3 2 0
                m_geno = (m_entry_i == 1) ? 0 : m_entry_i;
                m_hard_prob = value;
            }
            m_exp_value += m_geno * value;
            --m_entry_i;
        }
        // call if sample is missing
        void set_value(uint32_t, genfile::MissingValue value) {}

        void finalise() {}
        void sample_completed()
        {
            if (m_shift == 0)
                m_genotype[m_index] = 0; // match behaviour of binaryplink
            m_genotype[m_index] |= m_geno << m_shift;
            m_shift += 2;
            if (m_shift == BITCT) {
                m_index++;
                m_shift = 0;
            }
            rs.push(m_exp_value);
        }
        double info_score() const
        {
            double p = rs.mean() / 2.0;
            double p_all = 2.0 * p * (1.0 - p);
            return (rs.var() / p_all);
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
        size_t m_entry_i = 3;
        bool m_include_sample = false;
    };
};

#endif

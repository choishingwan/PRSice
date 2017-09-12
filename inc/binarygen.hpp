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
    BinaryGen(std::string prefix, std::string pheno_file, bool header,
              std::string remove_sample, std::string keep_sample,
              std::string extract_snp, std::string exclude_snp, bool ignore_fid,
              int num_auto = 22, bool no_x = false, bool no_y = false,
              bool no_xy = false, bool no_mt = false, bool keep_ambig = false,
              const size_t thread = 1, bool verbose = false);
    ~BinaryGen();

private:
    typedef std::vector<std::vector<double>> Data;
    Sample get_sample(std::vector<std::string>& token, bool ignore_fid,
                      bool has_sex, int sex_col, std::vector<int>& sex_info);
    std::vector<Sample> preload_samples(std::string pheno, bool has_header,
                                        bool ignore_fid);
    std::unordered_map<std::string, int> m_sample_index_check;

    std::vector<Sample> load_samples(bool ignore_fid);


    std::vector<SNP> load_snps();

    std::string m_cur_file;
    inline void load_raw(uintptr_t* genotype, const uint32_t snp_index,
                         const std::string& file_name)
    {
        if (m_cur_file.empty() || file_name.compare(m_cur_file) != 0) {
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
        m_bgen_file.seekg(snp_index, std::ios_base::beg);

        Data                         probability;
        ProbSetter                   setter(&probability);
        std::vector<genfile::byte_t> buffer1, buffer2;

        genfile::bgen::read_and_parse_genotype_data_block<ProbSetter>(
            m_bgen_file, m_bgen_info[file_name], setter, &buffer1, &buffer2);
        for (size_t i_sample = 0; i_sample < probability.size(); ++i_sample) {
            auto&& prob = probability[i_sample];
            if (prob.size() != 3) {
                // this is likely phased
                fprintf(stderr, "ERROR: Currently don't support phased data\n");
                fprintf(
                    stderr,
                    "       (It is because the lack of development time)\n");
                throw std::runtime_error("");
            }
            uintptr_t cur_geno = 1;
            for (size_t g = 0; g < prob.size(); ++g) {
                if (prob[g] >= filter.hard_threshold) {
                    cur_geno = (g == 2) ? 0 : 3 - g; // binary code for plink
                    break;
                }
            }
            // now genotype contain the genotype of this sample after filtering
            // need to bit shift here
            int shift = (i_sample % BITCT * 2);
            int index = (i_sample * 2) / BITCT;
            genotype[index] |= cur_geno << shift;
        }
    };

    inline void read_genotype(uintptr_t* genotype, const uint32_t snp_index,
                              const std::string& file_name)
    {
        // the bgen library seems over complicated
        // try to use PLINK one. The difficulty is ifstream vs FILE
        // this is for the LD calculation
        uintptr_t final_mask = get_final_mask(m_founder_ct);
        uintptr_t unfiltered_sample_ctl =
            BITCT_TO_WORDCT(m_unfiltered_sample_ct);
        std::fill(m_tmp_genotype, m_tmp_genotype + unfiltered_sample_ctl * 2,
                  0);
        // std::memset(m_tmp_genotype, 0x0, m_unfiltered_sample_ctl * 2 *
        // sizeof(uintptr_t));
        if (load_and_collapse_incl(snp_index, file_name, m_unfiltered_sample_ct,
                                   m_founder_ct, m_founder_info, final_mask,
                                   false, m_tmp_genotype, genotype))
        {
            throw std::runtime_error("ERROR: Cannot read the bed file!");
        }
    };

    // borrowed from plink
    uint32_t load_and_collapse_incl(const uint32_t     snp_index,
                                    const std::string& file_name,
                                    uint32_t           unfiltered_sample_ct,
                                    uint32_t           sample_ct,
                                    const uintptr_t* __restrict sample_include,
                                    uintptr_t final_mask, uint32_t do_reverse,
                                    uintptr_t* __restrict rawbuf,
                                    uintptr_t* __restrict mainbuf)
    {
        assert(unfiltered_sample_ct);
        if (unfiltered_sample_ct == sample_ct) {
            rawbuf = mainbuf;
        }
        load_raw(rawbuf, snp_index, file_name);

        if (unfiltered_sample_ct != sample_ct) {
            copy_quaterarr_nonempty_subset(rawbuf, sample_include,
                                           unfiltered_sample_ct, sample_ct,
                                           mainbuf);
        }
        else
        {
            mainbuf[(unfiltered_sample_ct - 1) / BITCT2] &= final_mask;
        }
        if (do_reverse) {
            reverse_loadbuf(sample_ct, (unsigned char*) mainbuf);
        }

        // mainbuf should contains the information
        return 0;
    }

    void read_score(misc::vec2d<Sample_lite>& current_prs_score,
                    size_t start_index, size_t end_bound);
    void hard_code_score(misc::vec2d<Sample_lite>& current_prs_score,
                         size_t start_index, size_t end_bound);
    void dosage_score(misc::vec2d<Sample_lite>& current_prs_score,
                      size_t start_index, size_t end_bound);
    std::unordered_map<std::string, genfile::bgen::Context> m_bgen_info;
    std::unordered_map<std::string, uint32_t>               m_offset_map;
    std::ifstream                                           m_bgen_file;
    /** DON'T TOUCH      */
    struct ProbSetter
    {
        ProbSetter(Data* result) : m_result(result), m_sample_i(0) {}

        // Called once allowing us to set storage.
        void initialise(std::size_t number_of_samples,
                        std::size_t number_of_alleles)
        {
            m_result->clear();
            m_result->resize(number_of_samples);
        }

        // If present with this signature, called once after initialise()
        // to set the minimum and maximum ploidy and numbers of probabilities
        // among samples in the data. This enables us to set up storage for the
        // data ahead of time.
        void set_min_max_ploidy(uint32_t min_ploidy, uint32_t max_ploidy,
                                uint32_t min_entries, uint32_t max_entries)
        {
            for (std::size_t i = 0; i < m_result->size(); ++i) {
                m_result->at(i).reserve(max_entries);
            }
        }

        // Called once per sample to determine whether we want data for this
        // sample
        bool set_sample(std::size_t i)
        {
            m_sample_i = i;
            // Yes, here we want info for all samples.
            return true;
        }

        // Called once per sample to set the number of probabilities that are
        // present.
        void set_number_of_entries(std::size_t        ploidy,
                                   std::size_t        number_of_entries,
                                   genfile::OrderType order_type,
                                   genfile::ValueType value_type)
        {
            assert(value_type == genfile::eProbability);
            m_result->at(m_sample_i).resize(number_of_entries);
            m_entry_i = 0;
        }

        // Called once for each genotype (or haplotype) probability per sample.
        void set_value(uint32_t, double value)
        {
            m_result->at(m_sample_i).at(m_entry_i++) = value;
        }

        // Ditto, but called if data is missing for this sample.
        void set_value(uint32_t, genfile::MissingValue value)
        {
            // Here we encode missing probabilities with -1
            m_result->at(m_sample_i).at(m_entry_i++) = -1;
        }

        // If present with this signature, called once after all data has been
        // set.
        void finalise()
        {
            // nothing to do in this implementation.
        }

    private:
        Data*       m_result;
        std::size_t m_sample_i;
        std::size_t m_entry_i;
    };
    void read_genotype_data_block(std::istream&                 aStream,
                                  genfile::bgen::Context const& context,
                                  std::vector<genfile::byte_t>* buffer)
    {
        uint32_t payload_size = 0;
        if ((context.flags & genfile::bgen::e_Layout)
                == genfile::bgen::e_Layout2
            || ((context.flags & genfile::bgen::e_CompressedSNPBlocks)
                != genfile::bgen::e_NoCompression))
        {
            genfile::bgen::read_little_endian_integer(aStream, &payload_size);
        }
        else
        {
            payload_size = 6 * context.number_of_samples;
        }
        buffer->resize(payload_size);
        aStream.read(reinterpret_cast<char*>(&(*buffer)[0]), payload_size);
    }
};

#endif


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

#include "genotype.hpp"

class BinaryPlink : public Genotype
{
public:
    BinaryPlink(std::string prefix, std::string remove_sample,
                std::string keep_sample, std::string extract_snp,
                std::string exclude_snp, std::string fam_name,
                std::string log_file, bool ignore_fid, bool nonfounder,
                int num_auto = 22, bool no_x = false, bool no_y = false,
                bool no_xy = false, bool no_mt = false, bool keep_ambig = false,
                const size_t thread = 1, bool verbose = false);
    ~BinaryPlink();

private:
    uintptr_t m_bed_offset = 3;
    std::vector<Sample> load_samples(bool ignore_fid);
    std::vector<SNP> load_snps();
    std::vector<size_t> m_num_snp_per_file; // for bed file size check
    std::string m_fam_name = "";

    void check_bed();

    inline void read_genotype(uintptr_t* genotype, const SNP& snp,
                              const std::string& file_name)
    {
        uintptr_t final_mask = get_final_mask(m_founder_ct);
        uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
        size_t snp_index = snp.snp_id();
        bool jump = !(snp_index - m_prev_index == 1);
        if (m_cur_file.empty() || m_cur_file.compare(file_name) != 0) {
            if (m_bed_file.is_open()) {
                m_bed_file.close();
            }
            std::string bedname = file_name + ".bed";
            m_bed_file.open(bedname.c_str(), std::ios::binary);
            jump = true;
        }
        // don't do jumping unless we have to
        if (jump) {
            if (!m_bed_file.seekg(
                    m_bed_offset
                        + (snp_index * ((uint64_t) unfiltered_sample_ct4)),
                    std::ios_base::beg))
            {
                throw std::runtime_error("ERROR: Cannot read the bed file!");
            }
        }
        m_prev_index = snp_index;
        // std::fill(m_tmp_genotype.begin(), m_tmp_genotype.end(), 0);
        // std::memset(m_tmp_genotype, 0x0, m_unfiltered_sample_ctl * 2 *
        // sizeof(uintptr_t));
        // this is for LD calculation, thus doesn't really need to worry about
        // the founder_count otherwise we need to use the sample_ct instead This
        // is also ok for the LD reference file as the only way that will affect
        // the sample inclusion / exclusion is already handled during file
        // loading
        if (load_and_collapse_incl(m_unfiltered_sample_ct, m_founder_ct,
                                   m_founder_info.data(), final_mask, false,
                                   m_bed_file, m_tmp_genotype.data(), genotype))
        {
            throw std::runtime_error("ERROR: Cannot read the bed file!");
        }
    };

    void read_score(misc::vec2d<Sample_lite>& current_prs_score,
                    size_t start_index, size_t end_bound);

    std::ifstream m_bed_file;
    std::string m_cur_file;
    int m_prev_index = -2;


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
        if (unfiltered_sample_ct == sample_ct) {
            rawbuf = mainbuf;
        }
        if (!m_bed_file.read((char*) rawbuf, unfiltered_sample_ct4)) {
            return RET_READ_FAIL;
        }
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
};

#endif

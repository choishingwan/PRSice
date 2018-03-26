
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

class BinaryPlink : public Genotype
{
public:
    BinaryPlink(const std::string& prefix, const std::string& sample_file,
                const size_t thread = 1, const bool ignore_fid = false,
                const bool keep_nonfounder = false,
                const bool keep_ambig = false);
    ~BinaryPlink();

private:
    uintptr_t m_bed_offset = 3;
    std::vector<Sample> gen_sample_vector();

    std::vector<SNP>
    gen_snp_vector(const double geno, const double maf, const double info,
                   const double hard_threshold, const bool hard_coded,
                   const std::string& out_prefix, Genotype* target = nullptr);

    void check_bed(const std::string& bed_name, size_t num_marker);

    // this is for ld calculation only
    inline void read_genotype(uintptr_t* genotype,
                              const std::streampos byte_pos,
                              const std::string& file_name)
    {
        uintptr_t final_mask = get_final_mask(m_founder_ct);
        uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
        std::streampos snp_index = byte_pos;
        if (m_cur_file.empty() || m_cur_file.compare(file_name) != 0) {
            if (m_bed_file.is_open()) {
                m_bed_file.close();
            }
            std::string bedname = file_name + ".bed";
            m_bed_file.open(bedname.c_str(), std::ios::binary);
            m_prev_loc = 0;
            m_cur_file = file_name;
        }
        if ((m_prev_loc != snp_index)
            && !m_bed_file.seekg(snp_index, std::ios_base::beg))
        {
            throw std::runtime_error("ERROR: Cannot read the bed file!");
        }
        // so that we don't jump if we don't need to
        m_prev_loc = snp_index + (std::streampos) unfiltered_sample_ct4;
        if (load_and_collapse_incl(m_unfiltered_sample_ct, m_founder_ct,
                                   m_founder_info.data(), final_mask, false,
                                   m_bed_file, m_tmp_genotype.data(), genotype))
        {
            throw std::runtime_error("ERROR: Cannot read the bed file!");
        }
    };

    void read_score(size_t start_index, size_t end_bound,
                    const size_t region_index);

    std::ifstream m_bed_file;
    std::string m_cur_file;
    std::streampos m_prev_loc = 0;

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
        if (!bedfile.read((char*) rawbuf, unfiltered_sample_ct4)) {
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

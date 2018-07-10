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

#include "snp.hpp"

SNP::~SNP() {}

std::vector<size_t> SNP::sort_by_p_chr(const std::vector<SNP>& input)
{
    std::vector<size_t> idx(input.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&input](size_t i1, size_t i2) {
        // plink do it with respect to the location instead of statistic
        // chr first such that SNPs within the same chromosome will
        // be processed together
        if (input[i1].m_chr == input[i2].m_chr) {
            if (input[i1].m_p_value == input[i2].m_p_value) {
                if (input[i1].m_loc == input[i2].m_loc) {
                    // assume it is impossible for non-duplicated SNPs to have
                    // same everything
                    return std::abs(input[i1].m_stat) > fabs(input[2].m_stat);
                }
                else
                    return input[i1].m_loc < input[i2].m_loc;
            }
            else
                return input[i1].m_p_value < input[i2].m_p_value;
        }
        else
            return input[i1].m_chr < input[i2].m_chr;
    });
    return idx;
}

void SNP::sort_snp_for_perm(std::vector<size_t>& index,
                            const std::vector<SNP>& input)
{
    // now index is sorted, we want to sort the top select_size entry of index
    // by the file info
    std::sort(index.begin(), index.end(), [&input](size_t i1, size_t i2) {
        if (input[i1].m_target_file == input[i2].m_target_file) {
            return input[i1].m_target_byte_pos < input[i2].m_target_byte_pos;
        }
        else
            return (input[i1].m_target_file == input[i2].m_target_file);
    });
}

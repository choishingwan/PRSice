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

#include "snp.hpp"

SNP::~SNP() {}

std::vector<size_t> SNP::sort_by_p_chr(const std::vector<SNP>& input)
{
    std::vector<size_t> idx(input.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&input](size_t i1, size_t i2) {
        // plink do it w.r.t the name of the RS ID (ignoring the string part)
        // which is slightly too complicated for us. Will simply use location
        // instead

        // chr first such that SNPs within the same chromosome will be
        // processed together
        if (input[i1].m_chr == input[i2].m_chr)
        {
            if (misc::logically_equal(input[i1].m_p_value, input[i2].m_p_value))
            {
                // in theory, we can also add in the stat and se,
                // but as they are double, there might be problem
                // (have tried to use stat and that cause seg fault)
                if (input[i1].m_loc == input[i2].m_loc)
                    return input[i1].m_rs < input[i2].m_rs;
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

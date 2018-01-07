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

SNP::SNP()
{
    basic.chr = -1;
    basic.loc = -1; // default is -1 to indicate that it is not provided
    basic.valid = true;
    statistic.flipped = false;
    statistic.stat = 0.0;
    statistic.se = 0.0;
    statistic.p_value = 2.0; // set this to 2 such that only SNPs in base file
                             // have valid P-value
    threshold.p_threshold = 0.0;
    threshold.category = 0;
    clump_info.clumped = false;
}


SNP::SNP(const std::string& rs_id, const int chr, const int loc,
         const std::string& ref_allele, const std::string& alt_allele,
         const std::string& file_name, const std::streampos byte_pos)
{
    basic.ref = ref_allele;
    basic.alt = alt_allele;
    basic.rs = rs_id;
    basic.chr = chr;
    basic.loc = loc;
    basic.valid = true;
    statistic.se = 0.0;
    statistic.p_value = 0.0;
    statistic.stat = 0.0;
    statistic.flipped = false;
    threshold.p_threshold = 0.0;
    threshold.category = 0;
    clump_info.clumped = false;
    file_info.file = file_name;
    file_info.byte_pos = byte_pos;
}
SNP::~SNP() {}

std::vector<size_t> SNP::sort_by_p(const std::vector<SNP>& input)
{
    std::vector<size_t> idx(input.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&input](size_t i1, size_t i2) {
        // plink do it with respect to the location instead of statistic
        if (input[i1].statistic.p_value == input[i2].statistic.p_value) {
            if (input[i1].basic.chr == input[i2].basic.chr) {
                if (input[i1].basic.loc == input[i2].basic.loc) {
                    if (fabs(input[i1].statistic.stat)
                        == fabs(input[i2].statistic.stat))
                    {
                        return input[i1].statistic.se < input[i2].statistic.se;
                    }
                    else
                        return fabs(input[i1].statistic.stat)
                               > fabs(input[2].statistic.stat);
                }
                else
                    return input[i1].basic.loc < input[i2].basic.loc;
            }
            else
                return input[i1].basic.chr < input[i2].basic.chr;
        }
        else
            return input[i1].statistic.p_value < input[i2].statistic.p_value;
    });
    return idx;
}

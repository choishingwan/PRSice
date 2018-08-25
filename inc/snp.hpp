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

#ifndef SNP_H
#define SNP_H

#include "commander.hpp"
#include "misc.hpp"
#include "plink_common.hpp"
#include "region.hpp"
#include "storage.hpp"
#include <algorithm>
#include <limits.h>
#include <numeric>
#include <stdexcept>
#include <string>

class Genotype;
class SNP
{
public:
    SNP(){};
    SNP(const std::string& rs_id, const intptr_t chr, const intptr_t loc,
        const std::string& ref_allele, const std::string& alt_allele,
        const std::string& file_name, const std::streampos byte_pos,
        const uint32_t homcom_ct, const uint32_t het_ct,
        const uint32_t homrar_ct, const uint32_t missing)
        : m_alt(alt_allele)
        , m_ref(ref_allele)
        , m_rs(rs_id)
        , m_target_file(file_name)
        , m_ref_file(file_name)
        , m_target_byte_pos(byte_pos)
        , m_ref_byte_pos(byte_pos)
        , m_chr(chr)
        , m_loc(loc)
        , m_homcom(homcom_ct)
        , m_het(het_ct)
        , m_homrar(homrar_ct)
        , m_missing(missing)
    {
        m_has_count = true;
    };
    SNP(const std::string& rs_id, const intptr_t chr, const intptr_t loc,
        const std::string& ref_allele, const std::string& alt_allele,
        const std::string& file_name, const std::streampos byte_pos)
        : m_alt(alt_allele)
        , m_ref(ref_allele)
        , m_rs(rs_id)
        , m_target_file(file_name)
        , m_ref_file(file_name)
        , m_target_byte_pos(byte_pos)
        , m_ref_byte_pos(byte_pos)
        , m_chr(chr)
        , m_loc(loc)
    {
        m_has_count = false;
    };
    virtual ~SNP();

    void set_statistic(const double stat, const double p_value,
                       const intptr_t category, const double p_threshold)
    {
        m_stat = stat;
        m_p_value = p_value;
        m_category = category;
        m_p_threshold = p_threshold;
    };
    void add_reference(const std::string& ref_file,
                       const std::streampos ref_byte_pos)
    {
        m_ref_file = ref_file;
        m_ref_byte_pos = ref_byte_pos;
    }

    void update_target(const std::string& ref_file,
                       const std::streampos ref_byte_pos)
    {
        m_target_file = ref_file;
        m_target_byte_pos = ref_byte_pos;
    }
    inline void set_flipped() { m_flipped = true; };
    std::string get_rs() const { return m_rs; };
    /*!
     * \brief Function to sort a vector of SNP by their chr then by their
     * p-value
     * \param input is the vector containing the SNPs
     * \return return a vector containing index to the sort order of the input
     */
    static std::vector<size_t> sort_by_p_chr(const std::vector<SNP>& input);
    static void sort_snp_for_perm(std::vector<size_t>& index,
                                  const std::vector<SNP>& input);
    static bool compare_snp(const SNP& a, const SNP& b)
    {
        if (a.m_chr == b.m_chr) {
            if (a.m_loc == b.m_loc) {
                if (a.m_target_file == b.m_target_file) {
                    return a.m_target_byte_pos < b.m_target_byte_pos;
                }
                return a.m_target_file < b.m_target_file;
            }
            return a.m_loc < b.m_loc;
        }
        return a.m_chr < b.m_chr;
    }

    inline bool matching(intptr_t chr, intptr_t loc, std::string& ref,
                         std::string& alt, bool& flipped)
    {
        // should be trimmed
        if (chr != -1 && m_chr != -1 && chr != m_chr) {
            return false;
        }
        if (loc != -1 && m_loc != -1 && loc != m_loc) {
            return false;
        }
        if (m_ref == ref) {
            if (!m_alt.empty() && !alt.empty()) {
                return (m_alt == alt);
            }
            else
                return true;
        }
        else if (complement(m_ref) == ref)
        {
            if (!m_alt.empty() && !alt.empty()) {
                return (complement(m_alt) == alt);
            }
            else
                return true;
        }
        else if (!m_alt.empty() && !alt.empty())
        {
            if ((m_ref == alt) && (m_alt == ref)) {
                flipped = true;
                return true;
            }
            if ((complement(m_ref) == alt) && (complement(m_alt) == ref)) {
                flipped = true;
                return true;
            }
            return false;
        }
        else
            return false; // cannot flip nor match
    };

    intptr_t chr() const { return m_chr; };
    intptr_t loc() const { return m_loc; };
    intptr_t category() const { return m_category; };
    double p_value() const { return m_p_value; };
    double stat() const { return m_stat; };
    double get_threshold() const { return m_p_threshold; };
    std::streampos byte_pos() const { return m_target_byte_pos; };
    std::streampos ref_byte_pos() const { return m_ref_byte_pos; };
    std::string file_name() const { return m_target_file; };
    std::string ref_file_name() const { return m_ref_file; };
    std::string rs() const { return m_rs; };
    std::string ref() const { return m_ref; };
    std::string alt() const { return m_alt; };
    bool is_flipped() { return m_flipped; };

    inline bool in(size_t i) const
    {
        if (i / BITCT >= m_max_flag_index)
            throw std::out_of_range("Out of range for flag");
        return ((m_flags[i / BITCT] >> (i % BITCT)) & 1);
    }
    void set_flag(Region& region)
    {
        m_max_flag_index = BITCT_TO_WORDCT(region.size());
        m_flags.resize(m_max_flag_index);
        region.update_flag(m_chr, m_rs, m_loc, m_flags);
    };

    void set_clumped() { m_clumped = true; };
    void clump(SNP& target, double r2, bool use_proxy, double proxy = 2)
    {
        if (target.clumped()) return;
        bool completed = false;
        // if we want to use proxy, and that our r2 is higher than
        // the proxy threshold, we will do clumping
        if (use_proxy && r2 > proxy) {
            // proxy clump
            for (size_t i_flag = 0; i_flag < m_max_flag_index; ++i_flag) {
                // two become one
                m_flags[i_flag] |= target.m_flags[i_flag];
            }
            completed = true;
        }
        else
        {
            // otherwise, we will just do noraml clumping
            for (size_t i_flag = 0; i_flag < m_max_flag_index; ++i_flag) {
                target.m_flags[i_flag] =
                    target.m_flags[i_flag]
                    ^ (m_flags[i_flag] & target.m_flags[i_flag]);
                completed = (target.m_flags[i_flag] == 0);
            }
        }
        if (completed) {
            target.set_clumped();
            target.m_remove = true;
        }
        m_clumped = true;
        // protect from other SNPs tempering its flags
    }

    bool remove() const { return m_remove; };
    bool clumped() const { return m_clumped; };
    bool valid() const { return m_valid; };
    void invalidate() { m_valid = false; };
    /*!
     * \brief Set the lower boundary (index of m_existed_snp) of this SNP if it
     * is used as the index
     * \param low the designated bound index
     */
    void set_low_bound(intptr_t low) { m_low_bound = low; }
    /*!
     * \brief Set the upper boundary (index of m_existed_snp) of this SNP if it
     * is used as the index
     * \param up the designated bound index
     */
    void set_up_bound(intptr_t up) { m_up_bound = up; }
    bool get_counts(uint32_t& homcom, uint32_t& het, uint32_t& homrar,
                    uint32_t& missing)
    {
        homcom = m_homcom;
        het = m_het;
        homrar = m_homrar;
        missing = m_missing;
        return m_has_count;
    }
    void set_counts(uint32_t& homcom, uint32_t& het, uint32_t& homrar,
                    uint32_t& missing)
    {
        m_homcom = homcom;
        m_het = het;
        m_homrar = homrar;
        m_missing = missing;
        m_has_count = true;
    }
    intptr_t up_bound() const { return m_up_bound; };
    intptr_t low_bound() const { return m_low_bound; };

private:
    // basic info
    // actually, the packing of the data is problematic and to enhance
    // performance we might want to organize the data into way where it is
    // easier to "cache" also use data types that are more friendly?
    std::vector<uintptr_t> m_flags;
    std::string m_alt;
    std::string m_ref;
    std::string m_rs;
    std::string m_target_file;
    std::string m_ref_file;
    std::streampos m_target_byte_pos;
    std::streampos m_ref_byte_pos;
    double m_stat = 0.0;
    double m_p_value = 2.0;
    double m_p_threshold = 0;
    intptr_t m_chr = -1;
    intptr_t m_category = -1;
    intptr_t m_loc = -1;
    intptr_t m_low_bound = 0;
    intptr_t m_up_bound = 0;
    uint32_t m_homcom = 0;
    uint32_t m_het = 0;
    uint32_t m_homrar = 0;
    uint32_t m_missing = 0;
    bool m_has_count = false;
    bool m_clumped = false;
    bool m_valid = true;
    bool m_remove = false;
    bool m_flipped = false;
    // This indicate where this SNP's bound is at
    // useful for PRSlice and also clumping
    // thinking about it. Even if the location isn't given for
    // PRSet or PRSlice, we can still use the coordinates from
    // the target / reference file
    // the bound is [ )
    // prset related
    size_t m_max_flag_index = 0;

    inline std::string complement(const std::string& allele) const
    {
        if (allele.compare("A") == 0 || allele.compare("a") == 0) return "T";
        if (allele.compare("T") == 0 || allele.compare("t") == 0) return "A";
        if (allele.compare("G") == 0 || allele.compare("g") == 0) return "C";
        if (allele.compare("C") == 0 || allele.compare("c") == 0)
            return "G";
        else
            return allele; // Cannot flip, so will just return it as is
    }
};

#endif // SNP_H

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


#ifndef SNP_H
#define SNP_H

#include "commander.hpp"
#include "misc.hpp"
#include "plink_common.hpp"
#include "storage.hpp"
#include <algorithm>
#include <limits.h>
#include <numeric>
#include <stdexcept>
#include <string>
static_assert(
    sizeof(std::streamsize) <= sizeof(unsigned long long),
    "streampos larger than unsigned long long, don't know how to proceed. "
    "Please use PRSice on another machine");
class Genotype;
class SNP
{
public:
    SNP() {}
    SNP(const std::string& rs_id, const size_t chr, const size_t loc,
        const std::string& ref_allele, const std::string& alt_allele,
        const double& stat, const double& p_value,
        const unsigned long long category, const double p_threshold)
        : m_alt(alt_allele)
        , m_ref(ref_allele)
        , m_rs(rs_id)
        , m_stat(stat)
        , m_p_value(p_value)
        , m_p_threshold(p_threshold)
        , m_chr(chr)
        , m_loc(loc)
        , m_category(category)
    {
    }
    SNP(const std::string& rs_id, const size_t chr, const size_t loc,
        const std::string& ref_allele, const std::string& alt_allele,
        const size_t& idx, const std::streampos byte_pos)
        : m_alt(alt_allele)
        , m_ref(ref_allele)
        , m_rs(rs_id)
        , m_chr(chr)
        , m_loc(loc)
    {
        m_target.name_idx = idx;
        m_target.byte_pos = byte_pos;
        m_reference.name_idx = idx;
        m_reference.byte_pos = byte_pos;
    }

    virtual ~SNP();

    void update_file(const size_t& idx, const std::streampos byte_pos,
                     const bool is_ref)
    {
        auto&& target = is_ref ? m_reference : m_target;
        target.name_idx = idx;
        target.byte_pos = byte_pos;
    }
    void add_snp_info(const SNP& src, const bool flipping, const bool is_ref)
    {
        if (!is_ref)
        {
            m_target.name_idx = src.get_file_idx();
            m_target.byte_pos = src.get_byte_pos();
            m_chr = src.chr();
            m_loc = src.loc();
            m_ref = src.ref();
            m_alt = src.alt();
            m_flipped = flipping;
        }
        else
        {
            m_ref_flipped = flipping;
        }
        m_reference.name_idx = src.get_file_idx();
        m_reference.byte_pos = src.get_byte_pos();
    }

    /*!
     * \brief Function to sort a vector of SNP by their chr then by their
     * p-value
     * \param input is the vector containing the SNPs
     * \return return a vector containing index to the sort order of the input
     */
    static std::vector<size_t> sort_by_p_chr(const std::vector<SNP>& input);

    /*!
     * \brief Compare the current SNP with another SNP
     * \param chr is the chromosome encoding of the other SNP
     * \param loc is the coordinate of the other SNP
     * \param ref is the reference allele of the other SNP
     * \param alt is the alternative allele of the other SNP
     * \param flipped is used as a return value. If flipping is required,
     * flipped = true
     * \return true if it is a match
     */
    inline bool matching(const SNP& i, bool& flipped)
    {
        return matching(i.chr(), i.loc(), i.ref(), i.alt(), flipped);
    }
    inline bool matching(const size_t chr, const size_t loc,
                         const std::string& ref, const std::string& alt,
                         bool& flipped)
    {
        // should be trimmed
        if (chr != ~size_t(0) && m_chr != ~size_t(0) && chr != m_chr)
        { return false; }
        if (loc != ~size_t(0) && m_loc != ~size_t(0) && loc != m_loc)
        { return false; }
        flipped = false;
        if (m_ref == ref)
        {
            if (!m_alt.empty() && !alt.empty()) { return (m_alt == alt); }
            else
                return true;
        }
        else if (complement(m_ref) == ref)
        {
            if (!m_alt.empty() && !alt.empty())
            { return (complement(m_alt) == alt); }
            else
                return true;
        }
        else if (!alt.empty())
        {
            // here, we already know the refs don't match so we want it to match
            // with alt
            if ((m_ref == alt) && (m_alt.empty() || m_alt == ref))
            {
                flipped = true;
                return true;
            }
            if ((complement(m_ref) == alt)
                && (m_alt.empty() || complement(m_alt) == ref))
            {
                flipped = true;
                return true;
            }
            return false;
        }
        else
            return false; // cannot flip nor match
    }

    size_t chr() const { return m_chr; }
    size_t loc() const { return m_loc; }
    unsigned long long category() const { return m_category; }

    void set_category(unsigned long long& cur_category, double& cur_p_start,
                      const double& upper, const double& inter, bool& warning)
    {
        warning = false;
        if (m_p_value <= cur_p_start + inter)
        { // do nothing
        }
        else if (m_p_value > upper)
        {
            if (!misc::logically_equal(cur_p_start, upper))
            {
                cur_p_start = upper;
                ++cur_category;
            }
        }
        else
        {
            // this is a new threshold
            ++cur_category;
            // there will be imprecision w.r.t new
            if ((m_p_value - cur_p_start) / inter
                > std::numeric_limits<unsigned long long>::max())
            { warning = true; }
            // use log to help with the numeric stability
            double interval =
                std::log(m_p_value - cur_p_start) - std::log(inter);
            interval = std::floor(std::exp(interval));
            cur_p_start += std::exp(std::log(interval) + std::log(inter));
        }
        m_category = cur_category;
        m_p_threshold = cur_p_start;
        return;
    }
    /*!
     * \brief Get the p-value of the SNP
     * \return the p-value of the SNP
     */
    double p_value() const { return m_p_value; }
    /*!
     * \brief Get the effect size of the SNP
     * \return the effect size of the SNP
     */
    double stat() const { return m_stat; }

    /*!
     * \brief Return the p-value threshold of which this SNP falls into
     * \return  the p-value threshold
     */
    double get_threshold() const { return m_p_threshold; }
    void get_file_info(size_t& idx, std::streampos& byte_pos,
                       bool is_ref = false) const
    {
        auto&& from = is_ref ? m_reference : m_target;
        idx = from.name_idx;
        byte_pos = from.byte_pos;
    }
    size_t get_file_idx(bool is_ref = false) const
    {
        return is_ref ? m_reference.name_idx : m_target.name_idx;
    }
    std::streampos get_byte_pos(bool is_ref = false) const
    {
        return is_ref ? m_reference.byte_pos : m_target.byte_pos;
    }

    std::string rs() const { return m_rs; }
    std::string& rs() { return m_rs; }
    std::string ref() const { return m_ref; }
    std::string& ref() { return m_ref; }
    std::string alt() const { return m_alt; }
    std::string& alt() { return m_alt; }
    bool is_flipped() const { return m_flipped; }
    bool is_ref_flipped() const { return m_ref_flipped; }

    /*!
     * \brief check if this SNP is within the i th region
     * \param i is the index of the region
     * \return true if this SNP falls within the i th region
     */
    inline bool in(size_t i) const
    {
        if (i / BITCT >= m_clump_info.flags.size())
            throw std::out_of_range("Out of range for flag");
        return (IS_SET(m_clump_info.flags.data(), i));
    }
    std::vector<uintptr_t>& get_flag() { return m_clump_info.flags; }

    /*!
     * \brief Set the SNP to be clumped such that it will no longer be
     * considered in clumping
     */
    void set_clumped() { m_clump_info.clumped = true; }

    /*!
     * \brief This is the clumping algorithm. The current SNP will remove
     * another SNP if their R2 is higher than a threshold
     * \param target is the target SNP
     * \param r2 is the observed R2
     * \param use_proxy indicate if we want to perform proxy clump
     * \param proxy is the threshold for proxy clumping
     */
    void clump(SNP& target, double r2, bool use_proxy, double proxy = 2)
    {
        // if the target is already clumped, we will do nothing
        if (target.clumped()) return;
        // we need to check if the target SNP is completely clumped (e.g. no
        // longer representing any set)
        bool target_clumped = true;
        // if we want to use proxy, and that our r2 is higher than
        // the proxy threshold, we will do the proxy clumping
        // and the index SNP will get all membership (or) from the clumped
        if (use_proxy && r2 > proxy)
        {
            for (size_t i_flag = 0; i_flag < target.m_clump_info.flags.size();
                 ++i_flag)
            {
                m_clump_info.flags[i_flag] |= target.m_clump_info.flags[i_flag];
            }
            target_clumped = true;
        }
        else
        {
            for (size_t i_flag = 0; i_flag < target.m_clump_info.flags.size();
                 ++i_flag)
            {
                // For normal clumping, we will remove set identity from the
                // target SNP whenever both SNPs are within the same set.
                // i.e. if flag of SNP A (current) is 11011 and SNP B (target)
                // is 11110, by the end of clumping, it will become SNP A
                // =11011, SNP B = 00100
                // bit operation meaning:
                // ~m_clump_info = not in index
                // target.m_clump_info & ~m_clump_info = retain bit that are not
                // found in index
                target.m_clump_info.flags[i_flag] =
                    target.m_clump_info.flags[i_flag]
                    & ~m_clump_info.flags[i_flag];
                // if all flags of the target SNP == 0, it means that it no
                // longer represent any gene set and is consided as "clumped"
                target_clumped &= (target.m_clump_info.flags[i_flag] == 0);
            }
        }
        if (target_clumped) { target.set_clumped(); }
    }

    /*!
     * \brief Indicate if this snp is clumped
     * \return  Return true if this is clumped
     */
    bool clumped() const { return m_clump_info.clumped; }
    /*!
     * \brief Set the lower boundary (index of m_existed_snp) of this SNP if it
     * is used as the index
     * \param low the designated bound index
     */
    void set_low_bound(size_t low) { m_clump_info.low_bound = low; }
    /*!
     * \brief Set the upper boundary (index of m_existed_snp) of this SNP if it
     * is used as the index
     * \param up the designated bound index
     */
    void set_up_bound(size_t up) { m_clump_info.up_bound = up; }
    /*!
     * \brief get_counts will return the current genotype count for this SNP.
     * Return true if this was previously calculated (and indicate the need of
     * calculation)
     *
     * \param homcom is the count of homozygous common allele
     * \param het is the count of heterozygous
     * \param homrar is the count of homozygous rare allele
     * \param missing is the number of missing genotypes
     * \return true if calculation is already done
     */
    bool get_counts(uint32_t& homcom, uint32_t& het, uint32_t& homrar,
                    uint32_t& missing, const bool use_ref_maf) const
    {
        auto&& from = use_ref_maf ? m_ref_count : m_target_count;
        homcom = from.homcom;
        het = from.het;
        homrar = from.homrar;
        missing = from.missing;
        return from.has_count;
    }
    /*!
     * \brief This function will set the genotype count for the current SNP, and
     * will set the has_count to true
     *
     * \param homcom is the count of homozygous common allele
     * \param het is the count of heterozygous
     * \param homrar is the count of homozygous rare allele
     * \param missing is the number of missing genotypes
     */
    void set_counts(uint32_t homcom, uint32_t het, uint32_t homrar,
                    uint32_t missing, bool is_ref)
    {
        auto&& target = is_ref ? m_ref_count : m_target_count;
        if (m_ref_flipped && is_ref) std::swap(homcom, homrar);
        target.homcom = homcom;
        target.het = het;
        target.homrar = homrar;
        target.missing = missing;
        target.has_count = true;
    }


    /*!
     * \brief Obtain the upper bound of the clump region correspond to this SNP
     * \return the upper bound of the region
     */
    size_t up_bound() const { return m_clump_info.up_bound; }
    /*!
     * \brief Obtain the lower bound of the clump region correspond to this SNP
     * \return the lower bound of the region
     */
    size_t low_bound() const { return m_clump_info.low_bound; }
    void set_expected(double expected, bool is_ref = false)
    {
        if (is_ref)
            m_ref_expected_value = expected;
        else
            m_expected_value = expected;
    }
    bool has_expected() const { return m_has_expected; }
    bool has_ref_expected() const { return m_has_ref_expected; }
    double get_expected(bool use_ref_maf) const
    {
        if (use_ref_maf) return m_ref_expected_value;
        return m_expected_value;
    }
    void invalid() { m_is_valid = false; }
    bool stored_genotype() const { return !m_genotype.empty(); }
    void assign_genotype(const std::vector<uintptr_t>& genotype)
    {
        m_genotype = genotype;
    }
    std::vector<uintptr_t> get_genotype() const { return m_genotype; }
    static std::string complement(const std::string& allele)
    {
        // assume capitalized
        if (allele == "A") return "T";
        if (allele == "T") return "A";
        if (allele == "G") return "C";
        if (allele == "C")
            return "G";
        else
            return allele; // Cannot flip, so will just return it as is
    }
    std::vector<uintptr_t>& mod_genotype() { return m_genotype; }

private:
    /*static std::string g_separator;
    static bool g_use_chr;
    */
    AlleleCounts m_ref_count;
    AlleleCounts m_target_count;
    FileInfo m_target;
    FileInfo m_reference;
    SNPClump m_clump_info;
    std::vector<uintptr_t> m_genotype;
    std::string m_alt;
    std::string m_ref;
    std::string m_rs;
    double m_stat = 0.0;
    double m_p_value = 2.0;
    double m_p_threshold = 0;
    double m_expected_value = 0.0;
    double m_ref_expected_value = 0.0;
    size_t m_chr = ~size_t(0);
    size_t m_loc = ~size_t(0);
    unsigned long long m_category = 0;
    bool m_has_expected = false;
    bool m_has_ref_expected = false;
    bool m_flipped = false;
    bool m_ref_flipped = false;
    bool m_is_valid = true;
};

#endif // SNP_H

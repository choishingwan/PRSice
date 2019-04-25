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
    SNP() {}
    SNP(const std::string& rs_id, const int chr, const int loc,
        const std::string& ref_allele, const std::string& alt_allele,
        const std::string& file_name, const std::streampos byte_pos,
        const size_t homcom_ct, const size_t het_ct, const size_t homrar_ct,
        const size_t missing)
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
    }
    SNP(const std::string& rs_id, const int chr, const int loc,
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
    }

    SNP(const std::string& rs_id, const int chr, const int loc,
        const std::string& ref_allele, const std::string& alt_allele,
        const double& stat, const double& p_value, const int category,
        const double p_threshold)
        : m_alt(alt_allele)
        , m_ref(ref_allele)
        , m_rs(rs_id)
        , m_chr(chr)
        , m_loc(loc)
    {
        m_has_count = false;
        m_stat = stat;
        m_p_value = p_value;
        assert(category < 0);
        m_category = category;
        m_p_threshold = p_threshold;
    }

    virtual ~SNP();
    /*!
     * \brief Add the statistic information for this SNP
     * \param stat is the effect size
     * \param se is the standard error of the effect size
     * \param p_value is the p-value
     * \param category is the category of this SNP
     * \param p_threshold is the p-value threshold this SNP fall into
     */
    void set_statistic(const double& stat, const double& p_value,
                       const double& se, const double& maf, const int category,
                       const double p_threshold)
    {
        m_stat = stat;
        m_p_value = p_value;
        // by our algorithm, we should always have category bigger than or equal
        // to 0
        assert(category < 0);
        m_category = category;
        m_p_threshold = p_threshold;
        m_standard_error = se;
        m_maf = maf;
    }
    /*!
     * \brief This is to change the m_ref_file and m_ref_byte_pos to account for
     * reference panel. Without reference panel, ref_file and byte_pos equals to
     * those observed in target
     * \param ref_file is the file name of the reference panel that contain this
     * SNP
     * \param ref_byte_pos is the location of this SNP on the reference panel
     * file
     */
    void add_reference(const std::string& ref_file,
                       const std::streampos ref_byte_pos, const bool flip)
    {
        m_ref_file = ref_file;
        m_ref_byte_pos = ref_byte_pos;
        m_ref_flipped = flip;
    }
    void update_reference(const std::string& ref_file,
                          const std::streampos ref_byte_pos)
    {
        m_ref_file = ref_file;
        m_ref_byte_pos = ref_byte_pos;
    }
    void update_target(const std::string& target_file,
                       const std::streampos byte_pos)
    {
        m_target_file = target_file;
        m_target_byte_pos = byte_pos;
    }
    void add_target(const std::string& target_file,
                    const std::streampos target_byte_pos, const int chr,
                    const int loc, const std::string& ref,
                    const std::string& alt, const bool flipping)
    {
        m_target_file = target_file;
        m_target_byte_pos = target_byte_pos;
        m_chr = chr;
        m_loc = loc;
        m_flipped = flipping;
        m_ref = ref;
        m_alt = alt;
    }
    void add_reference(const std::string& ref_file,
                       const std::streampos ref_byte_pos, const size_t homcom,
                       const size_t het, const size_t homrar,
                       const size_t missing)
    {
        m_ref_file = ref_file;
        m_ref_byte_pos = ref_byte_pos;
        m_homcom = homcom;
        m_ref_het = het;
        m_ref_homrar = homrar;
        m_ref_missing = missing;
    }


    inline void set_flipped() { m_flipped = true; }
    std::string get_rs() const { return m_rs; }
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
     * \param alt is the alternative allele of teh other SNP
     * \param flipped is used as a return value. If flipping is required,
     * flipped = true
     * \return true if it is a match
     */
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
    }

    intptr_t chr() const { return m_chr; }
    intptr_t loc() const { return m_loc; }
    intptr_t category() const { return m_category; }
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
     * \brief Get the SE of the SNP
     * \return the standard error of the SNP
     */
    double get_se() const { return m_standard_error; }
    /*!
     * \brief Return the MAF of the SNP
     * \return the MAF of the SNP based on Base data
     */
    double get_maf() const { return m_maf; }
    /*!
     * \brief Return the p-value threshold of which this SNP falls into
     * \return  the p-value threshold
     */
    double get_threshold() const { return m_p_threshold; }
    /*!
     * \brief Update the beta by minusing the null from it
     * \param null is the null beta
     */
    void update_stat(const double& null)
    {
        double temp = std::abs(m_stat) - null;
        // we multiply the sign of the statistic to the adjusted beta
        m_stat = ((m_stat > 0) - (m_stat < 0)) * std::max(temp, 0.0);
    }
    std::streampos byte_pos() const { return m_target_byte_pos; }
    std::streampos ref_byte_pos() const { return m_ref_byte_pos; }
    std::string file_name() const { return m_target_file; }
    std::string ref_file_name() const { return m_ref_file; }
    std::string rs() const { return m_rs; }
    std::string ref() const { return m_ref; }
    std::string alt() const { return m_alt; }
    bool is_flipped() { return m_flipped; }

    /*!
     * \brief check if this SNP is within the i th region
     * \param i is the index of the region
     * \return true if this SNP falls within the i th region
     */
    inline bool in(size_t i) const
    {
        if (i / BITCT >= m_max_flag_index)
            throw std::out_of_range("Out of range for flag");
        return (IS_SET(m_flags.data(), i));
    }

    void set_flag( const size_t num_region, const std::vector<uintptr_t>& flags){
        m_max_flag_index = BITCT_TO_WORDCT(num_region);
        m_flags = flags;
    }
    /*!
     * \brief Set the SNP to be clumped such that it will no longer be
     * considered in clumping
     */
    void set_clumped() { m_clumped = true; }

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
        bool completed = false;
        // if we want to use proxy, and that our r2 is higher than
        // the proxy threshold, we will do the proxy clumping
        if (use_proxy && r2 > proxy) {
            // If the observed R2 is higher than the proxy clumping threshold,
            // we will capture the flag of the target SNP (using |= )
            for (size_t i_flag = 0; i_flag < m_max_flag_index; ++i_flag) {
                // two become one
                m_flags[i_flag] |= target.m_flags[i_flag];
            }
            // for proxy clumping, the target SNP will always be removed after
            // the clump as the current SNP will represent all sets the target
            // SNP is a member of
            completed = true;
        }
        else
        {
            // otherwise, we will just do noraml clumping
            for (size_t i_flag = 0; i_flag < m_max_flag_index; ++i_flag) {
                // For normal clumping, we will remove set identity from the
                // target SNP whenever both SNPs are within the same set.
                // i.e. if flag of SNP A (current) is 11011 and SNP B (target)
                // is 11110, by the end of clumping, it will become SNP A
                // =11111, SNP B = 00100
                target.m_flags[i_flag] =
                    target.m_flags[i_flag]
                    ^ (m_flags[i_flag] & target.m_flags[i_flag]);
                // if all flags of the target SNP == 0, it means that it no
                // longer represent any gene set and is consided as "clumped"
                completed = (target.m_flags[i_flag] == 0);
            }
        }
        if (completed) {
            // if the target SNP no longer represent any gene set, it is
            // considered as clumped and can be removed
            target.set_clumped();
        }
        m_clumped = true;
        // protect from other SNPs tempering its flags
    }

    /*!
     * \brief Indicate if this snp is clumped
     * \return  Return true if this is clumped
     */
    bool clumped() const { return m_clumped; }
    /*!
     * \brief Set the lower boundary (index of m_existed_snp) of this SNP if it
     * is used as the index
     * \param low the designated bound index
     */
    void set_low_bound(int low) { m_low_bound = low; }
    /*!
     * \brief Set the upper boundary (index of m_existed_snp) of this SNP if it
     * is used as the index
     * \param up the designated bound index
     */
    void set_up_bound(int up) { m_up_bound = up; }
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
    bool get_counts(size_t& homcom, size_t& het, size_t& homrar,
                    size_t& missing) const
    {
        homcom = m_homcom;
        het = m_het;
        homrar = m_homrar;
        missing = m_missing;
        return m_has_count;
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
    void set_counts(size_t homcom, size_t het, size_t homrar, size_t missing)
    {
        m_homcom = homcom;
        m_het = het;
        m_homrar = homrar;
        m_missing = missing;
        m_has_count = true;
    }
    void set_ref_counts(size_t homcom, size_t het, size_t homrar,
                        size_t missing)
    {
        if (m_ref_flipped) {
            // we flip the count here so that the count will be
            // identical to the allele identity in target
            // we process the score, we only need to consider
            // flipping w.r.t target and base
            m_ref_homcom = homrar;
            m_ref_het = het;
            m_ref_homrar = homcom;
            m_ref_missing = missing;
        }
        else
        {
            m_ref_homcom = homcom;
            m_ref_het = het;
            m_ref_homrar = homrar;
            m_ref_missing = missing;
        }
        m_has_ref_count = true;
    }
    /*!
     * \brief Obtain the upper bound of the clump region correspond to this SNP
     * \return the upper bound of the region
     */
    intptr_t up_bound() const { return m_up_bound; }
    /*!
     * \brief Obtain the lower bound of the clump region correspond to this SNP
     * \return the lower bound of the region
     */
    intptr_t low_bound() const { return m_low_bound; }

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
    double m_maf = 0.0;
    double m_standard_error = 0.0;
    int m_chr = -1;
    int m_category = -1;
    int m_loc = -1;
    int m_low_bound = 0;
    int m_up_bound = 0;
    size_t m_homcom = 0;
    size_t m_het = 0;
    size_t m_homrar = 0;
    size_t m_missing = 0;
    size_t m_ref_homcom = 0;
    size_t m_ref_het = 0;
    size_t m_ref_homrar = 0;
    size_t m_ref_missing = 0;
    bool m_has_count = false;
    bool m_has_ref_count = false;
    bool m_clumped = false;
    bool m_flipped = false;
    bool m_ref_flipped = false;
    // prset related
    size_t m_max_flag_index = 0;

    inline std::string complement(const std::string& allele) const
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
};

#endif // SNP_H

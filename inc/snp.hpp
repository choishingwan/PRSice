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
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>


class Genotype;
class SNP
{
public:
    SNP();
    SNP(const std::string rs_id, const int chr, const int loc,
        const std::string ref_allele, const std::string alt_allele,
        const std::string file_name, const int num_line);
    virtual ~SNP();

    void set_statistic(const double stat, const double se, const double p_value,
                       const int category, const double p_threshold)
    {
        statistic.stat = stat;
        statistic.se = se;
        statistic.p_value = p_value;
        threshold.category = category;
        threshold.p_threshold = p_threshold;
    };
    void set_flipped() { statistic.flipped = true; };
    std::string get_rs() const { return basic.rs; };
    static std::vector<size_t> sort_by_p(const std::vector<SNP>& input);

    bool operator==(const SNP& Ref) const
    {
        if (basic.chr == Ref.basic.chr && basic.loc == Ref.basic.loc
            && basic.rs.compare(Ref.basic.rs) == 0)
        {
            if (basic.ref.compare(Ref.basic.ref) == 0) {
                if (!basic.alt.empty() && !Ref.basic.alt.empty()) {
                    return basic.alt.compare(Ref.basic.alt) == 0;
                }
                else
                    return true;
            }
            else if (complement(basic.ref).compare(Ref.basic.ref) == 0)
            {
                if (!basic.alt.empty() && !Ref.basic.alt.empty()) {
                    return complement(basic.alt).compare(Ref.basic.alt) == 0;
                }
                else
                    return true;
            }
            else if (!basic.alt.empty() && !Ref.basic.alt.empty())
            {
                if (basic.ref.compare(Ref.basic.alt) == 0
                    && basic.alt.compare(Ref.basic.ref) == 0)
                {
                    return true;
                }
                if (complement(basic.ref).compare(Ref.basic.alt) == 0
                    && complement(basic.alt).compare(Ref.basic.alt) == 0)
                {
                    return true;
                }
                return false;
            }
            else
                return false; // cannot flip nor match
        }
        else
        {
            return false;
        }
    };

    inline void flipped() { statistic.flipped = true; };
    inline void fill_info(int chr, int loc, std::string alt)
    {
        if (basic.chr == -1) basic.chr = chr;
        if (basic.loc == -1) basic.loc = loc;
        if (basic.alt.empty()) basic.alt = alt;
    };
    inline bool matching(int chr, int loc, std::string ref, std::string alt,
                         bool& flipped)
    {
        misc::trim(ref);
        misc::trim(alt);
        misc::trim(basic.ref);
        misc::trim(basic.alt);
        if (chr != -1 && basic.chr != -1 && chr != basic.chr) return false;
        if (loc != -1 && basic.loc != -1 && loc != basic.loc) return false;
        if (basic.ref.compare(ref) == 0) {
            if (!basic.alt.empty() && !alt.empty()) {
                return basic.alt.compare(alt) == 0;
            }
            else
                return true;
        }
        else if (complement(basic.ref).compare(ref) == 0)
        {
            if (!basic.alt.empty() && !alt.empty()) {
                return complement(basic.alt).compare(alt) == 0;
            }
            else
                return true;
        }
        else if (!basic.alt.empty() && !alt.empty())
        {
            if (basic.ref.compare(alt) == 0 && basic.alt.compare(ref) == 0) {
                flipped = true;
                return true;
            }
            if (complement(basic.ref).compare(alt) == 0
                && complement(basic.alt).compare(ref) == 0)
            {
                flipped = true;
                return true;
            }
            return false;
        }
        else
            return false; // cannot flip nor match
    };

    int chr() const { return basic.chr; };
    int loc() const { return basic.loc; };
    int snp_id() const { return file_info.id; };
    int category() const { return threshold.category; };
    double p_value() const { return statistic.p_value; };
    double stat() const { return statistic.stat; };
    double get_threshold() const { return threshold.p_threshold; };
    std::string file_name() const { return file_info.file; };
    std::string rs() const { return basic.rs; };
    std::string ref() const { return basic.ref; };
    std::string alt() const { return basic.alt; };
    bool is_flipped() { return statistic.flipped; };
    /*
    inline bool in(size_t i) const
    {
            if(i/m_bit_size >= m_flags.size()) throw std::out_of_range("Out of
    range for flag"); return (m_flags[i/m_bit_size] >> i%m_bit_size) & ONE; // 1
    = true, 0 = false
    }
     */
    inline bool in(size_t i) const
    {
        if (i / BITCT >= m_max_flag_index)
            throw std::out_of_range("Out of range for flag");
        return ((m_flags[i / BITCT] >> (i % BITCT)) & 1);
    }

    // void set_flag(std::vector<long_type> flag) { m_flags = flag; };

    void set_flag(Region& region)
    {
        m_max_flag_index = BITCT_TO_WORDCT(region.size());
        m_flags.resize(m_max_flag_index);
        // m_flag = new uintptr_t[m_max_flag_index];
        region.check(std::to_string(basic.chr), basic.loc, m_flags);
    };


    void add_clump(std::vector<size_t>& i)
    {
        clump_info.target.insert(clump_info.target.end(), i.begin(), i.end());
    };
    void add_clump_r2(std::vector<double>& i)
    {
        clump_info.r2.insert(clump_info.r2.end(), i.begin(), i.end());
    };

    void add_clump(size_t index, double r2)
    {
        std::lock_guard<std::mutex> self_lock(thread_safty_inspector);
        clump_info.target.push_back(index);
        clump_info.r2.push_back(r2);
    }
    void set_clumped() { clump_info.clumped = true; };
    void clump(std::vector<SNP>& snp_list, double proxy = 2);
    bool clumped() const { return clump_info.clumped; };

    uintptr_t* clump_geno() { return clump_info.geno1.data(); };
    ;
    uint32_t* tot() { return clump_info.tot.data(); };
    uint32_t get_tot(int i) { return clump_info.tot.at(i); };
    bool clump_missing() const { return clump_info.contain_missing; };
    void clean_clump()
    {
        clump_info.geno1 = std::vector<uintptr_t>();
        clump_info.tot = std::vector<uint32_t>();
    }

    // Now only genotype can access this pointer
    void geno_size(size_t size, Passkey<Genotype>)
    {
        clump_info.geno1.resize(size, 0);
    }
    bool has_geno() const { return !clump_info.geno1.empty(); }
    uintptr_t* geno(Passkey<Genotype>)
    {
        // pretty sure this is unsafe practice
        // as we expose the pointer to other stuff
        // not sure how to limit the access of this function
        return clump_info.geno1.data();
    };
    void load_tot(const uint32_t founder_ctv3)
    {
        clump_info.tot.resize(6, 0);
        clump_info.tot[0] =
            popcount_longs(clump_info.geno1.data(), founder_ctv3);
        clump_info.tot[1] = popcount_longs(
            &(clump_info.geno1.data()[founder_ctv3]), founder_ctv3);
        clump_info.tot[2] = popcount_longs(
            &(clump_info.geno1.data()[2 * founder_ctv3]), founder_ctv3);
    }
    void set_contain_missing(int contain_miss)
    {
        clump_info.contain_missing = (contain_miss == 3);
    }

    bool valid() const { return basic.valid; };
    void invalidate() { basic.valid = false; };

    // need this so that we can allow a member mutex class
    // Copy assignment
    SNP& operator=(const SNP& other)
    {
        if (this == &other) return *this;
        std::lock(thread_safty_inspector, other.thread_safty_inspector);
        std::lock_guard<std::mutex> self_lock(thread_safty_inspector,
                                              std::adopt_lock);
        std::lock_guard<std::mutex> other_lock(other.thread_safty_inspector,
                                               std::adopt_lock);
        // now we need to copy everything
        clump_info.clumped = other.clump_info.clumped;
        clump_info.contain_missing = other.clump_info.contain_missing;
        clump_info.contain_geno = other.clump_info.contain_geno;
        clump_info.r2 = other.clump_info.r2;
        clump_info.target = other.clump_info.target;
        clump_info.tot = other.clump_info.tot;
        clump_info.geno1 = other.clump_info.geno1;
        // on purposely ignore genotype pointer (we will remove it later)
        basic.ref = other.basic.ref;
        basic.alt = other.basic.alt;
        basic.rs = other.basic.rs;
        basic.chr = other.basic.chr;
        basic.loc = other.basic.loc;
        file_info.file = other.file_info.file;
        file_info.id = other.file_info.id;
        statistic.stat = other.statistic.stat;
        statistic.se = other.statistic.se;
        statistic.p_value = other.statistic.p_value;
        statistic.flipped = other.statistic.flipped;
        threshold.category = other.threshold.category;
        threshold.p_threshold = other.threshold.p_threshold;
        m_max_flag_index = other.m_max_flag_index;
        m_flags = other.m_flags;
        return *this;
    }
    // Copy initialization
    SNP(const SNP& other)
    {
        std::lock_guard<std::mutex> lock(other.thread_safty_inspector);
        clump_info.clumped = other.clump_info.clumped;
        clump_info.contain_missing = other.clump_info.contain_missing;
        clump_info.contain_geno = other.clump_info.contain_geno;
        clump_info.r2 = other.clump_info.r2;
        clump_info.target = other.clump_info.target;
        clump_info.tot = other.clump_info.tot;
        clump_info.geno1 = other.clump_info.geno1;
        // on purposely ignore genotype pointer (we will remove it later)
        basic.ref = other.basic.ref;
        basic.alt = other.basic.alt;
        basic.rs = other.basic.rs;
        basic.chr = other.basic.chr;
        basic.loc = other.basic.loc;
        file_info.file = other.file_info.file;
        file_info.id = other.file_info.id;
        statistic.stat = other.statistic.stat;
        statistic.se = other.statistic.se;
        statistic.p_value = other.statistic.p_value;
        statistic.flipped = other.statistic.flipped;
        threshold.category = other.threshold.category;
        threshold.p_threshold = other.threshold.p_threshold;
        m_max_flag_index = other.m_max_flag_index;
        m_flags = other.m_flags;
    }

    // move assignment
    SNP& operator=(SNP&& other)
    {
        if (this == &other) return *this;
        std::lock(thread_safty_inspector, other.thread_safty_inspector);
        std::lock_guard<std::mutex> self_lock(thread_safty_inspector,
                                              std::adopt_lock);
        std::lock_guard<std::mutex> other_lock(other.thread_safty_inspector,
                                               std::adopt_lock);
        clump_info.clumped = std::move(other.clump_info.clumped);
        other.clump_info.clumped = false;
        clump_info.contain_missing =
            std::move(other.clump_info.contain_missing);
        other.clump_info.contain_missing = false;
        clump_info.contain_geno = std::move(other.clump_info.contain_geno);
        other.clump_info.contain_geno = false;
        clump_info.r2 = std::move(other.clump_info.r2);
        other.clump_info.r2.clear();
        clump_info.target = std::move(other.clump_info.target);
        other.clump_info.target.clear();
        clump_info.tot = std::move(other.clump_info.tot);
        other.clump_info.tot.clear();
        clump_info.geno1 = std::move(other.clump_info.geno1);
        other.clump_info.geno1.clear();
        // on purposely ignore genotype pointer (we will remove it later)
        basic.ref = std::move(other.basic.ref);
        other.basic.ref = "";
        basic.alt = std::move(other.basic.alt);
        other.basic.alt = "";
        basic.rs = std::move(other.basic.rs);
        other.basic.rs = "";
        basic.chr = std::move(other.basic.chr);
        other.basic.chr = 0;
        basic.loc = std::move(other.basic.loc);
        other.basic.loc = 0;
        file_info.file = std::move(other.file_info.file);
        other.file_info.file = "";
        file_info.id = std::move(other.file_info.id);
        other.file_info.id = 0;
        statistic.stat = std::move(other.statistic.stat);
        other.statistic.stat = 0.0;
        statistic.se = std::move(other.statistic.se);
        other.statistic.se = 0.0;
        statistic.p_value = std::move(other.statistic.p_value);
        other.statistic.p_value = 0.0;
        statistic.flipped = std::move(other.statistic.flipped);
        threshold.category = std::move(other.threshold.category);
        other.threshold.category = 0.0;
        threshold.p_threshold = std::move(other.threshold.p_threshold);
        other.threshold.p_threshold = 0.0;
        m_max_flag_index = std::move(other.m_max_flag_index);
        other.m_max_flag_index = 0;
        m_flags = std::move(other.m_flags);
        other.m_flags.clear();
        return *this;
    }

private:
    // basic info
    mutable std::mutex thread_safty_inspector;
    struct
    {
        bool clumped;
        bool contain_missing;
        bool contain_geno;
        std::vector<double> r2;
        std::vector<size_t> target;
        std::vector<uint32_t> tot;
        std::vector<uintptr_t> geno1;
    } clump_info;

    struct
    {
        std::string ref;
        std::string alt;
        std::string rs;
        int chr;
        int loc;
        bool valid;
    } basic;

    struct
    {
        std::string file;
        int id;
    } file_info;

    struct
    {
        double stat;
        double se;
        double p_value;
        bool flipped;
    } statistic;

    struct
    {
        int category;
        double p_threshold;
    } threshold;

    // This indicate where this SNP's bound is at
    // useful for PRSlice and also clumping
    // thinking about it. Even if the location isn't given for
    // PRSet or PRSlice, we can still use the coordinates from
    // the target / reference file
    // the bound is [ )
    // prset related
    size_t m_max_flag_index = 0;
    std::vector<uintptr_t> m_flags;

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

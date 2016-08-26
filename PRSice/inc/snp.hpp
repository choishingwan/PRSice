#ifndef SNP_H
#define SNP_H

#include <string>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>

#include "commander.hpp"
#include "misc.hpp"

class SNP
{
    public:
        SNP();
#if defined(__LP64__) || defined(_WIN64)
        SNP(const std::string rs_id, const std::string chr, const size_t loc, const std::string ref_allele, const std::string alt_allele, const double statistic, const double se, const double p_value, uint64_t *flag);
#else
        SNP(const std::string rs_id, const std::string chr, const size_t loc, const std::string ref_allele, const std::string alt_allele, const double statistic, const double se, const double p_value, uint64_t *flag);
#endif
        virtual ~SNP();
        std::string get_ref_allele() { return m_ref_allele; }
        std::string get_alt_allele() { return m_alt_allele; }
        std::string get_rs_id() { return m_rs_id; }
        std::string get_chr() { return m_chr; }
        size_t get_loc() { return m_loc; }
        double set_stat() { return m_stat; }
        double get_p_value() { return m_p_value; }
        static std::vector<size_t> sort_by_p(const std::vector<SNP> &input);
        static std::vector<int> get_index(const Commander &c_commander, const std::string &c_input);
    protected:
    private:
        static size_t index_check(const std::string &c_in);
        static size_t index_check(const std::string &c_in, const std::vector<std::string> &c_header, const std::string &typeOfError);
        std::string m_ref_allele;
        std::string m_alt_allele;
        std::string m_rs_id;
        std::string m_chr;
        size_t m_loc;
        double m_stat;
        double m_standard_error;
        double m_p_value;
        std::vector<size_t> m_clump_target; // index of SNPs that are clumped under this SNP
#if defined(__LP64__) || defined(_WIN64)
        uint64_t *m_flags;
#else
        uint32_t *m_flags;
#endif
};

#endif // SNP_H

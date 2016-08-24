#ifndef SNP_H
#define SNP_H

#include <string>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <map>

#include "commander.hpp"
#include "misc.hpp"

class SNP
{
    public:
        SNP();
        SNP(const std::string rs_id, const std::string chr, const size_t loc, const std::string ref_allele, const std::string alt_allele, const double statistic, const double se, const double p_value);
        virtual ~SNP();
        std::string get_ref_allele() { return m_ref_allele; }
        std::string get_alt_allele() { return m_alt_allele; }
        std::string get_rs_id() { return m_rs_id; }
        std::string get_chr() { return m_chr; }
        size_t get_loc() { return m_loc; }
        double set_stat() { return m_stat; }
        double get_p_value() { return m_p_value; }
    protected:
    private:
        static std::vector<int> get_index(const Commander &c_commander, const std::string &c_input);
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
};

#endif // SNP_H
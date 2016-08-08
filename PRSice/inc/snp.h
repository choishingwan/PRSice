#ifndef SNP_H
#define SNP_H

#include <string>
#include <fstream>
#include <boost/ptr_container/ptr_vector.hpp>
#include <stdexcept>
#include "usefulTools.h"
#include "commander.h"

class SNP
{
    public:
        SNP();
        SNP(const std::string rs_id, const std::string chr, const size_t loc, const std::string ref_allele, const std::string alt_allele, const double statistic, const double p_value);
        virtual ~SNP();
        static void read_snp(const Commander &c_commander, boost::ptr_vector<SNP> &snp_list);
        std::string get_ref_allele() { return m_ref_allele; }
        std::string get_alt_allele() { return m_alt_allele; }
        std::string get_rs_id() { return m_rs_id; }
        std::string get_chr() { return m_chr; }
        size_t get_loc() { return m_loc; }
        double set_stat() { return m_stat; }
        double get_p_value() { return m_p_value; }
    protected:
    private:
        std::string m_ref_allele;
        std::string m_alt_allele;
        std::string m_rs_id;
        std::string m_chr;
        size_t m_loc;
        double m_stat;
        double m_p_value;
};

#endif // SNP_H

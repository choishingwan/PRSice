#ifndef PRSICE_H
#define PRSICE_H

#include <string>
#include <fstream>
#include <stdexcept>
#include <Eigen/Dense>
#include <boost/ptr_container/ptr_vector.hpp>
#include <vector>
#include <algorithm>
#include <map>
#include <stdio.h>
#include <thread>
#include <unordered_map>
#include "commander.hpp"
#include "misc.hpp"
#include "plink.hpp"
#include "snp.hpp"
#include "region.hpp"
#include "regression.hpp"
#include "storage.hpp"
//This should be the class to handle all the procedures
class PRSice
{
public:
    PRSice();
    PRSice(std::string base_name, int index): m_current_base(base_name), m_base_index(index)
    {
        if(index < 0)
        {
            throw std::out_of_range("Index cannot be less than 0");
        }
    };
    virtual ~PRSice();
    void get_snp(const Commander &c_commander, Region &region, const double &threshold);
    void clump(const Commander &c_commander, const Region &region);

protected:
private:
//      This is used for thead safety
    static std::mutex score_mutex;
    int m_base_index;
    std::string m_current_base;
    std::vector<p_partition> m_partition;
    std::vector<std::vector<PRSice_result> > m_prs_results;
    std::vector<PRSice_best> m_best_threshold;
    std::vector<std::vector<prs_score> > m_best_score;
    std::vector<size_t> m_num_snp_included;
    std::vector<std::vector<prs_score> > m_current_prs;
    boost::ptr_vector<SNP> m_snp_list;
    std::unordered_map<std::string, size_t> m_snp_index;
    std::vector<std::string> m_chr_list;


    void check_inclusion(std::unordered_map<std::string, size_t> &inclusion, const std::string &c_target_bim_name,
                         size_t &num_ambig, size_t &not_found, size_t &num_duplicate);
    void run_prs(const Commander &c_commander, const std::map<std::string, size_t> &inclusion,
                 const Region &c_region);
    void categorize(const Commander &c_commander, const std::map<std::string, size_t> &inclusion, bool &pre_run);
    Eigen::VectorXd gen_pheno_vec(const std::string &c_target, const std::string c_pheno, const int pheno_index,
                                  bool target_binary, std::map<std::string, size_t> &fam_index, std::vector<prs_score> &prs_fam);
    Eigen::MatrixXd gen_cov_matrix(const std::string &c_cov_file, const std::vector<std::string> &c_cov_header,
                                   std::map<std::string, size_t> &fam_index);



    bool get_prs_score(const std::string &target, size_t &cur_index);

    void thread_score( Eigen::MatrixXd &independent_variables, const Eigen::VectorXd &c_pheno,
                       const std::map<std::string, size_t> &c_fam_index, size_t region_start, size_t region_end,
                       bool target_binary, double threshold, size_t thread);
    void calculate_scores(const Commander &c_commander,  const Region &c_region, size_t cur_start_index,
                          Eigen::MatrixXd &independent_variables, Eigen::VectorXd &phenotype, const std::string &pheno_name,
                          const std::map<std::string, size_t> &fam_index);
    void prs_output(const Commander &c_commander, const Region &c_region, const double null_r2,
                    const std::string pheno_name) const;

};

#endif // PRSICE_H

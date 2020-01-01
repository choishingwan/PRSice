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


#ifndef PRSICE_INC_STORAGE_HPP_
#define PRSICE_INC_STORAGE_HPP_
#include "enumerators.h"
#include <Eigen/Dense>
#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <vector>
// From http://stackoverflow.com/a/12927952/1441789

struct Regress
{
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PQR;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType Pmat;
    Eigen::Index rank;
    Eigen::MatrixXd Rinv;
    Eigen::MatrixXd YCov;
    Eigen::VectorXd beta;
    Eigen::VectorXd se;
    Eigen::VectorXd effects;
    Eigen::VectorXd prs;
    Eigen::VectorXd fitted;
    Eigen::VectorXd resid;
    Eigen::VectorXd se_base;
};

struct PRS
{
    double prs;
    size_t num_snp;
    PRS() : prs(0.0), num_snp(0) {}
};

struct Sample_ID
{
    std::string FID;
    std::string IID;
    std::string pheno;
    bool founder;
    Sample_ID(const std::string& F, const std::string& I, const std::string& P,
              const bool& Founder)
        : FID(F), IID(I), pheno(P), founder(Founder)
    {
    }
    Sample_ID() : FID(""), IID(""), pheno(""), founder(false) {}
};

struct BaseFile
{
    std::vector<size_t> column_index =
        std::vector<size_t>(+BASE_INDEX::MAX + 1, 0);
    std::vector<std::string> column_name = {
        "CHR", "A2", "BP", "SE", "INFO,0.9", "", "", "A1", "SNP", "P", ""};
    // use int as vector<bool> is abnormal
    std::vector<int> has_column = std::vector<int>(+BASE_INDEX::MAX + 1, false);
    std::string file_name;
    int is_index = false;
    int is_beta = false;
    int is_or = false;
};

struct GenoFile
{
    std::string file_name;
    std::string file_list;
    std::string keep;
    std::string remove;
    std::string type = "bed";
    int num_autosome = 22;
    int hard_coded = false;
    int is_ref = false;
};

struct Phenotype
{
    std::vector<std::string> pheno_col;
    std::vector<std::string> cov_colname;
    std::vector<std::string> factor_cov;
    std::vector<double> prevalence;
    std::vector<size_t> col_index_of_cov;
    std::vector<size_t> col_index_of_factor_cov;
    std::vector<size_t> pheno_col_idx;
    std::vector<bool> binary;
    std::vector<bool> skip_pheno;
    std::string pheno_file;
    std::string cov_file;
    int ignore_fid = false;
};
struct FileInfo
{
    size_t name_idx = ~size_t(0);
    long long byte_pos = 0;
};
struct Permutations
{
    size_t num_permutation = 0;
    std::random_device::result_type seed = std::random_device()();
    int logit_perm = false;
    bool run_perm = false;
    bool run_set_perm = false;
};

// use size_t for low bound and up bound just in case if we encounter some huge
// chromosomes
struct SNPClump
{
    std::vector<uintptr_t> flags;
    size_t low_bound = ~size_t(0);
    size_t up_bound = ~size_t(0);
    size_t max_flag_idx = 0;
    bool clumped = false;
};

struct AlleleCounts
{
    size_t homcom = 0;
    size_t het = 0;
    size_t homrar = 0;
    size_t missing = 0;
    bool has_count = false;
};

struct GeneSets
{
    std::vector<std::string> msigdb;
    std::vector<std::string> bed;
    std::vector<std::string> snp;
    std::vector<std::string> feature;
    std::string background;
    std::string gtf;
    size_t wind_3 = 0;
    size_t wind_5 = 0;
    int full_as_background = false;
    bool run = false;
};
struct PThresholding
{
    std::vector<double> bar_levels;
    double lower = 5e-8;
    double inter = 0.00005;
    double upper = 0.5;
    int fastscore = false;
    int no_full = false;
    // indicate if we want to do thresholding with set based analysis
    bool set_threshold = false;
};

struct CalculatePRS
{
    MISSING_SCORE missing_score = MISSING_SCORE::MEAN_IMPUTE;
    SCORING scoring_method = SCORING::AVERAGE;
    MODEL genetic_model = MODEL::ADDITIVE;
    int thread = 1;
    int no_regress = false;
    int non_cumulate = false;
    int use_ref_maf = false;
};

struct QCFiltering
{
    double dose_threshold = 0.0;
    double geno = 1.0;
    double hard_threshold = 0.1;
    double maf = 0.0;
    double maf_case = 0.0;
    double info_score = 0.0;
    bool provided_hard_thres = false;
    bool provided_dose_thres = false;
};

struct Clumping
{
    double r2 = 0.1;
    double proxy = -1;
    double pvalue = 1;
    size_t distance = 250000;
    int no_clump = false;
    bool use_proxy = false;
    bool provided_distance = false;
};
// Passkey idiom, allow safer access to
// the raw pointer info in SNP
template <typename T>
class Passkey
{
private:
    friend T;
    Passkey() {}
    Passkey(const Passkey&) {}
    Passkey& operator=(const Passkey&) = delete;
};

#endif /* PRSICE_INC_STORAGE_HPP_ */

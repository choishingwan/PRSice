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
#include <cstdint>
#include <memory>
#include <string>
#include <vector>
// From http://stackoverflow.com/a/12927952/1441789

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

struct MAF_Store
{
    double maf;
    size_t index;
    int category;
};


struct BaseFile
{

    std::vector<size_t> column_index;
    std::vector<std::string> column_name = {
        "CHR", "A2", "BP", "SE", "INFO,0.9", "", "", "A1", "SNP", "P", ""};
    std::vector<int> has_column; // use int as vector<bool> is abnormal
    std::string file_name;
    int is_index = false;
    int is_beta = false;
    int is_or = false;
    BaseFile()
    {
        column_index.resize(+BASE_INDEX::MAX, 0);
        has_column.resize(+BASE_INDEX::MAX, false);
    }
};

struct GenoFile
{
    std::string file_name;
    std::string file_list;
    std::string keep;
    std::string remove;
    std::string type;
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
    std::vector<bool> binary;
    std::string pheno_file;
    std::string cov_file;
    int ignore_fid = false;
};

struct GeneSets
{
    std::vector<std::string> msigdb;
    std::vector<std::string> bed;
    std::vector<std::string> snp;
    std::vector<std::string> feature;
    std::string background;
    std::string gtf;
    size_t perm;
    unsigned long long wind_3 = 0;
    unsigned long long wind_5 = 0;
    int full_as_background = false;
    bool run = false;
    bool run_perm = false;
};
struct PThresholding
{
    std::vector<double> bar_levels;
    double lower = 5e-8;
    double inter = 0.00005;
    double upper = 0.5;
    int fastscore = false;
    int no_full = false;
    bool set_threshold = false;
};

struct CalculatePRS
{
    MISSING_SCORE missing_score = MISSING_SCORE::MEAN_IMPUTE;
    SCORING scoring_method = SCORING::AVERAGE;
    MODEL genetic_model = MODEL::ADDITIVE;
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
    unsigned long long distance = 250000;
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

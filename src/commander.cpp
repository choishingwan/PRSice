﻿// This file is part of PRSice-2, copyright (C) 2016-2019
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

#include "commander.hpp"

Commander::Commander()
{
    set_help_message();
    m_reference.is_ref = true;
}

bool Commander::process_command(int argc, char* argv[], Reporter& reporter)
{
    bool early_termination = false;
    try
    {
        bool error = init(argc, argv, early_termination, reporter);
        if (early_termination) return false;
        error |= validate_command(reporter);
        std::string message = get_program_header(argv[0]);
        for (auto&& com : m_parameter_log)
        { message.append(" \\\n    --" + com.first + " " + com.second); }
        message.append("\n");
        reporter.report(message, false);
        if (!m_error_message.empty()) reporter.report(m_error_message);
        if (error) throw std::runtime_error(m_error_message);
        return true;
    }
    catch (const std::runtime_error& er)
    {
        reporter.simple_report(er.what());
        throw er;
    }
}
bool Commander::init(int argc, char* argv[], bool& early_termination,
                     Reporter& reporter)
{
    // initialize get_opt. To be honest, with PRSice usage, this shouldn't be
    // required, but then it is required for our unit test where we repeatedly
    // test the get_opt
    optind = 0;
    if (argc <= 1)
    {
        // use default PRSice
        reporter.simple_report(m_help_message);
        throw std::runtime_error("Please provide the required parameters");
    }
    const char* optString = "b:B:c:C:f:F:g:h?i:k:l:L:m:n:o:p:s:t:u:v";

    const struct option longOpts[] = {
        // parameters with short flags
        {"base", required_argument, nullptr, 'b'},
        {"bed", required_argument, nullptr, 'B'},
        {"cov-col", required_argument, nullptr, 'c'},
        // Allow short form for cov-file and pheno-file
        {"cov", required_argument, nullptr, 'C'},
        {"cov-file", required_argument, nullptr, 'C'},
        {"pheno", required_argument, nullptr, 'f'},
        {"pheno-file", required_argument, nullptr, 'f'},
        {"pheno-col", required_argument, nullptr, 'F'},
        {"gtf", required_argument, nullptr, 'g'},
        {"help", no_argument, nullptr, 'h'},
        {"interval", required_argument, nullptr, 'i'},
        {"prevalence", required_argument, nullptr, 'k'},
        {"lower", required_argument, nullptr, 'l'},
        {"ld", required_argument, nullptr, 'L'},
        {"msigdb", required_argument, nullptr, 'm'},
        {"thread", required_argument, nullptr, 'n'},
        {"out", required_argument, nullptr, 'o'},
        {"pvalue", required_argument, nullptr, 'p'},
        {"seed", required_argument, nullptr, 's'},
        {"target", required_argument, nullptr, 't'},
        {"upper", required_argument, nullptr, 'u'},
        {"version", no_argument, nullptr, 'v'},
        // flags, only need to set them to true
        {"allow-inter", no_argument, &m_allow_inter, 1},
        {"all-score", no_argument, &m_print_all_scores, 1},
        {"beta", no_argument, &m_base_info.is_beta, 1},
        {"clump-only", no_argument, &m_clump_only, 1},
        {"fastscore", no_argument, &m_p_thresholds.fastscore, 1},
        {"full-back", no_argument, &m_prset.full_as_background, 1},
        {"hard", no_argument, &m_target.hard_coded, 1},
        {"ignore-fid", no_argument, &m_pheno_info.ignore_fid, 1},
        {"index", no_argument, &m_base_info.is_index, 1},
        {"keep-ambig", no_argument, &m_keep_ambig, 1},
        {"logit-perm", no_argument, &m_perm_info.logit_perm, 1},
        {"no-clump", no_argument, &m_clump_info.no_clump, 1},
        {"non-cumulate", no_argument, &m_prs_info.non_cumulate, 1},
        {"no-default", no_argument, &m_user_no_default, 1},
        {"no-full", no_argument, &m_p_thresholds.no_full, 1},
        {"no-regress", no_argument, &m_prs_info.no_regress, 1},
        {"nonfounders", no_argument, &m_include_nonfounders, 1},
        {"or", no_argument, &m_base_info.is_or, 1},
        {"print-snp", no_argument, &m_print_snp, 1},
        {"ultra", no_argument, &m_ultra_aggressive, 1},
        {"use-ref-maf", no_argument, &m_prs_info.use_ref_maf, 1},
        // long flags, need to work on them
        {"A1", required_argument, nullptr, 0},
        {"A2", required_argument, nullptr, 0},
        {"a1", required_argument, nullptr, 0},
        {"a2", required_argument, nullptr, 0},
        {"background", required_argument, nullptr, 0},
        {"bar-levels", required_argument, nullptr, 0},
        {"base-info", required_argument, nullptr, 0},
        {"base-maf", required_argument, nullptr, 0},
        {"binary-target", required_argument, nullptr, 0},
        {"bp", required_argument, nullptr, 0},
        {"chr", required_argument, nullptr, 0},
        {"chr-id", required_argument, nullptr, 0},
        {"clump-kb", required_argument, nullptr, 0},
        {"clump-p", required_argument, nullptr, 0},
        {"clump-r2", required_argument, nullptr, 0},
        {"cov-factor", required_argument, nullptr, 0},
        {"dose-thres", required_argument, nullptr, 0},
        {"exclude", required_argument, nullptr, 0},
        {"extract", required_argument, nullptr, 0},
        {"feature", required_argument, nullptr, 0},
        {"geno", required_argument, nullptr, 0},
        {"hard-thres", required_argument, nullptr, 0},
        {"id-delim", required_argument, nullptr, 0},
        {"info", required_argument, nullptr, 0},
        {"info-type", required_argument, nullptr, 0},
        {"keep", required_argument, nullptr, 0},
        {"ld-dose-thres", required_argument, nullptr, 0},
        {"ld-keep", required_argument, nullptr, 0},
        {"ld-list", required_argument, nullptr, 0},
        {"ld-type", required_argument, nullptr, 0},
        {"ld-remove", required_argument, nullptr, 0},
        {"ld-maf", required_argument, nullptr, 0},
        {"ld-geno", required_argument, nullptr, 0},
        {"ld-hard-thres", required_argument, nullptr, 0},
        {"ld-info", required_argument, nullptr, 0},
        {"maf", required_argument, nullptr, 0},
        {"memory", required_argument, nullptr, 0},
        {"missing", required_argument, nullptr, 0},
        {"model", required_argument, nullptr, 0},
        {"num-auto", required_argument, nullptr, 0},
        {"perm", required_argument, nullptr, 0},
        {"proxy", required_argument, nullptr, 0},
        {"remove", required_argument, nullptr, 0},
        {"score", required_argument, nullptr, 0},
        {"set-perm", required_argument, nullptr, 0},
        {"snp", required_argument, nullptr, 0},
        {"snp-set", required_argument, nullptr, 0},
        {"stat", required_argument, nullptr, 0},
        {"target-list", required_argument, nullptr, 0},
        {"type", required_argument, nullptr, 0},
        {"wind-5", required_argument, nullptr, 0},
        {"wind-3", required_argument, nullptr, 0},
        {"x-range", required_argument, nullptr, 0},
        {nullptr, 0, nullptr, 0}};
    try
    {
        return parse_command(argc, argv, optString, longOpts, early_termination,
                             reporter);
    }
    catch (const std::runtime_error& er)
    {
        throw er;
    }
}

bool Commander::parse_command(int argc, char* argv[], const char* optString,
                              const struct option longOpts[],
                              bool& early_termination, Reporter& reporter)
{
    int32_t max_threads = maximum_thread();
    int longIndex = 0;
    int opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    std::string command;
    bool error = false;
    // TODO: Standardize function return such that they will return error if
    // failed
    // opterr = 0;
    while (opt != -1)
    {
        switch (opt)
        {
        case 0:
            command = longOpts[longIndex].name;
            if (longOpts[longIndex].flag != nullptr) break;
            // reorganize all long ops according to alphabetical order
            else if (command == "A1" || command == "a1")
                set_string(optarg, "a1", +BASE_INDEX::EFFECT);
            else if (command == "A2" || command == "a2")
                set_string(optarg, "a2", +BASE_INDEX::NONEFFECT);
            else if (command == "background")
                set_string(optarg, command, m_prset.background);
            else if (command == "bar-levels")
            {
                error |= !load_numeric_vector<double>(
                    optarg, command, m_p_thresholds.bar_levels);
                m_p_thresholds.set_threshold = true;
            }
            else if (command == "base-info")
                set_string(optarg, command, +BASE_INDEX::INFO);
            else if (command == "base-maf")
                set_string(optarg, command, +BASE_INDEX::MAF);
            else if (command == "binary-target")
                error |=
                    !parse_binary_vector(optarg, command, m_pheno_info.binary);
            else if (command == "bp")
                set_string(optarg, command, +BASE_INDEX::BP);
            else if (command == "chr")
                set_string(optarg, command, +BASE_INDEX::CHR);
            else if (command == "chr-id")
                set_string(optarg, command, m_chr_id_formula);
            else if (command == "clump-kb")
            {
                error |= !parse_unit_value(optarg, command, 1,
                                           m_clump_info.distance);
                m_clump_info.provided_distance = true;
            }
            else if (command == "clump-p")
                error |=
                    !set_numeric<double>(optarg, command, m_clump_info.pvalue);
            else if (command == "clump-r2")
                error |= !set_numeric<double>(optarg, command, m_clump_info.r2);
            else if (command == "cov-factor")
                load_string_vector(optarg, command, m_pheno_info.factor_cov);
            else if (command == "dose-thres")
                error |= !set_numeric<double>(optarg, command,
                                              m_target_filter.dose_threshold);
            else if (command == "exclude")
                set_string(optarg, command, m_exclude_file);
            else if (command == "extract")
                set_string(optarg, command, m_extract_file);
            else if (command == "feature")
                load_string_vector(optarg, command, m_prset.feature);
            else if (command.compare("geno") == 0)
                error |=
                    !set_numeric<double>(optarg, command, m_target_filter.geno);
            else if (command == "hard-thres")
                error |= !set_numeric<double>(optarg, command,
                                              m_target_filter.hard_threshold);
            else if (command == "id-delim")
                set_string(optarg, command, m_id_delim, m_set_delim, true);
            else if (command == "info")
                error |= !set_numeric<double>(optarg, command,
                                              m_target_filter.info_score);
            else if (command == "info-type")
                error |= !set_info(optarg);
            else if (command == "keep")
                set_string(optarg, command, m_target.keep);
            else if (command == "ld-dose-thres")
                error |= !set_numeric<double>(optarg, command,
                                              m_ref_filter.dose_threshold);
            else if (command == "ld-geno")
                error |=
                    !set_numeric<double>(optarg, command, m_ref_filter.geno);
            else if (command == "ld-hard-thres")
                error |= !set_numeric<double>(optarg, command,
                                              m_ref_filter.hard_threshold);
            else if (command == "ld-info")
                error |= !set_numeric<double>(optarg, command,
                                              m_ref_filter.info_score);
            else if (command == "ld-keep")
                set_string(optarg, command, m_reference.keep);
            else if (command == "ld-list")
                set_string(optarg, command, m_reference.file_list);
            else if (command == "ld-maf")
                error |=
                    !set_numeric<double>(optarg, command, m_ref_filter.maf);
            else if (command == "ld-remove")
                set_string(optarg, command, m_reference.remove);
            else if (command == "ld-type")
                set_string(optarg, command, m_reference.type);
            else if (command == "maf")
                error |=
                    !set_numeric<double>(optarg, command, m_target_filter.maf);
            else if (command == "memory")
                error |= !set_memory(optarg);
            else if (command == "missing")
                error |= !set_missing(optarg);
            else if (command == "model")
                error |= !set_model(optarg);
            else if (command == "num-auto")
            {
                error |= !set_numeric<int>(optarg, "num-auto",
                                           m_target.num_autosome);
                m_reference.num_autosome = m_target.num_autosome;
            }
            else if (command == "pearson")
            {
                error = true;
                m_error_message.append("Error: This feature is "
                                       "disabled until further notice.\n");
            }
            else if (command == "perm")
            {
                error |= !set_numeric<size_t>(optarg, command,
                                              m_perm_info.num_permutation);
                m_perm_info.run_perm = true;
            }
            else if (command == "proxy")
                error |=
                    !set_numeric<double>(optarg, command, m_clump_info.proxy,
                                         m_clump_info.use_proxy);
            else if (command == "remove")
                set_string(optarg, command, m_target.remove);
            else if (command == "score")
                error |= !set_score(optarg);
            else if (command == "set-perm")
            {
                error |= !set_numeric<size_t>(optarg, command,
                                              m_perm_info.num_permutation);
                m_perm_info.run_set_perm = true;
            }
            else if (command == "snp")
                set_string(optarg, command, +BASE_INDEX::RS);
            else if (command == "snp-set")
            {
                load_string_vector(optarg, command, m_prset.snp);
                m_prset.run = true;
            }
            else if (command == "stat")
                set_string(optarg, command, +BASE_INDEX::STAT);
            else if (command == "target-list")
                set_string(optarg, command, m_target.file_list);
            else if (command == "type")
                set_string(optarg, command, m_target.type);
            else if (command == "wind-3")
                error |= !parse_unit_value(optarg, command, 0, m_prset.wind_3);
            else if (command == "wind-5")
                error |= !parse_unit_value(optarg, command, 0, m_prset.wind_5);
            else if (command.compare("x-range") == 0)
                set_string(optarg, command, m_exclusion_range);
            else
            {
                throw std::runtime_error(
                    "Error: Undefined operator: " + command
                    + ", please use --help for more information!");
            }
            break;
        case 'b': set_string(optarg, "base", m_base_info.file_name); break;
        case 'B':
            load_string_vector(optarg, "bed", m_prset.bed);
            m_prset.run = true;
            break;
        case 'c':
            load_string_vector(optarg, "cov-col", m_pheno_info.cov_colname);
            break;
        case 'C': set_string(optarg, "cov", m_pheno_info.cov_file); break;
        case 'f': set_string(optarg, "pheno", m_pheno_info.pheno_file); break;
        case 'F':
            load_string_vector(optarg, "pheno-col", m_pheno_info.pheno_col);
            break;
        case 'g': set_string(optarg, "gtf", m_prset.gtf, m_prset.run); break;
        case 'i':
            error |=
                !set_numeric<double>(optarg, "interval", m_p_thresholds.inter,
                                     m_p_thresholds.set_threshold);
            break;
        case 'k':
            error |= !load_numeric_vector<double>(optarg, "prevalence",
                                                  m_pheno_info.prevalence);
            break;
        case 'l':
            error |= !set_numeric<double>(optarg, "lower", m_p_thresholds.lower,
                                          m_p_thresholds.set_threshold);
            break;
        case 'L': set_string(optarg, "ld", m_reference.file_name); break;
        case 'm':
            load_string_vector(optarg, "msigdb", m_prset.msigdb);
            m_prset.run = true;
            break;
        case 'n':
            if (strcmp("max", optarg) == 0)
            {
                m_prs_info.thread = max_threads;
                m_parameter_log["thread"] = std::to_string(m_prs_info.thread);
            }
            else
            {
                error |= !set_numeric<int>(optarg, "thread", m_prs_info.thread);
                if (m_prs_info.thread > max_threads)
                { m_prs_info.thread = max_threads; }
                m_parameter_log["thread"] = std::to_string(m_prs_info.thread);
            }
            break;
        case 'o': set_string(optarg, "out", m_out_prefix); break;
        case 'p': set_string(optarg, "pvalue", +BASE_INDEX::P); break;
        case 's':
            error |= !set_numeric<std::random_device::result_type>(
                optarg, "seed", m_perm_info.seed);
            break;
        case 't': set_string(optarg, "target", m_target.file_name); break;
        case 'u':
            error |= !set_numeric<double>(optarg, "upper", m_p_thresholds.upper,
                                          m_p_thresholds.set_threshold);
            break;
        case 'h':
            reporter.simple_report(m_help_message);
            early_termination = true;
            return true;
        case 'v':
            reporter.simple_report(version + " (" + date + ") ");
            early_termination = true;
            return true;
        case '?': throw std::runtime_error("");
        default:
            throw std::runtime_error("Error: Undefined operator, please use "
                                     "--help for more "
                                     "information!");
        }
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    }
    if (m_allow_inter) m_parameter_log["allow-inter"] = "";
    if (m_p_thresholds.fastscore) m_parameter_log["fastscore"] = "";
    if (m_pheno_info.ignore_fid) m_parameter_log["ignore-fid"] = "";
    if (m_include_nonfounders) m_parameter_log["nonfounders"] = "";
    if (m_base_info.is_index) m_parameter_log["index"] = "";
    if (m_keep_ambig) m_parameter_log["keep-ambig"] = "";
    if (m_perm_info.logit_perm) m_parameter_log["logit-perm"] = "";
    if (m_clump_info.no_clump) m_parameter_log["no-clump"] = "";
    if (m_p_thresholds.no_full) m_parameter_log["no-full"] = "";
    if (m_prs_info.no_regress) m_parameter_log["no-regress"] = "";
    if (m_prs_info.non_cumulate) m_parameter_log["non-cumulate"] = "";
    if (m_print_all_scores) m_parameter_log["all-score"] = "";
    if (m_print_snp) m_parameter_log["print-snp"] = "";
    if (m_base_info.is_beta) m_parameter_log["beta"] = "";
    if (m_base_info.is_or) m_parameter_log["or"] = "";
    if (m_target.hard_coded) m_parameter_log["hard"] = "";
    if (m_ultra_aggressive) m_parameter_log["ultra"] = "";
    if (m_prs_info.use_ref_maf) m_parameter_log["use-ref-maf"] = "";
    if (m_user_no_default) m_parameter_log["no-default"] = "";
    return error;
}

bool Commander::validate_command(Reporter& reporter)
{
    bool error = !base_check();
    error |= !clump_check();
    error |= !covariate_check();
    error |= !filter_check();
    error |= !misc_check();
    error |= !ref_check();
    // pheno_check must come after base check because we want the beta / or
    // information for defining the default
    error |= !pheno_check();
    error |= !prset_check();
    error |= !prsice_check();
    error |= !target_check();
    // check all flags
    std::string log_name = m_out_prefix + ".log";
    try
    {
        reporter.initialize(log_name);
    }
    catch (const std::runtime_error& er)
    {
        std::string error_reason = er.what();
        if (error_reason.find('/') != std::string::npos
            || error_reason.find('\\') != std::string::npos)
            reporter.report(error_reason
                            + ". Maybe the path to file does not exists?");
        else
            return true;
    }
    return error;
}

std::string Commander::get_program_header(const std::string& name)
{
    std::chrono::time_point<std::chrono::system_clock> start;
    start = std::chrono::system_clock::now();
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    struct tm* timeinfo;
    char buffer[80];
    timeinfo = localtime(&start_time);
    strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", timeinfo);
    std::string message = "\n\nPRSice " + version + " (" + date + ") \n";
    message.append("https://github.com/choishingwan/PRSice\n");
    message.append("(C) 2016-2020 Shing Wan (Sam) Choi and Paul F. O'Reilly\n");
    message.append("GNU General Public License v3\n\n");
    message.append("If you use PRSice in any published work, please cite:\n");
    message.append("Choi SW, O'Reilly PF.\n");
    message.append(
        "PRSice-2: Polygenic Risk Score Software for Biobank-Scale Data.\n");
    message.append("GigaScience 8, no. 7 (July 1, 2019)\n\n");
    std::string time_str(buffer);
    message.append(time_str + "\n" + name);
    return message;
}
// Default constructor of Command
// Responsible for setting all the default values
// initialize the parameters, then call the
// parameter processing function
// Default destructor of Command, do nothing
Commander::~Commander() {}

// Function to set the help message
// avoid having large chunk of un-foldable code
void Commander::set_help_message()
{
    m_help_message =
        "usage: PRSice [options] <-b base_file> <-t target_file>\n"
        // Base file
        "\nBase File:\n"
        "    --a1                    Column header containing allele 1 "
        "(effective allele)\n"
        "                            Default: A1\n"
        "    --a2                    Column header containing allele 2 "
        "(non-effective allele)\n"
        "                            Default: A2\n"
        "    --base          | -b    Base association file\n"
        "    --base-info             Base INFO score filtering. Format should "
        "be\n"
        "                            <Column name>:<Threshold>. SNPs with info "
        "\n"
        "                            score less than <Threshold> will be "
        "ignored\n"
        "                            Column name default: INFO\n"
        "                            Threshold default: 0.9\n"
        "    --base-maf              Base MAF filtering. Format should be\n"
        "                            <Column name>:<Threshold>. SNPs with maf\n"
        "                            less than <Threshold> will be ignored. "
        "An\n"
        "                            additional column can also be added "
        "(e.g.\n"
        "                            also filter MAF for cases), using the\n"
        "                            following format:\n"
        "                            <Column name>:<Threshold>,<Column "
        "name>:<Threshold>\n"
        "    --beta                  Whether the test statistic is in the form "
        "of \n"
        "                            BETA or OR. If set, test statistic is "
        "assume\n"
        "                            to be in the form of BETA. Mutually "
        "exclusive\n"
        "                            from --or \n"
        "    --bp                    Column header containing the SNP "
        "coordinate\n"
        "                            Default: BP\n"
        "    --chr                   Column header containing the chromosome\n"
        "                            Default: CHR\n"
        "    --index                 If set, assume the INDEX instead of NAME "
        "for\n"
        "                            the corresponding columns are provided. "
        "Index\n"
        "                            should be 0-based (start counting from "
        "0)\n"
        "    --no-default            Remove all default options. If set, "
        "PRSice\n"
        "                            will not set any default column name and "
        "you\n"
        "                            must manually provide all required "
        "columns\n"
        "                            (--snp, --stat, --A1, --pvalue)\n"
        "    --or                    Whether the test statistic is in the form "
        "of \n"
        "                            BETA or OR. If set, test statistic is "
        "assume\n"
        "                            to be in the form of OR. Mutually "
        "exclusive \n"
        "                            from --beta \n"
        "    --pvalue        | -p    Column header containing the p-value\n"
        "                            Default: P\n"
        "    --snp                   Column header containing the SNP ID\n"
        "                            Default: SNP\n"
        "    --stat                  Column header containing the summary "
        "statistic\n"
        "                            If --beta is set, default as BETA. "
        "Otherwise,\n"
        "                            will search for OR or BETA from the "
        "header\n"
        "                            of the base file\n"
        // TARGET FILE
        "\nTarget File:\n"
        "    --binary-target         Indicate whether the target phenotype\n"
        "                            is binary or not. Either T or F should "
        "be\n"
        "                            provided where T represent a binary "
        "phenotype.\n"
        "                            For multiple phenotypes, the input should "
        "be\n"
        "                            separated by comma without space. \n"
        "                            Default: T if --beta and F if --beta is "
        "not\n"
        "    --geno                  Filter SNPs based on gentype missingness\n"
        "    --info                  Filter SNPs based on info score. Only "
        "used\n"
        "                            for imputed target\n"
        "    --keep                  File containing the sample(s) to be "
        "extracted from\n"
        "                            the target file. First column should be "
        "FID and\n"
        "                            the second column should be IID. If "
        "--ignore-fid is\n"
        "                            set, first column should be IID\n"
        "                            Mutually exclusive from --remove\n"
        "    --maf                   Filter SNPs based on minor allele "
        "frequency (MAF)\n"
        "    --nonfounders           Keep the nonfounders in the analysis\n"
        "                            Note: They will still be excluded from LD "
        "calculation\n"
        "    --pheno         | -f    Phenotype file containing the "
        "phenotype(s).\n"
        "                            First column must be FID of the samples "
        "and\n"
        "                            the second column must be IID of the "
        "samples.\n"
        "                            When --ignore-fid is set, first column "
        "must\n"
        "                            be the IID of the samples.\n"
        "                            Must contain a header if --pheno-col is\n"
        "                            specified\n"
        "    --pheno-col     | -F    Headers of phenotypes to be included from "
        "the\n"
        "                            phenotype file\n"
        "    --prevalence    | -k    Prevalence of all binary trait. If "
        "provided\n"
        "                            will adjust the ascertainment bias of the "
        "R2.\n"
        "                            Note that when multiple binary trait is "
        "found,\n"
        "                            prevalence information must be provided "
        "for\n"
        "                            all of them\n"
        "    --remove                File containing the sample(s) to be "
        "removed from\n"
        "                            the target file. First column should be "
        "FID and\n"
        "                            the second column should be IID. If "
        "--ignore-fid is\n"
        "                            set, first column should be IID\n"
        "                            Mutually exclusive from --keep\n"
        "    --target        | -t    Target genotype file. Currently support\n"
        "                            both BGEN and binary PLINK format. For \n"
        "                            multiple chromosome input, simply "
        "substitute\n"
        "                            the chromosome number with #. PRSice "
        "will\n"
        "                            automatically replace # with 1-22\n"
        "                            For binary plink format, you can also "
        "specify\n"
        "                            a seperate fam file by <prefix>,<fam "
        "file>\n"
        "    --target-list           File containing prefix of target "
        "genotype\n"
        "                            files. Similar to --target but allow more "
        "\n"
        "                            flexibility. Do not support external fam "
        "file\n"
        "                            at the moment\n"
        "    --type                  File type of the target file. Support bed "
        "\n"
        "                            (binary plink) and bgen format. Default: "
        "bed\n"
        // dosage
        "\nDosage:\n"
        "    --allow-inter           Allow the generate of intermediate file. "
        "This will\n"
        "                            speed up PRSice when using dosage data as "
        "clumping\n"
        "                            reference and for hard coding PRS "
        "calculation\n"
        "    --dose-thres            Translate any SNPs with highest genotype "
        "probability\n"
        "                            less than this threshold to missing call\n"
        "    --hard-thres            A hardcall is saved when the distance to "
        "the nearest\n"
        "                            hardcall is less than the hardcall "
        "threshold.\n"
        "                            Otherwise a missing code is saved\n"
        "                            Default is: "
        + misc::to_string(m_target_filter.hard_threshold)
        + "\n"
          "    --hard                  Use hard coding instead of dosage for "
          "PRS construction.\n"
          "                            Default is to use dosage instead of "
          "hard coding\n"
          // clumping
          "\nClumping:\n"
          "    --clump-kb              The distance for clumping in kb\n"
          "                            Default: "
        + misc::to_string(m_clump_info.distance / 1000)
        + "\n"
          "    --clump-r2              The R2 threshold for clumping\n"
          "                            Default: "
        + misc::to_string(m_clump_info.r2)
        + " (1mb for PRSet)\n"
          "    --clump-p               The p-value threshold use for "
          "clumping.\n"
          "                            Default: "
        + misc::to_string(m_clump_info.pvalue)
        + "\n"
          "    --ld            | -L    LD reference file. Use for LD "
          "calculation. If not\n"
          "                            provided, will use the post-filtered "
          "target genotype\n"
          "                            for LD calculation. Support multiple "
          "chromosome input\n"
          "                            Please see --target for more "
          "information\n"
          "    --ld-dose-thres         Translate any SNPs with highest "
          "genotype probability\n"
          "                            less than this threshold to missing "
          "call\n"
          "    --ld-geno               Filter SNPs based on genotype "
          "missingness\n"
          "    --ld-hard-thres         A hardcall is saved when the distance "
          "to the nearest\n"
          "                            hardcall is less than the hardcall "
          "threshold.\n"
          "                            Otherwise a missing code is saved\n"
          "                            Default is: "
        + misc::to_string(m_ref_filter.hard_threshold)
        + "\n"
          "    --ld-info               Filter SNPs based on info score. Only "
          "used\n"
          "                            for imputed LD reference\n"
          "    --ld-keep               File containing the sample(s) to be "
          "extracted from\n"
          "                            the LD reference file. First column "
          "should be FID and\n"
          "                            the second column should be IID. If "
          "--ignore-fid is\n"
          "                            set, first column should be IID\n"
          "                            Mutually exclusive from --ld-remove\n"
          "                            No effect if --ld was not provided\n"
          "    --ld-list               File containing prefix of LD reference "
          "files.\n"
          "                            Similar to --ld but allow more \n"
          "                            flexibility. Do not support external "
          "fam file\n"
          "                            at the moment\n"
          "    --ld-maf                Filter SNPs based on minor allele "
          "frequency\n"
          "    --ld-remove             File containing the sample(s) to be "
          "removed from\n"
          "                            the LD reference file. First column "
          "should be FID and\n"
          "                            the second column should be IID. If "
          "--ignore-fid is\n"
          "                            set, first column should be IID\n"
          "                            Mutually exclusive from --ld-keep\n"
          "    --ld-type               File type of the LD file. Support bed "
          "(binary plink)\n"
          "                            and bgen format. Default: bed\n"
          "    --no-clump              Stop PRSice from performing clumping\n"
          "    --proxy                 Proxy threshold for index SNP to be "
          "considered\n"
          "                            as part of the region represented by "
          "the clumped\n"
          "                            SNP(s). e.g. --proxy 0.8 means the "
          "index SNP will\n"
          "                            represent region of any clumped SNP(s) "
          "that has a\n"
          "                            R2>=0.8 even if the index SNP does not "
          "physically\n"
          "                            locate within the region\n"
          // Covariates
          "\nCovariate:\n"
          "    --cov           | -C    Covariate file. First column should be "
          "FID and \n"
          "                            the second column should be IID. If "
          "--ignore-fid\n"
          "                            is set, first column should be IID\n"
          "    --cov-col       | -c    Header of covariates. If not provided, "
          "will use\n"
          "                            all variables in the covariate file. By "
          "adding\n"
          "                            @ in front of the string, any numbers "
          "within [\n"
          "                            and ] will be parsed. E.g. @PC[1-3] "
          "will be\n"
          "                            read as PC1,PC2,PC3. Discontinuous "
          "input are also\n"
          "                            supported: @cov[1.3-5] will be parsed "
          "as \n"
          "                            cov1,cov3,cov4,cov5\n"
          "    --cov-factor            Header of categorical covariate(s). "
          "Dummy variable\n"
          "                            will be automatically generated. Any "
          "items in\n"
          "                            --cov-factor must also be found in "
          "--cov-col\n"
          "                            Also accept continuous input (start "
          "with @).\n"
          // PRSice
          "\nP-value Thresholding:\n"
          "    --bar-levels            Level of barchart to be plotted. When "
          "--fastscore\n"
          "                            is set, PRSice will only calculate the "
          "PRS for \n"
          "                            threshold within the bar level. Levels "
          "should be\n"
          "                            comma separated without space\n"
          "    --fastscore             Only calculate threshold stated in "
          "--bar-levels\n"
          "    --no-full               By default, PRSice will include the "
          "full model, \n"
          "                            i.e. p-value threshold = 1. Setting "
          "this flag will\n"
          "                            disable that behaviour\n"
          "    --interval      | -i    The step size of the threshold. "
          "Default: "
        + misc::to_string(m_p_thresholds.inter)
        + "\n"
          "    --lower         | -l    The starting p-value threshold. "
          "Default: "
        + misc::to_string(m_p_thresholds.lower)
        + "\n"
          "    --model                 Genetic model use for regression. The "
          "genetic\n"
          "                            encoding is based on the base data "
          "where the\n"
          "                            encoding represent number of the coding "
          "allele\n"
          "                            Available models include:\n"
          "                            add - Additive model, code as 0/1/2 "
          "(default)\n"
          "                            dom - Dominant model, code as 0/1/1\n"
          "                            rec - Recessive model, code as 0/0/1\n"
          "                            het - Heterozygous only model, code as "
          "0/1/0\n"
          "    --missing               Method to handle missing genotypes. By "
          "default, \n"
          "                            final scores are averages of valid "
          "per-allele \n"
          "                            scores with missing genotypes "
          "contribute an amount\n"
          "                            proportional to imputed allele "
          "frequency. To throw\n"
          "                            out missing observations instead "
          "(decreasing the\n"
          "                            denominator in the final average when "
          "this happens),\n"
          "                            use the 'SET_ZERO' modifier. "
          "Alternatively,\n"
          "                            you can use the 'CENTER' modifier to "
          "shift all scores\n"
          "                            to mean zero. \n"
          "    --no-regress            Do not perform the regression analysis "
          "and simply\n"
          "                            output all PRS.\n"
          "    --score                 Method to calculate the polygenic "
          "score.\n"
          "                            Available methods include:\n"
          "                            avg     - Take the average effect size "
          "(default)\n"
          "                            std     - Standardize the effect size \n"
          "                            con-std - Standardize the effect size "
          "using mean \n"
          "                                      and sd derived from control "
          "samples\n"
          "                            sum     - Direct summation of the "
          "effect size \n"
          "    --upper         | -u    The final p-value threshold. Default: "
        + misc::to_string(m_p_thresholds.upper)
        + "\n"
          "\nPRSet:\n"
          "    --background            String to indicate a background file. "
          "This string\n"
          "                            should have the format of Name:Type "
          "where type can be\n"
          "                            bed   - 0-based range with 3 column. "
          "Chr Start End\n"
          "                            range - 1-based range with 3 column. "
          "Chr Start End\n"
          "                            gene  - A file contain a column of gene "
          "name\n"
          "    --bed           | -B    Bed file containing the selected "
          "regions.\n"
          "                            Name of bed file will be used as the "
          "region\n"
          "                            identifier. WARNING: Bed file is "
          "0-based\n"
          "    --feature               Feature(s) to be included from the gtf "
          "file.\n"
          "                            Default: exon,CDS,gene,protein_coding.\n"
          "    --full-back             Use the whole genome as background for "
          "competitive\n"
          "                            p-value calculation\n"
          "    --gtf           | -g    GTF file containing gene boundaries. "
          "Required\n"
          "                            when --msigdb is used\n"
          "    --msigdb        | -m    MSIGDB file containing the pathway "
          "information.\n"
          "                            Require the gtf file\n"
          "    --snp-set               Provide a SNP set file containing the "
          "snp set(s).\n"
          "                            Two different file format is allowed:\n"
          "                            SNP list format - A file containing a "
          "single\n"
          "                                              column of SNP ID. "
          "Name of the\n"
          "                                              set will be the file "
          "name or\n"
          "                                              can be provided using "
          "\n"
          "                                              --snp-set File:Name\n"
          "                            MSigDB format   - Each row represent a "
          "single SNP \n"
          "                                              set with the first "
          "column \n"
          "                                              containing the name "
          "of the SNP\n"
          "                                              set.\n"
          "    --wind-3                Add N base(s) to the 3' region of each "
          "feature(s) \n"
          "    --wind-5                Add N base(s) to the 5' region of each "
          "feature(s) \n"
          // Misc
          "\nMisc:\n"
          "    --all-score             Output PRS for ALL threshold. WARNING: "
          "This\n"
          "                            will generate a huge file\n"
          "    --chr-id                Try to construct an RS ID for SNP based "
          "on its\n"
          "                            chromosome, coordinate, effective "
          "allele and \n"
          "                            non-effective allele.\n"
          "                            e.g. c:L-aBd is translated to: \n"
          "                            "
          "<chr>:<coordinate>-<effective><noneffective>d\n"
          "                            This is always true for target file, "
          "whereas for\n"
          "                            base file, this is only used if the RS "
          "ID \n"
          "                            wasn't provided\n"
          "    --exclude               File contains SNPs to be excluded from "
          "the\n"
          "                            analysis\n"
          "    --extract               File contains SNPs to be included in "
          "the \n"
          "                            analysis\n"
          "    --id-delim              This parameter causes sample IDs to be "
          "parsed as\n"
          "                            <FID><delimiter><IID>; the default "
          "delimiter\n"
          "                            is '_'. \n"
          "    --ignore-fid            Ignore FID for all input. When this is "
          "set,\n"
          "                            first column of all file will be assume "
          "to\n"
          "                            be IID instead of FID\n"
          "    --keep-ambig            Keep ambiguous SNPs. Only use this "
          "option\n"
          "                            if you are certain that the base and "
          "target\n"
          "                            has the same A1 and A2 alleles\n"
          "    --logit-perm            When performing permutation, still use "
          "logistic\n"
          "                            regression instead of linear "
          "regression. This\n"
          "                            will substantially slow down PRSice\n"
          "    --memory                Maximum memory usage allowed (in Mb). "
          "PRSice will try\n"
          "                            its best to honor this setting\n"
          "    --non-cumulate          Calculate non-cumulative PRS. PRS will "
          "be reset\n"
          "                            to 0 for each new P-value threshold "
          "instead of\n"
          "                            adding up\n"
          "    --out           | -o    Prefix for all file output\n"
          "    --perm                  Number of permutation to perform. This "
          "swill\n"
          "                            generate the empirical p-value. "
          "Recommend to\n"
          "                            use value larger than 10,000\n"
          "    --print-snp             Print all SNPs that remains in the "
          "analysis \n"
          "                            after clumping is performed. For PRSet, "
          "Y \n"
          "                            indicate the SNPs falls within the gene "
          "set \n"
          "                            of interest and N otherwise. If only "
          "PRSice \n"
          "                            is performed, a single \"gene set\" "
          "called \n"
          "                            \"Base\" will be presented with all "
          "entries\n"
          "                            marked as Y\n"
          "    --seed          | -s    Seed used for permutation. If not "
          "provided,\n"
          "                            system time will be used as seed. When "
          "same\n"
          "                            seed and same input is provided, same "
          "result\n"
          "                            can be generated\n"
          "    --thread        | -n    Number of thread use\n"
          "    --use-ref-maf           When specified, missingness imputation "
          "will be\n"
          "                            performed based on the reference "
          "samples\n"
          "    --ultra                 Ultra aggressive memory usage. When "
          "this is enabled\n"
          "                            PRSice and PRSet will try to load all "
          "genotypes into\n"
          "                            memory after clumping is performed. "
          "This should\n"
          "                            drastically speed up PRSice and PRSet "
          "at the expense\n"
          "                            of higher memory consumption.\n"
          "                            Has no effect for dosage score\n"
          "    --x-range               Range of SNPs to be excluded from the "
          "whole\n"
          "                            analysis. It can either be a single bed "
          "file\n"
          "                            or a comma seperated list of range. "
          "Range must\n"
          "                            be in the format of chr:start-end or "
          "chr:coordinate\n"
          "    --help          | -h    Display this help message\n";
}

std::vector<std::string> get_base_header(const std::string& file)
{
    if (file.empty())
    { throw std::runtime_error("Error: You must provide a base file\n"); }
    // get input header
    std::string header;
    const bool is_gz = misc::is_gz_file(file);
    if (is_gz)
    {
        GZSTREAM_NAMESPACE::igzstream in(file.c_str());
        if (!in.good())
        {
            throw std::runtime_error(
                "Error: Cannot open base file (gz) to read!\n");
        }
        std::getline(in, header);
        in.close();
    }
    else
    {
        std::ifstream base_test;
        base_test.open(file.c_str());
        if (!base_test.is_open())
        {
            throw std::runtime_error("Error: Cannot open base file: " + file
                                     + " to read!\n"
                                     + std::string(strerror(errno)) + "\n");
        }
        std::getline(base_test, header);
        base_test.close();
    }
    misc::trim(header);
    return misc::split(header);
}

// return false if invalid or if not found
bool Commander::get_statistic_column(
    const std::vector<std::string>& column_names)
{
    // don't allow both OR and BETA to be set
    if (m_base_info.is_or && m_base_info.is_beta) return false;
    const bool or_has_set = m_base_info.is_or;
    const bool beta_has_set = m_base_info.is_beta;

    // go through file and look for either OR or BETA
    m_base_info.column_name[+BASE_INDEX::STAT] = "OR";
    bool or_found = in_file(column_names, +BASE_INDEX::STAT, "Warning", false,
                            false, false);
    m_base_info.column_name[+BASE_INDEX::STAT] = "BETA";
    bool beta_found = in_file(column_names, +BASE_INDEX::STAT, "Warning", false,
                              false, false);
    if ((beta_has_set || or_has_set) && (beta_found || or_found))
    {
        if (beta_has_set && beta_found)
        {
            // don't need to do anything as beta_found come last
        }
        else if (or_has_set && or_found)
        {
            // we need to recall or to set the index
            m_base_info.column_name[+BASE_INDEX::STAT] = "OR";
            in_file(column_names, +BASE_INDEX::STAT, "Warning", false, false,
                    false);
        }
        else if (or_has_set && beta_found)
        {
        }
        else if (beta_has_set && or_found)
        {
            m_base_info.column_name[+BASE_INDEX::STAT] = "OR";
            in_file(column_names, +BASE_INDEX::STAT, "Warning", false, false,
                    false);
        }
        m_base_info.has_column[+BASE_INDEX::STAT] = true;
        return true;
    }
    else if (or_found && beta_found)
    {
        m_error_message.append("Error: Both OR and BETA "
                               "found in base file! We cannot determine "
                               "which statistic to use, please provide it "
                               "through --stat\n");
        return false;
    }
    else if (or_found || beta_found)
    {
        // don't need to update column_index as column_index is only set
        // when we find the item in the base file so this will only be set to
        // the correct one
        m_base_info.is_or = or_found;
        m_base_info.is_beta = beta_found;
        m_base_info.has_column[+BASE_INDEX::STAT] = true;
        return true;
    }
    else
    {
        // cannot find either
        m_error_message.append("Error: No statistic column in file!\n");
        return false;
    }
}

bool Commander::get_statistic_flag()
{
    std::string stat_temp = m_base_info.column_name[+BASE_INDEX::STAT];
    // guess flag based on stat provided
    misc::to_upper(stat_temp);
    if (stat_temp == "OR")
    {
        m_base_info.is_or = true;
        m_parameter_log["or"] = "";
    }
    else if (stat_temp == "BETA")
    {
        m_base_info.is_beta = true;
        m_parameter_log["beta"] = "";
    }
    else
    {
        m_error_message.append(
            "Error: Cannot determine if statistic is BETA or OR: "
            + m_base_info.column_name[+BASE_INDEX::STAT] + "\n");
        return false;
    }
    return true;
}

bool Commander::base_check()
{
    m_ran_base_check = true;
    std::vector<std::string> column_names =
        get_base_header(m_base_info.file_name);
    return base_column_check(column_names);
}

bool Commander::base_column_check(std::vector<std::string>& column_names)
{

    bool error = false;
    for (auto&& c : column_names) { misc::trim(c); }
    if (m_base_info.is_index)
    {
        // can't do much but to check the boundary
        for (size_t i = 0; i < column_names.size(); ++i)
        { column_names[i] = misc::to_string(i); }
    }
    if (!in_file(column_names, +BASE_INDEX::CHR, "Warning", m_user_no_default))
    { m_parameter_log.erase("chr"); }
    if (!in_file(column_names, +BASE_INDEX::NONEFFECT, "Warning",
                 m_user_no_default))
    { m_parameter_log.erase("A2"); }
    if (!in_file(column_names, +BASE_INDEX::BP, "Warning", m_user_no_default))
    { m_parameter_log.erase("bp"); }

    if (!in_file(column_names, +BASE_INDEX::EFFECT, "Error", m_user_no_default))
    {
        error = true;
        m_error_message.append(
            "Error: Column for the effective allele must be provided!\n");
    }
    if (!in_file(column_names, +BASE_INDEX::RS, "Warning", m_user_no_default))
    {
        if (m_chr_id_formula.empty())
        {
            error = true;
            m_error_message.append(
                "Error: Column for the SNP ID must be "
                "provided when --chr-id was not provided!\n");
        }
    }
    if (!in_file(column_names, +BASE_INDEX::P, "Error", m_user_no_default))
    {
        error = true;
        m_error_message.append(
            "Error: Column for the P-value must be provided!\n");
    }
    error |= !set_base_info_threshold(column_names);
    error |= !set_base_maf_filter(column_names);
    // now process the statistic column
    if (m_base_info.is_or && m_base_info.is_beta)
    {
        m_error_message.append("Error: Statistic cannot be both OR and beta\n");
        error = true;
    }
    bool user_provided_stat = m_base_info.has_column[+BASE_INDEX::STAT];
    bool has_col = in_file(column_names, +BASE_INDEX::STAT, "Error",
                           m_user_no_default, true, m_user_no_default);
    if (!has_col && user_provided_stat)
    {
        error = true;
        m_error_message.append("Error: Statistic column not found in file!\n");
    }
    else if (!m_user_no_default && !has_col)
    {
        // if allow default
        error |= !get_statistic_column(column_names);
    }
    // Statistic is ok, but beta and or not provided
    // use has_column vector instead of has_col as get_statistic_column might
    // have found the state column?
    if (m_base_info.has_column[+BASE_INDEX::STAT])
    {
        // flag not provided, need to guess
        if (!(m_base_info.is_or || m_base_info.is_beta))
        { error |= !get_statistic_flag(); }
    }
    m_base_info.column_index[+BASE_INDEX::MAX] = *max_element(
        m_base_info.column_index.begin(), m_base_info.column_index.end());
    add_base_to_log();
    return !error;
}

void Commander::add_base_to_log()
{
    for (size_t i = 0; i < +BASE_INDEX::MAX; ++i)
    {
        std::string flag = "";
        switch (i)
        {
        case +BASE_INDEX::P:
            if (m_base_info.has_column[i]) { flag = "pvalue"; }
            break;
        case +BASE_INDEX::BP:
            if (m_base_info.has_column[i]) { flag = "bp"; }
            break;
        case +BASE_INDEX::RS:
            if (m_base_info.has_column[i]) { flag = "snp"; }
            break;
        case +BASE_INDEX::CHR:
            if (m_base_info.has_column[i]) { flag = "chr"; }
            break;
        case +BASE_INDEX::INFO:
            if (m_base_info.has_column[i]) { flag = "base-info"; }
            break;
        case +BASE_INDEX::STAT:
            if (m_base_info.has_column[i]) { flag = "stat"; }
            break;
        case +BASE_INDEX::EFFECT:
            if (m_base_info.has_column[i]) { flag = "a1"; }
            break;
        case +BASE_INDEX::NONEFFECT:
            if (m_base_info.has_column[i]) { flag = "a2"; }
            break;
        }
        if (!flag.empty()) m_parameter_log[flag] = m_base_info.column_name[i];
    }
}
bool Commander::clump_check()
{
    bool error = false;
    if (m_clump_info.no_clump) return true;
    if (m_clump_info.use_proxy
        && !misc::within_bound<double>(m_clump_info.proxy, 0.0, 1.0))
    {
        error = true;
        m_error_message.append(
            "Error: Proxy threshold must be within 0 and 1!\n");
    }
    if (!misc::within_bound<double>(m_clump_info.pvalue, 0.0, 1.0))
    {
        error = true;
        m_error_message.append(
            "Error: P-value threshold must be within 0 and 1!\n");
    }
    if (!misc::within_bound<double>(m_clump_info.r2, 0.0, 1.0))
    {
        error = true;
        m_error_message.append("Error: R2 threshold must be within 0 and 1!\n");
    }
    m_parameter_log["clump-r2"] = std::to_string(m_clump_info.r2);
    m_parameter_log["clump-p"] = std::to_string(m_clump_info.pvalue);
    // we divided by 1000 here to make sure it is in KB (our preferred
    // format)
    if (!m_clump_info.provided_distance && m_prset.run)
    { m_clump_info.distance = 1000000; }
    m_parameter_log["clump-kb"] =
        std::to_string(m_clump_info.distance / 1000) + "kb";
    return !error;
}

bool Commander::need_target_as_reference()
{
    // default type is bed, won't be empty
    const bool not_provided = m_reference.type == "bed"
                              && m_reference.file_list.empty()
                              && m_reference.file_name.empty();
    if (!not_provided) return false;
    bool require_ref = false;
    // check maf
    require_ref |=
        !(misc::logically_equal(m_ref_filter.maf, 0.0)
          || misc::logically_equal(m_ref_filter.maf, m_target_filter.maf));
    // check geno
    require_ref |=
        !(misc::logically_equal(m_ref_filter.geno, 1.0)
          || misc::logically_equal(m_ref_filter.geno, m_target_filter.geno));
    // check sample selection
    require_ref |=
        !(m_reference.keep.empty() || m_reference.keep == m_target.keep);
    require_ref |=
        !(m_reference.remove.empty() || m_reference.remove == m_target.remove);
    // still don't need reference, but target isn't plink
    if (!require_ref && m_target.type != "bed")
    {
        // check info score
        require_ref |= !(misc::logically_equal(m_ref_filter.info_score, 0.0)
                         || misc::logically_equal(m_ref_filter.info_score,
                                                  m_target_filter.info_score));
        // check dose threshold
        require_ref |=
            !(misc::logically_equal(m_ref_filter.dose_threshold, 0.0)
              || misc::logically_equal(m_ref_filter.dose_threshold,
                                       m_target_filter.dose_threshold));
        // check hard threshold
        require_ref |=
            !(misc::logically_equal(m_ref_filter.hard_threshold, 0.1)
              || misc::logically_equal(m_ref_filter.hard_threshold,
                                       m_target_filter.hard_threshold));
        return require_ref;
    }
    else
        return require_ref;
}


bool Commander::ref_check()
{
    bool error = false;
    if (!m_reference.keep.empty() && !m_reference.remove.empty())
    {
        error = true;
        m_error_message.append("Error: Can only use either --keep or "
                               "--remove but not both\n");
    }

    if (need_target_as_reference())
    {
        // assign name over
        m_reference.file_list = m_target.file_list;
        if (!m_reference.file_list.empty())
        { m_parameter_log["ld-list"] = m_reference.file_list; }
        m_reference.file_name = m_target.file_name;
        if (!m_reference.file_name.empty())
        { m_parameter_log["ld"] = m_reference.file_name; }
        m_reference.type = m_target.type;
        m_parameter_log["ld-type"] = m_reference.type;
    }
    if (!m_reference.type.empty())
    {
        if (std::find(supported_types.begin(), supported_types.end(),
                      m_reference.type)
            == supported_types.end())
        {
            error = true;
            m_error_message.append(
                "Error: Unsupported LD format: " + m_reference.type + "\n");
        }
    }
    if (!misc::within_bound<double>(m_ref_filter.geno, 0.0, 1.0))
    {
        error = true;
        m_error_message.append("Error: LD genotype missingness threshold "
                               "must be larger than 0 and smaller than 1!\n");
    }
    // if the reference panel is bgen, or reference panel not provided
    // but the target is bgen then we will like to enforce hard
    // thresholding to the files for LD calculation
    if (!m_reference.file_name.empty() && !m_reference.file_list.empty())
    {
        error = true;
        m_error_message.append("Error: You can only use --ld or --ld-list "
                               "but not both\n");
    }

    if (m_reference.type == "bgen"
        || (m_reference.file_name.empty() && !m_reference.file_list.empty()
            && m_target.type == "bgen"))
    {
        if (!misc::within_bound<double>(m_ref_filter.hard_threshold, 0.0, 1.0))
        {
            error = true;
            m_error_message.append("Error: LD hard threshold must be larger "
                                   "than 0 and smaller than 1!\n");
        }
        if (!misc::within_bound<double>(m_ref_filter.dose_threshold, 0.0, 1.0))
        {
            error = true;
            m_error_message.append("Error: LD dosage threshold must be larger "
                                   "than 0 and smaller than 1!\n");
        }
        else if (!m_reference.file_name.empty()
                 || m_reference.file_name.empty())
        {
            m_parameter_log["ld-hard-thres"] =
                std::to_string(m_ref_filter.hard_threshold);
            m_parameter_log["ld-dose-thres"] =
                std::to_string(m_ref_filter.dose_threshold);
        }
        else
        {
            m_parameter_log["hard-thres"] =
                std::to_string(m_target_filter.hard_threshold);
            m_parameter_log["dose-thres"] =
                std::to_string(m_target_filter.dose_threshold);
        }
        if (!misc::within_bound<double>(m_ref_filter.info_score, 0.0, 1.0))
        {
            error = true;
            m_error_message.append("Error: LD INFO score threshold must be "
                                   "larger than 0 and smaller than 1!\n");
        }
    }
    if (!misc::within_bound<double>(m_ref_filter.maf, 0.0, 1.0))
    {
        error = true;
        m_error_message.append("Error: LD MAF threshold must be larger than "
                               "0 and smaller than 1!\n");
    }
    return !error;
}

size_t Commander::find_first_end(const std::string_view& cov, const size_t idx)
{
    if (cov.at(idx) != '[')
    {
        throw std::runtime_error(
            "Error: Invalid format. Expect string to start with [");
    }
    for (size_t i = idx + 1; i < cov.length(); ++i)
    {
        if (cov.at(i) == ']') return i;
        if (cov.at(i) == '[')
            throw std::runtime_error(
                "Error: Invalid format, we don't allow embedded range");
    }
    throw std::runtime_error(
        "Error: Invalid format, [ must accompany with a ]");
}
std::vector<size_t> Commander::parse_range(const std::string_view& cov)
{
    std::vector<size_t> res;
    if (cov.at(0) == '-' || cov.find("--") != std::string::npos)
    { throw std::runtime_error("Error: Do not accept negative ranges"); }
    std::vector<std::string_view> token = misc::tokenize(cov, "-");
    // from_chars will be faster, but less robust (1.5 will be converted to 1)
    if (token.size() == 1) { res = {misc::convert<size_t>(std::string(cov))}; }
    else
    {
        size_t start, end;
        start = misc::convert<size_t>(std::string(token.front()));
        end = misc::convert<size_t>(std::string(token.back()));
        if (start > end) { std::swap(start, end); }
        res.resize(end - start + 1, start);
        std::iota(res.begin(), res.end(), start);
    }
    return res;
}
std::vector<size_t> Commander::get_range(const std::string_view& cov,
                                         const size_t start, const size_t end)
{
    // need to remove []
    if (end >= cov.length() || start >= end)
    { throw std::runtime_error("Error: Wrong start and end format"); }
    if (!(cov.at(start) == '[' && cov.at(end) == ']'))
    {
        throw std::runtime_error("Error: Invalid input. Expect something "
                                 "starts with [ and end with ]");
    }
    std::vector<std::string_view> token =
        misc::tokenize(cov.substr(start + 1, end - start - 1), ".");
    std::vector<size_t> results, tmp;
    for (auto&& value : token)
    {
        // now try to account for -
        tmp = parse_range(value);
        results.insert(results.end(), tmp.begin(), tmp.end());
    }
    std::sort(results.begin(), results.end());
    results.erase(std::unique(results.begin(), results.end()), results.end());
    return results;
}

void Commander::update_covariate_range(const std::vector<size_t>& range,
                                       std::vector<std::string>& res)
{
    if (range.empty())
    {
        throw std::runtime_error(
            "Error: Invalid input. Something is wrong with Sam");
    }
    if (res.empty())
    {
        res.reserve(range.size());
        for (auto&& value : range) { res.push_back(std::to_string(value)); }
    }
    else
    {
        // there are content in res, so we will duplicate it w.r.t number in
        // range
        std::vector<std::string> tmp;
        tmp.reserve(res.size() * range.size());
        for (auto&& r : res)
        {
            for (auto&& value : range)
            { tmp.push_back(r + std::to_string(value)); }
        }
        res.clear();
        res = tmp;
    }
}

std::vector<std::string>
Commander::transform_covariate(const std::string& cov_in)
{
    // do not allow embedded range
    if (cov_in.empty() || cov_in.at(0) != '@')
    { return std::vector<std::string> {cov_in}; }
    // remove first @
    std::string cov = cov_in;
    cov.erase(0, 1);
    std::vector<std::string> result;
    std::vector<size_t> range;
    std::size_t prev = 0, pos;
    while ((pos = cov.find_first_of("[", prev)) != std::string::npos)
    {
        if (pos > prev)
        {
            std::string_view substring = cov.substr(prev, pos - prev);
            if (result.empty())
                result.emplace_back(substring);
            else
            {
                for (size_t i = 0; i < result.size(); ++i)
                { result[i] = result[i].append(substring); }
            }
            size_t end = find_first_end(cov, pos);
            update_covariate_range(get_range(cov, pos, end), result);
            pos = end;
        }
        prev = pos + 1;
    }
    if (prev < cov.length())
    {
        if (result.empty())
            result.emplace_back(cov.substr(prev, std::string::npos));
        else
        {
            std::string_view substring = cov.substr(prev, std::string::npos);
            for (auto&& c : result) { c.append(substring); }
        }
    }
    return result;
}
/*
std::vector<std::string>
Commander::transform_covariate(const std::string& cov_in)
{
    if (cov_in.empty() || cov_in.at(0) != '@')
        return std::vector<std::string> {cov_in};
    std::vector<std::string> final_covariates;
    std::vector<std::string> open;
    std::vector<std::string> close;
    std::vector<std::string> info;
    std::vector<std::string> individual;
    std::vector<std::string> range;
    std::string cov = cov_in;
    std::string prefix, suffix;
    // Remove the first @
    cov.erase(0, 1);
    // find the start of range by identifying [
    open = misc::split(cov, "[");
    prefix = open.front();
    if (open.size() != 2)
    {
        throw std::runtime_error("Error: Currently only support simple "
                                 "list (i.e. with one set of [])");
    }
    // check for the second set, we do allow XXX[123]YYY
    close = misc::split(open.back(), "]");
    if (close.size() == 2)
        suffix = close.back();
    else
        suffix = "";
    if (close.size() > 2)
    {
        throw std::runtime_error("Error: Currently only support simple "
                                 "list (i.e. with one set of [])");
    }
    individual = misc::split(close.front(), ".");
    for (auto&& ind : individual)
    {
        if (ind.find("-") != std::string::npos)
        {
            // This is list
            range = misc::split(ind, "-");
            if (range.size() != 2)
            {
                throw std::runtime_error("Error: Invalid range format, range "
                                         "must be in the form of start-end");
            }
            try
            {
                size_t start = misc::string_to_size_t(range[0].c_str());
                size_t end = misc::string_to_size_t(range[1].c_str());
                if (start > end) { std::swap(start, end); }
                for (; start <= end; ++start)
                {
                    final_covariates.push_back(prefix + misc::to_string(start)
                                               + suffix);
                }
            }
            catch (const std::runtime_error&)
            {
                std::string error_message = "Error: Invalid parameter: " + ind
                                            + ", only allow integer!";
                throw std::runtime_error(error_message);
            }
        }
        else
        {
            // this is single value
            try
            {
                misc::convert<int>(ind);
                final_covariates.push_back(prefix + ind + suffix);
            }
            catch (const std::runtime_error&)
            {
                std::string error_message = "Error: Invalid parameter: " + ind
                                            + ", only allow integer!";
                throw std::runtime_error(error_message);
            }
        }
    }
    return final_covariates;
}*/

std::unordered_set<std::string> Commander::get_cov_names()
{
    // first, transform all the covariates
    // the actual column name to be included (after parsing)
    std::unordered_set<std::string> included;
    // contain all the input within the --cov-col
    std::unordered_set<std::string> ori_input;
    // vector contain the transformed column names
    std::vector<std::string> transformed_cov;
    for (auto cov : m_pheno_info.cov_colname)
    {
        if (cov.empty()) continue;
        ori_input.insert(cov);
        transformed_cov = transform_covariate(cov);
        for (auto&& trans : transformed_cov) { included.insert(trans); }
    }
    return included;
}

std::vector<std::string> Commander::get_covariate_header()
{
    std::ifstream cov_file;
    cov_file.open(m_pheno_info.cov_file.c_str());
    if (!cov_file.is_open())
    {
        m_error_message.append("Error: Cannot open covariate file: "
                               + m_pheno_info.cov_file + "\n");
        throw std::runtime_error("Cannot open");
    }
    std::string line;
    std::getline(cov_file, line);
    cov_file.close();
    misc::trim(line);
    if (line.empty())
    {
        m_error_message.append(
            "Error: First line of covariate file is empty!\n");
        throw std::runtime_error("Empty line");
    }
    auto cov_header = misc::split(line);
    return cov_header;
}
size_t Commander::find_cov_idx(
    const std::unordered_set<std::string>& included,
    const std::unordered_map<std::string, size_t>& ref_index,
    std::string& missing)
{
    missing = "";
    m_pheno_info.col_index_of_cov.clear();
    size_t valid_cov = 0;
    std::string comma = "";
    for (auto&& cov : included)
    {
        // now for each covariate found in the covariate file, we add their
        // index to the storage
        auto idx = ref_index.find(cov);
        if (idx != ref_index.end())
        {
            m_pheno_info.col_index_of_cov.push_back(idx->second);
            ++valid_cov;
        }
        else
        {
            missing.append(comma + cov);
            comma = ",";
        }
    }
    return valid_cov;
}

void Commander::reorganize_cov_name(const std::vector<std::string>& cov_header)
{
    m_pheno_info.cov_colname.clear();
    std::sort(m_pheno_info.col_index_of_cov.begin(),
              m_pheno_info.col_index_of_cov.end());
    for (auto&& c : m_pheno_info.col_index_of_cov)
    { m_pheno_info.cov_colname.push_back(cov_header[c]); }
}
bool Commander::process_factor_cov(
    const std::unordered_set<std::string>& included,
    const std::unordered_map<std::string, size_t>& ref_index,
    const std::unordered_set<std::string> ori_input)
{
    // now start to process the factor covariates

    std::vector<std::string> transformed_cov;
    for (auto cov : m_pheno_info.factor_cov)
    {
        if (cov.empty()) continue;
        transformed_cov = transform_covariate(cov);
        for (auto&& trans : transformed_cov)
        {
            auto&& ref = ref_index.find(trans);
            if (included.find(trans) != included.end()
                && ref != ref_index.end())
            { m_pheno_info.col_index_of_factor_cov.push_back(ref->second); }
            else if (ori_input.find(trans) == ori_input.end())
            {
                // only complain if transform input isn't found in transformed
                // --cov-col
                // so if @PC[1.3.5] isn't found, and cov-col is @PC[1-10], then
                // we still allow such input
                m_error_message.append("Error: All factor covariates must be "
                                       "found in covariate list. "
                                       + trans
                                       + " not found in covariate list");
                return false;
            }
        }
    }
    std::sort(m_pheno_info.col_index_of_factor_cov.begin(),
              m_pheno_info.col_index_of_factor_cov.end());
    return true;
}
void Commander::prepare_header_cov_check(
    const std::vector<std::string>& cov_header,
    std::unordered_map<std::string, size_t>& ref_index,
    std::unordered_set<std::string>& included)
{
    for (size_t i = 0; i < cov_header.size(); ++i)
    { ref_index[cov_header[i]] = i; }
    if (m_pheno_info.cov_colname.size() == 0)
    {
        for (size_t i = (1 + !m_pheno_info.ignore_fid); i < cov_header.size();
             ++i)
        { included.insert(cov_header[i]); }
    }
}
bool Commander::covariate_check()
{
    // it is valid to have empty covariate file
    if (m_pheno_info.cov_file.empty()) return true;
    bool error = false;
    // now try to read the header of the covariate file
    // first, transform all the covariates
    // the actual column name to be included (after parsing)
    std::unordered_set<std::string> included = get_cov_names();
    std::unordered_set<std::string> ori_input = included;
    try
    {
        auto cov_header = get_covariate_header();
        std::unordered_map<std::string, size_t> ref_index;
        prepare_header_cov_check(cov_header, ref_index, included);
        std::string missing;
        size_t valid_cov = find_cov_idx(included, ref_index, missing);
        if (!missing.empty())
        {
            m_error_message.append(
                "Warning: Covariate(s) missing from file: " + missing + ".\n");
        }
        if (valid_cov == 0)
        {
            error = true;
            m_error_message.append("Error: No valid Covariate!\n");
        }
        reorganize_cov_name(cov_header);
        error |= !process_factor_cov(included, ref_index, ori_input);
    }
    catch (std::runtime_error&)
    {
        return false;
    }
    return !error;
}


bool Commander::filter_check()
{
    bool error = false;
    if (m_target.type != "bgen")
    {
        if (m_target_filter.provided_hard_thres
            || m_target_filter.provided_dose_thres)
        {
            m_error_message.append("Warning: Hard thresholding will only be "
                                   "performed for imputation input.\n");
        }
        if (m_target_filter.info_score > 0)
        {
            m_error_message.append(
                "Warning: INFO score can only be calculated for "
                "imputation input.\n");
        }
    }
    if (m_target.type == "bgen")
    {
        if (!misc::within_bound(m_target_filter.info_score, 0.0, 1.0))
        {
            error = true;
            m_error_message.append(
                "Error: INFO score threshold cannot be bigger than 1.0 "
                "or smaller than 0.0\n");
        }
        if (!misc::within_bound(m_target_filter.hard_threshold, 0.0, 1.0))
        {
            error = true;
            m_error_message.append(
                "Error: Hard threshold must be between 0 and 1!\n");
        }
        if (!misc::within_bound(m_target_filter.dose_threshold, 0.0, 1.0))
        {
            error = true;
            m_error_message.append(
                "Error: Dosage threshold must be between 0 and 1!\n");
        }
    }


    if (!misc::within_bound<double>(m_target_filter.geno, 0.0, 1.0))
    {
        error = true;
        m_error_message.append("Error: Genotype missingness threshold cannot "
                               "be bigger than 1.0 "
                               "or smaller than 0.0\n");
    }
    if (!misc::within_bound<double>(m_target_filter.maf, 0.0, 1.0))
    {
        error = true;
        m_error_message.append("Error: MAF threshold cannot be bigger than 1.0 "
                               "or smaller than 0.0\n");
    }
    return !error;
}

bool Commander::misc_check()
{
    bool error = false;
    if (!m_extract_file.empty() && !m_exclude_file.empty())
    {
        error = true;
        m_error_message.append(
            "Error: Can only use --extract or --exclude but not both\n");
    }

    m_parameter_log["seed"] = misc::to_string(m_perm_info.seed);
    if (m_prs_info.thread <= 0)
    {
        error = true;
        m_error_message.append(
            "Error: Number of thread must be larger than 1\n");
    }
    if (m_keep_ambig)
    {
        m_error_message.append(
            "Warning: By selecting --keep-ambig, PRSice assume the base and "
            "target are reporting alleles on the same strand and will "
            "therefore only perform dosage flip for the ambiguous SNPs. If you "
            "are unsure of what the strand is, then you should not select the "
            "--keep-ambig option\n");
    }
    if (!m_perm_info.run_perm && !m_perm_info.run_set_perm
        && m_perm_info.logit_perm)
    {
        m_error_message.append("Warning: Permutation not required, "
                               "--logit-perm has no effect\n");
    }
    // for no regress, we will alway print the scores (otherwise no point
    // running PRSice)
    if (m_prs_info.no_regress) m_print_all_scores = true;
    // Just in case thread wasn't provided, we will print the default number
    // of thread used
    m_parameter_log["thread"] = std::to_string(m_prs_info.thread);
    m_parameter_log["out"] = m_out_prefix;
    const bool use_reference =
        !(m_reference.file_list.empty() && m_reference.file_name.empty());
    if (m_prs_info.use_ref_maf && !use_reference)
    {
        error = true;
        // for now, force use_ref_maf to be used together with --ld.
        // might be able to do otherwise, but that will be too troublesome
        m_error_message.append(
            "Error: Cannot use reference MAF for missingness "
            "imputation if reference file isn't used\n");
    }
    if (m_allow_inter)
    {
        if ((m_target.type != "bgen"
             && m_reference.type != "bgen") // none are bgen
            || (use_reference && m_reference.type != "bgen"
                && !m_target.hard_coded) // reference file isn't bgen and not
                                         // require hard coding PRS
            || (!use_reference
                && m_target.type != "bgen") // Doesn't have a reference file and
                                            // target isn't bgen either
        )
        {
            m_allow_inter = false;
            m_error_message.append(
                "Warning: Intermediate not required. Will not "
                "generate intermediate file\n");
        }
    }
    if (m_prs_info.no_regress && m_pheno_info.pheno_col.size() > 1)
    {
        m_error_message.append(
            "Warning: --no-regress requested but multiple "
            "phenotype provided. As regression isn't performed, we will not "
            "utilize any of the phenotype information\n");
    }
    if (m_target.type == "bgen" && !m_target.hard_coded && m_ultra_aggressive)
    {
        m_error_message.append("Warning: --ultra does not work with none "
                               "hard-coded bgen file. Will disable it\n");
        m_ultra_aggressive = false;
    }
    return !error;
}

bool Commander::prset_check()
{
    if (!m_prset.run) return true;
    bool error = false;
    if (m_prset.gtf.empty() && !m_prset.msigdb.empty())
    {
        error = true;
        m_error_message.append(
            "Error: Must provide a gtf file if msigdb is specified\n");
    }

    if (m_prset.feature.empty() && !m_prset.gtf.empty())
    {
        m_prset.feature = {"exon", "gene", "protein_coding", "CDS"};
        m_parameter_log["feature"] = "exon,gene,protein_coding,CDS";
    }
    if (m_perm_info.run_perm && m_perm_info.run_set_perm)
    {
        error = true;
        m_error_message.append(
            "Error: Currently only support either set-base "
            "permutation (for competitive p-value) or PRSice "
            "base permutation (--perm)");
    }
    if (m_prset.gtf.empty() && m_prset.background.empty()
        && m_perm_info.run_set_perm)
    {
        // by default, if background is not provided, we will use the gtf as
        // the background, otherwise, we will use the whole genome as the
        // background
        m_error_message.append(
            "Warrning: Background file and gtf file not "
            "provided. Will use the whole genome as the "
            "background for competitive p-value calculation");
        m_prset.full_as_background = true;
        m_parameter_log["full-back"] = "";
    }
    return !error;
}

bool Commander::prsice_check()
{
    bool error = false;
    // First pass cleaning of barlevel
    std::sort(m_p_thresholds.bar_levels.begin(),
              m_p_thresholds.bar_levels.end());
    m_p_thresholds.bar_levels.erase(
        std::unique(m_p_thresholds.bar_levels.begin(),
                    m_p_thresholds.bar_levels.end()),
        m_p_thresholds.bar_levels.end());
    if (m_prset.run)
    {
        // if prset is performed
        if (m_p_thresholds.bar_levels.empty())
        {
            // two different default. If any thresholding related parameter
            // were used, use the default PRSice threshold. Otherwise, only use
            // 1
            if (!m_p_thresholds.set_threshold && !m_p_thresholds.fastscore)
            {
                m_parameter_log["bar-levels"] = 1;
                m_p_thresholds.fastscore = true;
                m_p_thresholds.bar_levels = {1};
            }
            else
            {
                m_p_thresholds.bar_levels = {0.001, 0.05, 0.1, 0.2,
                                             0.3,   0.4,  0.5};
            }
        }
    }
    else if (m_p_thresholds.bar_levels.empty())
    {
        m_p_thresholds.bar_levels = {0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
    }
    if (!m_p_thresholds.no_full) m_p_thresholds.bar_levels.push_back(1);
    // Second pass (mainly aiming for 1)
    std::sort(m_p_thresholds.bar_levels.begin(),
              m_p_thresholds.bar_levels.end());
    m_p_thresholds.bar_levels.erase(
        std::unique(m_p_thresholds.bar_levels.begin(),
                    m_p_thresholds.bar_levels.end()),
        m_p_thresholds.bar_levels.end());
    // now check if there are any negative / out bound in bar level
    auto max_threshold = *max_element(m_p_thresholds.bar_levels.begin(),
                                      m_p_thresholds.bar_levels.end());
    auto min_threshold = *min_element(m_p_thresholds.bar_levels.begin(),
                                      m_p_thresholds.bar_levels.end());
    if (max_threshold > 1.0)
    {
        error = true;
        m_error_message.append("Error: Cannot have p-value level > 1\n");
    }
    if (min_threshold < 0.0)
    {
        error = true;
        m_error_message.append(
            "Error: Cannot have p-value level less than 0\n");
    }
    std::string bar_message = "";
    for (auto&& b : m_p_thresholds.bar_levels)
    {
        if (bar_message.empty())
            bar_message.append(misc::to_string(b));
        else
            bar_message.append("," + misc::to_string(b));
    }
    m_parameter_log["bar-levels"] = bar_message;
    if (!m_p_thresholds.fastscore)
    {
        // perform high resolution scoring
        if (m_p_thresholds.inter <= 0)
        {
            error = true;
            m_error_message.append("Error: Cannot have negative interval!\n");
        }
        if (m_p_thresholds.upper < m_p_thresholds.lower)
        {
            error = true;
            m_error_message.append(
                "Error: Upper bound must be larger than lower bound!\n");
        }
        if (!misc::within_bound(m_p_thresholds.upper, 0.0, 1.0)
            || !misc::within_bound(m_p_thresholds.lower, 0.0, 1.0))
        {
            error = true;
            m_error_message.append(
                "Error: Invalid p-value threshold boundary!\n");
        }
        if (!misc::within_bound(m_p_thresholds.inter, 0.0, 1.0))
        {
            error = true;
            m_error_message.append("Error: Invalid p-value step-size!\n");
        }
        m_parameter_log["interval"] = misc::to_string(m_p_thresholds.inter);
        m_parameter_log["lower"] = misc::to_string(m_p_thresholds.lower);
        m_parameter_log["upper"] = misc::to_string(m_p_thresholds.upper);
    }
    return !error;
}

bool Commander::target_check()
{
    bool error = false;
    if (m_target.file_name.empty() && m_target.file_list.empty())
    {
        error = true;
        m_error_message.append(
            "Error: You must provide a target file or a file "
            "containing all target prefixs!\n");
    }
    if (!m_target.file_name.empty() && !m_target.file_list.empty())
    {
        error = true;
        m_error_message.append(
            "Error: You can only use --target or --target-list "
            "but not both\n");
    }
    if (!m_target.keep.empty() && !m_target.remove.empty())
    {
        error = true;
        m_error_message.append(
            "Error: Can use either --keep or --remove but not both\n");
    }
    if (std::find(supported_types.begin(), supported_types.end(), m_target.type)
        == supported_types.end())
    {
        error = true;
        m_error_message.append(
            "Error: Unsupported target format: " + m_target.type + "\n");
    }
    // if use target hard coding, make sure we do hard thresholding
    if (m_target.type == "bgen" && m_target.hard_coded)
    {
        m_parameter_log["hard-thres"] =
            std::to_string(m_target_filter.hard_threshold);
        m_parameter_log["dose-thres"] =
            std::to_string(m_target_filter.dose_threshold);
    }
    if (m_target.num_autosome < 0)
    {
        error = true;
        m_error_message.append(
            "Error: Currently do not support negative number of autosome");
    }
    else
    {
        m_parameter_log["num-auto"] = std::to_string(m_target.num_autosome);
    }
    return !error;
}

bool Commander::pheno_check()
{
    assert(m_ran_base_check);
    // pheno check must be performed after base check
    bool error = false;
    if (!m_pheno_info.pheno_col.empty() && m_pheno_info.pheno_file.empty())
    {
        error = true;
        m_error_message.append("Error: You must provide a phenotype file for "
                               "multiple phenotype analysis");
        return !error;
    }
    // check for duplicates
    if (!m_pheno_info.pheno_col.empty())
    {
        std::unordered_set<std::string> phenos(m_pheno_info.pheno_col.begin(),
                                               m_pheno_info.pheno_col.end());
        if (phenos.size() != m_pheno_info.pheno_col.size())
        {
            error = true;
            m_error_message.append(
                "Error: Duplicated phenotype column detected. Please make sure "
                "you have provided the correct input\n");
            return !error;
        }
    }
    if (m_pheno_info.binary.empty())
    {
        // add the default
        const size_t repeat =
            m_pheno_info.pheno_col.empty() ? 1 : m_pheno_info.pheno_col.size();
        if (m_base_info.is_beta)
        {
            m_parameter_log["binary-target"] = "F";
            if (repeat > 1)
            {
                m_parameter_log["binary-target"] = std::to_string(repeat) + "F";
            }
            m_pheno_info.binary.resize(repeat, false);
        }
        else
        {
            m_parameter_log["binary-target"] = "T";
            if (repeat > 1)
            {
                m_parameter_log["binary-target"] = std::to_string(repeat) + "T";
            }
            m_pheno_info.binary.resize(repeat, true);
        }
    }
    // now check if binary-target is sensible
    if (m_pheno_info.pheno_col.size() != m_pheno_info.binary.size())
    {
        if (m_pheno_info.pheno_col.empty() && m_pheno_info.binary.size() == 1)
        {
            // this is ok
            // now that we have always initialized m_is_binary
            // before the check, there shouldn't be a case where
            // m_is_binary < 1
        }
        else
        {
            error = true;
            m_error_message.append(
                "Error: Number of target phenotypes doesn't "
                "match information of binary target! You must "
                "indicate whether the phenotype is binary "
                "using --binary-target\n");
        }
    }
    // check if we have sufficient amount of prevalence info
    size_t num_bin = 0;
    for (auto binary : m_pheno_info.binary)
    {
        if (binary) ++num_bin;
    }

    if (!m_pheno_info.prevalence.empty())
    {
        if (num_bin
            != m_pheno_info.prevalence.size()) // need to be all or nothing
        {
            error = true;
            m_error_message.append(
                "Error: Number of target prevalence doesn't match "
                "number of binary traits. You must provide a prevalence for "
                "all "
                "binary trait(s) or not provide any prevalence (all or "
                "nothing)\n");
        }
        for (auto&& prev : m_pheno_info.prevalence)
        {
            if (!misc::within_bound(prev, 0.0, 1.0))
            {
                error = true;
                m_error_message.append(
                    "Error: Prevalence must be within 0 and 1\n");
                break;
            }
        }
    }
    return !error;
}

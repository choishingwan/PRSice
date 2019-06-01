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

#include "commander.hpp"

Commander::Commander()
{
    m_base_col_index.resize(+BASE_INDEX::MAX + 1, 0);
    m_base_has_col.resize(+BASE_INDEX::MAX + 1, false);
    set_help_message();
}

bool Commander::init(int argc, char* argv[], Reporter& reporter)
{
    if (argc <= 1) {
        usage();
        throw std::runtime_error("Please provide the required parameters");
    }
    const char* optString = "b:B:c:C:f:F:g:h?i:k:l:L:m:n:o:p:s:t:u:v";
    const struct option longOpts[] = {
        // parameters with short flags
        {"base", required_argument, nullptr, 'b'},
        {"bed", required_argument, nullptr, 'B'},
        {"cov-col", required_argument, nullptr, 'c'},
        // Add short form (because I am lazy)
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
        {"beta", no_argument, &m_stat_is_beta, 1},
        {"fastscore", no_argument, &m_fastscore, 1},
        {"full-back", required_argument, &m_full_background, 1},
        {"hard", no_argument, &m_target_is_hard_coded, 1},
        {"ignore-fid", no_argument, &m_ignore_fid, 1},
        {"index", no_argument, &m_input_is_index, 1},
        {"keep-ambig", no_argument, &m_keep_ambig, 1},
        {"logit-perm", no_argument, &m_logit_perm, 1},
        {"no-clump", no_argument, &m_no_clump, 1},
        {"non-cumulate", no_argument, &m_non_cumulate_prs, 1},
        {"no-default", no_argument, &m_user_no_default, 1},
        {"no-full", no_argument, &m_no_full, 1},
        {"no-regress", no_argument, &m_no_regress, 1},
        {"nonfounders", no_argument, &m_include_nonfounders, 1},
        {"or", no_argument, &m_stat_is_or, 1},
        {"pearson", no_argument, &m_pearson, 1},
        {"print-snp", no_argument, &m_print_snp, 1},
        {"use-ref-maf", no_argument, &m_use_ref_maf, 1},
        // long flags, need to work on them
        {"A1", required_argument, nullptr, 0},
        {"A2", required_argument, nullptr, 0},
        {"background", required_argument, nullptr, 0},
        {"bar-levels", required_argument, nullptr, 0},
        {"base-info", required_argument, nullptr, 0},
        {"base-maf", required_argument, nullptr, 0},
        {"binary-target", required_argument, nullptr, 0},
        {"bp", required_argument, nullptr, 0},
        {"chr", required_argument, nullptr, 0},
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
        {"perm", required_argument, nullptr, 0},
        {"proxy", required_argument, nullptr, 0},
        {"remove", required_argument, nullptr, 0},
        {"score", required_argument, nullptr, 0},
        {"set-perm", required_argument, nullptr, 0},
        {"shrink-perm", required_argument, nullptr, 0},
        {"snp", required_argument, nullptr, 0},
        {"snp-set", required_argument, nullptr, 0},
        {"stat", required_argument, nullptr, 0},
        {"target-list", required_argument, nullptr, 0},
        {"type", required_argument, nullptr, 0},
        {"wind-5", required_argument, nullptr, 0},
        {"wind-3", required_argument, nullptr, 0},
        {"x-range", required_argument, nullptr, 0},
        {nullptr, 0, nullptr, 0}};
    return parse_command(argc, argv, optString, longOpts, reporter);
}

bool Commander::parse_command(int argc, char* argv[], const char* optString,
                              const struct option longOpts[],
                              Reporter& reporter)
{

    int32_t max_threads = 1;
#if defined(WIN32) || defined(_WIN32) \
    || defined(__WIN32) && !defined(__CYGWIN__)
    // max thread estimation using windows
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    max_threads = sysinfo.dwNumberOfProcessors;
    int32_t known_procs = max_threads;
#else
    int32_t known_procs = static_cast<int32_t>(sysconf(_SC_NPROCESSORS_ONLN));
    max_threads = (known_procs == -1) ? 1 : known_procs;
#endif
    int longIndex = 0;
    int opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    // storing all the used parameters
    // this allow us to show the users all parameters in effect
    std::map<std::string, std::string> message_store;
    std::string command;
    std::string error_messages = "";
    bool dummy = false;
    bool error = false;
    while (opt != -1) {
        switch (opt)
        {
        case 0:
            command = longOpts[longIndex].name;
            if (longOpts[longIndex].flag != nullptr) break;
            // reorganize all long ops according to alphabetical order
            else if (command == "A1")
                set_string(optarg, message_store, m_effect_allele,
                           m_provided_effect_allele, command, error_messages);
            else if (command == "A2")
                set_string(optarg, message_store, m_non_effect_allele,
                           m_provided_non_effect_allele, command,
                           error_messages);
            else if (command == "background")
                set_string(optarg, message_store, m_background, dummy, command,
                           error_messages);
            else if (command == "bar-levels")
            {
                error |= !load_numeric_vector<double>(
                    optarg, message_store, error_messages, m_barlevel, command);
                m_set_use_thresholds = true;
            }
            else if (command == "base-info")
                set_string(optarg, message_store, m_info_col,
                           m_provided_info_threshold, command, error_messages);
            else if (command == "base-maf")
                set_string(optarg, message_store, m_maf_col,
                           m_perform_base_maf_control_filter, command,
                           error_messages);
            else if (command == "binary-target")
                error |=
                    !parse_binary_vector(optarg, message_store, error_messages,
                                         m_is_binary, command);
            else if (command == "bp")
                set_string(optarg, message_store, m_bp, m_provided_bp, command,
                           error_messages);
            else if (command == "chr")
                set_string(optarg, message_store, m_chr, m_provided_chr_col,
                           command, error_messages);
            else if (command == "clump-kb")
            {
                m_clump_distance =
                    set_distance(optarg, command, 1000, message_store, error,
                                 error_messages);
                m_provided_clump_dist = true;
            }
            else if (command == "clump-p")
                error |=
                    !set_numeric<double>(optarg, message_store, error_messages,
                                         m_clump_p, dummy, command);
            else if (command == "clump-r2")
                error |=
                    !set_numeric<double>(optarg, message_store, error_messages,
                                         m_clump_r2, dummy, command);
            else if (command == "cov-factor")
                load_string_vector(optarg, message_store, m_factor_cov, command,
                                   error_messages);
            else if (command == "dose-thres")
                error |=
                    !set_numeric<double>(optarg, message_store, error_messages,
                                         m_target_dose_thres, dummy, command);
            else if (command == "exclude")
                set_string(optarg, message_store, m_exclude_file, dummy,
                           command, error_messages);
            else if (command == "extract")
                set_string(optarg, message_store, m_extract_file, dummy,
                           command, error_messages);
            else if (command == "feature")
                load_string_vector(optarg, message_store, m_feature, command,
                                   error_messages);
            else if (command.compare("geno") == 0)
                error |= !set_numeric<double>(optarg, message_store,
                                              error_messages, m_target_geno,
                                              m_target_geno_filter, command);
            else if (command == "hard-thres")
                error |=
                    !set_numeric<double>(optarg, message_store, error_messages,
                                         m_target_hard_threshold,
                                         m_target_hard_thresholding, command);
            else if (command == "id-delim")
                set_string(optarg, message_store, m_id_delim, m_set_delim,
                           command, error_messages, true);
            else if (command == "info")
                error |= !set_numeric<double>(
                    optarg, message_store, error_messages, m_target_info_score,
                    m_target_info_filter, command);
            else if (command == "keep")
                set_string(optarg, message_store, m_target_keep, dummy, command,
                           error_messages);
            else if (command == "ld-dose-thres")
                error |=
                    !set_numeric<double>(optarg, message_store, error_messages,
                                         m_ref_dose_thres, dummy, command);
            else if (command == "ld-geno")
                error |= !set_numeric<double>(
                    optarg, message_store, error_messages, m_ref_geno,
                    m_perform_ref_geno_filter, command);
            else if (command == "ld-hard-thres")
                error |= !set_numeric<double>(
                    optarg, message_store, error_messages, m_ref_hard_threshold,
                    m_perform_ref_hard_thresholding, command);
            else if (command == "ld-info")
                error |= !set_numeric<double>(
                    optarg, message_store, error_messages, m_ref_info_score,
                    m_perform_ref_info_filter, command);
            else if (command == "ld-keep")
                set_string(optarg, message_store, m_ref_keep, dummy, command,
                           error_messages);
            else if (command == "ld-list")
                set_string(optarg, message_store, m_ref_list,
                           m_ref_list_provided, command, error_messages);
            else if (command == "ld-maf")
                error |= !set_numeric<double>(
                    optarg, message_store, error_messages, m_ref_maf,
                    m_perform_ref_maf_filter, command);
            else if (command == "ld-remove")
                set_string(optarg, message_store, m_ref_remove, dummy, command,
                           error_messages);
            else if (command == "ld-type")
                set_string(optarg, message_store, m_ref_type, dummy, command,
                           error_messages);
            else if (command == "maf")
                error |= !set_numeric<double>(optarg, message_store,
                                              error_messages, m_target_maf,
                                              m_target_maf_filter, command);
            else if (command == "memory")
                error |= !set_memory(optarg, message_store, error_messages);
            else if (command == "missing")
                error |= !set_missing(optarg, message_store, error_messages);
            else if (command == "model")
                error |= !set_model(optarg, message_store, error_messages);
            else if (command == "perm")
            {
                // use double to account for scientific?
                if (std::string(optarg).at(0) == '-') {
                    error = true;
                    error_messages.append(
                        "Error: Negative permutation number detected!\n");
                }
                else
                {
                    error |= !set_numeric<size_t>(optarg, message_store,
                                                  error_messages, m_permutation,
                                                  dummy, command);
                    m_perform_permutation = true;
                }
            }
            else if (command == "proxy")
                error |= !set_numeric<double>(optarg, message_store,
                                              error_messages, m_proxy_threshold,
                                              m_use_proxy_clump, command);
            else if (command == "remove")
                set_string(optarg, message_store, m_target_remove, dummy,
                           command, error_messages);
            else if (command == "score")
                error |= !set_score(optarg, message_store, error_messages);
            else if (command == "set-perm")
            {
                if (std::string(optarg).at(0) == '-') {
                    error = true;
                    error_messages.append("Error: Negative set based "
                                          "permutation number detected!\n");
                }
                else
                {
                    error |= !set_numeric<size_t>(optarg, message_store,
                                                  error_messages, m_set_perm,
                                                  dummy, command);
                    m_perform_set_perm = true;
                }
            }
            else if (command == "snp")
                set_string(optarg, message_store, m_snp, m_provided_snp_id,
                           command, error_messages);
            else if (command == "snp-set")
                set_string(optarg, message_store, m_snp_set, m_perform_prset,
                           command, error_messages);
            else if (command == "stat")
                set_string(optarg, message_store, m_statistic,
                           m_provided_statistic, command, error_messages);
            else if (command == "target-list")
                set_string(optarg, message_store, m_target_file_list,
                           m_use_target_list, command, error_messages);
            else if (command == "type")
                set_string(optarg, message_store, m_target_type, dummy, command,
                           error_messages);
            else if (command == "wind-3")
                m_window_3 = set_distance(optarg, command, 1, message_store,
                                          error, error_messages);
            else if (command == "wind-5")
                m_window_5 = set_distance(optarg, command, 1, message_store,
                                          error, error_messages);
            else if (command.compare("x-range") == 0)
                set_string(optarg, message_store, m_exclusion_range, dummy,
                           command, error_messages);
            else
            {
                std::string er = "Error: Undefined operator: " + command
                                 + ", please use --help for more information!";
                throw std::runtime_error(er);
            }
            break;
        case 'b':
            set_string(optarg, message_store, m_base_file, dummy, "base",
                       error_messages);
            break;
        case 'B':
            load_string_vector(optarg, message_store, m_bed_files, "bed",
                               error_messages);
            m_perform_prset = true;
            break;
        case 'c':
            load_string_vector(optarg, message_store, m_cov_colname, "cov-col",
                               error_messages);
            break;
        case 'C':
            set_string(optarg, message_store, m_cov_file, dummy, "cov-file",
                       error_messages);
            break;
        case 'f':
            set_string(optarg, message_store, m_pheno_file, dummy, "pheno-file",
                       error_messages);
            break;
        case 'F':
            load_string_vector(optarg, message_store, m_pheno_col, "pheno-col",
                               error_messages);
            break;
        case 'g':
            set_string(optarg, message_store, m_gtf, m_perform_prset, "gtf",
                       error_messages);
            break;
        case 'i':
            error |= !set_numeric<double>(optarg, message_store, error_messages,
                                          m_inter_threshold,
                                          m_set_use_thresholds, "interval");
            break;
        case 'k':
            error |= !load_numeric_vector<double>(optarg, message_store,
                                                  error_messages, m_prevalence,
                                                  "prevalence");
            break;
        case 'l':
            error |= !set_numeric<double>(optarg, message_store, error_messages,
                                          m_lower_threshold,
                                          m_set_use_thresholds, "lower");
            break;
        case 'L':
            set_string(optarg, message_store, m_ref_file, m_use_reference, "ld",
                       error_messages);
            break;
        case 'm':
            set_string(optarg, message_store, m_msigdb, m_perform_prset,
                       "msigdb", error_messages);
            break;
        case 'n':
            if (strcmp("max", optarg) == 0) {
                m_thread = max_threads;
                message_store["thread"] = std::to_string(m_thread);
            }
            else
            {
                error |=
                    !set_numeric<int>(optarg, message_store, error_messages,
                                      m_thread, dummy, "thread");
                if (m_thread > max_threads) {
                    m_thread = max_threads;
                    message_store["thread"] = std::to_string(m_thread);
                }
            }
            break;
        case 'o':
            set_string(optarg, message_store, m_out_prefix, dummy, "out",
                       error_messages);
            break;
        case 'p':
            set_string(optarg, message_store, m_p_value, m_provided_p_value,
                       "pvalue", error_messages);
            break;
        case 's':
            error |= !set_numeric<std::random_device::result_type>(
                optarg, message_store, error_messages, m_seed, m_provided_seed,
                "seed");
            break;
        case 't':
            set_string(optarg, message_store, m_target_file, dummy, "target",
                       error_messages);
            break;
        case 'u':
            error |= !set_numeric<double>(optarg, message_store, error_messages,
                                          m_upper_threshold,
                                          m_set_use_thresholds, "upper");
            break;
        case 'h':
        case '?': usage(); return false;
        case 'v':
            std::cerr << version << " (" << date << ") " << std::endl;
            return false;
        default:
            throw "Error: Undefined operator, please use --help for more "
                  "information!";
        }
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    }
    error |= !base_check(message_store, error_messages);
    error |= !clump_check(message_store, error_messages);
    error |= !covariate_check(error_messages);
    error |= !filter_check(error_messages);
    error |= !misc_check(message_store, error_messages);
    error |= !ref_check(message_store, error_messages);
    error |= !prset_check(message_store, error_messages);
    error |= !prsice_check(message_store, error_messages);
    error |= !target_check(message_store, error_messages);
    // check all flags
    std::string log_name = m_out_prefix + ".log";
    reporter.initiailize(log_name);


    if (m_print_all_scores) message_store["all-score"] = "";
    if (m_allow_inter) message_store["allow-inter"] = "";
    if (m_stat_is_or) message_store["or"] = "";
    if (m_stat_is_beta) message_store["beta"] = "";
    if (m_fastscore) message_store["fastscore"] = "";
    if (m_target_is_hard_coded) message_store["hard"] = "";
    if (m_ignore_fid) message_store["ignore-fid"] = "";
    if (m_input_is_index) message_store["index"] = "";
    if (m_keep_ambig) message_store["keep-ambig"] = "";
    if (m_logit_perm) message_store["logit-perm"] = "";
    if (m_no_clump) message_store["no-clump"] = "";
    if (m_user_no_default) message_store["no-default"] = "";
    if (m_no_full) message_store["no-full"] = "";
    if (m_no_regress) message_store["no-regress"] = "";
    if (m_include_nonfounders) message_store["nonfounders"] = "";
    if (m_pearson) message_store["pearson"] = "";
    if (m_print_snp) message_store["print-snp"] = "";
    if (m_use_ref_maf) message_store["use-ref-maf"] = "";
    std::chrono::time_point<std::chrono::system_clock> start;
    start = std::chrono::system_clock::now();
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    struct tm* timeinfo;
    char buffer[80];
    timeinfo = localtime(&start_time);

    strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", timeinfo);

    std::string message = "\nPRSice " + version + " (" + date + ") \n";
    message.append("https://github.com/choishingwan/PRSice\n");
    message.append("(C) 2016-2017 Shing Wan (Sam) Choi, Jack Euesden, Cathryn "
                   "M. Lewis, Paul F. O'Reilly\n");
    message.append("GNU General Public License v3\n\n");
    message.append("If you use PRSice in any published work, please cite:\n");
    message.append("Jack Euesden Cathryn M. Lewis Paul F. O'Reilly (2015)\n");
    message.append("PRSice: Polygenic Risk Score software.\n");
    message.append("Bioinformatics 31 (9): 1466-1468\n\n");


    std::string time_str(buffer);
    std::string prog_name = argv[0];
    message.append(time_str + "\n" + prog_name);

    for (auto&& com : message_store) {
        message.append(" \\\n    --" + com.first + " " + com.second);
    }
    message.append("\n");
    reporter.report(message, false);
    if (!error_messages.empty()) reporter.report(error_messages);
    if (error) throw std::runtime_error(error_messages);
    return true;
}


// Default constructor of Command
// Responsible for setting all the default values
// initialize the parameters, then call the
// parameter processing function
// Default destructor of Command, do nothing
Commander::~Commander()
{
    // dtor
}

// Function to set the help message
// avoid having large chunk of un-foldable code
void Commander::set_help_message()
{
    m_help_message =
        "usage: PRSice [options] <-b base_file> <-t target_file>\n"
        // Base file
        "\nBase File:\n"
        "    --A1                    Column header containing allele 1 "
        "(effective allele)\n"
        "                            Default: A1\n"
        "    --A2                    Column header containing allele 2 "
        "(non-effective allele)\n"
        "                            Default: A2\n"
        "    --base          | -b    Base association file\n"
        "    --base-info             Base INFO score filtering. Format should "
        "be\n"
        "                            <Column name>,<Threshold>. SNPs with info "
        "\n"
        "                            score less than <Threshold> will be "
        "ignored\n"
        "                            Column name default: INFO\n"
        "                            Threshold default: 0.9\n"
        "    --base-maf             Base MAF filtering. Format should be\n"
        "                            <Column name>,<Threshold>. SNPs with maf\n"
        "                            less than <Threshold> will be ignored. "
        "An\n"
        "                            additional column can also be added "
        "(e.g.\n"
        "                            also filter MAF for cases), using the\n"
        "                            following format:\n"
        "                            <Column name>,<Threshold>:<Column "
        "name>,<Threshold>\n"
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
        "    --index                 If set, assume the INDEX instead of NAME  "
        "for\n"
        "                            the corresponding columns are provided. "
        "Index\n"
        "                            should be 0-based (start counting from "
        "0)\n"
        "    --no-default            Remove all default options. If set, "
        "PRSice\n"
        "                            will not set any default column name and "
        "you\n"
        "                            will have to ensure all required columns "
        "are\n"
        "                            provided. (--snp, --stat, --A1, "
        "--pvalue)\n"
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
        "                            all of them (Either adjust all binary "
        "traits,\n"
        "                            or don't adjust at all)\n"
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
        + misc::to_string(m_target_hard_threshold)
        + "\n"
          "    --hard                  Use hard coding instead of dosage for "
          "PRS construction.\n"
          "                            Default is to use dosage instead of "
          "hard coding\n"
          // clumping
          "\nClumping:\n"
          "    --clump-kb              The distance for clumping in kb\n"
          "                            Default: "
        + misc::to_string(m_clump_distance / 1000)
        + "\n"
          "    --clump-r2              The R2 threshold for clumping\n"
          "                            Default: "
        + misc::to_string(m_clump_r2)
        + "\n"
          "    --clump-p               The p-value threshold use for "
          "clumping.\n"
          "                            Default: "
        + misc::to_string(m_clump_p)
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
        + misc::to_string(m_ref_hard_threshold)
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
        + misc::to_string(m_inter_threshold)
        + "\n"
          "    --lower         | -l    The starting p-value threshold. "
          "Default: "
        + misc::to_string(m_lower_threshold)
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
          "                            avg - Take the average effect size "
          "(default)\n"
          "                            std - Standardize the effect size \n"
          "                            sum - Direct summation of the effect "
          "size \n"
          "    --upper         | -u    The final p-value threshold. Default: "
        + misc::to_string(m_upper_threshold)
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
          "                                              containing the name of"
          " the SNP\n"
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
          "    --logit-perm            When performing permutation, still use "
          "logistic\n"
          "                            regression instead of linear "
          "regression. This\n"
          "                            will substantially slow down PRSice\n"
          "    --keep-ambig            Keep ambiguous SNPs. Only use this "
          "option\n"
          "                            if you are certain that the base and "
          "target\n"
          "                            has the same A1 and A2 alleles\n"
          "    --memory                Maximum memory usage allowed. PRSice "
          "will try\n"
          "                            its best to honor this setting\n"
          "    --non-cumulate          Calculate non-cumulative PRS. PRS will "
          "be reset\n"
          "                            to 0 for each new P-value threshold "
          "instead of\n"
          "                            adding up\n"
          "    --out           | -o    Prefix for all file output\n"
          "    --pearson               Use Pearson Correlation for LD "
          "calculation\n"
          "                            instead of the maximum likelihood "
          "haplotype\n"
          "                            frequency estimates. This will slightly "
          "\n"
          "                            decrease the accuracy of LD estimates, "
          "but\n"
          "                            should increase the speed of clumping\n"
          "    --perm                  Number of permutation to perform. This "
          "swill\n"
          "                            generate the empirical p-value. "
          "Recommend to\n"
          "                            use value larger than 10,000\n"
          "    --print-snp             Print all SNPs used to construct the "
          "best PRS\n"
          "    --seed          | -s    Seed used for permutation. If not "
          "provided,\n"
          "                            system time will be used as seed. When "
          "same\n"
          "                            seed and same input is provided, same "
          "result\n"
          "                            can be generated\n"
          "    --thread        | -n    Number of thread use\n"
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

// Print the help message
void Commander::usage() { fprintf(stderr, "%s\n", m_help_message.c_str()); }


bool Commander::base_check(std::map<std::string, std::string>& message,
                           std::string& error_message)
{
    bool error = false;
    bool has_col = false;
    size_t col_index;
    if (m_base_file.empty()) {
        error_message.append("Error: You must provide a base file\n");
        return false;
    }
    // get input header
    std::string header;
    if (m_base_file.substr(m_base_file.find_last_of(".") + 1).compare("gz")
        == 0)
    {
        GZSTREAM_NAMESPACE::igzstream in(m_base_file.c_str());
        if (!in.good()) {
            error_message.append(
                "Error: Cannot open base file (gz) to read!\n");
            return false;
        }
        std::getline(in, header);
        in.close();
    }
    else
    {
        std::ifstream base_test;
        base_test.open(m_base_file.c_str());
        if (!base_test.is_open()) {
            error_message.append("Error: Cannot open base file: " + m_base_file
                                 + " to read!\n");
            error_message.append("       " + std::string(strerror(errno))
                                 + "\n");
            return false;
        }
        std::getline(base_test, header);
        base_test.close();
    }
    misc::trim(header);
    std::vector<std::string> column_names = misc::split(header);
    if (m_user_no_default) {
        // remove all the default
        if (!m_provided_chr_col) m_chr = "";
        if (!m_provided_effect_allele) m_effect_allele = "";
        if (!m_provided_non_effect_allele) m_non_effect_allele = "";
        if (!m_provided_statistic) m_statistic = "";
        if (!m_provided_snp_id) m_snp = "";
        if (!m_provided_bp) m_bp = "";
        if (!m_provided_p_value) m_p_value = "";
        if (!m_provided_info_threshold) m_info_col = "";
    }
    if (m_input_is_index) {
        // can't do much but to check the boundary
        for (size_t i = 0; i < column_names.size(); ++i) {
            column_names[i] = std::to_string(i);
        }
    }
    has_col = index_check(m_chr, column_names, col_index);
    if (has_col) m_base_col_index[+BASE_INDEX::CHR] = col_index;
    m_base_has_col[+BASE_INDEX::CHR] = has_col;
    if (has_col)
        message["chr"] = m_chr;
    else if (m_provided_chr_col)
    {
        error_message.append("Warning: " + m_chr + " not found in base file\n");
        message.erase("chr");
    }
    has_col = index_check(m_effect_allele, column_names, col_index);
    if (has_col) m_base_col_index[+BASE_INDEX::REF] = col_index;
    m_base_has_col[+BASE_INDEX::REF] = has_col;
    if (has_col)
        message["A1"] = m_effect_allele;
    else if (m_provided_effect_allele)
    {
        error = true;
        error_message.append("Error: " + m_effect_allele
                             + " not found in base file\n");
    }
    has_col = index_check(m_non_effect_allele, column_names, col_index);
    if (has_col) m_base_col_index[+BASE_INDEX::ALT] = col_index;
    m_base_has_col[+BASE_INDEX::ALT] = has_col;
    if (has_col)
        message["A2"] = m_non_effect_allele;
    else if (m_provided_non_effect_allele)
    {
        error_message.append("Warning: " + m_non_effect_allele
                             + " not found in base file\n");
        message.erase("A2");
    }
    has_col = index_check(m_snp, column_names, col_index);
    if (has_col) m_base_col_index[+BASE_INDEX::RS] = col_index;
    m_base_has_col[+BASE_INDEX::RS] = has_col;
    if (has_col)
        message["snp"] = m_snp;
    else if (m_provided_snp_id)
    {
        error = true;
        error_message.append("Error: " + m_snp + " not found in base file\n");
    }
    has_col = index_check(m_bp, column_names, col_index);
    if (has_col) m_base_col_index[+BASE_INDEX::BP] = col_index;
    m_base_has_col[+BASE_INDEX::BP] = has_col;
    if (has_col)
        message["bp"] = m_bp;
    else if (m_provided_bp)
    {
        error_message.append("Warning: " + m_bp + " not found in base file\n");
        message.erase("bp");
    }

    has_col = index_check(m_p_value, column_names, col_index);
    if (has_col) m_base_col_index[+BASE_INDEX::P] = col_index;
    m_base_has_col[+BASE_INDEX::P] = has_col;
    if (has_col)
        message["pvalue"] = m_p_value;
    else if (m_provided_p_value)
    {
        error = true;
        error_message.append("Error: " + m_p_value
                             + " not found in base file\n");
    }
    std::string tmp_error_message;
    bool tmp_error;
    tmp_error = !set_base_info_threshold(column_names, tmp_error_message);
    if (m_provided_info_threshold) {
        // only need to provide the error message when user wants the filtering
        error |= tmp_error;
        error_message.append(tmp_error_message);
    }
    tmp_error_message.clear();
    tmp_error = !set_base_maf_filter(column_names, tmp_error_message);
    if (!m_maf_col.empty()) {
        error |= tmp_error;
        error_message.append(tmp_error_message);
    }
    // now process the statistic column
    if (m_stat_is_or && m_stat_is_beta) {
        error_message.append("Error: Statistic cannot be both OR and beta\n");
        error = true;
    }
    has_col = index_check(m_statistic, column_names, col_index);
    if (has_col) m_base_col_index[+BASE_INDEX::STAT] = col_index;
    m_base_has_col[+BASE_INDEX::STAT] = has_col;
    if (has_col)
        message["stat"] = m_statistic;
    else if (m_provided_statistic)
    {
        error_message.append("Error: " + m_statistic
                             + " not found in base file\n");
    }
    else if (!m_user_no_default)
    {
        // we can't find the default statistic column (BETA)
        // and user allow default option
        if (m_stat_is_or) {
            // search for OR
            has_col = index_check("OR", column_names, col_index);
            if (has_col) m_base_col_index[+BASE_INDEX::STAT] = col_index;
            m_base_has_col[+BASE_INDEX::STAT] = has_col;
            if (!has_col) {
                error = true;
                error_message.append("Error: Cannot find appropriate "
                                     "statistic column in base file!\n");
            }
            else
            {
                message["stat"] = "OR";
            }
        }
        else if (m_stat_is_beta)
        {
            // search for BETA
            has_col = index_check("BETA", column_names, col_index);
            if (has_col) m_base_col_index[+BASE_INDEX::STAT] = col_index;
            m_base_has_col[+BASE_INDEX::STAT] = has_col;
            if (!has_col) {
                error = true;
                error_message.append("Error: Cannot find appropriate "
                                     "statistic column in base file!\n");
            }
            else
            {
                message["stat"] = "BETA";
            }
        }
        else
        {
            // go through file and look for either OR or BETA
            bool or_found = false, beta_found = false;
            for (size_t i = 0; i < column_names.size(); ++i) {
                std::string temp = column_names[i];
                std::transform(temp.begin(), temp.end(), temp.begin(),
                               ::toupper);
                if (temp == "OR") {
                    m_stat_is_beta = false;
                    m_stat_is_or = true;
                    m_statistic = column_names[i];
                    message["stat"] = column_names[i];
                    or_found = true;
                    message["or"] = "";
                    message.erase("beta");
                    m_base_col_index[+BASE_INDEX::STAT] = i;
                    m_base_has_col[+BASE_INDEX::STAT] = true;
                }
                else if (temp == "BETA")
                {
                    m_stat_is_beta = true;
                    m_stat_is_or = false;
                    m_statistic = column_names[i];
                    message["stat"] = column_names[i];
                    message["beta"] = "";
                    beta_found = true;
                    message.erase("or");
                    m_base_col_index[+BASE_INDEX::STAT] = i;
                    m_base_has_col[+BASE_INDEX::STAT] = true;
                }
                if (beta_found && or_found) {
                    error_message.append(
                        "Error: Both OR and BETA "
                        "found in base file! We cannot determine "
                        "which statistic to use, please provide it "
                        "through --stat\n");
                    error = true;
                    break;
                }
            }
        }
    }

    // Statistic is ok, but beta or or not provided
    if (m_base_has_col[+BASE_INDEX::STAT]) {
        if (!m_stat_is_or && !m_stat_is_beta) {
            std::string stat_temp = m_statistic;
            std::transform(stat_temp.begin(), stat_temp.end(),
                           stat_temp.begin(), ::toupper);
            if (stat_temp == "OR") {
                m_stat_is_or = true;
                message["or"] = "";
            }
            else if (stat_temp == "BETA")
            {
                m_stat_is_beta = true;
                message["beta"] = "";
            }
        }
    }
    // now check all required columns are here
    if (!m_base_has_col[+BASE_INDEX::P]) {
        // we can actually losen this requirement if user doesn't
        // perform clumping and p-value thresholding
        error = true;
        error_message.append("Error: No p-value column (" + m_p_value
                             + ") in file!\n");
    }
    if (!m_base_has_col[+BASE_INDEX::STAT]) {
        error = true;
        error_message.append("Error: No statistic column (" + m_statistic
                             + ") in file!\n");
    }
    if (!m_base_has_col[+BASE_INDEX::RS]) {
        error = true;
        error_message.append("Error: No SNP name column (" + m_snp
                             + ") in file!\n");
    }
    if (!m_base_has_col[+BASE_INDEX::REF]) {
        error = true;
        error_message.append("Error: No Reference allele column ("
                             + m_effect_allele + ") in file!\n");
    }
    // we don't need bp and chr as we can always get those from the bim file
    // use a for loop as it is short enough and we only bother with those
    // we have index for
    size_t max_index = 0;
    for (size_t i = 0; i < m_base_col_index.size(); ++i) {
        if (m_base_has_col[i] && max_index < m_base_col_index[i]) {
            max_index = m_base_col_index[i];
        }
    }
    m_base_col_index[+BASE_INDEX::MAX] = max_index;
    return !error;
}

bool Commander::clump_check(std::map<std::string, std::string>& message,
                            std::string& error_message)
{
    bool error = false;
    if (!m_no_clump) {
        if (m_use_proxy_clump
            && (m_proxy_threshold < 0 || m_proxy_threshold > 1))
        {
            error = true;
            error_message.append(
                "Error: Proxy threshold must be within 0 and 1!\n");
        }
        if (m_clump_p < 0.0 || m_clump_p > 1.0) {
            error = true;
            error_message.append(
                "Error: P-value threshold must be within 0 and 1!\n");
        }
        if (m_clump_r2 < 0.0 || m_clump_r2 > 1.0) {
            error = true;
            error_message.append(
                "Error: R2 threshold must be within 0 and 1!\n");
        }
        if (!m_ref_type.empty()) {
            if (std::find(supported_types.begin(), supported_types.end(),
                          m_ref_type)
                == supported_types.end())
            {
                error = true;
                error_message.append(
                    "Error: Unsupported LD format: " + m_ref_type + "\n");
            }
        }
        if (m_clump_distance < 0.0) {
            error = true;
            error_message.append(
                "Error: Clumping distance must be positive!\n");
        }
        message["clump-r2"] = std::to_string(m_clump_r2);
        message["clump-p"] = std::to_string(m_clump_p);
        // we divided by 1000 here to make sure it is in KB (our preferred
        // format)
        if (!m_provided_clump_dist && m_perform_prset) {
            // change default based on PRSet or not
            m_clump_distance = 1000000;
        }
        message["clump-kb"] = std::to_string(m_clump_distance / 1000);
    }
    return !error;
}


bool Commander::ref_check(std::map<std::string, std::string>& message,
                          std::string& error_message)
{
    bool error = false;
    if (!m_ref_keep.empty() && !m_ref_remove.empty()) {
        error = true;
        error_message.append("Error: Can only use either --keep or "
                             "--remove but not both\n");
    }
    // only check the range is correct when it is needed
    if (m_perform_ref_geno_filter && (m_ref_geno < 0 || m_ref_geno > 1)) {
        error = true;
        error_message.append("Error: LD genotype missingness threshold "
                             "must be larger than 0 and smaller than 1!\n");
    }
    // if the reference panel is bgen, or reference panel not provided
    // but the target is bgen then we will like to enforce hard
    // thresholding to the files for LD calculation
    if (!m_ref_file.empty() && m_ref_list_provided) {
        error = true;
        error_message.append(
            "Error: You can only use --target or --target-list "
            "but not both\n");
    }
    if (m_ref_type == "bgen"
        || (m_ref_file.empty() && !m_ref_list_provided
            && m_target_type == "bgen"))
    {
        if (m_ref_hard_threshold > 1 || m_ref_hard_threshold < 0) {
            error = true;
            error_message.append("Error: LD hard threshold must be larger "
                                 "than 0 and smaller than 1!\n");
        }
        else if (!m_ref_file.empty() || m_ref_list_provided)
        {
            // reference file is provided, so hard thresholding is added
            // to ld-hard-thres (doens't matter if user has provided
            // this parameter or not)
            message["ld-hard-thres"] = std::to_string(m_ref_hard_threshold);
        }
        else
        {
            // reference file is not provided, so hard thresholding is
            // added to hard-thres as it is to be applied to the
            // target-ish (doens't matter if user has provided this
            // parameter or not)
            message["hard-thres"] = std::to_string(m_target_hard_threshold);
        }
    }
    if (m_perform_ref_maf_filter && (m_ref_maf > 1 || m_ref_maf < 0)) {
        error = true;
        error_message.append("Error: LD MAF threshold must be larger than "
                             "0 and smaller than 1!\n");
    }
    if (m_perform_ref_info_filter
        && (m_ref_info_score < 0 || m_ref_info_score > 1))
    {
        error = true;
        error_message.append("Error: LD INFO score threshold must be "
                             "larger than 0 and smaller than 1!\n");
    }
    return !error;
}
std::vector<std::string> Commander::transform_covariate(std::string& cov)
{
    std::vector<std::string> final_covariates;
    std::vector<std::string> open;
    std::vector<std::string> close;
    std::vector<std::string> info;
    std::vector<std::string> individual;
    std::vector<std::string> range;
    std::vector<int> numeric;
    std::vector<bool> list;
    if (cov.at(0) == '@') {
        // this needs transformation
        // remove the @ form the front of the string
        cov.erase(0, 1);
        // find the start of range by identifying [
        open = misc::split(cov, "[");
        for (auto&& o : open) {
            if (o.find("]") != std::string::npos) {
                // we also found close range in this
                close = misc::split(o, "]");
                // the first one will always be the list
                info.push_back(close[0]);
                list.push_back(true);
                // Nested List is not supported
                for (std::vector<std::string>::size_type cl = 1;
                     cl < close.size(); ++cl)
                {
                    info.push_back(close[cl]);
                    list.push_back(false);
                }
            }
            else
            {
                // TODO: this can only be a nested list, should issue a
                // warning
                info.push_back(o);
                list.push_back(false);
            }
        }

        for (std::vector<std::string>::size_type c = 0; c < info.size(); ++c) {
            if (list[c]) {
                individual = misc::split(info[c], ".");
                numeric.clear();
                for (auto&& ind : individual) {
                    if (ind.find("-") != std::string::npos) {
                        range = misc::split(ind, "-");
                        if (range.size() != 2) {
                            throw std::runtime_error(
                                "Error: Invalid range format, range "
                                "must be in the form of start-end");
                        }
                        try
                        {
                            size_t start = misc::convert<size_t>(range[0]);
                            size_t end = misc::convert<size_t>(range[1]);
                            if (start > end) {
                                std::swap(start, end);
                            }
                            for (std::vector<std::string>::size_type s = start;
                                 s <= end; ++s)
                            {
                                numeric.push_back(static_cast<int>(s));
                            }
                        }
                        catch (...)
                        {
                            std::string error_message =
                                "Error: Invalid parameter: " + range[0] + " or "
                                + range[1] + ", only allow integer!";
                            throw std::runtime_error(error_message);
                        }
                    }
                    else
                    {
                        try
                        {
                            int temp = misc::convert<int>(ind);
                            numeric.push_back(temp);
                        }
                        catch (...)
                        {
                            std::string error_message =
                                "Error: Invalid parameter: " + ind
                                + ", only allow integer!";
                            throw std::runtime_error(error_message);
                        }
                    }
                }

                // Now we have all the numeric parameters
                if (final_covariates.empty()) {
                    for (auto n : numeric) {
                        final_covariates.push_back(std::to_string(n));
                    }
                }
                else
                {
                    auto&& cur_size = final_covariates.size();
                    for (std::vector<std::string>::size_type final = 0;
                         final < cur_size; ++final)
                    {
                        std::string cur = final_covariates[final];
                        final_covariates[final].append(
                            std::to_string(numeric.front()));
                        for (std::vector<std::string>::size_type s = 1;
                             s < numeric.size(); ++s)
                        {
                            final_covariates.push_back(
                                cur + std::to_string(numeric[s]));
                        }
                    }
                }
            }
            else
            {
                for (std::vector<std::string>::size_type final = 0;
                     final < final_covariates.size(); ++final)
                {
                    final_covariates[final].append(info[c]);
                }
                if (final_covariates.empty())
                    final_covariates.push_back(info[c]);
            }
        }
    }
    else
    {
        // this doesn't need transformation, just return this
        final_covariates.push_back(cov);
    }
    return final_covariates;
}

bool Commander::covariate_check(std::string& error_message)
{
    // it is valid to have empty covariate file
    if (m_cov_file.empty()) return true;
    // first, transform all the covariates
    // the actual column name to be included (after parsing)
    std::unordered_set<std::string> included;
    // contain all the input within the --cov-col
    std::unordered_set<std::string> ori_input;
    // vector contain the transformed column names
    std::vector<std::string> transformed_cov;
    for (auto cov : m_cov_colname) {
        if (cov.empty()) continue;
        ori_input.insert(cov);
        // got annoyed with the input of PC.1 PC.2 PC.3, do this automatic
        // thingy to substitute them
        transformed_cov = transform_covariate(cov);
        for (auto&& trans : transformed_cov) {
            included.insert(trans);
        }
    }
    bool error = false;
    // now try to read the header of the covariate file
    std::ifstream cov_file;
    cov_file.open(m_cov_file.c_str());
    if (!cov_file.is_open()) {
        error_message.append("Error: Cannot open covariate file: " + m_cov_file
                             + "\n");
        // really can't do checking, will have to terminate now
        return false;
    }
    std::string line;
    std::getline(cov_file, line);
    if (line.empty()) {
        error_message.append("Error: First line of covariate file is empty!\n");
        return false;
    }
    // remove all the special characters
    misc::trim(line);
    cov_file.close();
    std::vector<std::string> cov_header = misc::split(line);
    std::string missing = "";
    std::unordered_map<std::string, uint32_t> ref_index;
    // now get the index for each column name in the covariate file
    for (std::vector<int>::size_type i = 0; i < cov_header.size(); ++i) {
        ref_index[cov_header[i]] = static_cast<uint32_t>(i);
    }

    // when user provide a covariate file but not the covariate name, we
    // will just read in every covariates
    if (m_cov_colname.size() == 0) {
        // add all headers to the covariate list
        for (std::vector<int>::size_type i = (1 + !m_ignore_fid);
             i < cov_header.size(); ++i)
        {
            included.insert(cov_header[i]);
        }
    }
    size_t valid_cov = 0;
    // covariate.covariates.clear();
    m_col_index_of_cov.clear();
    for (auto&& cov : included) {
        // now for each covariate found in the covariate file, we add their
        // index to the storage
        if (ref_index.find(cov) != ref_index.end()) {
            m_col_index_of_cov.push_back(ref_index[cov]);
            // covariate.covariates.push_back(cov);
            valid_cov++;
        }
        // store information of covariates not found in the covarite file
        else if (missing.empty())
        {
            missing = cov;
        }
        else
        {
            missing.append("," + cov);
        }
    }
    if (!missing.empty()) {
        error_message.append("Warning: Covariate(s) missing from file: "
                             + missing + ". Header of file is: " + line + "\n");
    }
    if (valid_cov == 0) {
        error = true;
        error_message.append("Error: No valid Covariate!\n");
    }
    // we will now push back the covariate name according to the order they
    // appeared in the file
    m_cov_colname.clear();
    std::sort(m_col_index_of_cov.begin(), m_col_index_of_cov.end());
    for (auto&& c : m_col_index_of_cov) {
        m_cov_colname.push_back(cov_header[c]);
    }
    // now start to process the factor covariates
    for (auto cov : m_factor_cov) {
        if (cov.empty()) continue;
        transformed_cov = transform_covariate(cov);
        for (auto&& trans : transformed_cov) {
            if (included.find(trans) != included.end()) {
                m_col_index_of_factor_cov.push_back(ref_index[trans]);
            }
            else if (ori_input.find(cov) == ori_input.end())
            {
                // only complain if untransform input isn't found in cov-col
                error = true;
                error_message.append("Error: All factor covariates must be "
                                     "found in covariate list. "
                                     + trans + " not found in covariate list");
            }
        }
    }
    std::sort(m_col_index_of_factor_cov.begin(),
              m_col_index_of_factor_cov.end());
    return !error;
}


bool Commander::filter_check(std::string& error_message)
{
    bool error = false;
    if (m_target_type != "bgen" && m_target_hard_thresholding) {
        error_message.append("Warning: Hard thresholding will only be "
                             "performed for imputation input.\n");
    }
    if (m_target_type != "bgen" && m_target_info_filter) {
        error_message.append("Warning: INFO score can only be calculated for "
                             "imputation input.\n");
    }
    if (m_target_type == "bgen" && m_target_hard_thresholding
        && (m_target_hard_threshold <= 0 || m_target_hard_threshold >= 1))
    {
        error = true;
        error_message.append(
            "Error: Hard threshold must be between 0 and 1!\n");
    }
    if (!m_extract_file.empty() && !m_exclude_file.empty()) {
        error = true;
        error_message.append(
            "Error: Can only use --extract or --exclude but not both\n");
    }

    if (m_target_info_filter
        && (m_target_info_score < 0 || m_target_info_score > 1))
    {
        error = true;
        error_message.append(
            "Error: INFO score threshold cannot be bigger than 1.0 "
            "or smaller than 0.0\n");
    }
    if (m_target_geno_filter && (m_target_geno < 0 || m_target_geno > 1)) {
        error = true;
        error_message.append("Error: Genotype missingness threshold cannot "
                             "be bigger than 1.0 "
                             "or smaller than 0.0\n");
    }
    if (m_target_maf_filter && (m_target_maf < 0 || m_target_maf > 1)) {
        error = true;
        error_message.append("Error: MAF threshold cannot be bigger than 1.0 "
                             "or smaller than 0.0\n");
    }
    return !error;
}

bool Commander::misc_check(std::map<std::string, std::string>& message,
                           std::string& error_message)
{
    bool error = false;
    if (!m_provided_seed) {
        message["seed"] = misc::to_string(m_seed);
    }
    if (m_thread <= 0) {
        error = true;
        error_message.append("Error: Number of thread must be larger than 1\n");
    }
    if (!m_perform_permutation && m_logit_perm) {
        error_message.append("Warning: Permutation not required, "
                             "--logit-perm has no effect\n");
    }
    // for no regress, we will alway print the scores (otherwise no point
    // running PRSice)
    if (m_no_regress) m_print_all_scores = true;
    // Just in case thread wasn't provided, we will print the default number
    // of thread used
    if (m_thread == 1) message["thread"] = "1";
    message["out"] = m_out_prefix;
    if (m_use_ref_maf && !m_use_reference) {
        error = true;
        // for now, force use_ref_maf to be used together with --ld.
        // might be able to do otherwise, but that will be too troublesome
        error_message.append("Error: Cannot use reference MAF for missingness "
                             "imputation if reference file isn't used\n");
    }
    if (m_allow_inter) {
        if ((m_target_type != "bgen" && m_ref_type != "bgen")
            || (m_use_reference && m_ref_type != "bgen"
                && !m_target_is_hard_coded)
            || (!m_use_reference && m_target_type != "bgen"))
        {
            m_allow_inter = false;
            error_message.append("Warning: Intermediate not required. Will not "
                                 "generate intermediate file\n");
        }
    }
    return !error;
}

bool Commander::prset_check(std::map<std::string, std::string>& message,
                            std::string& error_message)
{
    bool error = false;
    if (!m_perform_prset) return true;
    if (m_gtf.empty() && !m_msigdb.empty()) {
        error = true;
        error_message.append(
            "Error: Must provide a gtf file if msigdb is specified\n");
    }
    if (m_window_3 < 0 || m_window_5 < 0) {
        error = true;
        error_message.append(
            "Error: 5' and 3' extension must be larger than 0\n");
    }
    if (m_feature.empty()) {
        m_feature.push_back("exon");
        m_feature.push_back("gene");
        m_feature.push_back("protein_coding");
        m_feature.push_back("CDS");
        message["feature"] = "exon,gene,protein_coding,CDS";
    }
    if (m_perform_permutation && m_perform_set_perm) {
        error = true;
        error_message.append("Error: Currently only support either set-base "
                             "permutation (for competitive p-value) or PRSice "
                             "base permutation (--perm)");
    }
    if (m_gtf.empty() && m_background.empty() && m_perform_set_perm) {
        // by default, if background is not provided, we will use the gtf as
        // the background, otherwise, we will use the whole genome as the
        // background
        error_message.append("Warrning: Background file and gtf file not "
                             "provided. Will use the whole genome as the "
                             "background for competitive p-value calculation");
        m_full_background = true;
        message["full-back"] = "";
    }

    return !error;
}


bool Commander::prsice_check(std::map<std::string, std::string>& message,
                             std::string& error_message)
{
    bool error = false;
    // we should use assert, as there is no way the m_genetic_model is not
    // within the one we implemented here unless there're programming errors
    // add the message for model
    assert(m_genetic_model == MODEL::ADDITIVE
           || m_genetic_model == MODEL::DOMINANT
           || m_genetic_model == MODEL::RECESSIVE
           || m_genetic_model == MODEL::HETEROZYGOUS);
    switch (m_genetic_model)
    {
    case MODEL::ADDITIVE: message["model"] = "add"; break;
    case MODEL::DOMINANT: message["model"] = "dom"; break;
    case MODEL::RECESSIVE: message["model"] = "rec"; break;
    case MODEL::HETEROZYGOUS: message["model"] = "het"; break;
    }
    assert(m_missing_score == MISSING_SCORE::CENTER
           || m_missing_score == MISSING_SCORE::SET_ZERO
           || m_missing_score == MISSING_SCORE::MEAN_IMPUTE);
    switch (m_missing_score)
    {
    case MISSING_SCORE::CENTER: message["missing"] = "CENTER"; break;
    case MISSING_SCORE::SET_ZERO: message["missing"] = "SET_ZERO"; break;
    case MISSING_SCORE::MEAN_IMPUTE: message["missing"] = "MEAN_IMPUTE"; break;
    }
    assert(m_scoring_method == SCORING::SUM
           || m_scoring_method == SCORING::AVERAGE
           || m_scoring_method == SCORING::STANDARDIZE);
    switch (m_scoring_method)
    {
    case SCORING::SUM: message["score"] = "sum"; break;
    case SCORING::AVERAGE: message["score"] = "avg"; break;
    case SCORING::STANDARDIZE: message["score"] = "std"; break;
    }
    // First pass cleaning of barlevel
    std::sort(m_barlevel.begin(), m_barlevel.end());
    m_barlevel.erase(std::unique(m_barlevel.begin(), m_barlevel.end()),
                     m_barlevel.end());
    if (m_perform_prset) {
        // if prset is performed
        if (m_barlevel.empty()) {
            // two different default. If any thresholding related parameter
            // were used, use the default PRSice threshold. Otherwise, only use
            // 1
            if (!m_set_use_thresholds && !m_fastscore) {
                message["bar-levels"] = 1;
                m_fastscore = true;
                m_barlevel = {1};
            }
            else
            {
                m_barlevel = {0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
            }
        }
    }
    else if (m_barlevel.empty())
    {
        // if PRSice is used
        m_barlevel = {0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
    }
    if (!m_no_full) m_barlevel.push_back(1);
    // Second pass (mainly aiming for 1)
    std::sort(m_barlevel.begin(), m_barlevel.end());
    m_barlevel.erase(std::unique(m_barlevel.begin(), m_barlevel.end()),
                     m_barlevel.end());
    std::string bar_message = "";
    for (auto&& b : m_barlevel) {
        if (bar_message.empty())
            bar_message.append(misc::to_string(b));
        else
            bar_message.append("," + misc::to_string(b));
    }
    message["bar-levels"] = bar_message;
    if (!m_fastscore) {
        // perform high resolution scoring
        if (m_inter_threshold <= 0) {
            error = true;
            error_message.append("Error: Cannot have negative interval!\n");
        }
        if (m_upper_threshold < m_lower_threshold) {
            error = true;
            error_message.append(
                "Error: Upper bound must be larger than lower bound!\n");
        }
        if (m_upper_threshold < 0.0 || m_lower_threshold < 0.0) {
            error = true;
            error_message.append("Error: Cannot have negative bounds!\n");
        }

        message["interval"] = misc::to_string(m_inter_threshold);
        message["lower"] = misc::to_string(m_lower_threshold);
        message["upper"] = misc::to_string(m_upper_threshold);
    }
    return !error;
}

bool Commander::target_check(std::map<std::string, std::string>& message,
                             std::string& error_message)
{
    bool error = false;
    if (m_target_file.empty() && !m_use_target_list) {
        error = true;
        error_message.append("Error: You must provide a target file or a file "
                             "containing all target prefixs!\n");
    }
    if (!m_target_file.empty() && m_use_target_list) {
        error = true;
        error_message.append(
            "Error: You can only use --target or --target-list "
            "but not both\n");
    }
    if (!m_target_keep.empty() && !m_target_remove.empty()) {
        error = true;
        error_message.append(
            "Error: Can use either --keep or --remove but not both\n");
    }
    if (std::find(supported_types.begin(), supported_types.end(), m_target_type)
        == supported_types.end())
    {
        error = true;
        error_message.append(
            "Error: Unsupported target format: " + m_target_type + "\n");
    }
    // if use target hard coding, make sure we do hard thresholding
    if (m_target_type.compare("bgen") == 0 && m_target_is_hard_coded) {
        message["hard-thres"] = std::to_string(m_target_hard_threshold);
    }


    if (m_pheno_col.size() != 0 && m_pheno_file.empty()) {
        error = true;
        error_message.append("Error: You must provide a phenotype file for "
                             "multiple phenotype analysis");
    }
    if (m_is_binary.empty()) {
        // add the default
        if (m_stat_is_beta) {
            message["binary-target"] = "F";
            m_is_binary.push_back(false);
        }
        else
        {
            message["binary-target"] = "T";
            m_is_binary.push_back(true);
        }
    }
    // now check if the bar-level is sensible
    if (m_pheno_col.size() != m_is_binary.size()) {
        if (m_pheno_col.empty() && m_is_binary.size() == 1) {
            // this is ok
            // now that we have always initialized m_is_binary
            // before the check, there shouldn't be a case where
            // m_is_binary < 1
        }
        else
        {
            error = true;
            error_message.append("Error: Number of target phenotypes doesn't "
                                 "match information of binary target! You must "
                                 "indicate whether the phenotype is binary "
                                 "using --binary-target\n");
        }
    }
    // check if we have sufficient amount of prevalence info
    size_t num_bin = 0;
    for (auto binary : m_is_binary) {
        if (binary) num_bin++;
    }

    if (!m_prevalence.empty()) {
        if (num_bin != m_prevalence.size()) // need to be all or nothing
        {
            error = true;
            error_message.append(
                "Error: Number of target prevalence doesn't match "
                "number of binary traits. You must provide a prevalence for "
                "all "
                "binary trait(s) or not provide any prevalence (all or "
                "nothing)\n");
        }
        for (auto&& prev : m_prevalence) {
            if (prev > 1.0 || prev < 0.0) {
                error = true;
                error_message.append(
                    "Error: Prevalence cannot be bigger than 1.0 "
                    "or smaller than 0.0\n");
                break;
            }
        }
    }
    return !error;
}

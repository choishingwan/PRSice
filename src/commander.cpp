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

// function to process all parameter input
// return true when we need to continue the program (e.g. when
// --help isn't called)
Commander::Commander()
{
    base.name = "";
    base.chr = "CHR";
    base.effect_allele = "A1";
    base.non_effect_allele = "A2";
    base.statistic = "";
    base.snp = "SNP";
    base.bp = "BP";
    base.standard_error = "SE";
    base.p_value = "P";
    base.info_col = "INFO,0.9";
    base.maf_col = "";
    base.is_beta = false;
    base.is_index = false;
    base.no_default = false;
    base.info_score_threshold = 0;
    base.maf_control_threshold = 0;
    base.maf_case_threshold = 0;
    base.provided_chr = false;
    base.provided_effect_allele = false;
    base.provided_non_effect_allele = false;
    base.provided_statistic = false;
    base.provided_snp = false;
    base.provided_bp = false;
    base.provided_standard_error = false;
    base.provided_p_value = false;
    base.provided_info = false;
    base.col_index.resize(+BASE_INDEX::MAX + 1, -1);


    clumping.no_clump = false;
    clumping.proxy = -1;
    clumping.p_value = 1;
    clumping.r2 = 0.1;
    clumping.distance = 250;
    clumping.provided_proxy = false;

    covariate.file_name = "";

    misc.out = "PRSice";
    misc.non_cumulate = 0;
    misc.exclusion_range = "";
    misc.print_all_scores = false;
    misc.ignore_fid = false;
    misc.logit_perm = false;
    misc.memory = 0;
    misc.pearson = false;
    misc.permutation = 0;
    misc.print_snp = false;
    misc.provided_seed = false;
    misc.provided_memory = false;
    misc.thread = 1;
    misc.seed = 0;

    reference_panel.allow_inter = 0;
    reference_panel.file_name = "";
    reference_panel.multi_name = "";
    reference_panel.type = "bed";
    reference_panel.keep_file = "";
    reference_panel.remove_file = "";

    reference_snp_filtering.geno = 1.0;
    reference_snp_filtering.hard_threshold = 0.9;
    reference_snp_filtering.maf = 0;
    reference_snp_filtering.info_score = 0.0;

    p_thresholds.lower = 0.0001;
    p_thresholds.inter = 0.00005;
    p_thresholds.upper = 0.5;
    p_thresholds.fastscore = false;
    p_thresholds.no_full = false;
    p_thresholds.set_use_thresholds = false;

    prs_calculation.missing_score = "MEAN_IMPUTE";
    // Change this in next major release
    // prs_calculation.score_calculation = "sum";
    prs_calculation.score_calculation = "average";
    prs_calculation.model = MODEL::ADDITIVE;
    prs_calculation.no_regress = false;

    prs_snp_filtering.exclude_file = "";
    prs_snp_filtering.extract_file = "";
    prs_snp_filtering.geno = 1;
    prs_snp_filtering.hard_threshold = 0.9;
    prs_snp_filtering.maf = 0;
    prs_snp_filtering.info_score = 0;
    prs_snp_filtering.is_hard_coded = false;
    prs_snp_filtering.keep_ambig = false;
    prs_snp_filtering.predict_ambig = false;

    prset.gtf = "";
    prset.msigdb = "";
    prset.background = "";
    prset.single_snp_set = "";
    prset.multi_snp_sets = "";
    prset.perform_prset = false;
    prset.perform_set_perm = false;
    prset.set_perm = 10000;
    prset.window_5 = 0;
    prset.window_3 = 0;

    prslice.size = -1;
    prslice.provided = false;

    target.include_nonfounders = false;
    target.name = "";
    target.multi_name = "";
    target.pheno_file = "";
    target.type = "bed";
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
        {"base", required_argument, NULL, 'b'},
        {"bed", required_argument, NULL, 'B'},
        {"cov-col", required_argument, NULL, 'c'},
        // cov-header retain here for backward compatibility
        {"cov-header", required_argument, NULL, 'c'},
        {"cov-file", required_argument, NULL, 'C'},
        {"pheno-file", required_argument, NULL, 'f'},
        {"pheno-col", required_argument, NULL, 'F'},
        {"gtf", required_argument, NULL, 'g'},
        {"help", no_argument, NULL, 'h'},
        {"interval", required_argument, NULL, 'i'},
        {"prevalence", required_argument, NULL, 'k'},
        {"lower", required_argument, NULL, 'l'},
        {"ld", required_argument, NULL, 'L'},
        {"msigdb", required_argument, NULL, 'm'},
        {"thread", required_argument, NULL, 'n'},
        {"out", required_argument, NULL, 'o'},
        {"pvalue", required_argument, NULL, 'p'},
        {"seed", required_argument, NULL, 's'},
        {"target", required_argument, NULL, 't'},
        {"upper", required_argument, NULL, 'u'},
        {"version", no_argument, NULL, 'v'},
        // flags, only need to set them to true
        {"allow-inter", no_argument, &reference_panel.allow_inter, 1},
        {"all-score", no_argument, &misc.print_all_scores, 1},
        {"beta", no_argument, &base.is_beta, 1},
        {"hard", no_argument, &prs_snp_filtering.is_hard_coded, 1},
        {"ignore-fid", no_argument, &misc.ignore_fid, 1},
        {"index", no_argument, &base.is_index, 1},
        {"keep-ambig", no_argument, &prs_snp_filtering.keep_ambig, 1},
        {"logit-perm", no_argument, &misc.logit_perm, 1},
        {"no-clump", no_argument, &clumping.no_clump, 1},
        {"non-cumulate", no_argument, &misc.non_cumulate, 1},
        {"no-default", no_argument, &base.no_default, 1},
        {"no-full", no_argument, &p_thresholds.no_full, 1},
        {"no-regress", no_argument, &prs_calculation.no_regress, 1},
        {"nonfounders", no_argument, &target.include_nonfounders, 1},
        {"fastscore", no_argument, &p_thresholds.fastscore, 1},
        {"pearson", no_argument, &misc.pearson, 1},
        {"print-snp", no_argument, &misc.print_snp, 1},
        // long flags, need to work on them
        {"A1", required_argument, NULL, 0},
        {"A2", required_argument, NULL, 0},
        {"background", required_argument, NULL, 0},
        {"bar-levels", required_argument, NULL, 0},
        {"binary-target", required_argument, NULL, 0},
        {"bp", required_argument, NULL, 0},
        {"chr", required_argument, NULL, 0},
        {"clump-kb", required_argument, NULL, 0},
        {"clump-p", required_argument, NULL, 0},
        {"clump-r2", required_argument, NULL, 0},
        {"cov-factor", required_argument, NULL, 0},
        {"exclude", required_argument, NULL, 0},
        {"extract", required_argument, NULL, 0},
        {"feature", required_argument, NULL, 0},
        {"geno", required_argument, NULL, 0},
        {"hard-thres", required_argument, NULL, 0},
        {"info-base", required_argument, NULL, 0},
        {"info", required_argument, NULL, 0},
        {"keep", required_argument, NULL, 0},
        {"ld-keep", required_argument, NULL, 0},
        {"ld-list", required_argument, NULL, 0},
        {"ld-type", required_argument, NULL, 0},
        {"ld-remove", required_argument, NULL, 0},
        {"ld-maf", required_argument, NULL, 0},
        {"ld-geno", required_argument, NULL, 0},
        {"ld-hard-thres", required_argument, NULL, 0},
        {"ld-info", required_argument, NULL, 0},
        {"maf-base", required_argument, NULL, 0},
        {"maf", required_argument, NULL, 0},
        {"memory", required_argument, NULL, 0},
        {"missing", required_argument, NULL, 0},
        {"model", required_argument, NULL, 0},
        {"perm", required_argument, NULL, 0},
        {"proxy", required_argument, NULL, 0},
        {"prslice", required_argument, NULL, 0},
        {"remove", required_argument, NULL, 0},
        {"score", required_argument, NULL, 0},
        {"se", required_argument, NULL, 0},
        {"set-perm", required_argument, NULL, 0},
        {"snp", required_argument, NULL, 0},
        {"snp-set", required_argument, NULL, 0},
        {"snp-sets", required_argument, NULL, 0},
        {"stat", required_argument, NULL, 0},
        {"target-list", required_argument, NULL, 0},
        {"type", required_argument, NULL, 0},
        {"wind-5", required_argument, NULL, 0},
        {"wind-3", required_argument, NULL, 0},
        {"x-range", required_argument, NULL, 0},
        {NULL, 0, 0, 0}};
    return process(argc, argv, optString, longOpts, reporter);
}

bool Commander::process(int argc, char* argv[], const char* optString,
                        const struct option longOpts[], Reporter& reporter)
{

    int32_t max_threads = 1;
#ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    max_threads = sysinfo.dwNumberOfProcessors;
    int32_t known_procs = max_threads;
#else
    int32_t known_procs = sysconf(_SC_NPROCESSORS_ONLN);
    max_threads = (known_procs == -1) ? 1 : known_procs;
#endif


    int longIndex = 0;
    int opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    // storing all the used parameters
    // this allow us to show the users all parameters in effect
    std::map<std::string, std::string> message_store;
    std::string command;
    std::string error_messages = "";
    // the following two variables are used for scientific input
    // e.g. 1e6
    double dummy_double = 0.0;
    double intpart;
    bool dummy = false;
    bool error = false;
    while (opt != -1) {
        switch (opt)
        {
        case 0:
            command = longOpts[longIndex].name;
            if (longOpts[longIndex].flag != 0) break;
            // Long opts for base
            else if (command.compare("chr") == 0)
                set_string(optarg, message_store, base.chr, base.provided_chr,
                           command, error_messages);
            else if (command.compare("A1") == 0)
                set_string(optarg, message_store, base.effect_allele,
                           base.provided_effect_allele, command,
                           error_messages);
            else if (command.compare("A2") == 0)
                set_string(optarg, message_store, base.non_effect_allele,
                           base.provided_non_effect_allele, command,
                           error_messages);
            else if (command.compare("stat") == 0)
                set_string(optarg, message_store, base.statistic,
                           base.provided_statistic, command, error_messages);
            else if (command.compare("snp") == 0)
                set_string(optarg, message_store, base.snp, base.provided_snp,
                           command, error_messages);
            else if (command.compare("bp") == 0)
                set_string(optarg, message_store, base.bp, base.provided_bp,
                           command, error_messages);
            else if (command.compare("se") == 0)
                set_string(optarg, message_store, base.standard_error,
                           base.provided_standard_error, command,
                           error_messages);
            else if (command.compare("info-base") == 0)
                set_string(optarg, message_store, base.info_col,
                           base.provided_info, command, error_messages);
            else if (command.compare("maf-base") == 0)
                set_string(optarg, message_store, base.maf_col, dummy, command,
                           error_messages);
            // Long opts for clumping
            else if (command.compare("clump-p") == 0)
                set_numeric<double>(optarg, message_store, error_messages,
                                    clumping.p_value, dummy, error, command);
            else if (command.compare("clump-r2") == 0)
                set_numeric<double>(optarg, message_store, error_messages,
                                    clumping.r2, dummy, error, command);
            else if (command.compare("clump-kb") == 0)
                set_numeric<int>(optarg, message_store, error_messages,
                                 clumping.distance, dummy, error, command);
            else if (command.compare("proxy") == 0)
                set_numeric<double>(optarg, message_store, error_messages,
                                    clumping.proxy, clumping.provided_proxy,
                                    error, command);
            // Long opts for misc
            else if (command.compare("perm") == 0)
            {
                // use double to account for scientific?
                set_numeric<double>(optarg, message_store, error_messages,
                                    dummy_double, dummy, error, command);
                if (!error) {
                    std::modf(dummy_double, &intpart);
                    misc.permutation = intpart;
                }
            }
            else if (command.compare("x-range") == 0)
            {
                // Require additional processing
                set_string(optarg, message_store, misc.exclusion_range, dummy,
                           command, error_messages);
            }
            // Long opts for reference_panel
            else if (command.compare("ld-keep") == 0)
                set_string(optarg, message_store, reference_panel.keep_file,
                           dummy, command, error_messages);
            else if (command.compare("ld-list") == 0)
                set_string(optarg, message_store, reference_panel.multi_name,
                           dummy, command, error_messages);
            else if (command.compare("ld-remove") == 0)
                set_string(optarg, message_store, reference_panel.remove_file,
                           dummy, command, error_messages);
            else if (command.compare("ld-type") == 0)
                set_string(optarg, message_store, reference_panel.type, dummy,
                           command, error_messages);
            // Long opts for reference_snp_filtering
            else if (command.compare("ld-maf") == 0)
                set_numeric<double>(optarg, message_store, error_messages,
                                    reference_snp_filtering.maf, dummy, error,
                                    command);
            else if (command.compare("ld-geno") == 0)
                set_numeric<double>(optarg, message_store, error_messages,
                                    reference_snp_filtering.geno, dummy, error,
                                    command);
            else if (command.compare("ld-info") == 0)
                set_numeric<double>(optarg, message_store, error_messages,
                                    reference_snp_filtering.info_score, dummy,
                                    error, command);
            else if (command.compare("ld-hard-thres") == 0)
                set_numeric<double>(optarg, message_store, error_messages,
                                    reference_snp_filtering.hard_threshold,
                                    dummy, error, command);
            // Long opts for p_thresholds
            else if (command.compare("bar-levels") == 0)
            {
                load_numeric_vector<double>(
                    optarg, message_store, error_messages,
                    p_thresholds.barlevel, error, command);
                p_thresholds.set_use_thresholds = true;
            }
            // Long opts for prs_calculation
            else if (command.compare("model") == 0)
                set_model(optarg, message_store, error_messages, error);
            else if (command.compare("score") == 0)
                set_string(optarg, message_store,
                           prs_calculation.score_calculation, dummy, command,
                           error_messages);
            else if (command.compare("missing") == 0)
                set_string(optarg, message_store, prs_calculation.missing_score,
                           dummy, command, error_messages);
            // Long opts for prs_snp_filtering
            else if (command.compare("exclude") == 0)
                set_string(optarg, message_store,
                           prs_snp_filtering.exclude_file, dummy, command,
                           error_messages);
            else if (command.compare("extract") == 0)
                set_string(optarg, message_store,
                           prs_snp_filtering.extract_file, dummy, command,
                           error_messages);
            else if (command.compare("geno") == 0)
                set_numeric<double>(optarg, message_store, error_messages,
                                    prs_snp_filtering.geno, dummy, error,
                                    command);
            else if (command.compare("hard-thres") == 0)
                set_numeric<double>(optarg, message_store, error_messages,
                                    prs_snp_filtering.hard_threshold, dummy,
                                    error, command);
            else if (command.compare("maf") == 0)
                set_numeric<double>(optarg, message_store, error_messages,
                                    prs_snp_filtering.maf, dummy, error,
                                    command);
            else if (command.compare("info") == 0)
                set_numeric<double>(optarg, message_store, error_messages,
                                    prs_snp_filtering.info_score, dummy, error,
                                    command);
            // Long opts for PRSet
            else if (command.compare("feature") == 0)
                load_string_vector(optarg, message_store, prset.feature,
                                   command, error_messages);
            else if (command.compare("snp-set") == 0)
                set_string(optarg, message_store, prset.single_snp_set, dummy,
                           command, error_messages);
            else if (command.compare("snp-sets") == 0)
                set_string(optarg, message_store, prset.multi_snp_sets, dummy,
                           command, error_messages);
            else if (command.compare("set-perm") == 0)
                set_numeric<int>(optarg, message_store, error_messages,
                                 prset.set_perm, prset.perform_set_perm, error,
                                 command);
            else if (command.compare("background") == 0)
                set_string(optarg, message_store, prset.background, dummy,
                           command, error_messages);
            else if (command.compare("wind-5") == 0)
                set_numeric<int>(optarg, message_store, error_messages,
                                 prset.window_5, dummy, error, command);
            else if (command.compare("wind-3") == 0)
                set_numeric<int>(optarg, message_store, error_messages,
                                 prset.window_3, dummy, error, command);
            // Long opts for PRSlice
            else if (command.compare("prslice") == 0)
                set_numeric<int>(optarg, message_store, error_messages,
                                 prslice.size, prslice.provided, error,
                                 command);
            // Long opts for target
            else if (command.compare("keep") == 0)
                set_string(optarg, message_store, target.keep_file, dummy,
                           command, error_messages);
            else if (command.compare("remove") == 0)
                set_string(optarg, message_store, target.remove_file, dummy,
                           command, error_messages);
            else if (command.compare("type") == 0)
                set_string(optarg, message_store, target.type, dummy, command,
                           error_messages);
            else if (command.compare("target-list") == 0)
                set_string(optarg, message_store, target.multi_name, dummy,
                           command, error_messages);
            else if (command.compare("binary-target") == 0)
                load_binary_vector(optarg, message_store, error_messages,
                                   target.is_binary, error, command);
            else if (command.compare("memory") == 0)
                set_memory(optarg, message_store, error_messages, error);
            else if (command == "cov-factor")
            {
                load_string_vector(optarg, message_store,
                                   covariate.factor_covariates, "cov-factor",
                                   error_messages);
            }
            else
            {
                std::string er = "Error: Undefined operator: " + command
                                 + ", please use --help for more information!";
                throw std::runtime_error(er);
            }
            break;
        case 'b':
            set_string(optarg, message_store, base.name, dummy, "base",
                       error_messages);
            break;
        case 'B':
            load_string_vector(optarg, message_store, prset.bed, "bed",
                               error_messages);
            prset.perform_prset = true;
            break;
        case 'c':
            load_string_vector(optarg, message_store, covariate.covariates,
                               "cov-col", error_messages);
            break;
        case 'C':
            set_string(optarg, message_store, covariate.file_name, dummy,
                       "cov-file", error_messages);
            break;
        case 'f':
            set_string(optarg, message_store, target.pheno_file, dummy,
                       "pheno-file", error_messages);
            break;
        case 'F':
            load_string_vector(optarg, message_store, target.pheno_col,
                               "pheno-col", error_messages);
            break;
        case 'g':
            set_string(optarg, message_store, prset.gtf, prset.perform_prset,
                       "gtf", error_messages);
            break;
        case 'i':
            set_numeric<double>(
                optarg, message_store, error_messages, p_thresholds.inter,
                p_thresholds.set_use_thresholds, error, "interval");
            break;
        case 'k':
            load_numeric_vector<double>(optarg, message_store, error_messages,
                                        target.prevalence, error, "prevalence");
            break;
        case 'l':
            set_numeric<double>(
                optarg, message_store, error_messages, p_thresholds.lower,
                p_thresholds.set_use_thresholds, error, "lower");
            break;
        case 'L':
            set_string(optarg, message_store, reference_panel.file_name, dummy,
                       "ld", error_messages);
            break;
        case 'm':
            set_string(optarg, message_store, prset.msigdb, prset.perform_prset,
                       "msigdb", error_messages);
            break;
        case 'n':
            if (strcmp("max", optarg) == 0) {
                misc.thread = max_threads;
                message_store["thread"] = std::to_string(misc.thread);
            }
            else
            {
                set_numeric<int>(optarg, message_store, error_messages,
                                 misc.thread, dummy, error, "thread");
                if (misc.thread > max_threads) {
                    misc.thread = max_threads;
                    message_store["thread"] = std::to_string(misc.thread);
                }
            }
            break;
        case 'o':
            set_string(optarg, message_store, misc.out, dummy, "out",
                       error_messages);
            break;
        case 'p':
            set_string(optarg, message_store, base.p_value,
                       base.provided_p_value, "pvalue", error_messages);
            break;
        case 's':
            set_numeric<size_t>(optarg, message_store, error_messages,
                                misc.seed, misc.provided_seed, error, "seed");
            break;
        case 't':
            set_string(optarg, message_store, target.name, dummy, "target",
                       error_messages);
            break;
        case 'u':
            set_numeric<double>(
                optarg, message_store, error_messages, p_thresholds.upper,
                p_thresholds.set_use_thresholds, error, "upper");
            break;
        case 'h':
        case '?':
            usage();
            return false;
            break;
        case 'v':
            std::cerr << version << " (" << date << ") " << std::endl;
            return false;
            break;
        default:
            throw "Error: Undefined operator, please use --help for more "
                  "information!";
        }
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    }

    base_check(message_store, error, error_messages);
    clump_check(message_store, error, error_messages);
    covariate_check(error, error_messages);
    filter_check(error, error_messages);
    misc_check(message_store, error, error_messages);
    prset_check(message_store, error, error_messages);
    prsice_check(message_store, error, error_messages);
    prslice_check(error, error_messages);
    target_check(message_store, error, error_messages);
    if (prset.perform_prset && prslice.provided) {
        error = true;
        error_messages.append(
            "Error: PRSet and PRSlice cannot be performed together!\n");
    }
    // check all flags
    std::string log_name = misc.out + ".log";
    reporter.initiailize(log_name);


    if (base.is_beta) message_store["beta"] = "";
    if (base.is_index) message_store["index"] = "";
    if (base.no_default) message_store["no-default"] = "";
    if (clumping.no_clump) message_store["no-clump"] = "";
    if (prs_snp_filtering.is_hard_coded) message_store["hard"] = "";
    if (prs_snp_filtering.keep_ambig) message_store["keep-ambig"] = "";
    if (misc.print_all_scores) message_store["all-score"] = "";
    if (misc.ignore_fid) message_store["ignore-fid"] = "";
    if (misc.logit_perm) message_store["logit-perm"] = "";
    if (misc.pearson) message_store["pearson"] = "";
    if (misc.print_snp) message_store["print-snp"] = "";
    if (p_thresholds.fastscore) message_store["fastscore"] = "";
    if (p_thresholds.no_full) message_store["no-full"] = "";
    if (prs_calculation.no_regress) message_store["no-regress"] = "";
    if (target.include_nonfounders) message_store["nonfounders"] = "";
    if (reference_panel.allow_inter) message_store["allow-intermediate"] = "";
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
{help_message =
	       "usage: PRSice [options] <-b base_file> <-t target_file>\n"
	           // Base file
	       "\nBase File:\n"
	       "    --A1                    Column header containing allele 1 (effective allele)\n"
	       "                            Default: A1\n"
	       "    --A2                    Column header containing allele 2 (non-effective allele)\n"
	       "                            Default: A2\n"
	       "    --base          | -b    Base association file\n"
	       "    --beta                  Whether the test statistic is in the form of \n"
	       "                            BETA or OR. If set, test statistic is assume\n"
	       "                            to be in the form of BETA.\n"
	       "    --bp                    Column header containing the SNP coordinate\n"
	       "                            Default: BP\n"
	       "    --chr                   Column header containing the chromosome\n"
	       "                            Default: CHR\n"
	       "    --index                 If set, assume the INDEX instead of NAME  for\n"
	       "                            the corresponding columns are provided. Index\n"
	       "                            should be 0-based (start counting from 0)\n"
	       "    --info-base             Base INFO score filtering. Format should be\n"
	       "                            <Column name>,<Threshold>. SNPs with info \n"
	       "                            score less than <Threshold> will be ignored\n"
	       "                            Column name default: INFO\n"
	       "                            Threshold default: 0.9\n"
	       "    --maf-base              Base MAF filtering. Format should be\n"
	       "                            <Column name>,<Threshold>. SNPs with maf\n"
	       "                            less than <Threshold> will be ignored. An\n"
	       "                            additional column can also be added (e.g.\n"
	       "                            also filter MAF for cases), using the\n"
	       "                            following format:\n"
	       "                            <Column name>,<Threshold>:<Column name>,<Threshold>\n"
	       "    --no-default            Remove all default options. If set, PRSice\n"
	       "                            will not set any default column name and you\n"
	       "                            will have to ensure all required columns are\n"
	       "                            provided. (--snp, --stat, --A1, --pvalue)\n"
	       "    --pvalue        | -p    Column header containing the p-value\n"
	       "                            Default: P\n"
	       "    --se                    Column header containing the standard error\n"
	       "                            Default: SE\n"
	       "    --snp                   Column header containing the SNP ID\n"
	       "                            Default: SNP\n"
	       "    --stat                  Column header containing the summary statistic\n"
	       "                            If --beta is set, default as BETA. Otherwise,\n"
	       "                            will search for OR or BETA from the header\n"
	       "                            of the base file\n"
	           // TARGET FILE
	       "\nTarget File:\n"
	       "    --binary-target         Indicate whether the target phenotype\n"
	       "                            is binary or not. Either T or F should be\n"
	       "                            provided where T represent a binary phenotype.\n"
	       "                            For multiple phenotypes, the input should be\n"
	       "                            separated by comma without space. \n"
	       "                            Default: T if --beta and F if --beta is not\n"
	       "    --geno                  Filter SNPs based on gentype missingness\n"
	       "    --info                  Filter SNPs based on info score. Only used\n"
	       "                            for imputed target\n"
	       "    --keep                  File containing the sample(s) to be extracted from\n"
	       "                            the target file. First column should be FID and\n"
	       "                            the second column should be IID. If --ignore-fid is\n"
	       "                            set, first column should be IID\n"
	       "                            Mutually exclusive from --remove\n"
	       "    --maf                   Filter SNPs based on minor allele frequency (MAF)\n"
	       "    --nonfounders           Keep the nonfounders in the analysis\n"
	       "                            Note: They will still be excluded from LD calculation\n"
	       "    --pheno-col     | -F    Headers of phenotypes to be included from the\n"
	       "                            phenotype file\n"
	       "    --pheno-file    | -f    Phenotype file containing the phenotype(s).\n"
	       "                            First column must be FID of the samples and\n"
	       "                            the second column must be IID of the samples.\n"
	       "                            When --ignore-fid is set, first column must\n"
	       "                            be the IID of the samples.\n"
	       "                            Must contain a header if --pheno-col is\n"
	       "                            specified\n"
	       "    --prevalence    | -k    Prevalence of all binary trait. If provided\n"
	       "                            will adjust the ascertainment bias of the R2.\n"
	       "                            Note that when multiple binary trait is found,\n"
	       "                            prevalence information must be provided for\n"
	       "                            all of them (Either adjust all binary traits,\n"
	       "                            or don't adjust at all)\n"
	       "    --remove                File containing the sample(s) to be removed from\n"
	       "                            the target file. First column should be FID and\n"
	       "                            the second column should be IID. If --ignore-fid is\n"
	       "                            set, first column should be IID\n"
	       "                            Mutually exclusive from --keep\n"
	       "    --target        | -t    Target genotype file. Currently support\n"
	       "                            both BGEN and binary PLINK format. For \n"
	       "                            multiple chromosome input, simply substitute\n"
	       "                            the chromosome number with #. PRSice will\n"
	       "                            automatically replace # with 1-22\n"
	       "                            For binary plink format, you can also specify\n"
	       "                            a seperate fam file by <prefix>,<fam file>\n"
	       "    --target-list           File containing prefix of target genotype\n"
	       "                            files. Similar to --target but allow more \n"
	       "                            flexibility. Do not support external fam file\n"
	       "                            at the moment\n"
	       "    --type                  File type of the target file. Support bed \n"
	       "                            (binary plink) and bgen format. Default: bed\n"
	       //dosage
	       "\nDosage:\n"
	       "    --allow-inter           Allow the generate of intermediate file. This will\n"
	       "                            speed up PRSice when using dosage data as clumping\n"
	       "                            reference and for hard coding PRS calculation\n"
	       "    --hard-thres            Hard threshold for dosage data. Any call less than\n"
	       "                            this will be treated as missing. Note that if dosage\n"
	       "                            data is used as a LD reference, it will always be\n"
	       "                            hard coded to calculate the LD\n"
	       "    --hard                  Use hard coding instead of dosage for PRS construction.\n"
	       "                            Default is to use dosage instead of hard coding\n"
	       // clumping
	       "\nClumping:\n"
	       "    --clump-kb              The distance for clumping in kb\n"
	       "                            Default: "+ std::to_string(clumping.distance)+ "\n"
	       "    --clump-r2              The R2 threshold for clumping\n"
	       "                            Default: "+ std::to_string(clumping.r2)+ "\n"
	       "    --clump-p               The p-value threshold use for clumping.\n"
	       "                            Default: "+ std::to_string(clumping.p_value)+ "\n"
	       "    --ld            | -L    LD reference file. Use for LD calculation. If not\n"
	       "                            provided, will use the post-filtered target genotype\n"
	       "                            for LD calculation. Support multiple chromosome input\n"
	       "                            Please see --target for more information\n"
	       "    --ld-list               File containing prefix of LD reference files.\n"
	       "                            Similar to --ld but allow more \n"
	       "                            flexibility. Do not support external fam file\n"
	       "                            at the moment\n"
	       "    --ld-geno               Filter SNPs based on genotype missingness\n"
	       "    --ld-info               Filter SNPs based on info score. Only used\n"
	       "                            for imputed LD reference\n"
	       "    --ld-hard-thres         Hard threshold for dosage data. Any call less than\n"
	       "                            this will be treated as missing.\n"
	       "                            Default: "+std::to_string(reference_snp_filtering.hard_threshold)+"\n"
	       "    --ld-keep               File containing the sample(s) to be extracted from\n"
	       "                            the LD reference file. First column should be FID and\n"
	       "                            the second column should be IID. If --ignore-fid is\n"
	       "                            set, first column should be IID\n"
	       "                            Mutually exclusive from --ld-remove\n"
	       "                            No effect if --ld was not provided\n"
	       "    --ld-maf                Filter SNPs based on minor allele frequency\n"
	       "    --ld-remove             File containing the sample(s) to be removed from\n"
	       "                            the LD reference file. First column should be FID and\n"
	       "                            the second column should be IID. If --ignore-fid is\n"
	       "                            set, first column should be IID\n"
	       "                            Mutually exclusive from --ld-keep\n"
	       "    --ld-type               File type of the LD file. Support bed (binary plink)\n"
	       "                            and bgen format. Default: bed\n"
	       "    --no-clump              Stop PRSice from performing clumping\n"
	       "    --proxy                 Proxy threshold for index SNP to be considered\n"
	       "                            as part of the region represented by the clumped\n"
	       "                            SNP(s). e.g. --proxy 0.8 means the index SNP will\n"
	       "                            represent region of any clumped SNP(s) that has a\n"
	       "                            R2>=0.8 even if the index SNP does not physically\n"
	       "                            locate within the region\n"
	       // Covariates
	       "\nCovariate:\n"
	       "    --cov-col       | -c    Header of covariates. If not provided, will use\n"
	       "                            all variables in the covariate file. By adding\n"
	       "                            @ in front of the string, any numbers within [\n"
	       "                            and ] will be parsed. E.g. @PC[1-3] will be\n"
	       "                            read as PC1,PC2,PC3. Discontinuous input are also\n"
	       "                            supported: @cov[1.3-5] will be parsed as \n"
	       "                            cov1,cov3,cov4,cov5\n"
	       "    --cov-factor            Header of categorical covariate(s). Dummy variable\n"
	       "                            will be automatically generated. Any items in\n"
	       "                            --cov-factor must also be found in --cov-col\n"
	       "                            Also accept continuous input (start with @).\n"
	       "    --cov-file      | -C    Covariate file. First column should be FID and \n"
	       "                            the second column should be IID. If --ignore-fid\n"
	       "                            is set, first column should be IID\n"
	       //PRSice
	       "\nP-value Thresholding:\n"
	       "    --bar-levels            Level of barchart to be plotted. When --fastscore\n"
	       "                            is set, PRSice will only calculate the PRS for \n"
	       "                            threshold within the bar level. Levels should be\n"
	       "                            comma separated without space\n"
	       "    --fastscore             Only calculate threshold stated in --bar-levels\n"
	       "    --no-full               By default, PRSice will include the full model, \n"
	       "                            i.e. p-value threshold = 1. Setting this flag will\n"
	       "                            disable that behaviour\n"
	       "    --interval      | -i    The step size of the threshold. Default: "+ std::to_string(p_thresholds.inter)+ "\n"
	       "    --lower         | -l    The starting p-value threshold. Default: " + std::to_string(p_thresholds.lower)+ "\n"
	       "    --model                 Genetic model use for regression. The genetic\n"
	       "                            encoding is based on the base data where the\n"
	       "                            encoding represent number of the coding allele\n"
	       "                            Available models include:\n"
	       "                            add - Additive model, code as 0/1/2 (default)\n"
	       "                            dom - Dominant model, code as 0/1/1\n"
	       "                            rec - Recessive model, code as 0/0/1\n"
	       "                            het - Heterozygous only model, code as 0/1/0\n"
	       "    --missing               Method to handle missing genotypes. By default, \n"
	       "                            final scores are averages of valid per-allele \n"
	       "                            scores with missing genotypes contribute an amount\n"
	       "                            proportional to imputed allele frequency. To throw\n"
	       "                            out missing observations instead (decreasing the\n"
	       "                            denominator in the final average when this happens),\n"
	       "                            use the 'no_mean_imputation' modifier. Alternatively,\n"
	       "                            you can use the 'center' modifier to shift all scores\n"
	       "                            to mean zero. \n"
	       "    --no-regress            Do not perform the regression analysis and simply\n"
	       "                            output all PRS.\n"
	       "    --score                 Method to calculate the polygenic score.\n"
	       "                            Available methods include:\n"
	       "                            avg - Take the average effect size (default)\n"
	       "                            std - Standardize the effect size \n"
	       "                            sum - Direct summation of the effect size \n"
	       "    --upper         | -u    The final p-value threshold. Default: "+ std::to_string(p_thresholds.upper)+ "\n"
	       "\nPRSet:\n"
	       "    --bed           | -B    Bed file containing the selected regions.\n"
	       "                            Name of bed file will be used as the region\n"
	       "                            identifier. WARNING: Bed file is 0-based\n"
	       "    --feature               Feature(s) to be included from the gtf file.\n"
	       "                            Default: exon,CDS,gene,protein_coding.\n"
	       "    --gtf           | -g    GTF file containing gene boundaries. Required\n"
	       "                            when --msigdb is used\n"
	       "    --msigdb        | -m    MSIGDB file containing the pathway information.\n"
	       "                            Require the gtf file\n"
	       "    --snp-set               Provide a SNP set file containing a single snp set.\n"
	       "                            Name of SNP set file will be used as the region\n"
	       "                            identifier. This file should contain only one column.\n"
	       "    --snp-sets              Provide a SNP set file containing multiple snp sets.\n"
	       "                            Each row represent a single SNP set with the first\n"
	       "                            column containing name of the SNP set.\n"
	       //PRSlice
	       "\nPRSlice:\n"
	       "    --prslice               Perform PRSlice where the whole genome is first cut\n"
	       "                            into bin size specified by this option. PRSice will\n"
	       "                            then be performed on each bin. Bins are then sorted\n"
	       "                            according to the their R2. PRSice is then performed\n"
	       "                            again to find the best bin combination.\n"
	       "                            This cannot be performed together with PRSet\n"
	       "                            (Currently not implemented)\n"
	       //Misc
	       "\nMisc:\n"
	       "    --all-score             Output PRS for ALL threshold. WARNING: This\n"
	       "                            will generate a huge file\n"
	       "    --non-cumulate          Calculate non-cumulative PRS. PRS will be reset\n"
	       "                            to 0 for each new P-value threshold instead of\n"
	       "                            adding up\n"
	       "    --exclude               File contains SNPs to be excluded from the\n"
	       "                            analysis\n"
	       "    --extract               File contains SNPs to be included in the \n"
	       "                            analysis\n"
	       "    --ignore-fid            Ignore FID for all input. When this is set,\n"
	       "                            first column of all file will be assume to\n"
	       "                            be IID instead of FID\n"
	       "    --logit-perm            When performing permutation, still use logistic\n"
	       "                            regression instead of linear regression. This\n"
	       "                            will substantially slow down PRSice\n"
	       "    --keep-ambig            Keep ambiguous SNPs. Only use this option\n"
	       "                            if you are certain that the base and target\n"
	       "                            has the same A1 and A2 alleles\n"
	       "    --out           | -o    Prefix for all file output\n"
	       "    --pearson               Use Pearson Correlation for LD calculation\n"
	       "                            instead of the maximum likelihood haplotype\n"
	       "                            frequency estimates. This will slightly \n"
	       "                            decrease the accuracy of LD estimates, but\n"
	       "                            should increase the speed of clumping\n"
	       "    --perm                  Number of permutation to perform. This swill\n"
	       "                            generate the empirical p-value. Recommend to\n"
	       "                            use value larger than 10,000\n"
	       "    --print-snp             Print all SNPs used to construct the best PRS\n"
	       "    --seed          | -s    Seed used for permutation. If not provided,\n"
	       "                            system time will be used as seed. When same\n"
	       "                            seed and same input is provided, same result\n"
	       "                            can be generated\n"
	       "    --thread        | -n    Number of thread use\n"
	       "    --x-range               Range of SNPs to be excluded from the whole\n"
	       "                            analysis. It can either be a single bed file\n"
	       "                            or a comma seperated list of range. Range must\n"
	       "                            be in the format of chr:start-end or chr:coordinate\n"
	       "    --help          | -h    Display this help message\n";
}

// Print the help message
void Commander::usage() { fprintf(stderr, "%s\n", help_message.c_str()); }


void Commander::base_check(std::map<std::string, std::string>& message,
                           bool& error, std::string& error_message)
{
    if (base.name.empty()) {
        error = true;
        error_message.append("Error: You must provide a base file\n");
    }
    else
    {
        // check the base file and get the corresponding index

        std::string line;
        if (base.name.substr(base.name.find_last_of(".") + 1).compare("gz")
            == 0)
        {
            GZSTREAM_NAMESPACE::igzstream in(base.name.c_str());
            if (!in.good()) {
                error = true;
                error_message.append(
                    "Error: Cannot open base file (gz) to read!\n");
                return;
            }
            std::getline(in, line);
            in.close();
        }
        else
        {
            std::ifstream base_test;
            base_test.open(base.name.c_str());
            if (!base_test.is_open()) {
                error = true;
                error_message.append("Error: Cannot open base file to read!\n");
                return;
            }
            std::getline(base_test, line);
            base_test.close();
        }
        // check the base file header is correct
        std::vector<std::string> token = misc::split(line);
        int max_size = token.size();
        if (base.no_default) {
            // remove all the default
            if (!base.provided_chr) base.chr = "";
            if (!base.provided_effect_allele) base.effect_allele = "";
            if (!base.provided_non_effect_allele) base.non_effect_allele = "";
            if (!base.provided_statistic) base.statistic = "";
            if (!base.provided_snp) base.snp = "";
            if (!base.provided_bp) base.bp = "";
            if (!base.provided_standard_error) base.standard_error = "";
            if (!base.provided_p_value) base.p_value = "";
            if (!base.provided_info) base.info_col = "";
        }
        if (!base.is_index) {
            if (!base.no_default) {
                if (!base.statistic.empty()) {
                    // if statistics is provided, we can guess if it
                    // is beta or not
                    if (base.statistic.length() == 2
                        && toupper(base.statistic[0]) == 'O'
                        && toupper(base.statistic[1]) == 'R')
                    {
                        base.is_beta = false;
                    }
                    else if (base.statistic.length() == 4
                             && toupper(base.statistic[0]) == 'B'
                             && toupper(base.statistic[1]) == 'E'
                             && toupper(base.statistic[2]) == 'T'
                             && toupper(base.statistic[3]) == 'A')
                    {
                        // although user cannot do --no-beta, it is a crazy
                        // use case where BETA != beta, right?
                        // TODO: add no-beta flag...
                        base.is_beta = true;
                        message["beta"] = "";
                    }
                }
                else if (base.statistic.empty() && base.is_beta)
                {
                    base.statistic = "BETA";
                    message["stat"] = "BETA";
                }
                else if (base.statistic.empty())
                {
                    bool found_or = false, found_beta = false;
                    for (size_t i = 0; i < token.size(); ++i) {
                        if (token[i].length() == 2
                            && toupper(token[i][0]) == 'O'
                            && toupper(token[i][1] == 'R'))
                        {
                            base.is_beta = false;
                            base.statistic = token[i];
                            message["stat"] = token[i];
                            found_or = true;
                        }
                        else if (token[i].length() == 4
                                 && toupper(token[i][0]) == 'B'
                                 && toupper(token[i][1]) == 'E'
                                 && toupper(token[i][2]) == 'T'
                                 && toupper(token[i][3]) == 'A')
                        {
                            base.is_beta = true;
                            base.statistic = token[i];
                            // Again, this will be problematic if the BETA
                            // is actually OR...
                            message["stat"] = token[i];
                            message["beta"] = "";
                            found_beta = true;
                        }
                        if (found_beta && found_or) {
                            error = true;
                            error_message.append(
                                "Error: Both OR and BETA "
                                "found in base file! We cannot determine "
                                "which statistic to use, please provide it "
                                "through --stat\n");
                            break;
                        }
                    }
                }
            }
            base.col_index[+BASE_INDEX::CHR] = index_check(base.chr, token);
            if (base.col_index[+BASE_INDEX::CHR] != -1)
                message["chr"] = base.chr;
            base.col_index[+BASE_INDEX::REF] =
                index_check(base.effect_allele, token);
            if (base.col_index[+BASE_INDEX::REF] != -1)
                message["A1"] = base.effect_allele;
            base.col_index[+BASE_INDEX::ALT] =
                index_check(base.non_effect_allele, token);
            if (base.col_index[+BASE_INDEX::ALT] != -1)
                message["A2"] = base.non_effect_allele;
            base.col_index[+BASE_INDEX::STAT] =
                index_check(base.statistic, token);
            if (base.col_index[+BASE_INDEX::STAT] != -1)
                message["stat"] = base.statistic;
            base.col_index[+BASE_INDEX::RS] = index_check(base.snp, token);
            if (base.col_index[+BASE_INDEX::RS] != -1)
                message["snp"] = base.snp;
            base.col_index[+BASE_INDEX::BP] = index_check(base.bp, token);
            if (base.col_index[+BASE_INDEX::BP] != -1) message["bp"] = base.bp;
            base.col_index[+BASE_INDEX::SE] =
                index_check(base.standard_error, token);
            if (base.col_index[+BASE_INDEX::SE] != -1)
                message["se"] = base.standard_error;
            base.col_index[+BASE_INDEX::P] = index_check(base.p_value, token);
            if (base.col_index[+BASE_INDEX::P] != -1)
                message["pvalue"] = base.p_value;

            if (!base.info_col.empty()) {
                std::vector<std::string> info = misc::split(base.info_col, ",");
                base.col_index[+BASE_INDEX::INFO] = index_check(info[0], token);
                if (info.size() != 2) {
                    error = true;
                    error_message.append("Error: Invalid format of "
                                         "--info-base. Should be "
                                         "ColName,Threshold.\n");
                }
                else
                {
                    try
                    {
                        base.info_score_threshold =
                            misc::convert<double>(info[1]);
                        if (base.info_score_threshold < 0
                            || base.info_score_threshold > 1)
                        {
                            error = true;
                            error_message.append(
                                "Error: Base INFO threshold must "
                                "be within 0 and 1!\n");
                        }
                        else
                        {
                            message["info-base"] = base.info_col;
                        }
                    }
                    catch (const std::runtime_error& er)
                    {
                        error = true;
                        error_message.append(
                            "Error: Invalid argument "
                            "passed to --info-base: "
                            + base.info_col
                            + "! Second argument must be numeric\n");
                    }
                }
            }
            // comma separate
            if (!base.maf_col.empty()) {
                std::string maf_error =
                    "Error: Invalid format of --maf-base. "
                    "Should be ColName,Threshold."
                    "or ColName,Threshold:ColName,Threshold.\n";
                std::vector<std::string> maf_type =
                    misc::split(base.maf_col, ":");
                if (maf_type.size() == 0 || maf_type.size() > 2) {
                    error = true;
                    error_message.append("Error: Currently only support at "
                                         "most 2 MAF filtering for base");
                }
                else
                {
                    std::vector<std::string> maf =
                        misc::split(maf_type[0], ",");
                    if (maf.size() != 2) {
                        error = true;
                        error_message.append(maf_error);
                    }
                    base.col_index[+BASE_INDEX::MAF] =
                        index_check(maf[0], token);
                    try
                    {
                        base.maf_control_threshold =
                            misc::convert<double>(maf[1]);
                        if (base.maf_control_threshold < 0
                            || base.maf_control_threshold > 1)
                        {
                            error = true;
                            error_message.append(
                                "Error: Base MAF threshold must "
                                "be within 0 and 1!\n");
                        }
                        message["maf-base"] = base.maf_col;
                    }
                    catch (const std::runtime_error& er)
                    {
                        error = true;
                        error_message.append(
                            "Error: Invalid argument passed to --maf-base: "
                            + base.maf_col + "! Threshold must be numeric\n");
                    }
                    if (maf_type.size() > 1) {
                        maf = misc::split(maf_type[1], ",");
                        if (maf.size() != 2) {
                            error = true;
                            error_message.append(maf_error);
                        }
                        base.col_index[+BASE_INDEX::MAF_CASE] =
                            index_check(maf[0], token);
                        try
                        {
                            base.maf_case_threshold =
                                misc::convert<double>(maf[1]);
                            if (base.maf_case_threshold < 0
                                || base.maf_case_threshold > 1)
                            {
                                error = true;
                                error_message.append(
                                    "Error: Base MAF threshold must "
                                    "be within 0 and 1!\n");
                            }
                            message["maf-base"] = base.maf_col;
                        }
                        catch (const std::runtime_error& er)
                        {
                            error = true;
                            error_message.append(
                                "Error: Invalid argument "
                                "passed to --maf-base: "
                                + base.maf_col
                                + "! Threshold must be numeric\n");
                        }
                    }
                }
            }
            // no default for MAF as there can be many different MAF
            // headers
        }
        else
        { // only required for index, as the defaults are in string
            if (base.provided_chr) {
                base.col_index[+BASE_INDEX::CHR] = index_check(
                    base.chr, max_size, error, error_message, "CHR");
            }
            if (base.provided_effect_allele) {
                base.col_index[+BASE_INDEX::REF] = index_check(
                    base.effect_allele, max_size, error, error_message, "REF");
            }
            if (base.provided_non_effect_allele) {
                base.col_index[+BASE_INDEX::ALT] =
                    index_check(base.non_effect_allele, max_size, error,
                                error_message, "ALT");
            }
            if (base.provided_bp) {
                base.col_index[+BASE_INDEX::BP] =
                    index_check(base.bp, max_size, error, error_message, "BP");
            }
            if (base.provided_standard_error) {
                base.col_index[+BASE_INDEX::SE] = index_check(
                    base.standard_error, max_size, error, error_message, "SE");
            }
            if (base.provided_info) {
                std::vector<std::string> info = misc::split(base.info_col, ",");
                base.col_index[+BASE_INDEX::INFO] = index_check(
                    info[0], max_size, error, error_message, "INFO");
                if (info.size() != 2) {
                    error = true;
                    error_message.append("Error: Invalid format of "
                                         "--info-base. Should be "
                                         "ColName,Threshold.\n");
                }
                try
                {
                    base.info_score_threshold = misc::convert<double>(info[1]);
                    if (base.info_score_threshold < 0
                        || base.info_score_threshold > 1)
                    {
                        error = true;
                        error_message.append("Error: Base INFO threshold "
                                             "must be within 0 and 1!\n");
                    }
                }
                catch (const std::runtime_error& er)
                {
                    error = true;
                    error_message.append(
                        "Error: Invalid argument passed to --info-base: "
                        + base.info_col
                        + "! Second argument must be numeric\n");
                }
            }
            if (!base.maf_col.empty()) {
                std::vector<std::string> maf = misc::split(base.maf_col, ",");
                base.col_index[+BASE_INDEX::MAF] =
                    index_check(maf[0], max_size, error, error_message, "MAF");
                std::string maf_error =
                    "Error: Invalid format of --maf-base. "
                    "Should be ColName,Threshold."
                    "or ColName,Threshold:ColName,Threshold.\n";
                std::vector<std::string> maf_type =
                    misc::split(base.maf_col, ":");
                if (maf_type.size() == 0 || maf_type.size() > 2) {
                    error = true;
                    error_message.append("Error: Currently only support at "
                                         "most 2 MAF filtering for base");
                }
                else
                {
                    std::vector<std::string> maf =
                        misc::split(maf_type[0], ",");
                    if (maf.size() != 2) {
                        error = true;
                        error_message.append(maf_error);
                    }
                    base.col_index[+BASE_INDEX::MAF] =
                        index_check(maf[0], token);
                    try
                    {
                        base.maf_control_threshold =
                            misc::convert<double>(maf[1]);
                        if (base.maf_control_threshold < 0
                            || base.maf_control_threshold > 1)
                        {
                            error = true;
                            error_message.append(
                                "Error: Base MAF threshold must "
                                "be within 0 and 1!\n");
                        }
                        message["maf-base"] = base.maf_col;
                    }
                    catch (const std::runtime_error& er)
                    {
                        error = true;
                        error_message.append(
                            "Error: Invalid argument passed to --maf-base: "
                            + base.maf_col + "! Threshold must be numeric\n");
                    }
                    if (maf_type.size() > 1) {
                        maf = misc::split(maf_type[1], ",");
                        if (maf.size() != 2) {
                            error = true;
                            error_message.append(maf_error);
                        }
                        base.col_index[+BASE_INDEX::MAF_CASE] =
                            index_check(maf[0], token);
                        try
                        {
                            base.maf_case_threshold =
                                misc::convert<double>(maf[1]);
                            if (base.maf_case_threshold < 0
                                || base.maf_case_threshold > 1)
                            {
                                error = true;
                                error_message.append(
                                    "Error: Base MAF threshold must "
                                    "be within 0 and 1!\n");
                            }
                            message["maf-base"] = base.maf_col;
                        }
                        catch (const std::runtime_error& er)
                        {
                            error = true;
                            error_message.append(
                                "Error: Invalid argument "
                                "passed to --maf-base: "
                                + base.maf_col
                                + "! Threshold must be numeric\n");
                        }
                    }
                }
            }
            base.col_index[+BASE_INDEX::P] =
                index_check(base.p_value, max_size, error, error_message, "P");
            base.col_index[+BASE_INDEX::STAT] = index_check(
                base.statistic, max_size, error, error_message, "STAT");
            base.col_index[+BASE_INDEX::RS] =
                index_check(base.snp, max_size, error, error_message, "RS");
        }

        // now check all required columns are here
        if (base.col_index[+BASE_INDEX::P] == -1) {
            // we can actually losen this requirement if user doesn't
            // perform clumping
            error = true;
            error_message.append("Error: No p-value column (" + base.p_value
                                 + ") in file!\n");
        }
        if (base.col_index[+BASE_INDEX::STAT] == -1) {
            error = true;
            error_message.append("Error: No statistic column (" + base.statistic
                                 + ") in file!\n");
        }
        if (base.col_index[+BASE_INDEX::RS] == -1) {
            error = true;
            error_message.append("Error: No SNP name column (" + base.snp
                                 + ") in file!\n");
        }
        if (base.col_index[+BASE_INDEX::REF] == -1) {
            error = true;
            error_message.append("Error: No Reference allele column ("
                                 + base.effect_allele + ") in file!\n");
        }

        double max_index =
            *max_element(base.col_index.begin(), base.col_index.end());
        base.col_index[+BASE_INDEX::MAX] = max_index;
    }
}

void Commander::clump_check(std::map<std::string, std::string>& message,
                            bool& error, std::string& error_message)
{
    if (!clumping.no_clump) {
        if (!reference_panel.keep_file.empty()
            && !reference_panel.remove_file.empty())
        {
            error = true;
            error_message.append(
                "Error: Can only use either --keep or --remove but not both\n");
        }
        // require clumping
        if (clumping.provided_proxy
            && (clumping.proxy < 0 || clumping.proxy > 1))
        {
            error = true;
            error_message.append(
                "Error: Proxy threshold must be within 0 and 1!\n");
        }
        if (clumping.p_value < 0.0 || clumping.p_value > 1.0) {
            error = true;
            error_message.append(
                "Error: P-value threshold must be within 0 and 1!\n");
        }
        if (clumping.r2 < 0.0 || clumping.r2 > 1.0) {
            error = true;
            error_message.append(
                "Error: R2 threshold must be within 0 and 1!\n");
        }
        if (!reference_panel.type.empty()) {
            bool alright = false;
            for (auto&& type : supported_types) {
                if (reference_panel.type.compare(type) == 0) {
                    alright = true;
                    break;
                }
            }
            if (!alright) {
                error = true;
                error_message.append("Error: Unsupported LD format: "
                                     + reference_panel.type + "\n");
            }
        }
        if (clumping.distance < 0.0) {
            error = true;
            error_message.append(
                "Error: Clumping distance must be positive!\n");
        }
        // if user does provide a input, this will just overwrite with the same
        // value This allow us to not keep the bool information
        message["clump-r2"] = std::to_string(clumping.r2);
        message["clump-p"] = std::to_string(clumping.p_value);
        message["clump-kb"] = std::to_string(clumping.distance);
        // now check the snp filtering
        // we automatically ignore any geno that are larger than 1
        // also output an error message

        if (reference_snp_filtering.geno != 0
            && (reference_snp_filtering.geno < 0
                || reference_snp_filtering.geno > 1))
        {
            error = true;
            error_message.append("Error: LD genotype missingness threshold "
                                 "must be larger than 0 and smaller than 1!\n");
        }
        if (reference_panel.type.compare("bgen") == 0
            || (reference_panel.file_name.empty()
                && reference_panel.multi_name.empty()
                && target.type.compare("bgen") == 0))
        {
            if (reference_snp_filtering.hard_threshold > 1
                || reference_snp_filtering.hard_threshold < 0)
            {
                error = true;
                error_message.append("Error: LD hard threshold must be larger "
                                     "than 0 and smaller than 1!\n");
            }
            else if (!reference_panel.file_name.empty()
                     || !reference_panel.multi_name.empty())
            {
                // this is to be consistent where ld parameter won't
                // apply to target file
                message["ld-hard-thres"] =
                    std::to_string(reference_snp_filtering.hard_threshold);
            }
            else
            {
                // we use the prs_snp_filtering.hard_threshold
                // such that if user has provided this function, this function
                // will still provide the correct output
                message["hard-thres"] =
                    std::to_string(prs_snp_filtering.hard_threshold);
            }
        }
        if (reference_snp_filtering.maf != 0
            && (reference_snp_filtering.maf > 1
                || reference_snp_filtering.maf < 0))
        {
            error = true;
            error_message.append("Error: LD MAF threshold must be larger than "
                                 "0 and smaller than 1!\n");
        }
        if (reference_snp_filtering.info_score < 0
            || reference_snp_filtering.info_score > 1)
        {
            error = true;
            error_message.append("Error: LD INFO score threshold must be "
                                 "larger than 0 and smaller than 1!\n");
        }
    }
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
        cov.erase(0, 1);
        open = misc::split(cov, "[");
        for (auto o : open) {
            if (o.find("]") != std::string::npos) {
                close = misc::split(o, "]");
                // the first one will always be the list
                info.push_back(close[0]);
                list.push_back(true);
                // Nested List is not supported
                for (size_t cl = 1; cl < close.size(); ++cl) {
                    info.push_back(close[cl]);
                    list.push_back(false);
                }
            }
            else
            {
                info.push_back(o);
                list.push_back(false);
            }
        }

        for (size_t c = 0; c < info.size(); ++c) {
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
                            for (size_t s = start; s <= end; ++s) {
                                numeric.push_back(s);
                            }
                        }
                        catch (const std::runtime_error& error)
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
                        catch (const std::runtime_error& error)
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
                    size_t cur_size = final_covariates.size();
                    for (size_t final = 0; final < cur_size; ++final) {
                        std::string cur = final_covariates[final];
                        final_covariates[final].append(
                            std::to_string(numeric.front()));
                        for (size_t s = 1; s < numeric.size(); ++s) {
                            final_covariates.push_back(
                                cur + std::to_string(numeric[s]));
                        }
                    }
                }
            }
            else
            {
                for (size_t final = 0; final < final_covariates.size(); ++final)
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
        final_covariates.push_back(cov);
    }
    return final_covariates;
}

void Commander::covariate_check(bool& error, std::string& error_message)
{
    if (covariate.file_name.empty()) return;
    // first, transform all the covariates
    std::unordered_set<std::string> included;
    std::vector<std::string> transformed_cov;
    for (auto cov : covariate.covariates) {
        if (cov.empty()) continue;
        // got annoyed with the input of PC.1 PC.2 PC.3, do this automatic
        // thingy to substitute them
        transformed_cov = transform_covariate(cov);
        for (auto&& trans : transformed_cov) {
            included.insert(trans);
        }
    }
    std::ifstream cov_file;
    cov_file.open(covariate.file_name.c_str());
    if (!cov_file.is_open()) {
        error = true;
        error_message.append(
            "Error: Cannot open covariate file: " + covariate.file_name + "\n");
        return;
    }
    std::string line;
    std::getline(cov_file, line);
    if (line.empty()) {
        error = true;
        error_message.append("Error: First line of covariate file is empty!\n");
        return;
    }
    cov_file.close();
    std::vector<std::string> cov_header = misc::split(line);
    std::string missing = "";
    std::unordered_map<std::string, uint32_t> ref_index;
    for (size_t i = 0; i < cov_header.size(); ++i) {
        ref_index[cov_header[i]] = i;
    }
    if (covariate.covariates.size() == 0) {
        // add all headers to the covariate list
        for (size_t i = (1 + !misc.ignore_fid); i < cov_header.size(); ++i) {
            included.insert(cov_header[i]);
        }
    }
    size_t valid_cov = 0;
    for (auto&& cov : included) {
        if (ref_index.find(cov) != ref_index.end()) {
            covariate.covariate_index.push_back(ref_index[cov]);
            covariate.covariates.push_back(cov);
            valid_cov++;
        }
        else if (missing.empty())
            missing = cov;
        else
            missing.append("," + cov);
    }
    if (!missing.empty()) {
        error_message.append("Warning: Covariate(s) missing from file: "
                             + missing + ". Header of file is: " + line + "\n");
    }
    if (valid_cov == 0) {
        error = true;
        error_message.append("Error: No valid Covariate!\n");
    }
    covariate.covariates.clear();
    std::sort(covariate.covariate_index.begin(),
              covariate.covariate_index.end());
    for (auto&& c : covariate.covariate_index) {
        covariate.covariates.push_back(cov_header[c]);
    }
    for (auto cov : covariate.factor_covariates) {
        if (cov.empty()) continue;
        transformed_cov = transform_covariate(cov);
        for (auto&& trans : transformed_cov) {
            if (included.find(trans) != included.end()) {
                covariate.factor_index.push_back(ref_index[trans]);
            }
            else
            {
                error = true;
                error_message.append("Error: All factor covariates must be "
                                     "found in covariate list. "
                                     + trans + " not found in covariate list");
            }
        }
    }
    std::sort(covariate.factor_index.begin(), covariate.factor_index.end());
}


void Commander::filter_check(bool& error, std::string& error_message)
{
    if (target.type.compare("bgen") == 0
        && (prs_snp_filtering.hard_threshold <= 0
            || prs_snp_filtering.hard_threshold >= 1))
    {
        error = true;
        error_message.append(
            "Error: Hard threshold must be between 0 and 1!\n");
    }
    if (!prs_snp_filtering.extract_file.empty()
        && !prs_snp_filtering.exclude_file.empty())
    {
        error = true;
        error_message.append(
            "Error: Can only use --extract or --exclude but not both\n");
    }

    if ((prs_snp_filtering.info_score < 0 || prs_snp_filtering.info_score > 1))
    {
        error = true;
        error_message.append(
            "Error: INFO score threshold cannot be bigger than 1.0 "
            "or smaller than 0.0\n");
    }
    if ((prs_snp_filtering.geno < 0 || prs_snp_filtering.geno > 1)) {
        error = true;
        error_message.append(
            "Error: Genotype missingness threshold cannot be bigger than 1.0 "
            "or smaller than 0.0\n");
    }
    if ((prs_snp_filtering.maf < 0 || prs_snp_filtering.maf > 1)) {
        error = true;
        error_message.append("Error: MAF threshold cannot be bigger than 1.0 "
                             "or smaller than 0.0\n");
    }
}

void Commander::misc_check(std::map<std::string, std::string>& message,
                           bool& error, std::string& error_message)
{
    if (misc.permutation < 0) {
        error = true;
        error_message.append("Error: Negative number of permutation!\n");
    }
    if (!misc.provided_seed) {
        misc.seed = std::random_device()();
        message["seed"] = std::to_string(misc.seed);
    }
    if (misc.thread <= 0) {
        error = true;
        error_message.append("Error: Number of thread must be larger than 1\n");
    }
    if (misc.permutation <= 0 && misc.logit_perm) {
        error_message.append(
            "Warning: Permutation not required, --logit-perm has no effect\n");
    }
    if (prs_calculation.no_regress) misc.print_all_scores = true;
    if (misc.thread == 1) message["thread"] = "1";
    message["out"] = misc.out;
}

void Commander::prset_check(std::map<std::string, std::string>& message,
                            bool& error, std::string& error_message)
{
    if (!prset.perform_prset) return;
    if (prset.gtf.empty() && !prset.msigdb.empty()) {
        error = true;
        error_message.append(
            "Error: Must provide a gtf file if msigdb is specified\n");
    }
    if (prset.feature.empty()) {
        prset.feature.push_back("exon");
        prset.feature.push_back("gene");
        prset.feature.push_back("protein_coding");
        prset.feature.push_back("CDS");
        message["feature"] = "exon,gene,protein_coding,CDS";
    }
    // don't check file exist here?
    if (prset.set_perm <= 0) {
        error = true;
        error_message.append(
            "Error: Negative number of set permutation provided!");
    }
    if (prset.set_perm != 0 && misc.permutation > 0) {
        error = true;
        error_message.append("Error: Currently only support either set-base "
                             "permutation (for competitive p-value) or PRSice "
                             "base permutation (--perm)");
    }
}


void Commander::prsice_check(std::map<std::string, std::string>& message,
                             bool& error, std::string& error_message)
{

    switch (prs_calculation.model)
    {
    case MODEL::ADDITIVE: message["model"] = "add"; break;
    case MODEL::DOMINANT: message["model"] = "dom"; break;
    case MODEL::RECESSIVE: message["model"] = "rec"; break;
    case MODEL::HETEROZYGOUS: message["model"] = "het"; break;
    default: error = true; error_message.append("Error: Unrecognized model!");
    }
    if (p_thresholds.barlevel.size() == 0 && !prset.perform_prset) {
        // always output the bar_message so that we can tell R what to do next
        p_thresholds.barlevel = {0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
        if (!p_thresholds.no_full) p_thresholds.barlevel.push_back(1);
    }
    if (prset.perform_prset) {
        if (!p_thresholds.set_use_thresholds && !p_thresholds.fastscore) {
            // if user use fastscore or provided any threshold, then we will
            // not kick in this default behaviour
            message["bar-levels"] = 1;
            p_thresholds.fastscore = true;
            p_thresholds.barlevel = {1};
        }
        else if (p_thresholds.barlevel.size() == 0)
        {
            // return to default of PRSice otherwise
            p_thresholds.barlevel = {0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
            if (!p_thresholds.no_full) p_thresholds.barlevel.push_back(1);
        }
    }
    else if (!p_thresholds.fastscore)
    {
        // deleted the no-regress check. If someone want to output
        // a large file, then so be it
        if (p_thresholds.inter <= 0) {
            error = true;
            error_message.append("Error: Cannot have negative interval!\n");
        }
        if (p_thresholds.upper < p_thresholds.lower) {
            error = true;
            error_message.append(
                "Error: Upper bound must be larger than lower bound!\n");
        }
        if (p_thresholds.upper < 0.0 || p_thresholds.lower < 0.0) {
            error = true;
            error_message.append("Error: Cannot have negative bounds!\n");
        }

        message["interval"] = std::to_string(p_thresholds.inter);
        message["lower"] = std::to_string(p_thresholds.lower);
        message["upper"] = std::to_string(p_thresholds.upper);
    }
    if (!p_thresholds.no_full) p_thresholds.barlevel.push_back(1);
    std::sort(p_thresholds.barlevel.begin(), p_thresholds.barlevel.end());
    p_thresholds.barlevel.erase(
        std::unique(p_thresholds.barlevel.begin(), p_thresholds.barlevel.end()),
        p_thresholds.barlevel.end());
    std::string bar_message = "";
    for (auto&& b : p_thresholds.barlevel) {
        if (bar_message.empty())
            bar_message.append(misc::to_string(b));
        else
            bar_message.append("," + misc::to_string(b));
    }
    message["bar-levels"] = bar_message;
}

void Commander::prslice_check(bool& error, std::string& error_message)
{
    if (prslice.provided) {
        if (misc.print_all_scores) {
            error = true;
            error_message.append("Error: Cannot output PRS for all threshold "
                                 "when using PRSlice!\n");
        }
        if (prslice.size <= 0) {
            error = true;
            error_message.append(
                "Error: PRSlice size cannot be less than 1!\n");
        }
    }
}

void Commander::target_check(std::map<std::string, std::string>& message,
                             bool& error, std::string& error_message)
{
    if (target.name.empty() && target.multi_name.empty()) {
        error = true;
        error_message.append("Error: You must provide a target file or a file "
                             "containing all target prefixs!\n");
    }
    if (!target.keep_file.empty() && !target.remove_file.empty()) {
        error = true;
        error_message.append(
            "Error: Can only use either --keep or --remove but not both\n");
    }
    bool alright = false;
    for (auto&& type : supported_types) {
        if (target.type.compare(type) == 0) {
            alright = true;
            break;
        }
    }
    if (!alright) {
        error = true;
        error_message.append("Error: Unsupported target format: " + target.type
                             + "\n");
    }
    if (target.type.compare("bgen") == 0 && prs_snp_filtering.is_hard_coded) {
        message["hard-thres"] =
            std::to_string(prs_snp_filtering.hard_threshold);
    }


    if (target.pheno_col.size() != 0 && target.pheno_file.empty()) {
        error = true;
        error_message.append("Error: You must provide a phenotype file for "
                             "multiple phenotype analysis");
    }
    if (target.pheno_file.empty() && target.is_binary.empty()) {
        if (base.is_beta) {
            message["binary-target"] = "F";
            target.is_binary.push_back(false);
        }
        else
        {
            message["binary-target"] = "T";
            target.is_binary.push_back(true);
        }
    }
    else
    {
        if (target.pheno_col.empty() && target.is_binary.size() == 1) {
            // this is ok
        }
        else if (target.pheno_col.empty() && target.is_binary.empty())
        {
            if (base.is_beta) {
                message["binary-target"] = "F";
                target.is_binary.push_back(false);
            }
            else
            {
                message["binary-target"] = "T";
                target.is_binary.push_back(true);
            }
        }
        else if (target.pheno_col.size() <= 1 && target.is_binary.empty())
        {
            if (base.is_beta) {
                message["binary-target"] = "F";
                target.is_binary.push_back(false);
            }
            else
            {
                message["binary-target"] = "T";
                target.is_binary.push_back(true);
            }
        }
        else if (target.pheno_col.size() != target.is_binary.size())
        {
            error = true;
            error_message.append("Error: Number of target phenotypes doesn't "
                                 "match information of binary target! You must "
                                 "indicate whether the phenotype is binary "
                                 "using --binary-target\n");
        }
    }

    size_t num_bin = 0;
    for (auto&& binary : target.is_binary) {
        if (binary) num_bin++;
    }
    if (!target.prevalence.empty()
        && num_bin > target.prevalence.size()) // need to be all or nothing
    {
        error = true;
        error_message.append(
            "Error: Number of target prevalence doesn't match "
            "number of binary traits. You must provide a prevalence for all "
            "binary trait(s) or not provide any prevalence (all or nothing)\n");
    }
    for (auto&& prev : target.prevalence) {
        if (prev > 1.0 || prev < 0.0) {
            error = true;
            error_message.append("Error: Prevalence cannot be bigger than 1.0 "
                                 "or smaller than 0.0\n");
            break;
        }
    }
}

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

bool Commander::process(int argc, char* argv[], const char* optString,
                        const struct option longOpts[], Reporter& reporter)
{
    int longIndex = 0;
    int opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    std::string command;
    std::map<std::string, std::string> message_store;
    std::string error_messages = "";
    std::string temp_string = "";
    size_t temp_int = 0;
    bool dummy = false;
    bool error = false;
    while (opt != -1) {
        switch (opt)
        {
        case 0:
            command = longOpts[longIndex].name;
            if (longOpts[longIndex].flag != 0)
                break;
            else if (command.compare("chr") == 0)
                set_string(optarg, message_store, base.chr, base.provided_chr,
                           command, error_messages);
            else if (command.compare("A1") == 0)
                set_string(optarg, message_store, base.ref_allele,
                           base.provided_ref, command, error_messages);
            else if (command.compare("A2") == 0)
                set_string(optarg, message_store, base.alt_allele,
                           base.provided_alt, command, error_messages);
            else if (command.compare("stat") == 0)
                set_string(optarg, message_store, base.statistic,
                           base.provided_stat, command, error_messages);
            else if (command.compare("snp") == 0)
                set_string(optarg, message_store, base.snp, base.provided_snp,
                           command, error_messages);
            else if (command.compare("bp") == 0)
                set_string(optarg, message_store, base.bp, base.provided_bp,
                           command, error_messages);
            else if (command.compare("se") == 0)
                set_string(optarg, message_store, base.standard_error,
                           base.provided_se, command, error_messages);
            else if (command.compare("model") == 0)
                set_model(optarg, message_store, error_messages, error);
            else if (command.compare("cov-header")
                     == 0) // cerr for backward compatibility
                load_string_vector(optarg, message_store, covariate.covariates,
                                   "cov-col", error_messages);
            else if (command.compare("keep") == 0)
                set_string(optarg, message_store, target.keep_file,
                           target.keep_sample, command, error_messages);
            else if (command.compare("exclude") == 0)
                set_string(optarg, message_store, filter.exclude_file,
                           filter.exclude, command, error_messages);
            else if (command.compare("extract") == 0)
                set_string(optarg, message_store, filter.extract_file,
                           filter.extract, command, error_messages);
            else if (command.compare("ld-keep") == 0)
                set_string(optarg, message_store, clumping.keep_file,
                           clumping.keep_sample, command, error_messages);
            else if (command.compare("ld-remove") == 0)
                set_string(optarg, message_store, clumping.remove_file,
                           clumping.remove_sample, command, error_messages);
            else if (command.compare("remove") == 0)
                set_string(optarg, message_store, target.remove_file,
                           target.remove_sample, command, error_messages);
            else if (command.compare("ld-type") == 0)
                set_string(optarg, message_store, clumping.type,
                           clumping.use_type, command, error_messages);
            else if (command.compare("maf-base") == 0)
                set_string(optarg, message_store, base.maf, base.provided_maf,
                           command, error_messages);
            else if (command.compare("type") == 0)
                set_string(optarg, message_store, target.type, target.use_type,
                           command, error_messages);
            else if (command.compare("score") == 0)
                set_string(optarg, message_store, prsice.missing_score, dummy,
                           command, error_messages);
            else if (command.compare("hard-thres") == 0)
                set_numeric<double>(optarg, message_store, error_messages,
                                    filter.hard_threshold,
                                    filter.use_hard_thres, error, command);
            else if (command.compare("info") == 0)
                set_numeric<double>(optarg, message_store, error_messages,
                                    filter.info_score, filter.info_filtering,
                                    error, command);
            else if (command.compare("clump-p") == 0)
                set_numeric<double>(optarg, message_store, error_messages,
                                    clumping.p_value, clumping.provide_p, error,
                                    command);
            else if (command.compare("clump-r2") == 0)
                set_numeric<double>(optarg, message_store, error_messages,
                                    clumping.r2, clumping.provide_r2, error,
                                    command);
            else if (command.compare("clump-kb") == 0)
                set_numeric<int>(optarg, message_store, error_messages,
                                 clumping.distance, clumping.provide_distance,
                                 error, command);
            else if (command.compare("prslice") == 0)
                set_numeric<int>(optarg, message_store, error_messages,
                                 prslice.size, prslice.provided, error,
                                 command);
            else if (command.compare("proxy") == 0)
                set_numeric<double>(optarg, message_store, error_messages,
                                    clumping.proxy, clumping.provide_proxy,
                                    error, command);
            else if (command.compare("perm") == 0)
                set_numeric<int>(optarg, message_store, error_messages,
                                 misc.permutation, misc.provided_permutation,
                                 error, command);
            else if (command.compare("binary-target") == 0)
                load_binary_vector(optarg, message_store, error_messages,
                                   target.is_binary, error, command);
            else if (command.compare("pheno-col") == 0)
                load_string_vector(optarg, message_store, target.pheno_col,
                                   command, error_messages);
            else if (command.compare("feature") == 0)
                load_string_vector(optarg, message_store, prset.feature,
                                   command, error_messages);
            else if (command.compare("bar-levels") == 0)
                load_numeric_vector<double>(optarg, message_store,
                                            error_messages, prsice.barlevel,
                                            error, command);
            else if (command.compare("memory") == 0)
                set_numeric<int>(optarg, message_store, error_messages,
                                 misc.memory, misc.provided_memory, error,
                                 command);
            else if (command.compare("info-base") == 0)
                set_string(optarg, message_store, base.info_col, base.use_info,
                           command, error_messages);

            else
            {
                std::string er = "Undefined operator: " + command
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
            set_string(optarg, message_store, covariate.name, dummy, "cov-file",
                       error_messages);
            break;
        case 'f':
            set_string(optarg, message_store, target.pheno_file, dummy,
                       "pheno-file", error_messages);
            break;
        case 'g':
            set_string(optarg, message_store, prset.gtf, prset.perform_prset,
                       "gtf", error_messages);
            break;
        case 'i':
            set_numeric<double>(optarg, message_store, error_messages,
                                prsice.inter, prsice.provide_inter, error,
                                "interval");
            break;
        case 'k':
            load_numeric_vector<double>(optarg, message_store, error_messages,
                                        target.prevalence, error, "prevalence");
            break;
        case 'l':
            set_numeric<double>(optarg, message_store, error_messages,
                                prsice.lower, prsice.provide_lower, error,
                                "lower");
            break;
        case 'L':
            set_string(optarg, message_store, clumping.ld, dummy, "ld",
                       error_messages);
            break;
        case 'm':
            set_string(optarg, message_store, prset.msigdb, prset.perform_prset,
                       "msigdb", error_messages);
            break;
        case 'n':
            temp_string = optarg;
            temp_int = std::thread::hardware_concurrency();
            if (temp_string.compare("max") == 0) {
                misc.thread = temp_int;
                misc.provide_thread = true;
            }
            else
            {
                temp_int = 1;
                try
                {
                    temp_int = misc::convert<int>(optarg);
                }
                catch (const std::runtime_error& er)
                {
                    error_messages.append(
                        "ERROR: Non numeric argument passed to thread: "
                        + std::string(optarg) + "!\n");
                }
                if (temp_int > std::thread::hardware_concurrency()) {
                    misc.thread = std::thread::hardware_concurrency();
                    misc.provide_thread = true;
                }
                else
                {
                    misc.thread = temp_int;
                    misc.provide_thread = true;
                }
            }

            if (message_store.find("thread") != message_store.end()) {
                error_messages.append(
                    "Warning: Duplicated argument --thread\n");
            }
            message_store["thread"] = std::to_string(misc.thread);
            break;
        case 'o':
            set_string(optarg, message_store, misc.out, misc.provided_output,
                       "out", error_messages);
            break;
        case 'p':
            set_string(optarg, message_store, base.p_value, base.provided_p,
                       "pvalue", error_messages);
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
            set_numeric<double>(optarg, message_store, error_messages,
                                prsice.upper, prsice.provide_upper, error,
                                "upper");
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
            throw "Undefined operator, please use --help for more information!";
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
            "ERROR: PRSet and PRSlice cannot be performed together!\n");
    }
    // check all flags
    std::string log_name = misc.out + ".log";
    reporter.initiailize(log_name);
    if (base.beta) message_store["beta"] = "";
    if (base.index) message_store["index"] = "";
    if (clumping.no_clump) message_store["no-clump"] = "";
    if (filter.hard_coding) message_store["hard"] = "";
    if (filter.keep_ambig) message_store["keep-ambig"] = "";
    if (misc.all) message_store["all"] = "";
    if (misc.ignore_fid) message_store["ignore-fid"] = "";
    if (misc.logit_perm) message_store["logit-perm"] = "";
    if (misc.print_snp) message_store["print-snp"] = "";
    if (prsice.fastscore) message_store["fastscore"] = "";
    if (prsice.full) message_store["full"] = "";
    if (prsice.no_regress) message_store["no-regress"] = "";
    if (target.nonfounders) message_store["nonfounders"] = "";
    if ((clumping.ld.empty() && target.type.compare("bgen") == 0)
        || clumping.type.compare("bgen") == 0)
    {
        if (!filter.use_hard_thres) {
            message_store["hard-thres"] = std::to_string(filter.hard_threshold);
        }
    }
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
    message.append("If you use PRSice in any publised work, please cite:\n");
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


Commander::Commander()
{

    base.beta = false;
    base.name = "";
    base.chr = "CHR";
    base.ref_allele = "A1";
    base.alt_allele = "A2";
    base.statistic = "OR";
    base.snp = "SNP";
    base.bp = "BP";
    base.standard_error = "SE";
    base.p_value = "P";
    base.info_col = "INFO,0.9";
    base.maf = "";
    base.info_score = 0.9;
    base.index = false;
    base.provided_chr = false;
    base.provided_ref = false;
    base.provided_alt = false;
    base.provided_maf = false;
    base.provided_stat = false;
    base.provided_snp = false;
    base.provided_bp = false;
    base.provided_se = false;
    base.provided_p = false;
    base.use_info = false;
    base.col_index.resize(+BASE_INDEX::MAX + 1, -1);


    clumping.distance = 250;
    clumping.keep_sample = false;
    clumping.ld = "";
    clumping.no_clump = false;
    clumping.provide_proxy = false;
    clumping.proxy = -1.0;
    clumping.p_value = 1.0;
    clumping.r2 = 0.1;
    clumping.remove_sample = false;
    clumping.type = "bed";
    clumping.provide_p = false;
    clumping.provide_r2 = false;
    clumping.provide_distance = false;

    covariate.name = "";
    covariate.ancestry_dim = "MDS";

    filter.exclude = false;
    filter.extract = false;
    filter.geno = 0.0;
    filter.mind = 0.0;
    filter.maf = 0.01;
    filter.hard_coding = false;
    filter.hard_threshold = 0.9;
    filter.info_filtering = false;
    filter.info_score = 0.9;
    filter.keep_ambig = false;
    filter.use_prob = false;
    filter.use_maf = false;
    filter.use_mind = false;
    filter.use_hard_thres = false;
    filter.use_geno = false;

    misc.all = false;
    misc.ignore_fid = false;
    misc.logit_perm = false;
    misc.memory = 10000;
    misc.out = "PRSice";
    misc.permutation = 10000;
    misc.print_snp = false;
    misc.print_all_samples = false;
    misc.provided_permutation = false;
    misc.provided_output = false;
    misc.provided_seed = false;
    misc.seed = 0;
    misc.thread = 1;
    misc.provide_thread = false;

    prset.gtf = "";
    prset.msigdb = "";
    prset.perform_prset = false;

    prsice.missing_score = "";
    prsice.model = +MODEL::ADDITIVE;
    prsice.lower = 0.0001;
    prsice.upper = 0.5;
    prsice.inter = 0.00005;
    prsice.fastscore = false;
    prsice.no_regress = false;
    prsice.full = false;
    prsice.provide_lower = false;
    prsice.provide_upper = false;
    prsice.provide_inter = false;

    prslice.size = -1;
    prslice.provided = false;

    species.num_auto = 22;
    species.no_x = false;
    species.no_y = false;
    species.no_xy = false;
    species.no_mt = false;
    species.double_set = false;

    target.remove_sample = false;
    target.keep_sample = false;
    target.nonfounders = false;
    target.name = "";
    target.pheno_file = "";
    target.type = "bed";
    info();
}

bool Commander::init(int argc, char* argv[], Reporter& reporter)
{
    if (argc <= 1) {
        usage();
        throw std::runtime_error("Please provide the required parameters");
    }
    const char* optString = "a:b:B:c:C:f:g:i:l:L:m:n:o:p:t:u:h?";
    const struct option longOpts[] = {
        // parameters with short flags
        {"ancestry", required_argument, NULL, 'a'},
        {"base", required_argument, NULL, 'b'},
        {"bed", required_argument, NULL, 'B'},
        {"cov-col", required_argument, NULL, 'c'},
        {"cov-header", required_argument, NULL,
         0}, // retain here for backward compatibility
        {"cov-file", required_argument, NULL, 'C'},
        {"pheno-file", required_argument, NULL, 'f'},
        {"gtf", required_argument, NULL, 'g'},
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
        // flags, only need to set them to true
        {"all", no_argument, &misc.all, 1},
        {"beta", no_argument, &base.beta, 1},
        {"full", no_argument, &prsice.full, 1},
        {"hard", no_argument, &filter.hard_coding, 1},
        {"ignore-fid", no_argument, &misc.ignore_fid, 1},
        {"index", no_argument, &base.index, 1},
        {"keep-ambig", no_argument, &filter.keep_ambig, 0},
        {"logit-perm", no_argument, &misc.logit_perm, 1},
        {"no-clump", no_argument, &clumping.no_clump, 1},
        {"no-regression", no_argument, &prsice.no_regress, 1},
        {"no-x", no_argument, &species.no_x, 1},
        {"no-y", no_argument, &species.no_y, 1},
        {"no-xy", no_argument, &species.no_xy, 1},
        {"no-mt", no_argument, &species.no_mt, 1},
        {"nonfounders", no_argument, &target.nonfounders, 1},
        {"fastscore", no_argument, &prsice.fastscore, 1},
        {"print-snp", no_argument, &misc.print_snp, 1},
        // long flags, need to work on them
        {"A1", required_argument, NULL, 0},
        {"A2", required_argument, NULL, 0},
        {"bar-levels", required_argument, NULL, 0},
        {"binary-target", required_argument, NULL, 0},
        {"bp", required_argument, NULL, 0},
        {"chr", required_argument, NULL, 0},
        {"clump-kb", required_argument, NULL, 0},
        {"clump-p", required_argument, NULL, 0},
        {"clump-r2", required_argument, NULL, 0},
        {"exclude", required_argument, NULL, 0},
        {"extract", required_argument, NULL, 0},
        {"feature", required_argument, NULL, 0},
        {"hard-thres", required_argument, NULL, 0},
        {"info-base", required_argument, NULL, 0},
        {"info", required_argument, NULL, 0},
        {"keep", required_argument, NULL, 0},
        {"ld-keep", required_argument, NULL, 0},
        {"ld-type", required_argument, NULL, 0},
        {"ld-remove", required_argument, NULL, 0},
        {"maf-base", required_argument, NULL, 0},
        {"memory", required_argument, NULL, 0},
        {"model", required_argument, NULL, 0},
        {"num-auto", required_argument, NULL, 0},
        {"perm", required_argument, NULL, 0},
        {"pheno-col", required_argument, NULL, 0},
        {"proxy", required_argument, NULL, 0},
        {"prslice", required_argument, NULL, 0},
        {"remove", required_argument, NULL, 0},
        {"score", required_argument, NULL, 0},
        {"se", required_argument, NULL, 0},
        {"snp", required_argument, NULL, 0},
        {"stat", required_argument, NULL, 0},
        {"type", required_argument, NULL, 0},
        // species flag
        {"cow", no_argument, NULL, 0},
        {"dog", no_argument, NULL, 0},
        {"horse", no_argument, NULL, 0},
        {"mouse", no_argument, NULL, 0},
        {"rice", no_argument, NULL, 0},
        {"sheep", no_argument, NULL, 0},
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'v'},
        {NULL, 0, 0, 0}};
    return process(argc, argv, optString, longOpts, reporter);
}


Commander::~Commander()
{
    // dtor
}

void Commander::info()
{
    help_message =
        "usage: PRSice [options] <-b base_file> <-t target_file>\n"
        "\nBase File:\n"
        "    --base          | -b    Base association file\n"
        "    --beta                  Whether the test statistic is in the form "
        "of \n"
        "                            BETA or OR. If set, test statistic is "
        "assume\n"
        "                            to be in the form of BETA.\n"
        "    --A1                    Column header containing allele 1 (effect "
        "allele)\n"
        "                            Default: A1\n"
        "    --A2                    Column header containing allele 2 "
        "(reference allele)\n"
        "                            Default: A2\n"
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
        "    --info-base             Base INFO score filtering. Format should "
        "be\n"
        "                            <Column name>,<Threshold>. SNPs with info "
        "\n"
        "                            score less than <Threshold> will be "
        "ignored\n"
        "                            Column name default: INFO\n"
        "                            Threshold default: 0.9\n"
        "    --maf-base              Base MAF filtering. Format should be\n"
        "                            <Column name>,<Threshold>. SNPs with maf "
        "\n"
        "                            less than <Threshold> will be ignored\n"
        "    --pvalue        | -p    Column header containing the p-value\n"
        "                            Default: P\n"
        "    --se                    Column header containing the standard "
        "error\n"
        "                            Default: SE\n"
        "    --snp                   Column header containing the SNP ID\n"
        "                            Default: SNP\n"
        "    --stat                  Column header containing the summary "
        "statistic\n"
        "                            If --beta is set, default as BETA. "
        "Otherwise,\n"
        "                            will search for OR or BETA from the "
        "header\n"
        "                            of the base file\n"
        "\nClumping:\n"
        "    --clump-kb              The distance for clumping in kb\n"
        "                            Default: "
        + std::to_string(clumping.distance)
        + "\n"
          "    --clump-r2              The R2 threshold for clumping\n"
          "                            Default: "
        + std::to_string(clumping.r2)
        + "\n"
          "    --clump-p               The p-value threshold use for "
          "clumping.\n"
          "                            Default: "
        + std::to_string(clumping.p_value)
        + "\n"
          "    --ld            | -L    LD reference file. Use for LD "
          "calculation. If not\n"
          "                            provided, will use the post-filtered "
          "target genotype\n"
          "                            for LD calculation. Support multiple "
          "chromosome input\n"
          "                            Please see --target for more "
          "information\n"
          "    --ld-keep               File containing the sample(s) to be "
          "extracted from\n"
          "                            the LD reference file. First column "
          "should be FID and\n"
          "                            the second column should be IID. If "
          "--ignore-fid is\n"
          "                            set, first column should be IID\n"
          "                            Mutually exclusive from --ld-remove\n"
          "                            No effect if --ld was not provided\n"
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
          "\nCovariate:\n"
          "    --cov-file      | -C    Covariate file. First column should be "
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
          "\nDosage:\n"
          "    --hard-thres            Hard threshold for dosage data. Any "
          "call less than\n"
          "                            this will be treated as missing. Note "
          "that if dosage\n"
          "                            data is used as a LD reference, it will "
          "always be\n"
          "                            hard coded to calculate the LD\n"
          "                            Default: "
        + std::to_string(filter.hard_threshold)
        + "\n"
          "    --hard                  Use hard coding instead of dosage for "
          "PRS construction.\n"
          "                            Default is to use dosage instead of "
          "hard coding\n"
          "\nPRSet:\n"
          "    --bed           | -B    Bed file containing the selected "
          "regions.\n"
          "                            Name of bed file will be used as the "
          "region\n"
          "                            identifier. WARNING: Bed file is "
          "0-based\n"
          "    --feature               Feature(s) to be included from the gtf "
          "file.\n"
          "                            Default: exon,CDS,gene,protein_coding.\n"
          "    --gtf           | -g    GTF file containing gene boundaries. "
          "Required\n"
          "                            when --msigdb is used\n"
          "    --msigdb        | -m    MSIGDB file containing the pathway "
          "information.\n"
          "                            Require the gtf file\n"
          "\nPRSice:\n"
          "    --bar-levels            Level of barchart to be plotted. When "
          "--fastscore\n"
          "                            is set, PRSice will only calculate the "
          "PRS for \n"
          "                            threshold within the bar level. Levels "
          "should be\n"
          "                            comma separated without space\n"
          "    --fastscore             Only calculate threshold stated in "
          "--bar-levels\n"
          "    --full                  Include the full model in the analysis\n"
          "    --interval      | -i    The step size of the threshold. "
          "Default: "
        + std::to_string(prsice.inter)
        + "\n"
          "    --lower         | -l    The starting p-value threshold. "
          "Default: "
        + std::to_string(prsice.lower)
        + "\n"
          "    --model                 Genetic model use for regression. The "
          "genetic\n"
          "                            encoding is based on the base data "
          "where the\n"
          "                            encoding represent number of the "
          "effective allele\n"
          "Available"
          "                            models include:\n"
          "                            add - Additive model, code as 0/1/2 "
          "(default)\n"
          "                            dom - Dominant model, code as 0/1/1\n"
          "                            rec - Recessive model, code as 0/0/1\n"
          "                            het - Heterozygous only model, code as "
          "0/1/0\n"
          "    --no-regress            Do not perform the regression analysis "
          "and simply\n"
          "                            output all PRS.\n"
          "    --score                 Method to handle missing genotypes. By "
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
          "                            use the 'no_mean_imputation' modifier. "
          "Alternatively,\n"
          "                            you can use the 'center' modifier to "
          "shift all scores\n"
          "                            to mean zero. \n"
          "    --upper         | -u    The final p-value threshold. Default: "
        + std::to_string(prsice.upper)
        + "\n"
          "\nPRSlice:\n"
          "    --prslice               Perform PRSlice where the whole genome "
          "is first cut\n"
          "                            into bin size specified by this option. "
          "PRSice will\n"
          "                            then be performed on each bin. Bins are "
          "then sorted\n"
          "                            according to the their R2. PRSice is "
          "then performed\n"
          "                            again to find the best bin "
          "combination.\n"
          "                            This cannot be performed together with "
          "PRSet\n"
          "                            (Currently not implemented)\n"
          "\nTarget File:\n"
          "    --binary-target         Indicate whether the target phenotype\n"
          "                            is binary or not. Either T or F should "
          "be\n"
          "                            provided where T represent a binary "
          "phenotype.\n"
          "                            For multiple phenotypes, the input "
          "should be\n"
          "                            separated by comma without space. \n"
          "                            Default: T if --beta and F if --beta is "
          "not\n"
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
          "    --remove                File containing the sample(s) to be "
          "removed from\n"
          "                            the target file. First column should be "
          "FID and\n"
          "                            the second column should be IID. If "
          "--ignore-fid is\n"
          "                            set, first column should be IID\n"
          "                            Mutually exclusive from --keep\n"
          "    --pheno-file    | -f    Phenotype file containing the "
          "phenotype(s).\n"
          "                            First column must be FID of the samples "
          "and\n"
          "                            the second column must be IID of the "
          "samples.\n"
          "                            When --ignore-fid is set, first column "
          "must\n"
          "                            be the IID of the samples.\n"
          "                            Must contain a header if --pheno-col "
          "is\n"
          "                            specified\n"
          "    --pheno-col             Headers of phenotypes to be included "
          "from the\n"
          "                            phenotype file\n"
          "    --prevalence    | -k    Prevalence of all binary trait. If "
          "provided\n"
          "                            will adjust the ascertainment bias of "
          "the R2.\n"
          "                            Note that when multiple binary trait is "
          "found,\n"
          "                            prevalence information must be provided "
          "for\n"
          "                            all of them (Either adjust all binary "
          "traits,\n"
          "                            or don't adjust at all)\n"
          "    --nonfounders           Keep the nonfounders in the analysis\n"
          "                            Note: They will still be excluded from "
          "LD calculation\n"
          "    --target        | -t    Target genotype file. Currently "
          "support\n"
          "                            both BGEN and binary PLINK format. For "
          "\n"
          "                            multiple chromosome input, simply "
          "substitute\n"
          "                            the chromosome number with #. PRSice "
          "will\n"
          "                            automatically replace # with 1-22\n"
          "                            For binary plink format, you can also "
          "specify\n"
          "                            a seperate fam file by <prefix>,<fam "
          "file>\n"
          "    --type                  File type of the target file. Support "
          "bed \n"
          "                            (binary plink) and bgen format. "
          "Default: bed\n"
          "\nMisc:\n"
          "    --all                   Output PRS for ALL threshold. WARNING: "
          "This\n"
          "                            will generate a huge file\n"
          "    --exclude               File contains SNPs to be excluded from "
          "\n"
          "                            analysis\n"
          "    --extract               File contains SNPs to be included in "
          "the \n"
          "                            analysis\n"
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
          "    --out           | -o    Prefix for all file output\n"
          "    --perm                  Number of permutation to perform. This "
          "will\n"
          "                            generate the empirical p-value. "
          "Recommend to\n"
          "                            use value larger than 10,000\n"
          "    --seed          | -s    Seed used for permutation. If not "
          "provided,\n"
          "                            system time will be used as seed. When "
          "same\n"
          "                            seed and same input is provided, same "
          "result\n"
          "                            can be generated\n"
          "    --print-snp             Print all SNPs used to construct the "
          "best PRS\n"
          "    --thread        | -n    Number of thread use\n"
          "    --help          | -h    Display this help message\n";
}

void Commander::usage() { fprintf(stderr, "%s\n", help_message.c_str()); }


void Commander::base_check(std::map<std::string, std::string>& message,
                           bool& error, std::string& error_message)
{
    if (base.name.empty()) {
        error = true;
        error_message.append("ERROR: You must provide a base file\n");
    }
    else
    {
        // check the base file and get the corresponding index
        std::ifstream base_test;
        base_test.open(base.name.c_str());
        if (!base_test.is_open()) {
            error = true;
            error_message.append("ERROR: Cannot open base file to read!\n");
        }
        else
        {
            // check the base file header is correct
            std::string line;
            std::getline(base_test, line);
            base_test.close();
            std::vector<std::string> token = misc::split(line);
            int max_size = token.size();
            if (!base.index) {
                if (base.provided_stat) {
                    // if statistics is provided, we can guess if it
                    // is beta or not
                    if (base.statistic.length() == 2
                        && toupper(base.statistic[0]) == 'O'
                        && toupper(base.statistic[1]) == 'R')
                    {
                        base.beta = false;
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
                        base.beta = true;
                        message["beta"] = "";
                    }
                }
                else if (!base.provided_stat && base.beta)
                {
                    base.provided_stat = true;
                    base.statistic = "BETA";
                    message["stat"] = "BETA";
                }
                else if (!base.provided_stat)
                {
                    for (size_t i = 0; i < token.size(); ++i) {
                        if (token[i].length() == 2
                            && toupper(token[i][0]) == 'O'
                            && toupper(token[i][1] == 'R'))
                        {
                            base.provided_stat = true;
                            base.beta = false;
                            base.statistic = token[i];
                            message["stat"] = "OR";
                            break;
                        }
                        else if (token[i].length() == 4
                                 && toupper(token[i][0]) == 'B'
                                 && toupper(token[i][1]) == 'E'
                                 && toupper(token[i][2]) == 'T'
                                 && toupper(token[i][3]) == 'A')
                        {
                            base.provided_stat = true;
                            base.beta = true;
                            base.statistic = token[i];
                            // Again, this will be problematic if the BETA
                            // is actually OR... (Special crazy user)
                            message["stat"] = "BETA";
                            message["beta"] = "";
                            break;
                        }
                    }
                }
                base.col_index[+BASE_INDEX::CHR] = index_check(base.chr, token);
                if (!base.provided_chr
                    && base.col_index[+BASE_INDEX::CHR] != -1)
                    message["chr"] = base.chr;
                base.col_index[+BASE_INDEX::REF] =
                    index_check(base.ref_allele, token);
                if (!base.provided_ref
                    && base.col_index[+BASE_INDEX::REF] != -1)
                    message["A1"] =
                        base.ref_allele; // actually the alternative allele
                base.col_index[+BASE_INDEX::ALT] =
                    index_check(base.alt_allele, token);
                if (!base.provided_alt
                    && base.col_index[+BASE_INDEX::ALT] != -1)
                    message["A2"] = base.alt_allele; // the effective allele
                base.col_index[+BASE_INDEX::STAT] =
                    index_check(base.statistic, token);
                if (!base.provided_stat
                    && base.col_index[+BASE_INDEX::STAT] != -1)
                    message["stat"] = base.statistic;
                base.col_index[+BASE_INDEX::RS] = index_check(base.snp, token);
                if (!base.provided_snp && base.col_index[+BASE_INDEX::RS] != -1)
                    message["snp"] = base.snp;
                base.col_index[+BASE_INDEX::BP] = index_check(base.bp, token);
                if (!base.provided_bp && base.col_index[+BASE_INDEX::BP] != -1)
                    message["bp"] = base.bp;
                base.col_index[+BASE_INDEX::SE] =
                    index_check(base.standard_error, token);
                if (!base.provided_se && base.col_index[+BASE_INDEX::SE] != -1)
                    message["se"] = base.standard_error;
                base.col_index[+BASE_INDEX::P] =
                    index_check(base.p_value, token);
                if (!base.provided_p && base.col_index[+BASE_INDEX::P] != -1)
                    message["pvalue"] = base.p_value;


                std::vector<std::string> info = misc::split(base.info_col, ",");
                base.col_index[+BASE_INDEX::INFO] = index_check(info[0], token);
                if (info.size() != 2) {
                    error = true;
                    error_message.append(
                        "ERROR: Invalid format of --info-base.\n");
                    error_message.append(
                        "       Should be ColName,Threshold.\n");
                }
                try
                {
                    base.info_score = misc::convert<double>(info[1]);
                    if (base.info_score < 0 || base.info_score > 1) {
                        error = true;
                        error_message.append("ERROR: Base INFO threshold must "
                                             "be within 0 and 1!\n");
                        base.use_info = false;
                    }
                }
                catch (const std::runtime_error& er)
                {
                    error = true;
                    error_message.append(
                        "ERROR: Invalid argument passed to --info-base: "
                        + base.info_col + "!\n");
                    error_message.append(
                        "       Second argument must be numeric\n");
                }

                // Found INFO even if user didn't define INFO
                if (!base.use_info && base.col_index[+BASE_INDEX::INFO] != -1) {
                    // as default will always be of the correct format,
                    // we don't need to worry about the error messages above
                    message["info-base"] = base.info_col;
                }

                // comma separate
                if (base.provided_maf) {
                    // won't do it unless we have something provided
                    // as no default, will cause problem
                    std::vector<std::string> maf = misc::split(base.maf, ",");
                    if (maf.size() != 2) {
                        error = true;
                        error_message.append(
                            "ERROR: Invalid format of --maf-base.\n");
                        error_message.append(
                            "       Should be ColName,Threshold.\n");
                    }
                    base.col_index[+BASE_INDEX::MAF] =
                        index_check(maf[0], token);
                    try
                    {
                        base.maf_threshold = misc::convert<double>(maf[1]);
                        if (base.maf_threshold < 0 || base.maf_threshold > 1) {
                            error = true;
                            error_message.append(
                                "ERROR: Base MAF threshold must "
                                "be within 0 and 1!\n");
                        }
                        message["maf-base"] = base.maf;
                    }
                    catch (const std::runtime_error& er)
                    {
                        error = true;
                        error_message.append(
                            "ERROR: Invalid argument passed to --maf-base: "
                            + base.maf + "!\n");
                        error_message.append(
                            "       Second argument must be numeric\n");
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
                if (base.provided_ref) {
                    base.col_index[+BASE_INDEX::REF] = index_check(
                        base.ref_allele, max_size, error, error_message, "REF");
                }
                if (base.provided_alt) {
                    base.col_index[+BASE_INDEX::ALT] = index_check(
                        base.alt_allele, max_size, error, error_message, "ALT");
                }
                if (base.provided_bp) {
                    base.col_index[+BASE_INDEX::BP] = index_check(
                        base.bp, max_size, error, error_message, "BP");
                }
                if (base.provided_se) {
                    base.col_index[+BASE_INDEX::SE] =
                        index_check(base.standard_error, max_size, error,
                                    error_message, "SE");
                }
                if (base.use_info) {
                    std::vector<std::string> info =
                        misc::split(base.info_col, ",");
                    base.col_index[+BASE_INDEX::INFO] = index_check(
                        info[0], max_size, error, error_message, "INFO");
                    if (info.size() != 2) {
                        error = true;
                        error_message.append(
                            "ERROR: Invalid format of --info-base.\n");
                        error_message.append(
                            "       Should be ColName,Threshold.\n");
                    }
                    try
                    {
                        base.info_score = misc::convert<double>(info[1]);
                        if (base.info_score < 0 || base.info_score > 1) {
                            error = true;
                            error_message.append("ERROR: Base INFO threshold "
                                                 "must be within 0 and 1!\n");
                        }
                    }
                    catch (const std::runtime_error& er)
                    {
                        error = true;
                        error_message.append(
                            "ERROR: Invalid argument passed to --info-base: "
                            + base.info_col + "!\n");
                        error_message.append(
                            "       Second argument must be numeric\n");
                    }
                }
                if (base.provided_maf) {
                    std::vector<std::string> maf = misc::split(base.maf, ",");
                    base.col_index[+BASE_INDEX::MAF] = index_check(
                        maf[0], max_size, error, error_message, "MAF");
                    if (maf.size() != 2) {
                        error = true;
                        error_message.append(
                            "ERROR: Invalid format of --maf-base.\n");
                        error_message.append(
                            "       Should be ColName,Threshold.\n");
                    }
                    try
                    {
                        base.maf_threshold = misc::convert<double>(maf[1]);
                        if (base.maf_threshold < 0 || base.maf_threshold > 1) {
                            error = true;
                            error_message.append("ERROR: Base MAF threshold "
                                                 "must be within 0 and 1!\n");
                        }
                    }
                    catch (const std::runtime_error& er)
                    {
                        error = true;
                        error_message.append(
                            "ERROR: Invalid argument passed to --maf-base: "
                            + base.maf + "!\n");
                        error_message.append(
                            "       Second argument must be numeric\n");
                    }
                }
                base.col_index[+BASE_INDEX::P] = index_check(
                    base.p_value, max_size, error, error_message, "P");
                base.col_index[+BASE_INDEX::STAT] = index_check(
                    base.statistic, max_size, error, error_message, "STAT");
                base.col_index[+BASE_INDEX::RS] =
                    index_check(base.snp, max_size, error, error_message, "RS");
            }
            if (base.col_index[+BASE_INDEX::P] == -1) {
                error = true;
                error_message.append("ERROR: No p-value column (" + base.p_value
                                     + ") in file!\n");
            }
            else
                base.provided_p = true;
            if (base.col_index[+BASE_INDEX::STAT] == -1) {
                error = true;
                error_message.append("ERROR: No statistic column ("
                                     + base.statistic + ") in file!\n");
            }
            else
                base.provided_stat = true;
            if (base.col_index[+BASE_INDEX::RS] == -1) {
                error = true;
                error_message.append("ERROR: No SNP name column (" + base.snp
                                     + ") in file!\n");
            }
            else
                base.provided_snp = true;
            if (base.col_index[+BASE_INDEX::REF] == -1) {
                error = true;
                error_message.append("ERROR: No Reference allele column ("
                                     + base.ref_allele + ") in file!\n");
            }
            else
                base.provided_ref = true;

            double max_index =
                *max_element(base.col_index.begin(), base.col_index.end());
            base.col_index[+BASE_INDEX::MAX] = max_index;
        }
    }
}

void Commander::clump_check(std::map<std::string, std::string>& message,
                            bool& error, std::string& error_message)
{
    if (!clumping.no_clump) {
        if (clumping.keep_sample && clumping.remove_sample) {
            error = true;
            error_message.append(
                "ERROR: Can only use either --keep or --remove but not both\n");
        }
        // require clumping
        if (clumping.provide_proxy && clumping.proxy <= 0) {
            error = true;
            error_message.append(
                "ERROR: Proxy threshold cannot be negative!\n");
        }
        if (clumping.p_value < 0.0 || clumping.p_value > 1.0) {
            error = true;
            error_message.append(
                "ERROR: P-value threshold must be within 0 and 1!\n");
        }
        if (clumping.r2 < 0.0 || clumping.r2 > 1.0) {
            error = true;
            error_message.append(
                "ERROR: R2 threshold must be within 0 and 1!\n");
        }
        if (!clumping.type.empty()) {
            bool alright = false;
            for (auto&& type : supported_types) {
                if (clumping.type.compare(type) == 0) {
                    alright = true;
                    break;
                }
            }
            if (!alright) {
                error = true;
                error_message.append(
                    "ERROR: Unsupported LD format: " + clumping.type + "\n");
            }
        }
        if (!clumping.provide_r2)
            message["clump-r2"] = std::to_string(clumping.r2);
        if (!clumping.provide_p)

            message["clump-p"] = std::to_string(clumping.p_value);
        if (!clumping.provide_distance)
            message["clump-kb"] = std::to_string(clumping.distance);
        if (clumping.distance < 0.0) {
            error = true;
            error_message.append(
                "ERROR: Clumping distance must be positive!\n");
        }
    }
}


void Commander::covariate_check(bool& error, std::string& error_message)
{
    if (covariate.name.empty() || covariate.covariates.size() == 0) return;
    std::ifstream cov_file;
    cov_file.open(covariate.name.c_str());
    if (!cov_file.is_open()) {
        error = true;
        error_message.append(
            "ERROR: Cannot open covariate file: " + covariate.name + "\n");
        return;
    }
    std::string line;
    std::getline(cov_file, line);
    if (line.empty()) {
        error = true;
        error_message.append("ERROR: First line of covariate file is empty!\n");
        return;
    }
    cov_file.close();
    // obtain the header information
    std::unordered_set<std::string> included;
    for (auto cov : covariate.covariates) {
        if (cov.empty()) continue;
        if (included.find(cov)
            == included.end()) // to avoid duplicated covariance headers
        {
            // got annoyed with the input of PC.1 PC.2 PC.3, do this automatic
            // thingy to substitute them
            if (cov.at(0) == '@') {
                cov.erase(0, 1);
                std::vector<std::string> open = misc::split(cov, "[");
                std::vector<std::string> info;
                std::vector<bool> list;
                for (auto o : open) {
                    if (o.find("]") != std::string::npos) {
                        std::vector<std::string> close = misc::split(o, "]");
                        // the first one will always be the list
                        info.push_back(close[0]);
                        list.push_back(true);
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
                std::vector<std::string> final_covariates;
                for (size_t c = 0; c < info.size(); ++c) {
                    if (list[c]) {
                        std::vector<std::string> individual =
                            misc::split(info[c], ".");
                        std::vector<int> numeric;
                        for (auto&& ind : individual) {
                            if (ind.find("-") != std::string::npos) {
                                std::vector<std::string> range =
                                    misc::split(ind, "-");
                                if (range.size() != 2) {
                                    throw std::runtime_error(
                                        "ERROR: Invalid range format, range "
                                        "must be in the form of start-end");
                                }
                                try
                                {
                                    size_t start =
                                        misc::convert<size_t>(range[0]);
                                    size_t end =
                                        misc::convert<size_t>(range[1]);
                                    if (start > end) {
                                        int temp = end;
                                        end = start;
                                        start = temp;
                                    }
                                    for (size_t s = start; s <= end; ++s) {
                                        numeric.push_back(s);
                                    }
                                }
                                catch (const std::runtime_error& error)
                                {
                                    std::string error_message =
                                        "ERROR: Invalid parameter: " + range[0]
                                        + " or " + range[1]
                                        + ", only allow integer!";
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
                                        "ERROR: Invalid parameter: " + ind
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
                        for (size_t final = 0; final < final_covariates.size();
                             ++final)
                        {
                            final_covariates[final].append(info[c]);
                        }
                        if (final_covariates.empty())
                            final_covariates.push_back(info[c]);
                    }
                }
                for (auto res : final_covariates) {
                    if (included.find(res) == included.end()) {
                        included.insert(res);
                    }
                }
            }
            else
                included.insert(cov);
        }
    }

    std::vector<std::string> token = misc::split(line);
    std::string missing = "";
    std::unordered_set<std::string> ref;
    for (auto&& head : token) {
        ref.insert(head);
    }
    size_t valid_cov = 0;
    std::vector<std::string> final_cov;
    for (auto&& cov : included) {
        if (ref.find(cov) != ref.end()) {
            final_cov.push_back(cov);
            valid_cov++;
        }
        else if (missing.empty())
            missing = cov;
        else
            missing.append("," + cov);
    }
    if (!missing.empty())
        error_message.append(
            "WARNING: Covariate(s) missing from file: " + missing + "\n");
    if (valid_cov == 0) {
        error = true;
        error_message.append("ERROR: No valid Covariate!\n");
    }
    covariate.covariates = final_cov;
}


void Commander::filter_check(bool& error, std::string& error_message)
{
    if (filter.use_hard_thres
        && (filter.hard_threshold < 0 || filter.hard_threshold > 1))
    {
        error = true;
        error_message.append("ERROR: Negative number of permutation!\n");
    }
    if (filter.extract && filter.exclude) {
        error = true;
        error_message.append(
            "ERROR: Hard threshold should be between 0 and 1\n");
    }

    if (filter.info_filtering
        && (filter.info_score < 0 || filter.info_score > 1))
    {
        error = true;
        error_message.append("ERROR: INFO score cannot be bigger than 1.0 "
                             "or smaller than 0.0\n");
        break;
    }
}

void Commander::misc_check(std::map<std::string, std::string>& message,
                           bool& error, std::string& error_message)
{
    if (misc.provided_permutation && misc.permutation < 0) {
        error = true;
        error_message.append("ERROR: Negative number of permutation!\n");
    }
    if (misc.provided_seed && misc.seed < 0) {
        error = true;
        error_message.append("ERROR: Negative seed!\n");
    }
    else if (!misc.provided_seed)
    {
        misc.seed = std::random_device()();
        message["seed"] = std::to_string(misc.seed);
    }
    if (misc.thread <= 0) {
        error = true;
        error_message.append("ERROR: Number of thread must be larger than 1\n");
    }
    if (!misc.provided_permutation && misc.logit_perm) {
        error_message.append(
            "WARNING: Permutation not required, --logit-perm has no effect\n");
    }
    if (prsice.no_regress) misc.all = true;
    if (!misc.provide_thread) message["thread"] = "1";
    if (!misc.provided_output) message["out"] = misc.out;
}

void Commander::prset_check(std::map<std::string, std::string>& message,
                            bool& error, std::string& error_message)
{
    if (!prset.perform_prset) return;
    if (!prset.gtf.empty() && prset.msigdb.empty()) {
        error = true;
        error_message.append(
            "ERROR: Must provide a gtf file if msigdb is specified\n");
    }
    if (prset.feature.empty()) {
        prset.feature.push_back("exon");
        prset.feature.push_back("gene");
        prset.feature.push_back("protein_coding");
        prset.feature.push_back("CDS");
        message["feature"] = "exon,gene,protein_coding,CDS";
    }
}


void Commander::prsice_check(std::map<std::string, std::string>& message,
                             bool& error, std::string& error_message)
{
    if (!prsice.provided_model) {
        message["model"] = "add";
    }
    if (prsice.fastscore && prsice.barlevel.size() == 0 && !prset.perform_prset)
    {
        message["bar-levels"] = "0.001,0.05,0.1,0.2,0.3,0.4,0.5";
        prsice.barlevel = {0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
    }
    std::sort(prsice.barlevel.begin(), prsice.barlevel.end());
    prsice.barlevel.erase(
        std::unique(prsice.barlevel.begin(), prsice.barlevel.end()),
        prsice.barlevel.end());
    if (prset.perform_prset) {
        if (!prsice.provide_inter && !prsice.provide_upper
            && !prsice.provide_lower && !prsice.fastscore)
        {
            message["bar-levels"] = 1;
            prsice.fastscore = true;
            prsice.barlevel = {1};
        }
    }
    else if (!prsice.fastscore)
    {
        if (prsice.no_regress) {
            error = true;
            error_message.append(
                "ERROR: no-regress can only be used with fastscore!\n");
        }
        if (prsice.inter <= 0) {
            error = true;
            error_message.append("ERROR: Cannot have negative interval!\n");
        }
        if (prsice.upper < prsice.lower) {
            error = true;
            error_message.append(
                "ERROR: Upper bound must be larger than lower bound!\n");
        }
        if (prsice.upper < 0.0 || prsice.lower < 0.0) {
            error = true;
            error_message.append("ERROR: CAnnot have negative bounds!\n");
        }
        if (!prsice.provide_inter)
            message["interval"] = std::to_string(prsice.inter);
        if (!prsice.provide_lower)
            message["lower"] = std::to_string(prsice.lower);
        if (!prsice.provide_upper)
            message["upper"] = std::to_string(prsice.upper);
    }
}

void Commander::prslice_check(bool& error, std::string& error_message)
{
    if (prslice.provided) {
        if (misc.all) {
            error = true;
            error_message.append("ERROR: Cannot output PRS for all threshold "
                                 "when using PRSlice!\n");
        }
        if (prslice.size <= 0) {
            error = true;
            error_message.append(
                "ERROR: PRSlice size cannot be less than 1!\n");
        }
    }
}

void Commander::target_check(std::map<std::string, std::string>& message,
                             bool& error, std::string& error_message)
{
    if (target.name.empty()) {
        error = true;
        error_message.append("ERROR: You must provide a target file!\n");
    }
    if (target.keep_sample && target.remove_sample) {
        error = true;
        error_message.append(
            "ERROR: Can only use either --keep or --remove but not both\n");
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
        error_message.append("ERROR: Unsupported target format: " + target.type
                             + "\n");
    }
    if (target.type.compare("bgen") == 0 && filter.hard_coding
        && !filter.use_hard_thres)
    {
        message["hard-thres"] = std::to_string(filter.hard_threshold);
    }
    if (target.pheno_col.size() != 0 && target.pheno_file.empty()) {
        error = true;
        error_message.append("ERROR: You must provide a phenotype file for "
                             "multiple phenotype analysis");
    }
    if (target.pheno_file.empty() && target.is_binary.empty()) {
        message["binary-target"] = "T";
        target.is_binary.push_back(true);
    }
    else
    {
        if (target.pheno_col.empty() && target.is_binary.size() == 1) {
            // this is ok
        }
        else if (target.pheno_col.empty() && target.is_binary.empty())
        {
            if (base.beta) {
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
            if (base.beta) {
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
            error_message.append("ERROR: Number of target phenotypes doesn't "
                                 "match information of binary\n");
            error_message.append("       target! You must indicate whether the "
                                 "phenotype is binary using\n");
            error_message.append("       --binary-target\n");
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
        error_message.append("ERROR: Number of target prevalence doesn't match "
                             "number of binary traits\n");
        error_message.append("       You must provide a prevalence for all "
                             "binary trait(s) or not \n");
        error_message.append(
            "       provide any prevalence (all or nothing)\n");
    }
    for (auto&& prev : target.prevalence) {
        if (prev > 1.0 || prev < 0.0) {
            error = true;
            error_message.append("ERROR: Prevalence cannot be bigger than 1.0 "
                                 "or smaller than 0.0\n");
            break;
        }
    }
}

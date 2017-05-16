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

bool Commander::initialize(int argc, char *argv[])
{
	program_info();
	info();
    if(argc<=1)
    {
        usage();
        throw std::runtime_error("Please provide the required parameters");
    }
    static const char *optString = "a:b:B:c:C:f:g:i:l:L:m:n:o:p:t:u:h?";
    static const struct option longOpts[]=
    {
        // parameters with short flags
        {"ancestry",required_argument,NULL,'a'},
        {"base",required_argument,NULL,'b'},
        {"bed",required_argument,NULL,'B'},
        {"cov-header",required_argument,NULL,'c'},
        {"cov-file",required_argument,NULL,'C'},
        {"pheno-file",required_argument,NULL,'f'},
        {"gtf",required_argument,NULL,'g'},
        {"interval",required_argument,NULL,'i'},
        {"prevalence",required_argument,NULL,'k'},
        {"lower",required_argument,NULL,'l'},
        {"ld",required_argument,NULL,'L'},
        {"msigdb",required_argument,NULL,'m'},
        {"thread",required_argument,NULL,'n'},
        {"out",required_argument,NULL,'o'},
        {"pvalue",required_argument,NULL,'p'},
        {"seed",required_argument,NULL,'s'},
        {"target",required_argument,NULL,'t'},
        {"upper",required_argument,NULL,'u'},
		// flags, only need to set them to true
        {"all",no_argument,&misc.all,1},
        {"beta",required_argument,&base.beta,1},
        {"full",no_argument,&prsice.full,1},
        {"ignore-fid",no_argument,&misc.ignore_fid,1},
        {"index",no_argument,&base.index,1},
		{"no-clump",no_argument,&clumping.no_clump,1},
		{"no-regression",no_argument,&prsice.no_regress,1},
		{"no-x",no_argument,&species.no_x,1},
		{"no-y",no_argument,&species.no_y,1},
		{"no-xy",no_argument,&species.no_xy,1},
		{"no-mt",no_argument,&species.no_mt,1},
		{"fastscore",no_argument,&prsice.fastscore,1},
        {"print-snp",no_argument,&misc.print_snp,1},
		// long flags, need to work on them
        {"A1",required_argument,NULL,0},
        {"A2",required_argument,NULL,0},
        {"bar-levels",required_argument,NULL,0},
        {"binary-target",required_argument,NULL,0},
        {"bp",required_argument,NULL,0},
        {"chr",required_argument,NULL,0},
        {"clump-kb",required_argument,NULL,0},
        {"clump-p",required_argument,NULL,0},
        {"clump-r2",required_argument,NULL,0},
        {"feature",required_argument,NULL,0},
        {"info",required_argument,NULL,0},
        {"keep",required_argument,NULL,0},
		{"ld-type", required_argument, NULL, 0},
		{"num-auto", required_argument, NULL, 0},
        {"perm",required_argument,NULL,0},
        {"pheno-col",required_argument,NULL,0},
        {"proxy",required_argument,NULL,0},
        {"prslice",required_argument,NULL,0},
        {"remove",required_argument,NULL,0},
		{"score", required_argument, NULL, 0},
        {"se",required_argument,NULL,0},
        {"snp",required_argument,NULL,0},
        {"stat",required_argument,NULL,0},
		{"type", required_argument, NULL, 0},
		//species flag
        {"cow",no_argument,NULL,0},
        {"dog",no_argument,NULL,0},
        {"horse",no_argument,NULL,0},
        {"mouse",no_argument,NULL,0},
        {"rice",no_argument,NULL,0},
        {"sheep",no_argument,NULL,0},
        {"help",no_argument,NULL,'h'},
        {NULL, 0, 0, 0}
    };
    bool error = false;
    bool species_error = false;
    int longIndex=0;
    int opt = 0;

    std::string command ="";
    opt=getopt_long(argc, argv, optString, longOpts, &longIndex);
    std::string error_message = "";
    //Start reading all the parameters and perform the qc at the same time
    while(opt!=-1)
    {
        switch(opt)
        {
        case 0:
            command = longOpts[longIndex].name;
            if (longOpts[longIndex].flag != 0) break;
            else if(command.compare("chr")==0)  set_chr(optarg);
            else if(command.compare("A1")==0) set_ref(optarg);
            else if(command.compare("A2")==0) set_alt(optarg);
            else if(command.compare("stat")==0) set_stat(optarg);
            else if(command.compare("snp")==0) set_snp(optarg);
            else if(command.compare("bp")==0) set_bp(optarg);
            else if(command.compare("se")==0) set_se(optarg);
            else if(command.compare("keep")==0) set_keep(optarg);
            else if(command.compare("remove")==0) set_remove(optarg);
            else if(command.compare("ld-type")==0) clumping.type = optarg;
            else if(command.compare("type")==0) target.type = optarg;
            else if(command.compare("score")==0) prsice.missing_score = optarg;
            else if(command.compare("cow")==0) set_species(29, false, false, true, false,
            		error, error_message, species_error);
            else if(command.compare("dog")==0) set_species(38, false, false, false, false,
            		error, error_message, species_error);
            else if(command.compare("horse")==0) set_species(31, false, false, true, true,
            		error, error_message, species_error);
            else if(command.compare("mouse")==0) set_species(19, false, false, true, true,
            		error, error_message, species_error);
            else if(command.compare("rice")==0) set_species(-12, false, false, false, false,
            		error, error_message, species_error);
            else if(command.compare("sheep")==0) set_species(26, false, false, true, true,
            		error, error_message, species_error);
            else if(command.compare("num-auto")==0)
            {
                try{
                    species.num_auto = misc::convert<int>(optarg);
                }
                catch(const std::runtime_error &er)
                {
                    error_message.append("ERROR: Number of autosomal chromosome provided isn't numeric\n");
                    error=true;
                }
            }
            else if(command.compare("info")==0)
            {
                try{
                    filter.info_score = misc::convert<double>(optarg);
                    filter.use_info = true;
                }
                catch(const std::runtime_error &er)
                {
                    error_message.append("ERROR: Info score isn't numeric!\n");
                    error = true;
                }
            }
            else if(command.compare("clump-p")==0)
            {
                try{
                    clumping.p_value = misc::convert<double>(optarg);
                }
                catch(const std::runtime_error &er)
                {
                    error_message.append("ERROR: Clump P-value threshold provided isn't numeric\n");
                    error=true;
                }
            }
            else if(command.compare("clump-r2")==0)
            {
                try{
                    clumping.r2 = misc::convert<double>(optarg);
                }
                catch(const std::runtime_error &er)
                {
                    error_message.append("ERROR: Clump R2 threshold provided isn't numeric\n");
                    error=true;
                }
            }
            else if(command.compare("clump-kb")==0)
            {
                try{
                    clumping.distance = misc::convert<int>(optarg)*1000; //kb
                }
                catch(const std::runtime_error &er)
                {
                    error_message.append("ERROR: Clump Distance provided isn't numeric\n");
                    error=true;
                }
            }
            else if(command.compare("prslice")==0)
            {
                try{
                    set_prslice(misc::convert<int>(optarg));
                }
                catch(const std::runtime_error &er)
                {
                    error_message.append("ERROR: PRSlice size provided isn't numeric\n");
                    error=true;
                }
            }
            else if(command.compare("binary-target")==0)
            {
            		std::vector<std::string> token = misc::split(optarg, ",");
            	    for(auto &&bin : token) target.is_binary.push_back(misc::to_bool(bin));
            }
            else if(command.compare("bar-levels")==0)
            {
                std::vector<std::string> token = misc::split(optarg, ",");
                try
                {
                    for(auto &&bar : token) prsice.barlevel.push_back(misc::convert<double>(bar));
                }
                catch(const std::runtime_error &er)
                {
                    error_message.append("ERROR: None numeric barchart level\n");
                    error=true;
                }
            }
            else if(command.compare("proxy")==0)
            {
                try
                {
                    clumping.proxy = misc::convert<double>(optarg);
                    clumping.provide_proxy = true;
                }
                catch(const std::runtime_error &er)
                {
                    error_message.append("ERROR: Proxy provided isn't numeric\n");
                    error=true;
                }
            }
            else if(command.compare("pheno-col")==0)
            {
                std::vector<std::string> token = misc::split(optarg, ",");
                target.pheno_col.insert(target.pheno_col.end(), token.begin(), token.end());
            }
            else if(command.compare("perm")==0)
            {
                try{
                    set_permutation(optarg);
                }
                catch(const std::runtime_error &er)
                {
                    error_message.append("ERROR: Number of Permutation provided isn't numeric\n");
                    error=true;
                }
            }
            else if(command.compare("feature")==0)
            {
                std::vector<std::string> token = misc::split(optarg, ",");
                prset.feature.insert(prset.feature.end(), token.begin(), token.end());
            }
            else
            {
                std::string er = "Undefined operator: "+command+", please use --help for more information!";
                throw std::runtime_error(er);
            }
        break;

        case 'a':
            covariate.ancestry_dim = optarg;
            fprintf(stderr, "Currently we have not implement this function\n");
        break;
        case 'b':
            base.name = optarg;
        break;
        case 'B':
        {
            std::vector<std::string> token = misc::split(optarg, ",");
            prset.bed.insert(prset.bed.end(), token.begin(), token.end());
        }
        break;
        case 'c':
        {
            std::vector<std::string> token= misc::split(optarg, ",");
            covariate.covariates.insert(covariate.covariates.end(), token.begin(), token.end());
        }
        break;
        case 'C':
        	    covariate.name = optarg;
        	    break;
        case 'f':
            target.pheno_file = optarg;
        	break;
        case 'g':
            prset.gtf = optarg;
        	break;
        case 'i':
            try{
                prsice.inter = misc::convert<double>(optarg);
            }
            catch(const std::runtime_error &er)
            {
                error_message.append("ERROR: Interval provided isn't numeric\n");
                error=true;
            }
            break;
        case 'k':
            try{
                std::vector<std::string> token = misc::split(optarg, ",");
                for(auto &&prev : token) target.prevalence.push_back(misc::convert<double>(prev));
            }
            catch(const std::runtime_error &er)
            {
                error_message.append("ERROR: Prevalence provided isn't numeric\n");
                error=true;
            }
            break;
        case 'l':
            try{
                prsice.lower = misc::convert<double>(optarg);
            }
            catch(const std::runtime_error &er)
            {
                error_message.append("ERROR: Lower bound provided isn't numeric\n");
                error=true;
            }
        	break;
        case 'L':
            clumping.ld = optarg;
        	break;
        case 'm':
            prset.msigdb = optarg;
        	break;
        case 'n':
            try{
                misc.thread = misc::convert<int>(optarg);
            }
            catch(const std::runtime_error &er)
            {
                error_message.append("ERROR: Number of thread provided isn't numeric\n");
                error=true;
            }
        	break;
        case 'o':
            misc.out = optarg;
        	break;
        case 'p':
            set_p(optarg);
            break;
        case 's':
            set_seed(optarg);
            break;
        case 't':
            target.name = optarg;
        	break;
        case 'u':
            try{
                prsice.upper = misc::convert<double>(optarg);
            }
            catch(const std::runtime_error &er)
            {
                error_message.append("ERROR: Upper bound provided isn't numeric\n");
                error=true;
            }
        	break;
        case 'h':
        case '?':
            usage();
            return false;
            break;
        default:
            throw "Undefined operator, please use --help for more information!";
        }
        opt=getopt_long(argc, argv, optString, longOpts, &longIndex);
    }

    base_check(error, error_message);
    clump_check(error, error_message);
    misc_check(error, error_message);
    prset_check(error, error_message);
    prsice_check(error, error_message);
    prslice_check(error, error_message);
    target_check(error, error_message);
    if((prset.bed.size() != 0 || !prset.gtf.empty()) && prslice.provided)
    {
        error = true;
        error_message.append("ERROR: PRSet and PRSlice cannot be performed together!\n");
    }
    if(error) throw std::runtime_error(error_message);
    fprintf(stderr, "\n");
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
	base.index = false;
    base.provided_chr = false;
    base.provided_ref = false;
    base.provided_alt = false;
    base.provided_stat = false;
    base.provided_snp = false;
    base.provided_bp = false;
    base.provided_se = false;
    base.provided_p = false;
    base.col_index.resize(+BASE_INDEX::MAX+1, -1);


    clumping.distance = 250000;
    clumping.keep_sample = false;
	clumping.ld = "";
	clumping.no_clump = false;
	clumping.provide_proxy = false;
	clumping.proxy =-1.0;
	clumping.p_value =1.0;
	clumping.r2 = 0.1;
    clumping.remove_sample = false;
    clumping.type = "bed";

	covariate.name = "";
	covariate.ancestry_dim="MDS";

	filter.geno = 0.0;
	filter.mind = 0.0;
	filter.info_score = 0.8;
	filter.maf = 0.01;
	filter.use_maf = false;
	filter.use_mind = false;
	filter.use_info = false;
	filter.use_geno = false;

	misc.all = false;
	misc.out = "PRSice";
	misc.print_snp = false;
	misc.ignore_fid = false;
	misc.permutation = 0;
	misc.provided_permutation = false;
	misc.provided_seed =false;
	misc.seed = 0;
	misc.thread = 1;

	prset.gtf = "";
	prset.msigdb = "";

	prsice.missing_score = "";
	prsice.lower = 0.0001;
	prsice.upper = 0.5;
	prsice.inter = 0.00005;
	prsice.fastscore = false;
	prsice.no_regress = false;
	prsice.full = true;

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
	target.name = "";
	target.pheno_file = "";
	target.type = "bed";
}

Commander::~Commander()
{
    //dtor
}

void Commander::info()
{
    help_message =
            "usage: PRSice [options] <-b base_file> <-t target_file>\n"
            "\nBase File:\n"
            "    --base          | -b    Base association file\n"
            "    --beta                  Whether the test statistic is in the form of \n"
            "                            BETA or OR. If set, test statistic is assume\n"
            "                            to be in the form of BETA.\n"
            "    --A1                    Column header containing the reference allele\n"
            "                            Default: A1\n"
            "    --A2                    Column header containing the alternative allele\n"
            "                            Default: A2\n"
            "    --bp                    Column header containing the SNP coordinate\n"
            "                            Default: BP\n"
            "    --chr                   Column header containing the chromosome\n"
            "                            Default: CHR\n"
            "    --index                 If set, assume the INDEX instead of NAME of\n"
            "                            the corresponding columns are provided. Index\n"
            "                            should be 0-based (start counting from 0)\n"
            "    --pvalue        | -p    Column header containing the p-value\n"
            "                            Default: P\n"
            "    --se                    Column header containing the standard error\n"
            "                            Default: SE\n"
            "    --snp                   Column header containing the SNP ID\n"
            "                            Default: SNP\n"
            "    --stat                  Column header containing the summary statistic\n"
            "                            If --beta is set, default as BETA. Otherwise,\n"
            "                            try and search for OR or BETA from the header\n"
            "                            of the base file\n"
            "\nClumping:\n"
            "    --clump-kb              The distance for clumping in kb\n"
            "                            Default: "+std::to_string(clumping.distance/1000)+"\n"
            "    --clump-r2              The R2 threshold for clumping\n"
            "                            Default: "+std::to_string(clumping.r2)+"\n"
            "    --clump-p               The p-value threshold use for clumping.\n"
            "                            Default: "+std::to_string(clumping.p_value)+"\n"
            "    --ld            | -L    LD reference file. Use for LD calculation. If not\n"
            "                            provided, will use the post-filtered target genotype\n"
            "                            for LD calculation. Support multiple chromosome input\n"
            "                            Please see --target for more information\n"
            "    --ld-keep               File containing the sample(s) to be extracted from\n"
            "                            the LD reference file. First column should be FID and\n"
            "                            the second column should be IID. If --ignore-fid is\n"
            "                            set, first column should be IID\n"
            "                            Mutually exclusive from --ld-remove\n"
            "    --ld-remove             File containing the sample(s) to be removed from\n"
            "                            the LD reference file. First column should be FID and\n"
            "                            the second column should be IID. If --ignore-fid is\n"
            "                            set, first column should be IID\n"
            "                            Mutually exclusive from --ld-keep\n"
            "    --ld-type               File type of the LD file. Support bed (binary plink)\n"
            "                            and bgen format. Default: bed\n"
            "    --no-clump              Avoid performing clumping\n"
            "    --proxy                 Proxy threshold for index SNP to be considered\n"
            "                            as part of the region represented by the clumped\n"
            "                            SNP(s). e.g. --proxy 0.8 means the index SNP will\n"
            "                            represent region of any clumped SNP(s) that has a\n"
            "                            R2>=0.8 even if the index SNP does not physically\n"
            "                            locate within the region\n"
            "\nCovariate:\n"
            "    --cov-file      | -C    Covariate file. First column should be FID and \n"
            "                            the second column should be IID. If --ignore-fid\n"
            "                            is set, first column should be IID\n"
            "    --cov-header    | -c    Header of covariates. If not provided, will use\n"
            "                            all variables in the covariate file\n"
            "\nDosage:\n"
            "    --info                  Info threshold for dosage data. Any call less than\n"
            "                            this will be treated as missing. Default: "+std::to_string(filter.info_score)+
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
            "\nPRSice:\n"
            "    --bar-levels            Level of barchart to be plotted. When --fastscore\n"
            "                            is set, PRSice will only calculate the PRS for \n"
            "                            threshold within the bar level. Levels should be\n"
            "                            comma separated without space\n"
            "    --fastscore             Only calculate threshold stated in --bar-levels\n"
            "    --full                  Include the full model in the analysis\n"
            "    --interval      | -i    The step size of the threshold. Default: "+std::to_string(prsice.inter)+"\n"
            "    --lower         | -l    The starting p-value threshold. Default: "+std::to_string(prsice.lower)+"\n"
            "    --no-regress            Do not perform the regression analysis and simply\n"
            "                            output all PRS.\n"
            "    --score                 Method to handle missing genotypes. By default, \n"
            "                            final scores are averages of valid per-allele \n"
            "                            scores with missing genotypes contribute an amount\n"
            "                            proportional to imputed allele frequency. To throw\n"
            "                            out missing observations instead (decreasing the\n"
            "                            denominator in the final average when this happens),\n"
            "                            use the 'no_mean_imputation' modifier. Alternatively,\n"
            "                            you can use the 'center' modifier to shift all scores\n"
            "                            to mean zero. \n"
            "    --upper         | -u    The final p-value threshold. Default: "+std::to_string(prsice.upper)+"\n"
            "\nPRSlice:\n"
            "    --prslice               Perform PRSlice where the whole genome is first cut\n"
            "                            into bin size specified by this option. PRSice will\n"
            "                            then be performed on each bin. Bins are then sorted\n"
            "                            according to the their R2. PRSice is then performed\n"
            "                            again to find the best bin combination.\n"
            "                            This cannot be performed together with PRSet"
            "\nTarget File:\n"
            "    --binary-target         Indicate whether the target phenotype\n"
            "                            is binary or not. Either T or F should be\n"
            "                            provided where T represent a binary phenotype.\n"
            "                            For multiple phenotypes, the input should be\n"
            "                            separated by comma without space. Default: T\n"
            "    --keep                  File containing the sample(s) to be extracted from\n"
            "                            the target file. First column should be FID and\n"
            "                            the second column should be IID. If --ignore-fid is\n"
            "                            set, first column should be IID\n"
            "                            Mutually exclusive from --remove\n"
            "    --pheno-file    | -f    Phenotype file containing the phenotype(s).\n"
            "                            First column must be FID of the samples and\n"
            "                            the second column must be IID of the samples.\n"
            "                            When --ignore-fid is set, first column must\n"
            "                            be the IID of the samples.\n"
            "                            Must contain a header if --pheno-col is\n"
            "                            specified\n"
            "    --pheno-col             Headers of phenotypes to be included from the\n"
            "                            phenotype file\n"
            "    --prevalence    | -k    Prevalence of all binary trait. If provided\n"
            "    --remove                will adjust the ascertainment bias of the R2.\n"
            "                            Note that when multiple binary trait is found,\n"
            "                            you must provide prevalence information for\n"
            "                            all of them.\n"
            "    --target        | -t    Target genotype file. Currently support\n"
            "                            both BGEN and binary PLINK format. For \n"
            "                            multiple chromosome input, simply substitute\n"
            "                            the chromosome number with #. PRSice will\n"
            "                            automatically replace # with 1-22\n"
            "    --type                  File type of the target file. Support bed \n"
            "                            (binary plink) and bgen format. Default: bed\n"
            "\nMisc:\n"
            "    --all                   Output PRS for ALL threshold. WARNING: This\n"
            "                            will generate a huge file\n"
            "    --ignore-fid            Ignore FID for all input. When this is set,\n"
            "                            first column of most file will be assume to\n"
            "                            be IID instead of FID\n"
            "    --out           | -o    Prefix for all file output\n"
            "    --perm                  Number of permutation to perform. This will\n"
            "                            generate the empirical p-value for the BEST\n"
            "                            threshold\n"
            "    --seed          | -s    Seed used for permutation. If not provided,\n"
            "    --print-snp             system time will be used as seed. When same\n"
            "                            seed and same input is provided, same result\n"
            "                            should be generated\n"
            "    --thread        | -n    Number of thread use\n"
            "    --help          | -h    Display this help message\n";

}

void Commander::usage()
{
    fprintf(stderr, "%s\n", help_message.c_str());
	/*
    fprintf(stderr, "Usage: PRSice [Options] \n\n");
    std::string category = "";
    size_t width = 80;
    size_t max_prev = 0;
    for(auto message : m_help_messages)
    {
    		max_prev=(max_prev > std::get<help_index::LONG>(message).length()) ? max_prev: std::get<help_index::LONG>(message).length();
    }
    max_prev+=13;
    int remain = width-max_prev;
    std::string pad="";
    for(size_t i = 0; i < max_prev; ++i) pad+=" ";
    for(auto message : m_help_messages)
    {
    		if(category.empty() || std::get<help_index::CATEGORY>(message)!=category)
    		{
    			fprintf(stderr, "\n%s:\n",  std::get<help_index::CATEGORY>(message).c_str());
    			category =  std::get<help_index::CATEGORY>(message);
    		}
    		char short_flag = std::get<help_index::SHORT>(message);
    		std::string cur_message = "    ";
    		if(short_flag=='\0') cur_message+="     --"+std::get<help_index::LONG>(message);
    		else cur_message+="-"+std::string(1,short_flag)+" | --"+std::get<help_index::LONG>(message);
    		int diff = max_prev-cur_message.length();
    		if(diff > 0){
    			for(size_t i = 0; i < diff; ++i) cur_message.append(" ");
    		}

    		std::string description = std::get<help_index::DESCRIPTION>(message);
    		std::cerr << cur_message;
    		size_t cur_index = 0;
    		size_t last_space = 0;
    		for(size_t i = 0; i < description.length(); ++i)
    		{
    			if(description[i]==' ') last_space=i;
    			if(i-cur_index>=remain)
    			{
    				if(i != last_space)
    				{
    					if(cur_index != 0) std::cerr << pad;
    					std::cerr << description.substr(cur_index, last_space-cur_index) << std::endl;
    					i=last_space+1;
    					cur_index=i;
    				}
    			}
    			else if(i==description.length()-1)
    			{
				if(cur_index != 0) std::cerr << pad;
    				std::cerr << description.substr(cur_index, i-cur_index+1) << std::endl;
    			}
    		}
    }
    */
}

void Commander::program_info()
{
	std::cerr << std::endl;
	std::cerr << "PRSice "<< version << " (" << date << ") https://github.com/choishingwan/PRSice"<< std::endl;
	std::cerr << "(C) 2016-2017 Jack Euesden, Cathryn M. Lewis, Paul F. O'Reilly, Sam Choi" << std::endl;
	std::cerr << "GNU General Public License v3" << std::endl << std::endl;
}

void Commander::user_input() const{
	fprintf(stderr, "Base file: ");
	if (!base.beta) {
		fprintf(stderr, "%s(OR) ", base.name.c_str());
	} else {
		fprintf(stderr, "%s(Beta) ", base.name.c_str());
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "\nUser Defined Column Headers\n");
	fprintf(stderr, "==============================\n");
	if (base.col_index[+BASE_INDEX::CHR]!=-1)
		fprintf(stderr, "Chr            : %s\n", base.chr.c_str());
	fprintf(stderr, "SNP            : %s\n",base.snp.c_str());
	if (base.col_index[+BASE_INDEX::BP]!=-1)
		fprintf(stderr, "BP             : %s\n", base.bp.c_str());
	fprintf(stderr, "Ref Allele     : %s\n", base.ref_allele.c_str());
	if (base.col_index[+BASE_INDEX::ALT]!=-1)
		fprintf(stderr, "Alt Allele     : %s\n", base.alt_allele.c_str());
	fprintf(stderr, "Statistic      : %s\n", base.statistic.c_str());
	if (base.col_index[+BASE_INDEX::SE]!=-1)
		fprintf(stderr, "Standard Error : %s\n", base.standard_error.c_str());
	fprintf(stderr, "P-value        : %s\n", base.p_value.c_str());
	fprintf(stderr, "\nClumping Parameters: \n");
	fprintf(stderr, "==============================\n");
	fprintf(stderr, "P-Threshold  : %f\n", clumping.p_value);
	fprintf(stderr, "R2-Threshold : %f\n", clumping.r2);
	fprintf(stderr, "Window Size  : %d\n", clumping.distance);
	if(!prsice.no_regress)
	{
		fprintf(stderr, "\nThreshold Selected: \n");
		fprintf(stderr, "==============================\n");
		if(prsice.fastscore)
		{
			fprintf(stderr, "Fastscore    : TRUE\n");
			fprintf(stderr, "Barlevels    : ");
			if(prsice.barlevel.empty())
			{
				throw std::runtime_error("Cannot perform fastscore without bar level information!");
			}
			fprintf(stderr, "%f", prsice.barlevel.front());
			for(size_t i_bar = 1; i_bar < prsice.barlevel.size(); ++i_bar)
			{
				fprintf(stderr, ", %f", prsice.barlevel[i_bar]);
			}
			fprintf(stderr, "\n");
		}
		else
		{
			fprintf(stderr, "Fastscore    : FALSE\n");
			fprintf(stderr, "Lower Bound  : %f\n", prsice.lower);
			fprintf(stderr, "Upper Bound  : %f\n", prsice.upper);
			fprintf(stderr, "Intervals    : %f\n", prsice.inter);
		}
	}
}



void Commander::base_check(bool &error, std::string &error_message)
{

	if(base.name.empty())
	{
		error = true;
		error_message.append("ERROR: You must provide a base file\n");
	}
	else
	{

		// check the base file and get the corresponding index
		std::ifstream base_test;
		base_test.open(base.name.c_str());
		if(!base_test.is_open())
		{
			error =true;
			error_message.append("ERROR: Cannot open base file to read!\n");
		}
		else
		{
		    std::string line;
			std::getline(base_test, line);
			base_test.close();
			std::vector<std::string> token = misc::split(line);
			int max_size = token.size();
			if(!base.index)
			{
			    if(base.provided_stat){
	                if(base.statistic.length() ==2 && toupper(base.statistic[0])=='O'
	                        && toupper(base.statistic[1])=='R')
	                {
                        base.beta = false;
                        fprintf(stderr, "Base assumed be OR\n");
	                }
	                else if(base.statistic.length()==4 && toupper(base.statistic[0])=='B'
	                        && toupper(base.statistic[1])=='E' && toupper(base.statistic[2])=='T'
	                                && toupper(base.statistic[3])=='A')
	                {
                        base.beta = true;
                        fprintf(stderr, "Base assumed be BETA\n");
	                }
	            }
	            else if(!base.provided_stat && base.beta)
	            {
	                base.provided_stat = true;
	                base.statistic = "BETA";
	                fprintf(stderr, "Base statistic not provided, assumed to be %s\n", base.statistic.c_str());
	            }
	            else if(!base.provided_stat)
	            {
	                for(size_t i = 0; i < token.size(); ++i)
	                {
	                    if(token[i].length()==2 && toupper(token[i][0])=='O'
	                        && toupper(token[i][1]=='R'))
	                    {
	                        base.provided_stat = true;
	                        base.beta = false;
	                        base.statistic = token[i];
	                        fprintf(stderr, "Base statistic guessed to be %s (%s)\n",
	                                token[i].c_str(),"OR");
	                        break;
	                    }
	                    else if(token[i].length()==4 && toupper(token[i][0])=='B'
                            && toupper(token[i][1])=='E' && toupper(token[i][2])=='T'
                                    && toupper(token[i][3])=='A')
	                    {
                            base.provided_stat = true;
                            base.beta = true;
                            base.statistic = token[i];
                            fprintf(stderr, "Base statistic guessed to be %s (%s)\n",
                                    token[i].c_str(),"BETA");
                            break;
	                    }
	                }
	            }

			    base.col_index[+BASE_INDEX::CHR] = index_check(base.chr, token);
				base.col_index[+BASE_INDEX::REF] = index_check(base.ref_allele, token);
				base.col_index[+BASE_INDEX::ALT] = index_check(base.alt_allele, token);
				base.col_index[+BASE_INDEX::STAT] = index_check(base.statistic, token);
				base.col_index[+BASE_INDEX::RS] = index_check(base.snp, token);
				base.col_index[+BASE_INDEX::BP] = index_check(base.bp, token);
				base.col_index[+BASE_INDEX::SE] = index_check(base.standard_error, token);
				base.col_index[+BASE_INDEX::P] = index_check(base.p_value, token);
			}
			else
			{ // only required for index, as the defaults are in string
				if(base.provided_chr)
				{
					base.col_index[+BASE_INDEX::CHR] = index_check(base.chr, max_size, error, error_message, "CHR");
				}
				if(base.provided_ref)
				{
					base.col_index[+BASE_INDEX::REF] = index_check(base.ref_allele, max_size, error, error_message, "REF");
				}
				if(base.provided_alt)
				{
					base.col_index[+BASE_INDEX::ALT] = index_check(base.alt_allele, max_size, error, error_message, "ALT");
				}
				if(base.provided_bp)
				{
					base.col_index[+BASE_INDEX::BP] = index_check(base.bp, max_size, error, error_message, "BP");
				}
				if(base.provided_se)
				{
					base.col_index[+BASE_INDEX::SE] = index_check(base.standard_error, max_size, error, error_message, "SE");
				}

				base.col_index[+BASE_INDEX::P] = index_check(base.p_value, max_size, error, error_message, "P");
				base.col_index[+BASE_INDEX::STAT] = index_check(base.statistic, max_size, error, error_message, "STAT");
				base.col_index[+BASE_INDEX::RS] = index_check(base.snp, max_size, error, error_message, "RS");

			}
			if(base.col_index[+BASE_INDEX::P]==-1)
			{
				error = true;
				error_message.append("ERROR: No p-value column ("+base.p_value+") in file!\n");
			}
			else base.provided_p = true;
			if(base.col_index[+BASE_INDEX::STAT]==-1)
			{
				error = true;
				error_message.append("ERROR: No statistic column ("+base.statistic+") in file!\n");
			}
			else base.provided_stat = true;
			if(base.col_index[+BASE_INDEX::RS]==-1)
			{
				error = true;
				error_message.append("ERROR: No SNP name column ("+base.snp+") in file!\n");
			}
			else base.provided_snp = true;
			if(base.col_index[+BASE_INDEX::REF]==-1)
			{
				error = true;
				error_message.append("ERROR: No Reference allele column ("+base.ref_allele+") in file!\n");
			}else base.provided_ref = true;

			double max_index = *max_element(base.col_index.begin(), base.col_index.end());
			base.col_index[+BASE_INDEX::MAX] = max_index;
		}
	}
}

void Commander::clump_check(bool &error, std::string &error_message)
{
	if(!clumping.no_clump)
	{
		// require clumping
		if(clumping.provide_proxy && clumping.proxy <= 0)
		{
			error = true;
			error_message.append("ERROR: Proxy threshold cannot be negative!\n");
		}
		if(clumping.p_value < 0.0  || clumping.p_value > 1.0)
		{
			error = true;
			error_message.append("ERROR: P-value threshold must be within 0 and 1!\n");
		}
		if(clumping.r2 < 0.0  || clumping.r2 > 1.0)
		{
			error = true;
			error_message.append("ERROR: R2 threshold must be within 0 and 1!\n");
		}
		if(clumping.distance < 0.0)
		{
			error = true;
			error_message.append("ERROR: Clumping distance must be positive!\n");
		}
		if(!clumping.type.empty())
		{
			bool alright = false;
			for(auto &&type: supported_types)
			{
				if(clumping.type.compare(type)==0)
				{
					alright = true;
					break;
				}
			}
			if(!alright)
			{
				error =true;
				error_message.append("ERROR: Unsupported LD format: "+clumping.type+"\n");
			}
		}
	}
}

void Commander::misc_check(bool &error, std::string &error_message)
{
	if(misc.provided_permutation && misc.permutation < 0)
	{
		error =true;
		error_message.append("ERROR: Negative number of permutation!\n");
	}
	if(misc.thread <= 0)
	{
		error = true;
		error_message.append("ERROR: Number of thread must be larger than 1\n");
	}
	if(filter.use_info && (filter.info_score < 0 || filter.info_score > 1))
	{
	    error = true;
	    error_message.append("ERROR: Info score should be between 0 and 1\n");
	}
}

void Commander::prset_check(bool &error, std::string &error_message)
{
	if(!prset.gtf.empty() && prset.msigdb.empty())
	{
		error = true;
		error_message.append("ERROR: Must provide a gtf file if msigdb is specified\n");
	}
	if(prset.feature.empty())
	{
		prset.feature.push_back("exon");
		prset.feature.push_back("gene");
		prset.feature.push_back("protein_coding");
		prset.feature.push_back("CDS");
	}
}


void Commander::prsice_check(bool &error, std::string &error_message)
{
	if(prsice.fastscore && prsice.barlevel.size()==0)
	{
		fprintf(stderr, "barlevel not provided.\n");
		fprintf(stderr, "Will set to default: 0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5\n");
		prsice.barlevel = {0.001,0.05,0.1,0.2,0.3,0.4,0.5};
	}
	std::sort(prsice.barlevel.begin(), prsice.barlevel.end());
	prsice.barlevel.erase( std::unique(prsice.barlevel.begin(), prsice.barlevel.end()),
			prsice.barlevel.end());
	if(!prsice.fastscore)
	{
		if(prsice.no_regress)
		{
			error = true;
			error_message.append("ERROR: no-regress can only be used with fastscore!\n");
		}
		if(prsice.inter <= 0)
		{
			error = true;
			error_message.append("ERROR: Cannot have negative interval!\n");
		}
		if(prsice.upper < prsice.lower)
		{
			error = true;
			error_message.append("ERROR: Upper bound must be larger than lower bound!\n");
		}
	}
}

void Commander::prslice_check(bool &error, std::string &error_message)
{
	if(prslice.provided)
	{
		if(misc.all)
		{
			error = true;
			error_message.append("ERROR: Cannot output PRS for all threshold when using PRSlice!\n");
		}
		if(prslice.size <=0)
		{
			error = true;
			error_message.append("ERROR: PRSlice size cannot be less than 1!\n");
		}
	}
}

void Commander::target_check(bool &error, std::string &error_message)
{
	if(target.name.empty())
	{
		error = true;
		error_message.append("ERROR: You must provide a target file!\n");
	}

	bool alright = false;
	for(auto &&type: supported_types)
	{
		if(target.type.compare(type)==0)
		{
			alright = true;
			break;
		}
	}
	if(!alright)
	{
		error =true;
		error_message.append("ERROR: Unsupported target format: "+target.type+"\n");
	}
	if(target.pheno_file.empty() && target.is_binary.empty())
	{
		fprintf(stderr, "Target assumed to be binary\n");
		target.is_binary.push_back(true);
	}
	else
	{
		if(target.pheno_col.empty() && target.is_binary.size()==1)
		{
			// this is ok
		}
		else if(target.pheno_col.empty() && target.is_binary.empty())
		{
		    fprintf(stderr, "Phenotype assumed to be binary\n");
		    target.is_binary.push_back(true);
		}
		else if(target.pheno_col.size() <= 1 && target.is_binary.empty())
        {
            fprintf(stderr, "%s assumed to be binary\n", target.pheno_col.front().c_str());
            target.is_binary.push_back(true);
        }
		else if(target.pheno_col.size() != target.is_binary.size())
		{
			error = true;
	        error_message.append("ERROR: Number of target phenotypes doesn't match information of binary\n");
	        error_message.append("       target! You must indicate whether the phenotype is binary using\n");
	        error_message.append("       --binary-target\n");
		}
	}

	size_t num_bin=0;
	for(auto &&binary : target.is_binary)
	{
	    if(binary) num_bin++;
	}
	if(!target.prevalence.empty() && num_bin > target.prevalence.size()) // need to be all or nothing
	{
	    error = true;
	    error_message.append("ERROR: Number of target prevalence doesn't match number of binary traits\n");
	    error_message.append("       You must provide a prevalence for all binary trait(s) or not \n");
	    error_message.append("       provide any prevalence (all or nothing)\n");
	}
	for(auto &&prev : target.prevalence)
	{
	    if(prev > 1.0 || prev < 0.0)
	    {
	        error = true;
	        error_message.append("ERROR: Prevalence cannot be bigger than 1.0 or smaller than 0.0\n");
	        break;
	    }
	}

}

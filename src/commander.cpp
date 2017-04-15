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
        {"lower",required_argument,NULL,'l'},
        {"ld",required_argument,NULL,'L'},
        {"msigdb",required_argument,NULL,'m'},
        {"thread",required_argument,NULL,'n'},
        {"out",required_argument,NULL,'o'},
        {"pvalue",required_argument,NULL,'p'},
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
            }else if(command.compare("proxy")==0)
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
            		misc.permutation = misc::convert<int>(optarg);
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
    base.provided_beta = false;
    base.col_index.resize(+BASE_INDEX::MAX+1, -1);

	clumping.ld = "";
	clumping.type = "bed";
	clumping.no_clump = false;
	clumping.provide_proxy = false;
	clumping.proxy =-1.0;
	clumping.p_value =1.0;
	clumping.r2 = 0.1;
	clumping.distance = 250000;

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
	/*
	m_help_messages.push_back(help("Required", 'b', "base", "Base association files. "
			"User can provide multiple base files"));
	m_help_messages.push_back(help("Required", 't', "target", "Plink binary file prefix for target files. "
			"Currently only support plink binary inputs. For multiple target phenotypes, user should use the "
			"--pheno_file option together with the pheno_col option. For multiple chromosome input, one "
			"should substitute the chromosome number with #. For example, if you files are presented as "
			"genotype_chr1_test, genotype_chr2_test, genotype_chr3_test, then you can use: genotype_chr#_test."
			"Please note that the substitute is based on your base file. So if your base file code chromosome "
			"with chr, e.g. chr1 chr2 etc, then in our example case, you should code your plink file as "
			"genotype_#_test"));
	m_help_messages.push_back(help("Required", '\0', "binary-target", "Indicate whether the target sample has "
			"binary phenotype or not. For each phenotype, user need to provide either T or F where T"
			" means the phenotype is binary"));
	m_help_messages.push_back(help("Required", '\0', "beta", "Indicate whether the test statistic is beta "
			"instead of OR. Must be of the same length as base"));
	m_help_messages.push_back(help("Options", 'f', "pheno-file", "Phenotype file containing the target "
			"phenotype(s). If provided, the fam file of the target is ignored. First column must be "
			"IID of the samples. Must contain a header if pheno_col is specified"));
	m_help_messages.push_back(help("Options", '\0',"pheno-col", "Headers of pheenotypes from phenotype file"));
	m_help_messages.push_back(help("Options", 'L',"ld", "Plink binary file prefix for the reference file "
			"used for LD calculation. If not provided, will use the target genotype for the LD calculation. "
			"Can also use multiple chromosome plink file. Please see --target for more information."));
	m_help_messages.push_back(help("Options", 'c', "cov-header", "Header of covariates. If not provided, "
			"will use all variable in the covariate file as the covarite."));
	m_help_messages.push_back(help("Options", 'C', "cov-file", "Covarite file. Formate should be: ID Cov1 Cov2"));
	m_help_messages.push_back(help("Options", '\0', "full", "Also include the full model in the PRSice output"));
	m_help_messages.push_back(help("Options", '\0', "all", "Output PRS for ALL threshold. Can only be used together "
			"with fastscore to avoid huge output files."));
	m_help_messages.push_back(help("Options",'\0', "no-regress", "Do not perform the regression analysis and "
			"simply output all PRS. Can only be used together with fastscore to avoid huge output files. If "
			"you must, you can modify bar_levels to obtain the fine scale PRS outputs"));
	m_help_messages.push_back(help("Options",'o', "out", "Prefix of all output. Default: "+m_out));
	m_help_messages.push_back(help("Scoring options",'l', "lower", "The starting p-value threshold. "
			"Default: "+std::to_string(m_lower)));
	m_help_messages.push_back(help("Scoring options",'u', "upper", "The final p-value threshold. "
			"Default: "+std::to_string(m_upper)));
	m_help_messages.push_back(help("Scoring options",'i', "interval", "The step size of the threshold. "
			"Default: "+std::to_string(m_inter)));
	m_help_messages.push_back(help("Scoring options",'\0', "fastscore", "Calculate the minimum amount of threshold as "
			"required by the bar_level option"));
	m_help_messages.push_back(help("Scoring options",'\0', "score", "Method to handle missing genotypes. By default, "
			"final scores are averages of valid per-allele scores with missing genotypes contribute an amount "
			"proportional to imputed allele frequency. To throw out missing observations instead (decreasing "
			"the denominator in the final average when this happens), use the 'no_mean_imputation' modifier."
			" Alternatively, you can use the 'center' modifier to shift all scores to mean zero. "));
	m_help_messages.push_back(help("File Headers", '\0', "chr", "Column header of Chromosome <Required>"));
	m_help_messages.push_back(help("File Headers", '\0', "A1", "Column header of Reference Allele <Required>"));
	m_help_messages.push_back(help("File Headers", '\0', "A2", "Column header of Alternaative Allele"));
	m_help_messages.push_back(help("File Headers", '\0', "stat", "Column header of test statistic <Required>"));
	m_help_messages.push_back(help("File Headers", '\0', "snp", "Column header of SNP id"));
	m_help_messages.push_back(help("File Headers", '\0', "bp", "Column header of SNP location"));
	m_help_messages.push_back(help("File Headers", '\0', "se", "Column header of Standard Error"));
	m_help_messages.push_back(help("File Headers", 'p', "pvalue", "Column header of p-value <Required> "));
	m_help_messages.push_back(help("File Headers", '\0', "index", "Indicate all the above options are providing "
			"the INDEX of the corresponding column. (Index should be 0-based). Useful when your base file "
			"each have a different header but the column index remains the same"));
	m_help_messages.push_back(help("Clumping", '\0', "clump-p", "The p-value threshold use for clumping."
			"Default: "+std::to_string(m_clump)));
	m_help_messages.push_back(help("Clumping", '\0', "clump-r2", "The R2 threshold for clumping. Please note that "
			"as we did not implement the maximum likelihood R2 calculation, the clumping result can differ "
			"slightly from plink. Default: "+std::to_string(m_clump_r2)));
	m_help_messages.push_back(help("Clumping", '\0', "clump-kb", "The distance for clumping in kb."
			"Default: "+std::to_string(m_clump_kb/1000)));
	m_help_messages.push_back(help("PRSet", 'B', "bed", "Bed file containing the selected regions. "
			"Name of bed file will be used as the region identifier."));
	m_help_messages.push_back(help("PRSet", 'g', "gtf", "GTF file containing gene boundaries. Required "
			"when --msigdb is set."));
	m_help_messages.push_back(help("PRSet", 'm', "msigdb", "MSIGDB file containing the pathway information "
			"require the gtf file."));
	m_help_messages.push_back(help("PRSet", '\0', "print-all", "Print the detail report for all sets"));
	m_help_messages.push_back(help("PRSet", '\0', "proxy", "Proxy threshold for index SNP to be considered "
			"as part of the region represented by the clumped SNPs. e.g. --proxy 0.8 means the index SNP will "
			"represent the region of any clumped SNPs that has a R2 >= 0.8 with it even if it is not physically "
			"within these regions"));
	m_help_messages.push_back(help("PRSet", '\0', "feature", "Features to be included from the gtf file. Default "
			"is exon, CDS, gene and protein_coding. If this parameter is provided, all default will be ignored."));
	m_help_messages.push_back(help("PRSlice", '\0', "prslice", "Perform PRSlice where the whole genome is first "
			"cut into bin size specified by this option. PRSice will then be performed on each bin. Bins are "
			"then sorted according to the their R2. PRSice is then performed again to find the best bin combination."
			" This cannot be performed together with PRSet"));
	m_help_messages.push_back(help("Plotting", '\0', "bar-levels", "Level of barchart to be plotted. When fastscore "
			"is set, PRSice will only calculate the PRS for threshold within the bar level"));
	m_help_messages.push_back(help("Misc", '\0', "ignore-fid", "Ignore the FID field for covariate and phenotype "
			"matching. When set, assume first column of phenotype file as IID, otherwise, assume first column "
			"of phenotype file as FID and the second column as IID."));
	m_help_messages.push_back(help("Misc", '\0', "print-snp", "Print out the SNP(s) used for constructing "
			"the best PRS score"));
	m_help_messages.push_back(help("Misc", '\0', "perm", "Number of permutation to perform. When this parameter is provided,"
			" permutation will be performed to obtain an empirical P-value. This will significantly increases the run time "
			"of PRSice."));
	m_help_messages.push_back(help("Misc", 'T', "thread", "Number of thread use"));
	m_help_messages.push_back(help("Misc", 'h', "help", "Display this help message"));

	*/
}
void Commander::usage()
{
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
	std::cerr << "GNU General Public License v3" << std::endl;
}

void Commander::user_input() const{
	fprintf(stderr, "\nUser Input\n");
	fprintf(stderr, "==============================\n");
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
	if(misc.permutation < 0)
	{
		error =true;
		error_message.append("ERROR: Negative number of permutation!\n");
	}
	if(misc.thread <= 0)
	{
		error = true;
		error_message.append("ERROR: Number of thread must be larger than 1\n");
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
		if(target.pheno_col.size() <= 1 && target.is_binary.empty())
		{
			fprintf(stderr, "%s assumed to be binary\n", target.pheno_col.front().c_str());
			target.is_binary.push_back(true);
		}
		else if(target.pheno_col.empty() && target.is_binary.size()==1)
		{
			// this is ok
		}
		else if(target.pheno_col.size() != target.is_binary.size())
		{
			error = true;
	        error_message.append("ERROR: Number of target phenotypes doesn't match information of binary\n");
	        error_message.append("       target! You must indicate whether the phenotype is binary using\n");
	        error_message.append("       --binary-target\n");
		}
	}

}

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
    static const char *optString = "b:t:c:C:a:f:L:p:T:u:l:i:B:g:m:o:h?";
    static const struct option longOpts[]=
    {
        {"base",required_argument,NULL,'b'},
        {"target",required_argument,NULL,'t'},
        {"covar_header",required_argument,NULL,'c'},
        {"covar_file",required_argument,NULL,'C'},
        {"ancestry",required_argument,NULL,'a'},
        {"pheno_file",required_argument,NULL,'f'},
        {"pheno_col",required_argument,NULL,0},
        {"ld",required_argument,NULL,'L'},
        {"pvalue",required_argument,NULL,'p'},
        {"thread",required_argument,NULL,'T'},
        {"upper",required_argument,NULL,'u'},
        {"lower",required_argument,NULL,'l'},
        {"interval",required_argument,NULL,'i'},
        {"bed",required_argument,NULL,'B'},
        {"gtf",required_argument,NULL,'g'},
        {"msigdb",required_argument,NULL,'m'},
        {"out",required_argument,NULL,'o'},
        {"beta",required_argument,NULL,0},
        {"chr",required_argument,NULL,0},
        {"A1",required_argument,NULL,0},
        {"A2",required_argument,NULL,0},
        {"stat",required_argument,NULL,0},
        {"snp",required_argument,NULL,0},
        {"bp",required_argument,NULL,0},
        {"se",required_argument,NULL,0},
        {"clump_p",required_argument,NULL,0},
        {"clump_r2",required_argument,NULL,0},
        {"clump_kb",required_argument,NULL,0},
        {"binary_target",required_argument,NULL,0},
        {"bar_levels",required_argument,NULL,0},
        {"gen_bed",no_argument,NULL,0},
        {"index",no_argument,NULL,0},
        {"all",no_argument,NULL,0},
        {"full",no_argument,NULL,0},
        {"print_all",no_argument,NULL,0},
        {"no_regression",no_argument,NULL,0},
        {"fastscore",no_argument,NULL,0},
        {"proxy",required_argument,NULL,0},
        {"feature",required_argument,NULL,0},
        {"perm",required_argument,NULL,0},
        {"prslice",required_argument,NULL,0},
        {"no_clump",no_argument,NULL,0},
        {"help",no_argument,NULL,'h'},
        {NULL, 0, 0, 0}
    };

    bool error = false;
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
            if(command.compare("chr")==0) m_chr=optarg;
            else if(command.compare("A1")==0) m_ref_allele = optarg;
            else if(command.compare("A2")==0) m_alt_allele=optarg;
            else if(command.compare("stat")==0)
            {
            	m_statistic = optarg;
            	m_stat_provided = true;
            }
            else if(command.compare("snp")==0) m_snp=optarg;
            else if(command.compare("bp")==0) m_bp=optarg;
            else if(command.compare("se")==0) m_standard_error = optarg;
            else if(command.compare("index")==0) m_index = true;
            else if(command.compare("all")==0) m_all = true;
            else if(command.compare("full")==0) m_full = true;
            else if(command.compare("print_all")==0) m_print_all = true;
            else if(command.compare("no_clump")==0) m_no_clump = true;
            else if(command.compare("clump_p")==0)
            {
                double temp = atof(optarg);
                if(temp < 0.0 || temp > 1.0)
                {
                    error = true;
                    error_message.append("Clumping p-values must be >=0 and <= 1.0\n");
                }
                else m_clump = temp;
            }
            else if(command.compare("clump_r2")==0)
            {
                double temp = atof(optarg);
                if(temp < 0.0 || temp > 1.0)
                {
                    error = true;
                    error_message.append("Clumping R2 must be >=0 and <= 1.0\n");
                }
                if(temp == 0.0)
                {
                    std::cerr << "WARNING: As clumping R2==0, no clumping will be performed" << std::endl;
                }
                else m_clump_r2 = temp;
            }
            else if(command.compare("clump_kb")==0)
            {
                int temp = atoi(optarg);
                if(temp <= 0.0)
                {
                    error = true;
                    error_message.append("Clumping window size must be larger than 0kb\n");
                }
                else m_clump_kb = temp*1000; //change it to kb, might want to allow different units
            }
            else if(command.compare("prslice")==0)
            {
                int temp = atoi(optarg);
                if(temp <= 0.0)
                {
                    error = true;
                    error_message.append("PRSlice bin size must be larger than 0\n");
                }
                else m_prslice_size = temp*1000; //change it to kb, might want to allow different units
            }
            else if(command.compare("binary_target")==0)
            {
            		std::vector<std::string> token = misc::split(optarg, ", ");
            	    for(size_t i = 0; i < token.size(); ++i) m_target_is_binary.push_back(misc::to_bool(token[i]));
            }
            else if(command.compare("bar_levels")==0)
            {
                std::vector<std::string> token = misc::split(optarg, ", ");
                try
                {
                    for(size_t i = 0; i < token.size(); ++i)
                    {
                        double temp = misc::convert<double>(token[i]);
                        if(temp < 0)
                        {
                            error_message.append("ERROR: Negative barchart level\n");
                            error=true;
                            break;
                        }
                        m_barlevel.push_back(temp);
                    }
                }
                catch(const std::runtime_error &er)
                {
                    error_message.append("ERROR: None numeric barchart level\n");
                    error=true;
                }
            }
            else if(command.compare("beta")==0)
            {
                std::vector<std::string> token = misc::split(optarg, ", ");
                for(size_t i = 0; i < token.size(); ++i) m_use_beta.push_back(misc::to_bool(token[i]));
            }
            else if(command.compare("gen_bed")==0) m_gen_bed = true;
            else if(command.compare("no_regression")==0) m_no_regress = true;
            else if(command.compare("fastscore")==0) m_fastscore = true;
            else if(command.compare("proxy")==0)
            {
                try
                {
                    m_proxy = misc::convert<double>(optarg);
                    if(m_proxy<=0.0)
                    {
                        error_message.append("ERROR: Proxy must be bigger than 0.0\n");
                        error=true;
                    }
                }
                catch(const std::runtime_error &er)
                {
                    error_message.append("ERROR: Proxy provided isn't numeric\n");
                    error=true;
                }
            }
            else if(command.compare("pheno_col")==0)
            {
                std::vector<std::string> token = misc::split(optarg, ", ");
                m_pheno_col.insert(m_pheno_col.end(), token.begin(), token.end());
            }
            else if(command.compare("perm")==0)
            {
            	int temp = atoi(optarg);
            	if(temp < 0.0)
            	{
            		error = true;
            		error_message.append("Number of permutation must be bigger than 0\n");
            	}
            	else m_permutation = temp;

            }
            else if(command.compare("feature")==0)
            {
            	std::vector<std::string> token = misc::split(optarg, ", ");
            	m_feature.insert(m_feature.end(), token.begin(), token.end());
            }
            else
            {
                std::string er = "Undefined operator: "+command+", please use --help for more information!";
                throw std::runtime_error(er);
            }
            break;
        case 'b':
        {
            std::vector<std::string> token= misc::split(optarg, ", ");
            m_base.insert(m_base.end(), token.begin(), token.end());
            if(m_base.size() ==0)
            {
                error = true;
                error_message.append("You must provide at least one valid base file name\n");
            }
        }
        break;
        case 't':
            m_target = optarg;
            break;
        case 'c':
        {
            std::vector<std::string> token= misc::split(optarg, ", ");
            m_covariates.insert(m_covariates.end(), token.begin(), token.end());
        }
        break;
        case 'C':
            m_covariate_file = optarg;
            break;
        case 'a':
            m_ancestry_dim = optarg;
            if(m_ancestry_dim.compare("MDS") != 0 && m_ancestry_dim.compare("mds") != 0 &&
                    m_ancestry_dim.compare("PCA") != 0 && m_ancestry_dim.compare("pca") != 0 )
            {
                error = true;
                error_message.append("Only support PCA and MDS for the calculation of ancestry information\n");
            }
            fprintf(stderr, "Currently we have not implement this function\n");
            break;
        case 'f':
            m_pheno_file = optarg;
            break;
        case 'p': // the index/header of p-value in the file
            m_p_value = optarg;
            break;
        case 'L':
            m_ld_prefix=optarg;
            break;
        case 'T':
        {
            int temp = atoi(optarg);
            if(temp<=0)
            {
                std::cerr << "Number of thread cannot be less than 1" << std::endl;
                std::cerr << "Will set to 1 instead" << std::endl;
                m_thread = 1;
            }
            else m_thread = temp;
        }
        break;
        case 'u':
        {
            double temp = atof(optarg);
            if(temp <= 0.0 || temp > 1.0)
            {
                error = true;
                error_message.append("Upper interval must be > 0 and <= 1.0\n");
            }
            else m_upper = temp;
        }
        break;
        case 'l':
        {
            double temp = atof(optarg);
            if(temp < 0.0 || temp >1.0)
            {
                error = true;
                error_message.append("Lower interval must be >= 0 and < 1.0\n");
            }
            else m_lower = temp;
        }
        break;
        case 'i':
        {
            double temp = atof(optarg);
            if(temp <= 0.0 || temp > 1.0)
            {
                error = true;
                error_message.append("Interval must be >=0 and <= 1.0\n");
            }
            else m_inter = temp;
        }
        break;
        case 'B':
        {
            std::vector<std::string> token = misc::split(optarg, ", ");
            m_bed_list.insert(m_bed_list.end(), token.begin(), token.end());
        }
        break;
        case 'g':
            m_gtf = optarg;
            break;
        case 'm':
            m_msigdb= optarg;
            break;
        case 'o':
            m_out = optarg;
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
    if(!m_stat_provided && m_beta_provided)
    {
    	m_statistic = (m_use_beta.front())? "BETA" : "OR";
    }
    if(m_base.size()==0)
    {
        error=true;
        error_message.append("There are no base file to run\n");
    }
    if(m_target.size()==0)
    {
        error=true;
        error_message.append("There are no target file to run\n");
    }
    // Start performing the check on the inputs

    if(!m_msigdb.empty() && m_gtf.empty())
    {
        error = true;
        error_message.append("Must provide the GTF file when only MSIGDB file is provided\n");
    }
    if(m_gen_bed && m_gtf.empty())
    {
        fprintf(stderr, "ERROR: Cannot generate gene bed file without given the gtf file!\n");
        fprintf(stderr, "       Will not generate the gene bed file\n");
    }
    if(m_out.empty())
    {
        fprintf(stderr, "WARNING: Output prefix is empty, will set it to PRSice\n");
        m_out = "PRSice";
    }
    // add default binary
    if((!m_msigdb.empty() || m_bed_list.size() != 0) && m_chr.empty() && m_bp.empty())
    {
        fprintf(stderr, "WARNING: For pathway/region PRSice to work, you must provide\n");
        fprintf(stderr, "         the chromosome and bp information\n");
        fprintf(stderr, "         As chromosome / bp information were not provided, we\n");
        fprintf(stderr, "         will disable the pathwya/region PRSice analysis\n");
        m_msigdb = "";
        m_bed_list.clear();
    }
    if(m_use_beta.size()==0)
    {
        for(size_t i = 0; i < m_base.size(); ++i)
        {
            m_use_beta.push_back(false); // default is binary
        }
        m_beta_provided = true;
    }
    else if(m_use_beta.size() != m_base.size())
    {
    		error=true;
        error_message.append("ERROR: Number of beta doesn't match number of base file!\n");
        error_message.append("       Default value only work when all base file are using OR and\n");
        error_message.append("       when --beta is not used\n");
    }
    if(m_pheno_col.size()!=0 && m_pheno_col.size()!=m_target_is_binary.size())
    {
    		error=true;
        error_message.append("ERROR: Number of target phenotypes doesn't match information of binary\n");
        error_message.append("       target! You must indicate whether the phenotype is binary using\n");
        error_message.append("       --binary_target\n");

    }
    if(m_no_regress && !m_fastscore)
    {
        fprintf(stderr, "WARNING: To limit the amount of output,\n");
        fprintf(stderr, "         no-regress can only be used with\n");
        fprintf(stderr, "         fastscore. Will use fastscore\n");
        m_fastscore=true;
    }
    if(m_no_regress) m_all=true;
    if(m_fastscore && m_barlevel.size()==0)
    {
        fprintf(stderr, "barlevel not provided. Will set to default: 0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5\n");
        m_barlevel = {0.001,0.05,0.1,0.2,0.3,0.4,0.5};
    }
    std::sort(m_barlevel.begin(), m_barlevel.end());
    if(m_feature.empty())
    {
    	m_feature.push_back("exon");
    	m_feature.push_back("gene");
    	m_feature.push_back("protein_coding");
    	m_feature.push_back("CDS");
    }
    if(m_inter <=0.0 && !m_fastscore) // double comparison error, need to optimize it
    {
    	error = true;
    	error_message.append("Error: Interval cannot be 0!\n");
    }
    if(error) throw std::runtime_error(error_message);
    return true;
}


Commander::Commander()
{
    // should gives the default here
    m_target="";
    m_pheno_file="";
    m_covariate_file="";
    m_ancestry_dim="MDS";
    m_chr = "CHR";
    m_ref_allele="A1";
    m_alt_allele="A2";
    m_statistic = "OR";
    m_snp="SNP";
    m_bp="BP";
    m_standard_error = "SE";
    m_p_value = "P";
    m_ld_prefix="";
    m_gtf="";
    m_msigdb="";
    m_out = "PRSice";
    m_fastscore =false;
    m_index =false;
    m_gen_bed = false;
    m_no_regress = false;
    m_all = false;
    m_full = false;
    m_beta_provided = false;
    m_stat_provided = false;
    m_proxy = -1.0;
    m_no_clump = false;
    m_clump = 1.0;
    m_clump_r2 = 0.1;
    m_clump_kb = 250000;
    m_permutation = 0;
    m_lower = 0.0001;
    m_upper = 0.5;
    m_inter = 0.00005;
    m_prslice_size = -1.0;
    m_thread=1;
    m_print_all= false;
}

Commander::~Commander()
{
    //dtor
}

void Commander::info()
{
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
	m_help_messages.push_back(help("Required", '\0', "binary_target", "Indicate whether the target sample has "
			"binary phenotype or not. For each phenotype, user need to provide either T or F where T"
			" means the phenotype is binary"));
	m_help_messages.push_back(help("Required", '\0', "beta", "Indicate whether the test statistic is beta "
			"instead of OR. Must be of the same length as base"));
	m_help_messages.push_back(help("Options", 'f', "pheno_file", "Phenotype file containing the target "
			"phenotype(s). If provided, the fam file of the target is ignored. First column must be "
			"IID of the samples. Must contain a header if pheno_col is specified"));
	m_help_messages.push_back(help("Options", '\0',"pheno_col", "Headers of pheenotypes from phenotype file"));
	m_help_messages.push_back(help("Options", 'L',"ld", "Plink binary file prefix for the reference file "
			"used for LD calculation. If not provided, will use the target genotype for the LD calculation. "
			"Can also use multiple chromosome plink file. Please see --target for more information."));
	m_help_messages.push_back(help("Options", 'c', "covar_header", "Header of covariates. If not provided, "
			"will use all variable in the covariate file as the covarite."));
	m_help_messages.push_back(help("Options", 'C', "covar_file", "Covarite file. Formate should be: ID Cov1 Cov2"));
	m_help_messages.push_back(help("Options", '\0', "full", "Also include the full model in the PRSice output"));
	m_help_messages.push_back(help("Options", '\0', "all", "Output PRS for ALL threshold. Can only be used together "
			"with fastscore to avoid huge output files."));
	m_help_messages.push_back(help("Options",'\0', "no_regress", "Do not perform the regression analysis and "
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
	m_help_messages.push_back(help("Clumping", '\0', "clump_p", "The p-value threshold use for clumping."
			"Default: "+std::to_string(m_clump)));
	m_help_messages.push_back(help("Clumping", '\0', "clump_r2", "The R2 threshold for clumping. Please note that "
			"as we did not implement the maximum likelihood R2 calculation, the clumping result can differ "
			"slightly from plink. Default: "+std::to_string(m_clump_r2)));
	m_help_messages.push_back(help("Clumping", '\0', "clump_kb", "The distance for clumping in kb."
			"Default: "+std::to_string(m_clump_kb/1000)));
	m_help_messages.push_back(help("PRSet", 'B', "bed", "Bed file containing the selected regions. "
			"Name of bed file will be used as the region identifier."));
	m_help_messages.push_back(help("PRSet", 'g', "gtf", "GTF file containing gene boundaries. Required "
			"when --msigdb is set."));
	m_help_messages.push_back(help("PRSet", 'm', "msigdb", "MSIGDB file containing the pathway information "
			"require the gtf file."));
	m_help_messages.push_back(help("PRSet", '\0', "gen_bed", "Generate bed file of gene regions from "
			" the gtf file."));
	m_help_messages.push_back(help("PRSet", '\0', "print_all", "Print the detail report for all sets"));
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
	m_help_messages.push_back(help("Plotting", '\0', "bar_levels", "Level of barchart to be plotted. When fastscore "
			"is set, PRSice will only calculate the PRS for threshold within the bar level"));
	m_help_messages.push_back(help("Misc", '\0', "perm", "Number of permutation to perform. When this parameter is provided,"
			" permutation will be performed to obtain an empirical P-value. This will significantly increases the run time "
			"of PRSice."));
	m_help_messages.push_back(help("Misc", 'T', "thread", "Number of thread use"));
	m_help_messages.push_back(help("Misc", 'h', "help", "Display this help message"));
}
void Commander::usage()
{
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
}

void Commander::program_info()
{
	std::cerr << "####################################################################" << std::endl;
	std::cerr << "#                                                                  #" << std::endl;
	std::cerr << "#                                                                  #" << std::endl;
	std::cerr << "# PRSice: Polygenic Risk Score software                            #" << std::endl;
	std::cerr << "#                                                                  #" << std::endl;
	std::cerr << "# Jack Euesden, Cathryn M. Lewis, Paul F. O'Reilly 2014            #" << std::endl;
	std::cerr << "#                                                                  #" << std::endl;
	std::cerr << "#                                                                  #" << std::endl;
	std::cerr << "# If you use PRSice in published work, please cite:                #" << std::endl;
	std::cerr << "#                                                                  #" << std::endl;
	std::cerr << "#                                                                  #" << std::endl;
	std::cerr << "# \"PRSice: Polygenic Risk Score software\"                          #" << std::endl;
	std::cerr << "# Euesden, Lewis, O'Reilly, Bioinformatics (2015) 31 (9):1466-1468 #" << std::endl;
	std::cerr << "#                                                                  #" << std::endl;
	std::cerr << "#                                                                  #" << std::endl;
	std::cerr << "####################################################################" << std::endl;


}

void Commander::user_input() const{
	fprintf(stderr, "\nUser Input\n");
	fprintf(stderr, "==============================\n");
	fprintf(stderr, "Base files: ");
	for (size_t i_base = 0; i_base < m_base.size(); ++i_base) {
			if (!m_use_beta.at(i_base)) {
				fprintf(stderr, "%s(OR) ", m_base[i_base].c_str());
			} else {
				fprintf(stderr, "%s(Beta) ", m_base[i_base].c_str());
			}
		}
		fprintf(stderr, "\n");
		fprintf(stderr, "\nUser Defined Column Headers\n");
		fprintf(stderr, "==============================\n");
		if (!m_chr.empty())
			fprintf(stderr, "Chr            : %s\n", m_chr.c_str());
		fprintf(stderr, "SNP            : %s\n", m_snp.c_str());
		if (!m_bp.empty())
			fprintf(stderr, "BP             : %s\n", m_bp.c_str());
		fprintf(stderr, "Ref Allele     : %s\n", m_ref_allele.c_str());
		if (!m_alt_allele.empty())
			fprintf(stderr, "Alt Allele     : %s\n", m_alt_allele.c_str());
		if (!m_statistic.empty())
			fprintf(stderr, "Statistic      : %s\n", m_statistic.c_str());
		if (!m_standard_error.empty())
			fprintf(stderr, "Standard Error : %s\n", m_standard_error.c_str());
		fprintf(stderr, "P-value        : %s\n", m_p_value.c_str());
		fprintf(stderr, "\nClumping Parameters: \n");
		fprintf(stderr, "==============================\n");
		fprintf(stderr, "P-Threshold  : %f\n", m_clump);
		fprintf(stderr, "R2-Threshold : %f\n", m_clump_r2);
		fprintf(stderr, "Window Size  : %zu\n", m_clump_kb);
}

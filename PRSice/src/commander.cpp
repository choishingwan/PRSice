#include "../inc/commander.hpp"

bool Commander::initialize(int argc, char *argv[]){
    if(argc<=1){
        usage();
        throw std::runtime_error("Please provide the required parameters");
    }
    static const char *optString = "b:t:c:C:a:f:L:p:T:u:l:i:B:g:m:o:h?";
    static const struct option longOpts[]={
        {"base",required_argument,NULL,'b'},
        {"target",required_argument,NULL,'t'},
        {"covar_header",required_argument,NULL,'c'},
        {"covar_file",required_argument,NULL,'C'},
        {"ancestry",required_argument,NULL,'a'},
        {"pheno_file",required_argument,NULL,'f'},
        {"ld",required_argument,NULL,'L'},
        {"beta",required_argument,NULL,0},
        {"chr",required_argument,NULL,0},
        {"A1",required_argument,NULL,0},
        {"A2",required_argument,NULL,0},
        {"stat",required_argument,NULL,0},
        {"snp",required_argument,NULL,0},
        {"bp",required_argument,NULL,0},
        {"se",required_argument,NULL,0},
        {"pvalue",required_argument,NULL,'p'},
        {"upper",required_argument,NULL,'u'},
        {"lower",required_argument,NULL,'l'},
        {"interval",required_argument,NULL,'i'},
        {"clump-p",required_argument,NULL,0},
//        {"clump_p2",required_argument,NULL,0},
        {"clump_r2",required_argument,NULL,0},
        {"clump_kb",required_argument,NULL,0},
        {"binary_target",required_argument,NULL,0},
        {"thread",required_argument,NULL,'T'},
        {"bed",required_argument,NULL,'B'},
        {"gtf",required_argument,NULL,'g'},
        {"gen_bed",no_argument,NULL,0},
        {"msigdb",required_argument,NULL,'m'},
        {"out",required_argument,NULL,'o'},
        {"index",no_argument,NULL,0},
        {"no-regression",no_argument,NULL,0},
        {"fastscore",no_argument,NULL,0},
        {"proxy",no_argument,NULL,0},
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
    while(opt!=-1){
        switch(opt){
            case 0:
                command = longOpts[longIndex].name;
                if(command.compare("chr")==0) m_chr=optarg;
                else if(command.compare("A1")==0) m_ref_allele = optarg;
                else if(command.compare("A2")==0) m_alt_allele=optarg;
                else if(command.compare("stat")==0) m_statistic = optarg;
                else if(command.compare("snp")==0) m_snp=optarg;
                else if(command.compare("bp")==0) m_bp=optarg;
                else if(command.compare("se")==0) m_standard_error = optarg;
                else if(command.compare("index")==0) m_index = true;
                else if(command.compare("clump-p")==0){
                    double temp = atof(optarg);
                    if(temp < 0.0 || temp > 1.0){
                        error = true;
                        error_message.append("Clumping p-values must be >=0 and <= 1.0\n");
                    }
                    else m_clump = temp;
                }
                else if(command.compare("clump_r2")==0){
                    double temp = atof(optarg);
                    if(temp < 0.0 || temp > 1.0){
                        error = true;
                        error_message.append("Clumping R2 must be >=0 and <= 1.0\n");
                    }
                    if(temp == 0.0){
                        std::cerr << "WARNING: As clumping R2==0, no clumping will be performed" << std::endl;
                    }
                    else m_clump_r2 = temp;
                }
                else if(command.compare("clump_kb")==0){
                    int temp = atoi(optarg);
                    if(temp <= 0.0){
                        error = true;
                        error_message.append("Clumping window size must be larger than 0kb\n");
                    }
                    else m_clump_kb = temp;
                }
                else if(command.compare("binary_target")==0){
                		std::vector<std::string> token = misc::split(optarg, ", ");
                		for(size_t i = 0; i < token.size(); ++i) m_target_is_binary.push_back(misc::to_bool(token[i]));
                }
                else if(command.compare("beta")==0){
                		std::vector<std::string> token = misc::split(optarg, ", ");
                		for(size_t i = 0; i < token.size(); ++i) m_use_beta.push_back(misc::to_bool(token[i]));
                	}
                else if(command.compare("gen_bed")==0) m_gen_bed = true;
                else if(command.compare("no-regression")==0) m_no_regress = true;
                else if(command.compare("fastscore")==0) m_fastscore = true;
                else if(command.compare("proxy")==0) m_proxy = true;
                else{
                		std::string er = "Undefined operator: "+command+", please use --help for more information!";
                		throw std::runtime_error(er);
                }
                break;
            case 'b':
                m_base= misc::split(optarg, ", ");
                if(m_base.size() ==0){
                    error = true;
                    error_message.append("You must provide at least one valid base file name\n");
                }
                break;
            case 't':
                m_target = misc::split(optarg, ", ");
                if(m_target.size() ==0){
                    error = true;
                    error_message.append("You must provide at least one valid target file name\n");
                }
                break;
            case 'c':
                m_covariates = misc::split(optarg, ", ");
                break;
            case 'C':
                m_covariate_file = optarg;
                break;
            case 'a':
                m_ancestry_dim = optarg;
                if(m_ancestry_dim.compare("MDS") != 0 && m_ancestry_dim.compare("mds") != 0 &&
                   m_ancestry_dim.compare("PCA") != 0 && m_ancestry_dim.compare("pca") != 0 ){
                    error = true;
                    error_message.append("Only support PCA and MDS for the calculation of ancestry information\n");
                }
                break;
            case 'f':
                m_pheno_file = optarg;
                if(m_pheno_file.empty()){
                    error = true;
                    error_message.append("No parameter is given for phenotype file\n");
                }
                break;
            case 'p':
                m_p_value = optarg;
                break;
            case 'L':
                m_ld_prefix=optarg;
                break;
            case 'T':
            {
                int temp = atoi(optarg);
                if(temp<=0){
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
                if(temp <= 0.0 || temp > 1.0){
                    error = true;
                    error_message.append("Upper interval must be > 0 and <= 1.0\n");
                }
                else m_upper = temp;
            }
                break;
            case 'l':
            {
                double temp = atof(optarg);
                if(temp < 0.0 || temp >1.0){
                    error = true;
                    error_message.append("Lower interval must be >= 0 and < 1.0\n");
                }
                else m_lower = temp;
            }
                break;
            case 'i':
            {
                double temp = atof(optarg);
                if(temp <= 0.0 || temp > 1.0){
                    error = true;
                    error_message.append("Interval must be >=0 and <= 1.0\n");
                }
                else m_inter = temp;
            }
                break;
            case 'B':
                m_bed_list= misc::split(optarg, ", ");
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
    if(m_target.size() != 1 && m_target.size() != m_target_is_binary.size()){
        error=true;
        error_message.append("Length of binary target list does not match number of target\n");
    }
    if(!m_msigdb.empty() && m_gtf.empty()){
        error = true;
        error_message.append("Must provide the GTF file when only MSIGDB file is provided\n");
    }
    if(m_gen_bed && m_gtf.empty()){
        fprintf(stderr, "ERROR: Cannot generate gene bed file without given the gtf file!\n");
        fprintf(stderr, "       Will not generate the gene bed file\n");
    }
    if(m_out.empty()){
        fprintf(stderr, "ERROR: Output prefix is empty, will set it to default\n");
        m_out = "PRSice";
    }
    // add default binary
    if(m_target_is_binary.size()<m_target.size() && m_target.size()==1){
    		 m_target_is_binary.push_back(true); // default is binary
    }
    else if(m_target_is_binary.size() != m_target.size()){
		error_message.append("ERROR: Number of binary target doesn't match number of target file!\n");
		error_message.append("       Default value only work when all target file are binary and\n");
		error_message.append("       when --binary_target is not used\n");
    }
    if(m_use_beta.size()<m_base.size() && m_base.size() ==1){
    		for(size_t i = m_use_beta.size(); i < m_base.size(); ++i){
    			m_use_beta.push_back(false); // default is binary
    		}
    }
    else if(m_use_beta.size() != m_base.size()){
		error_message.append("ERROR: Number of beta doesn't match number of base file!\n");
		error_message.append("       Default value only work when all base file are using OR and\n");
		error_message.append("       when --beta is not used\n");
    }
    if(m_no_regress && !m_fastscore){
    		fprintf(stderr, "WARNING: To limit the amount of output,\n");
    		fprintf(stderr, "         no-regress can only be used with\n");
    		fprintf(stderr, "         fastscore. Will use fastscore\n");
    		m_fastscore=true;
    }
    if(error) throw std::runtime_error(error_message);
    return true;
}


Commander::Commander()
{
    //ctor
}

Commander::~Commander()
{
    //dtor
}

void Commander::usage(){
    fprintf(stderr, "Usage: PRSice [Options] \n\n");
    fprintf(stderr, "Required Inputs:\n");
    fprintf(stderr, "         -b | --base         A list of base file. Should be comma separated\n");
    fprintf(stderr, "         -t | --target       A list of target file. Should be comma separated\n");
    fprintf(stderr, "                             Currently only support the PLINK binary input\n");
    fprintf(stderr, "         --binary_target     Indication of whether if the target is binary\n");
    fprintf(stderr, "                             or not. default is T. Should be of the same \n");
    fprintf(stderr, "                             length as target \n");
    fprintf(stderr, "         --beta              Indication of whether if the test statistic\n");
    fprintf(stderr, "                             is beta instead of OR. default is F. Should\n");
    fprintf(stderr, "                             be of the same length as base \n");
    fprintf(stderr, "\nOptions\n");
    fprintf(stderr, "         -f | --pheno_file   The phenotype file containing the target\n");
    fprintf(stderr, "                             phenotypes. If provided, the fam file of the\n");
    fprintf(stderr, "                             target is ignored\n");
    fprintf(stderr, "         -L | --ld           PLINK binary input for LD calculation. If not\n");
    fprintf(stderr, "                             provided, will use the target genotypes for the \n");
    fprintf(stderr, "                             calculation of the LD during clumping\n");
    fprintf(stderr, "         -c | --covar_header Header of covariates, if not provided, will \n");
    fprintf(stderr, "                             use all variable in the covariate file as the\n");
    fprintf(stderr, "                             covariate\n");
    fprintf(stderr, "         -C | --covar_file   Covariate file. Format should be:\n");
    fprintf(stderr, "                             ID Cov1 Cov2\n");
    fprintf(stderr, "                             Must contain a header\n");
    fprintf(stderr, "         -a | --ancestry     NOT DEVELOPED YET\n");
    fprintf(stderr, "         -o | --out          The prefix of all output. Default %s\n", m_out.c_str());
    fprintf(stderr, "\nScoring options:\n");
    fprintf(stderr, "         -l | --lower        The starting p-value threshold. default: %f\n", m_lower);
    fprintf(stderr, "         -u | --upper        The final p-value threshold. default: %f\n", m_upper);
    fprintf(stderr, "         -i | --interval     The step size of the threshold. default: %f\n", m_inter);
    fprintf(stderr, "\nFile Headers:\n");
    fprintf(stderr, "              --chr          Column header of Chromosome <Required>\n");
    fprintf(stderr, "              --A1           Column header of Reference Allele <Required>\n");
    fprintf(stderr, "              --A2           Column header of Alternaative Allele \n");
    fprintf(stderr, "              --stat         Column header of test statistic, either BETA or OR\n");
    fprintf(stderr, "              --snp          Column header of SNP id\n");
    fprintf(stderr, "              --bp           Column header of SNP location\n");
    fprintf(stderr, "              --se           Column header of Standard Error\n");
    fprintf(stderr, "              --pvalue       Column head of p-value <Required>\n");
    fprintf(stderr, "              --index        If the base file doesn't contain a header, you can\n");
    fprintf(stderr, "                             use this option, which essentially state that all \n");
    fprintf(stderr, "                             the above options are providing the INDEX of the\n");
    fprintf(stderr, "                             corresponding column. (Index should be 0-based)\n");
    fprintf(stderr, "\nClumping:\n");
    fprintf(stderr, "              --clump-p      The p-value threshold use for clumping. \n");
    fprintf(stderr, "                             Default is %f \n", m_clump);
//    fprintf(stderr, "         --clump_p2          \n");
    fprintf(stderr, "              --clump_r2     The R2 threshold for clumping.\n");
    fprintf(stderr, "                             Please note that as we did not implement\n");
    fprintf(stderr, "                             the maximum likelihood R2 calculation, the\n");
    fprintf(stderr, "                             clumping result can differ slightly from plink\n");
    fprintf(stderr, "              --clump_kb     The distance for clumping in kb\n");
    fprintf(stderr, "\nSelections:\n");
    fprintf(stderr, "         -B | --bed          Bed file containing the selected regions. \n");
    fprintf(stderr, "                             Name of bed file will be used as the region\n");
    fprintf(stderr, "                             identifier \n");
    fprintf(stderr, "         -g | --gtf          GTF file containing gene boundaries. Required\n");
    fprintf(stderr, "                             when --msigdb is set \n");
    fprintf(stderr, "         -m | --msigdb       MSIGDB file containing the pathway information\n");
    fprintf(stderr, "                             require the gtf file \n");
    fprintf(stderr, "              --gen_bed      Generate bed file of gene regions from \n");
    fprintf(stderr, "                             the gtf file \n");
    fprintf(stderr, "                             Default: false \n");
    fprintf(stderr, "              --proxy        Perform clumping with proxy SNPs. Any index SNP\n");
    fprintf(stderr, "                             are considered as part of the region of the\n");
    fprintf(stderr, "                             clumped SNPs\n");
    fprintf(stderr, "                             Default: false \n");
    fprintf(stderr, "\nMisc:\n");
    fprintf(stderr, "         -T | --thread       Number of thread used\n");
    fprintf(stderr, "         -h | --help         Display this help message\n");


}

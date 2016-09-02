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
        {"clump",required_argument,NULL,0},
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
                else if(command.compare("clump")==0){
                    double temp = atof(optarg);
                    if(temp < 0.0 || temp > 1.0){
                        error = true;
                        error_message.append("Clumping p-values must be >=0 and <= 1.0\n");
                    }
                    else m_clump = temp;
                }
//                else if(command.compare("clump_p2")==0){
//                    double temp = atof(optarg);
//                    if(temp < 0.0 || temp > 1.0){
//                        error = true;
//                        error_message.append("Clumping p-values must be >=0 and <= 1.0\n");
//                    }
//                    else m_clump_p2 = temp;
//                }
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
                if(m_covariate_files.size() != 0){
                    error = true;
                    error_message.append("Covariate header and covariate file option is mutually exclusive\n");
                    error_message.append("Please only use one of the options (You can submit multiple covariate files\n");
                }
                break;
            case 'C':
                m_covariate_files = misc::split(optarg, ", ");
                if(m_covariates.size() != 0){
                    error = true;
                    error_message.append("Covariate header and covariate file option is mutually exclusive\n");
                    error_message.append("Please only use one of the options (You can submit multiple covariate files\n");
                }
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
    fprintf(stderr, "         -b | --base         \n");
    fprintf(stderr, "         -t | --target       \n");
    fprintf(stderr, "         --binary_target     Indication of whether if the target is binary\n");
    fprintf(stderr, "                             or not. default is T. Should be of the same \n");
    fprintf(stderr, "                             length as target \n");
    fprintf(stderr, "         --beta              Indication of whether if the test statistic\n");
    fprintf(stderr, "                             is beta instead of OR. default is F. Should\n");
    fprintf(stderr, "                             be of the same length as base \n");
    fprintf(stderr, "\nOptions\n");
    fprintf(stderr, "         -f | --pheno_file   \n");
    fprintf(stderr, "         -L | --ld           \n");
    fprintf(stderr, "         -c | --covar_header \n");
    fprintf(stderr, "         -C | --covar_file   \n");
    fprintf(stderr, "         -a | --ancestry     \n");
    fprintf(stderr, "         -o | --out          The prefix of all output. Default %s\n", m_out.c_str());
    fprintf(stderr, "\nScoring options:\n");
    fprintf(stderr, "         -l | --lower        \n");
    fprintf(stderr, "         -u | --upper        \n");
    fprintf(stderr, "         -i | --interval     \n");
    fprintf(stderr, "\nFile Headers:\n");
    fprintf(stderr, "              --chr          \n");
    fprintf(stderr, "              --A1           \n");
    fprintf(stderr, "              --A2           \n");
    fprintf(stderr, "              --stat         \n");
    fprintf(stderr, "              --snp          \n");
    fprintf(stderr, "              --bp           \n");
    fprintf(stderr, "              --se           \n");
    fprintf(stderr, "              --pvalue       \n");
    fprintf(stderr, "              --index        \n");
    fprintf(stderr, "\nClumping:\n");
    fprintf(stderr, "              --clump        \n");
//    fprintf(stderr, "         --clump_p2          \n");
    fprintf(stderr, "              --clump_r2     \n");
    fprintf(stderr, "              --clump_kb     \n");
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
    fprintf(stderr, "\nMisc:\n");
    fprintf(stderr, "         -T | --thread       \n");
    fprintf(stderr, "         -h | --help         \n");


}

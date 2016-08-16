#include "commander.h"

Commander::Commander(const int argc, const char *argv[]){
    if(argc<=1){
        usage();
        throw "Please provide the required parameters";
    }
    static const char *optString = "b:t:c:C:a:f:L:p:T:u:l:i:h?";
    static const struct option longOpts[]={
        {"base",required_argument,NULL,'b'},
        {"target",required_argument,NULL,'t'},
        {"covar_header",required_argument,NULL,'c'},
        {"covar_file",required_argument,NULL,'C'},
        {"ancestry",required_argument,NULL,'a'},
        {"pheno_file",required_argument,NULL,'f'},
        {"ld",required_argument,NULL,'L'},
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
        {"clump_p1",required_argument,NULL,0},
        {"clump_p2",required_argument,NULL,0},
        {"clump_r2",required_argument,NULL,0},
        {"clump_kb",required_argument,NULL,0},
        {"thread",required_argument,NULL,'T'},
		{"help",no_argument,NULL,'h'},
		{NULL, 0, 0, 0}
	};
    
    bool error = false;
    int longIndex=0;
    int opt = 0;
    opt=getopt_long(argc, argv, optString, longOpts, &longIndex);
    std::vector<char> type;
    std::vector<std::string> columnName;
    std::map<std::string, char> duplication;
    std::string error_message = "";
    //Start reading all the parameters and perform the qc at the same time
    while(opt!=-1){
        switch(opt){
            case 0:
                std::string command = longOpts[longIndex].name;
                if(command.compare("chr")==0) m_chr=optarg;
                else if(command.compare("A1")==0) m_ref_allele = optarg;
                else if(command.compare("A2")==0) m_alt_allele=optarg;
                else if(command.compare("stat")==0) m_statistic = optarg;
                else if(command.compare("snp")==0) m_snp=optarg;
                else if(command.compare("bp")==0) m_bp=optarg;
                else if(command.compare("se")==0) m_standard_error = optarg;
                else if(command.compare("clump_p1")==0){
                    double temp = atof(optarg);
                    if(temp < 0.0 || temp > 1.0){
                        error = true;
                        error_message.append("Clumping p-values must be >=0 and <= 1.0\n");
                    }
                    else m_clump_p1 = temp;
                }
                else if(command.compare("clump_p2")==0){
                    double temp = atof(optarg);
                    if(temp < 0.0 || temp > 1.0){
                        error = true;
                        error_message.append("Clumping p-values must be >=0 and <= 1.0\n");
                    }
                    else m_clump_p2 = temp;
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
                break;
            case 'b':
                m_base= usefulTools::split(optarg, ", ");
                if(m_base.size() ==0){
                    error = true;
                    error_message.append("You must provide at least one valid base file name\n");
                }
                break;
            case 't':
                m_target = usefulTools::split(optarg, ", ");
                if(m_target.size() ==0){
                    error = true;
                    error_message.append("You must provide at least one valid target file name\n");
                }
                break;
            case 'c':
                m_covariates = usefulTools::split(optarg, ", ");
                if(m_covariate_files.size() != 0){
                    error = true;
                    error_message.append("Covariate header and covariate file option is mutually exclusive\n");
                    error_message.append("Please only use one of the options (You can submit multiple covariate files\n");
                }
                break;
            case 'C':
                m_covariate_files = usefulTools::split(optarg, ", ");
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
            case 'h':
            case '?':
                
                return;
                break;
            default:
                throw "Undefined operator, please use --help for more information!";
        }
        opt=getopt_long(argc, argv, optString, longOpts, &longIndex);
    }
    
    
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
}

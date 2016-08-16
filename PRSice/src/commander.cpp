#include "commander.h"

Commander::Commander(const int argc, const char *argv[]){
    if(argc<=1){
        usage();
        throw "Please provide the required parameters";
    }
    static const char *optString = "b:t:p:mh?";
    static const struct option longOpts[]={
	    //Qt specific parameter
        {"base",required_argument,NULL,'b'},
        {"target",required_argument,NULL,'t'},
        {"covar_header",required_argument,NULL,'c'},
        {"covar_file",required_argument,NULL,'C'},
        {"ancestry",required_argument,NULL,'a'},
        {"pheno_file",required_argument,NULL,'f'},
        {"chr",required_argument,NULL,0},
        {"A1",required_argument,NULL,0},
        {"A2",required_argument,NULL,0},
        {"stat",required_argument,NULL,0},
        {"snp",required_argument,NULL,0},
        {"bp",required_argument,NULL,0},
        {"se",required_argument,NULL,0},
        {"p",required_argument,NULL,'p'},
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
            case 'h':
            case '?':
                if(mode.compare("quant")==0) qtUsage();
                else if(mode.compare("binary")==0) btUsage();
                return false; //This is not an error, just tell the programme to end
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

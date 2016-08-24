#include "../inc/prsice.hpp"

void PRSice::process(const std::string &c_input, const Commander &c_commander){
	// As we don't use the ptr_vector, we need to read the SNP here
	// otherwise the SNP might go out of scope
	std::vector<SNP> snp_list;
    std::map<std::string, size_t> snp_index;
    std::ifstream snp_file;
    snp_file.open(c_input.c_str());
    if(!snp_file.is_open()){
    		std::string error_message = "Cannot open base file: "+c_input;
    		throw std::runtime_error(error_message);
    }
    std::string line;
    std::vector<int> index = SNP::get_index(c_commander, c_input);
    int max_index = index.back();
    size_t num_duplicated = 0;
    size_t num_stat_not_convertible = 0;
	size_t num_p_not_convertible = 0;
	bool se_error = false;
    while(std::getline(snp_file, line)){
    		misc::trim(line);
    		if(!line.empty()){
    			std::vector<std::string> token = misc::split(line);
    		   	if(token.size() <= max_index) throw std::runtime_error("More index than column in data");
    		   	if(snp_index.find(token[index[4]])==snp_index.end()) num_duplicated++;
    		   	else{
    		   		std::string rs_id = token[index[4]];
    		   		std::string chr = "";
    		   		if(index[0] >= 0) chr = token[index[0]];
    		   		std::string ref_allele = "";
    		   		if(index[1] >= 0) ref_allele = token[index[1]];
    		   		std::string alt_allele = "";
    		   		if(index[2] >= 0) alt_allele = token[index[2]];
    		   		double stat = 0.0;
    		   		if(index[3] >= 0){
    		   			//Check if it is double
    		   			try{
    		   				stat = misc::convert<double>(token[index[3]]);
    		   			}
    		   			catch(const std::runtime_error& error){ //we know only runtime error is throw
    		   				num_stat_not_convertible++;
    		        		}
    		   		}
    		   		double se = 0.0;
    		    		if(index[6] >= 0){
    		    			try{
    		    				se = misc::convert<double>(token[index[6]]);
    		    			}
    		    			catch(const std::runtime_error &error){
    		    				se_error = true;
    		    			}
    		    		}
    		    		double pvalue = 0.0;
    		    		if(index[7] >= 0){
    		    			try{
    		    				pvalue = misc::convert<double>(token[index[7]]);
    		    			}
    		    			catch(const std::runtime_error &error){
    		    				num_p_not_convertible++;
    		    			}
    		    		}
    		    		size_t loc = 0;
    		    		if(index[5]>=0){
    		    			int temp = atoi(token[index[5]].c_str());
    		    			if(temp <0) fprintf(stderr, "ERROR: %s has negative loci\n", rs_id);
    		    			else loc = temp;
    		    		}
    		      	snp_list.push_back(SNP(rs_id, chr, loc, ref_allele, alt_allele, stat, se, pvalue ));
    		   	}
    		}
    }
    snp_file.close();
    fprintf(stderr, "Number of duplicated SNPs : %zu\n", num_duplicated);
    fprintf(stderr, "Number of SNPs included   : %zu\n", snp_list.size());
	if(num_stat_not_convertible!=0) fprintf(stderr, "Failed to convert %zu OR/beta\n", num_stat_not_convertible);
	if(num_p_not_convertible!=0) fprintf(stderr, "Failed to convert %zu p-value\n", num_p_not_convertible);

	// Next, start generating the selection matrix
	// This should be the matrix deciding whether if the SNP will be included in the
	// specific threshold
	// When no inclusion criteria is given, we will generate a vector of 1, indicating
	// that we will include all SNPs in the analysis
	// This should be similar to the concept in SHREK


	// Maybe perform clumping?
	std::vector<std::string> target = c_commander.get_target();
	// We need to do this for each target if the LD is not provided
	if(!c_commander.ld_prefix().empty()){
		PLINK geno = PLINK(c_commander.ld_prefix());
		// Do clumping and at the same time generate the results?
	}
	else{
		for(size_t i_target = 0; i_target < target.size(); ++i_target){
			PLINK geno = PLINK(target[i_target]);
		}
	}

}

void PRSice::run(const Commander &c_commander){
    std::vector<std::string> base = c_commander.get_base();
    int num_base = base.size();
    if(num_base == 0) throw std::runtime_error("There is no base case to run");
    if(num_base < 0) throw std::runtime_error("Negative number of base");
    else{
        for(size_t i = 0; i < num_base; ++i){
            //We process each base case independently
            process(base[i], c_commander);
        }
    }
}

PRSice::PRSice()
{
    //ctor
}

PRSice::~PRSice()
{
    //dtor
}

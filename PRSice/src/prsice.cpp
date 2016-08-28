#include "prsice.hpp"

void PRSice::process(const std::string &c_input, const Commander &c_commander, Region &region){
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
	// only A1 is required
	// so ambiguous check will be when reading the LD file
	// same for flipping
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
    		    			if(temp <0) fprintf(stderr, "ERROR: %s has negative loci\n", rs_id.c_str());
    		    			else loc = temp;
    		    		}
#if defined(__LP64__) || defined(_WIND64)
    		    		uint64_t* flag = region.check(chr, loc);
#endif
    		      	snp_list.push_back(SNP(rs_id, chr, loc, ref_allele, alt_allele, stat, se, pvalue, region.check(chr, loc), region.size()));
    		   	}
    		}
    }
    snp_file.close();
    fprintf(stderr, "Number of duplicated SNPs : %zu\n", num_duplicated);
    fprintf(stderr, "Number of SNPs included   : %zu\n", snp_list.size());
	if(num_stat_not_convertible!=0) fprintf(stderr, "Failed to convert %zu OR/beta\n", num_stat_not_convertible);
	if(num_p_not_convertible!=0) fprintf(stderr, "Failed to convert %zu p-value\n", num_p_not_convertible);

    // Will do the selection on the fly?
    // Note: PLINK Clumping will discard any SNPs that doesn't pass clump-p1
	std::vector<size_t> p_sort_order = SNP::sort_by_p(snp_list);
    // Read target file first, only include SNPs that are also in the target
    // can also perform the ambiguous SNP removal at this point
    // Therefore, anything happened from this point onward should be target
    // specific
    std::vector<std::string> target = c_commander.get_target();
    for(size_t i_target = 0; i_target < target.size(); ++i_target){
        std::map<std::string, bool> inclusion;
        std::string target_bim_name = target[i]+".bim";
        get_inclusion(inclusion, target_bim_name, snp_list, snp_index);
        // Then read in the LD file, that can either be the target file or an
        // external reference
        // This should perform the clumping, which will produce a list of SNPs
        // that are supposedly included in the final PRS
        if(c_commander.ld_prefix()..empty()){
            // we will perform clumping using the target file
            // Clumping will update the m_clump_target of the SNP class
            // And should update the inclusion index we have
            // The region flag should also be updated such that
            // the clump index SNP should represent the region of all the
            // clumped SNPs
        }
        else{
            // we will perform clumping using the LD file
        }
        // So technically, from here, we just need to perform the PRS with
        // the inclusion map


    }

}

void PRSice::get_inclusion(std::map<std::string, bool> &inclusion, const std::string &target_bim_name,
                       std::vector<SNP> &snp_list, const std::map<std::string, size_t> &snp_index){
    std::ifstream target_file;
    target_file.open(target_bim_name.c_str());
    if(!target_file.is_open()){
        std::string error_message = "Cannot open target bim file: "+target_bim_name;
        throw std::runtime_error(error_message);
    }
    size_t num_ambig=0;
    while(std::getline(target_file, line)){
        misc::trim(line);
        if(!line.empty()){
            std::vector<std::string> token = misc::split(line);
            if(token.size() < 6) throw std::runtime_error("Malformed bim file. Should contain at least 6 column");
            std::string chr = token[0];
            std::string rsid = token[1];
            size_t loc = 0;
            int temp = 0;
            try{
                temp =misc::convert<int>(token[3]);
                if(temp < 0){
                    std::string error_message = "Negative coordinate of SNP in "+target_bim_name;
                    throw std::runtime_error(error_message);
                }
                loc = temp;
            }
            catch(std::runtime_error &error){
                std::string error_message = "Non-numeric coordinate of SNP in "+target_bim_name;
                throw std::runtime_error(error_message);
            }
            std::string ref_allele = token[4];
            std::string alt_allele = token[5];
            if(snp_index.find(rsid)==snp_index.end()){
                // No summary statistic for this SNP, so we will ignore it
                inclusion[rsid] = false;
            }
            else{
                // will do some soft checking, will issue warning if there are any problem
                // first check if ambiguous
                if( (ref_allele.compare("A") && alt_allele.compare("T")) ||
                    (ref_allele.compare("a") && alt_allele.compare("t")) ||
                    (ref_allele.compare("T") && alt_allele.compare("A")) ||
                    (ref_allele.compare("t") && alt_allele.compare("a")) ||
                    (ref_allele.compare("G") && alt_allele.compare("C")) ||
                    (ref_allele.compare("g") && alt_allele.compare("c")) ||
                    (ref_allele.compare("C") && alt_allele.compare("G")) ||
                    (ref_allele.compare("c") && alt_allele.compare("g")))
                {
                    num_ambig++;
                    inclusion[rsid] = false;
                }
                else{
                    // not ambiguous, now do soft checking
                    size_t index = snp_index[rsid];
                    bool same = snp_list[index].check_loc(chr, loc, ref_allele, alt_allele);
                    if(!same){
                        fprintf(stderr, "WARNING: %s differ between target and base file\n", rsid.c_str());
                        fprintf(stderr, "         It is advised that you check the files are \n");
                        fprintf(stderr, "         From the same genome build\n");
                    }
                    inclusion[rsid] = true;
                }
            }
        }
    }
    target_file.close();
}

// This will update the score for each individual
void PRSice::score(const std::map<std::string, bool> &inclusion, const std::string target_name,
                   const std::map<std::string, size_t> &snp_index, std::vector<SNP> &snp_list){

}

void PRSice::run(const Commander &c_commander, Region &region){
    std::vector<std::string> base = c_commander.get_base();
    int num_base = base.size();
    if(num_base == 0) throw std::runtime_error("There is no base case to run");
    if(num_base < 0) throw std::runtime_error("Negative number of base");
    else{
        for(size_t i = 0; i < num_base; ++i){
            //We process each base case independently
        		region.reset();
            process(base[i], c_commander, region);
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

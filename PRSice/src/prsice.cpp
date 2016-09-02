#include "prsice.hpp"

void PRSice::process(const std::string &c_input, bool beta, const Commander &c_commander, Region &region){
	// As we don't use the ptr_vector, we need to read the SNP here
	// otherwise the SNP might go out of scope

    std::vector<int> index = SNP::get_index(c_commander, c_input);
	std::vector<SNP> snp_list;
    std::map<std::string, size_t> snp_index;
    std::ifstream snp_file;
    snp_file.open(c_input.c_str());
    if(!snp_file.is_open()){
    		std::string error_message = "Cannot open base file: "+c_input;
    		throw std::runtime_error(error_message);
    }
    std::string line;
    int max_index = index.back();
    size_t num_duplicated = 0;
    size_t num_stat_not_convertible = 0;
	size_t num_p_not_convertible = 0;
	// only A1 is required
	// so ambiguous check will be when reading the LD file
	// same for flipping
	bool se_error = false;
	if(!c_commander.index()) std::getline(snp_file, line);
    while(std::getline(snp_file, line)){
    		misc::trim(line);
    		if(!line.empty()){
    			std::vector<std::string> token = misc::split(line);
    		   	if(token.size() <= max_index) throw std::runtime_error("More index than column in data");
    		   	if(snp_index.find(token[index[4]])!=snp_index.end()) num_duplicated++;
    		   	else{
    		   		std::string rs_id = token[index[4]];
    		   		std::string chr = "";
    		   		if(index[0] >= 0) chr = token[index[0]];
    		   		std::string ref_allele = "";
    		   		if(index[1] >= 0) ref_allele = token[index[1]];
    		   		std::string alt_allele = "";
    		   		if(index[2] >= 0) alt_allele = token[index[2]];
		    		double pvalue = 0.0;
		    		if(index[7] >= 0){
		    			try{
		    				pvalue = misc::convert<double>(token[index[7]]);
		    			}
		    			catch(const std::runtime_error &error){
		    				num_p_not_convertible++;
		    			}
		    		}
    		   		double stat = 0.0;
    		   		if(index[3] >= 0){
    		   			//Check if it is double
    		   			try{
    		   				stat = misc::convert<double>(token[index[3]]);
    		   				if(!beta) stat = log(stat);
    		   			}
    		   			catch(const std::runtime_error& error){ //we know only runtime error is throw
    		   				num_stat_not_convertible++;
    		        		}
    		   		}else stat = pvalue;
    		   		double se = 0.0;
    		    		if(index[6] >= 0){
    		    			try{
    		    				se = misc::convert<double>(token[index[6]]);
    		    			}
    		    			catch(const std::runtime_error &error){
    		    				se_error = true;
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
    		    		snp_index[rs_id] =snp_list.size();
    		      	snp_list.push_back(SNP(rs_id, chr, loc, ref_allele, alt_allele, stat, se, pvalue, region.check(chr, loc), region.size()));
    		   	}
    		}
    }
    snp_file.close();
    fprintf(stderr, "Number of duplicated SNPs : %zu\n", num_duplicated);
    fprintf(stderr, "Number of SNPs from base  : %zu\n", snp_list.size());
	if(num_stat_not_convertible!=0) fprintf(stderr, "Failed to convert %zu OR/beta\n", num_stat_not_convertible);
	if(num_p_not_convertible!=0) fprintf(stderr, "Failed to convert %zu p-value\n", num_p_not_convertible);

    // Read target file first, only include SNPs that are also in the target
    // can also perform the ambiguous SNP removal at this point
    // Therefore, anything happened from this point onward should be target
    // specific
    std::vector<std::string> target = c_commander.get_target();
	fprintf(stderr,"\nClumping Parameters: \n");
    fprintf(stderr,"==============================\n");
    fprintf(stderr,"P-Threshold  : %f\n", c_commander.get_clump_p());
    fprintf(stderr,"R2-Threshold : %f\n", c_commander.get_clump_r2());
    fprintf(stderr,"Window Size  : %zu\n", c_commander.get_clump_kb());

    for(size_t i_target = 0; i_target < target.size(); ++i_target){
    		fprintf(stderr,"\nStart processing: %s\n", target[i_target].c_str());
    	    fprintf(stderr,"==============================\n");
        std::map<std::string, size_t> inclusion;
        std::string target_bim_name = target[i_target]+".bim";
//        get_inclusion(inclusion, target_bim_name, snp_list, snp_index);
        // Then read in the LD file, that can either be the target file or an
        // external reference
        // This should perform the clumping, which will produce a list of SNPs
        // that are supposedly included in the final PRS
        bool has_ld = !c_commander.ld_prefix().empty();
        std::string ld_file = (has_ld)? c_commander.ld_prefix(): target[i_target];
        // Clumping will update the m_clump_target of the SNP class
        // And should update the inclusion index we have
        // The region flag should also be updated such that
        // the clump index SNP should represent the region of all the
        // clumped SNPs
        PLINK clump(ld_file, c_commander.get_thread());
        if(has_ld) fprintf(stderr,"\nIn LD Reference %s\n", ld_file.c_str());
        else fprintf(stderr,"\nStart performing clumping\n");
        clump.initialize(inclusion, snp_list, snp_index);
        if(has_ld){
        		// perform additional filtering
            fprintf(stderr,"\nIn target %s\n", ld_file.c_str());
        		update_inclusion(inclusion, target_bim_name, snp_list, snp_index);
        }
        clump.clumping(inclusion, snp_list, snp_index, c_commander.get_clump_p(), c_commander.get_clump_r2(), c_commander.get_clump_kb());
        // From here, inclusion include the SNPs that we want to use for PRS

        double bound_start = c_commander.get_lower();
        double bound_end = c_commander.get_upper();
        double bound_inter = c_commander.get_inter();
        std::vector<double> prs_score;
        // Each time, only read SNPs under the boundary and perform the analysis
        double current_upper = bound_start;
        double current_lower = 0.0;
        size_t num_snp = 0; // The number of SNPs included so far
        for(;current_upper<bound_end; current_upper+=bound_inter){
            fprintf(stderr,"\rCalculating cutoff %f", current_upper);
        		// will add anything from (lower, upper]
        		bool reg = score(inclusion, snp_list, target[i_target], prs_score, current_lower, current_upper, num_snp);
        		current_lower = current_upper;
        		// This should update the score
        		//TODO: PRSice can also calculate the PCA / MDS and use as covariate in its analysis
        		if(reg){
					std::ofstream test;
					std::string test_name = c_commander.get_out()+std::to_string(current_upper)+".debug";
					// TODO: The rest should be here
					test.open(test_name.c_str());
					for(size_t i = 0; i < prs_score.size(); ++i){
						test << prs_score[i] << std::endl;
					}
					test.close();
        		}
        }
        if(current_upper != bound_end){
            fprintf(stderr,"\rCalculating cutoff %f", bound_end);
        		bool reg = score(inclusion, snp_list, target[i_target], prs_score, current_lower, bound_end, num_snp);
        }
        fprintf(stderr, "\n");
    }
}

// basically update the score vector to contain the new polygenic score
bool PRSice::score(const std::map<std::string, size_t> &inclusion, std::vector<SNP> &snp_list,
		const std::string &target, std::vector<double> &prs_score,
		double threshold_lower, double threshold_upper, size_t &num_snp){
	// Here we will make a different inclusion for the inclusion
	if(threshold_lower==threshold_upper) return false; //nothing to do here
	std::map<std::string, size_t>::const_iterator i_inclusion;
	std::map<std::string, size_t> cur_inclusion;
	for(i_inclusion = inclusion.begin(); i_inclusion!=inclusion.end(); ++i_inclusion){
		double cur_p = snp_list[i_inclusion->second].get_p_value();
		if(cur_p < threshold_upper && cur_p >= threshold_lower){
			cur_inclusion[i_inclusion->first]=i_inclusion->second;
			num_snp++;
		}
	}
	if(cur_inclusion.size()==0) return false;
	PLINK prs(target);
	prs.initialize();
	prs.get_score(cur_inclusion, snp_list, prs_score);
	return true;
}

void PRSice::update_inclusion(std::map<std::string, size_t> &inclusion, const std::string &target_bim_name,
                       std::vector<SNP> &snp_list, const std::map<std::string, size_t> &snp_index){
    std::ifstream target_file;
    target_file.open(target_bim_name.c_str());
    if(!target_file.is_open()){
        std::string error_message = "Cannot open target bim file: "+target_bim_name;
        throw std::runtime_error(error_message);
    }
    size_t num_ambig=0, not_found=0;
    std::string line;
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
            if(inclusion.find(rsid) != inclusion.end() && snp_index.find(rsid)!=snp_index.end()){
                // will do some soft checking, will issue warning if there are any problem
                // first check if ambiguous
                if( (ref_allele.compare("A")==0 && alt_allele.compare("T")==0) ||
                    (ref_allele.compare("a")==0 && alt_allele.compare("t")==0) ||
                    (ref_allele.compare("T")==0 && alt_allele.compare("A")==0) ||
                    (ref_allele.compare("t")==0 && alt_allele.compare("a")==0) ||
                    (ref_allele.compare("G")==0 && alt_allele.compare("C")==0) ||
                    (ref_allele.compare("g")==0 && alt_allele.compare("c")==0) ||
                    (ref_allele.compare("C")==0 && alt_allele.compare("G")==0) ||
                    (ref_allele.compare("c")==0 && alt_allele.compare("g")==0))
                {
                    num_ambig++;
                    if(inclusion.find(rsid)!=inclusion.end()) inclusion.erase(rsid);
                }
                else{
                    // not ambiguous, now do soft checking
                    size_t index = snp_index.at(rsid);
                    bool same = snp_list[index].check_loc(chr, loc, ref_allele, alt_allele);
                    if(!same){
                        fprintf(stderr, "WARNING: %s differ between target and base file\n", rsid.c_str());
                        fprintf(stderr, "         It is advised that you check the files are \n");
                        fprintf(stderr, "         From the same genome build\n");
                    }
                    inclusion[rsid] = snp_index.at(rsid);
                }
            }
            else {
            		not_found++;
              	if(inclusion.find(rsid)!=inclusion.end()) inclusion.erase(rsid);
            }
        }
    }
    target_file.close();
    if(num_ambig != 0)	fprintf(stderr, "Number of ambiguous SNPs  : %zu\n", num_ambig);
    if(not_found != 0)	fprintf(stderr, "Number of SNPs not found  : %zu\n", not_found);
    fprintf(stderr, "Final number of SNPs      : %zu\n", inclusion.size());
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
            process(base[i], c_commander.get_base_binary(i), c_commander, region);
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

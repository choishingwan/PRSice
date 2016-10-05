#include "prsice.hpp"

std::mutex PRSice::score_mutex;

// Should be alright without problem
void PRSice::run(const Commander &c_commander, Region &region)
{
	// First, wrap the whole process with a for loop
	// to loop through all base phenotype
    std::vector<std::string> base = c_commander.get_base();
    int num_base = base.size();
    if(num_base == 0) throw std::runtime_error("There is no base case to run");
    if(num_base < 0) throw std::runtime_error("Negative number of base");
    else{
    		if(num_base > 1) fprintf(stderr, "Multiple base phenotype detected. You might want to run separate instance of PRSice to speed up the process\n");
        for(size_t i = 0; i < num_base; ++i){
    			fprintf(stderr,"\nStart processing: %s\n", base[i].c_str());
    			fprintf(stderr,"==============================\n");
        		m_current_base=base[i];
            process(base[i], c_commander.get_base_binary(i), c_commander, region);
            fprintf(stderr, "\n");
        }
    }
}

// Seems to be working alright (hope I know how to write test cases...)
void PRSice::process(const std::string &c_input, bool beta, const Commander &c_commander, Region &region)
{
	// we need to reset the region flag for region inclusion information
	region.reset();
	// snp_list should contain all SNP information from the base file
	boost::ptr_vector<SNP> snp_list;
	// and snp_index is used for quickly finding the index of specific SNP in snp_list
    std::map<std::string, size_t> snp_index;
    // Now obtaining the SNP information from the base file
    get_snp(snp_list, snp_index, c_input, beta, c_commander, region);
    // Then we will perform the rest of the process for each individual
    // target file.
    std::vector<std::string> target = c_commander.get_target();
	fprintf(stderr,"\nClumping Parameters: \n");
    fprintf(stderr,"==============================\n");
    fprintf(stderr,"P-Threshold  : %f\n", c_commander.get_clump_p());
    fprintf(stderr,"R2-Threshold : %f\n", c_commander.get_clump_r2());
    fprintf(stderr,"Window Size  : %zu\n", c_commander.get_clump_kb());

    for(size_t i_target = 0; i_target < target.size(); ++i_target){
    		fprintf(stderr,"\nStart processing: %s\n", target[i_target].c_str());
    	    fprintf(stderr,"==============================\n");
    	    // The inclusion map should contain the SNPs that are supposed to be used for the
    	    // calculation of PRS
        std::map<std::string, size_t> inclusion;
        bool has_ld = !c_commander.ld_prefix().empty();
        std::string ld_file = (has_ld)? c_commander.ld_prefix(): target[i_target];
        // create the plink class for clumping
        PLINK clump(ld_file, c_commander.get_thread());
        if(has_ld) fprintf(stderr,"\nIn LD Reference %s\n", ld_file.c_str());
        else fprintf(stderr,"\nStart performing clumping\n");
        // now initialize the plink class. This should also update inclusion such
        // that it will only include SNPs also found in the plink file
        clump.initialize(inclusion, snp_list, snp_index);
        if(has_ld){
        		// because we have independent LD file, we want to make sure
        		// only SNPs that are also found in the target file are used
        		// for clumping
        		std::string target_bim_name = target[i_target]+".bim";
            fprintf(stderr,"\nIn target %s\n", target_bim_name.c_str());
            // This will provide us the final inclusion map, which contains the
            // SNPs that are supposed to be used for clumping
        		update_inclusion(inclusion, target_bim_name, snp_list, snp_index);
        }
        // This will perform clumping. When completed, inclusion will contain the index SNPs
        // However, this should be changed in later version such that we can also
        // handle different regions
        clump.start_clumping(inclusion, snp_list, snp_index, c_commander.get_clump_p(),
        		c_commander.get_clump_r2(), c_commander.get_clump_kb(), c_commander.get_proxy());
        // Now begin the calculation of the PRS
        calculate_score(c_commander, c_commander.get_target_binary(i_target),
        		i_target, inclusion, snp_list, region);
    }
}

// Seems alright now
void PRSice::get_snp(boost::ptr_vector<SNP> &snp_list,
		std::map<std::string, size_t> &snp_index, const std::string &c_input, bool beta,
		const Commander &c_commander, Region &region)
{
	// First, we need to obtain the index of different columns from the base files
	// Might also want to include some warnings regarding OR or beta.
	// NOTE: -1 means missing and index is hard coded such that each index will represent
	//       specific header

	// just issue the warning. would terminate though
	if(beta && c_commander.statistic().compare("OR")==0) fprintf(stderr, "WARNING: OR detected but user suggest the input is beta!\n");
	std::vector<int> index = SNP::get_index(c_commander, c_input);
	// Open the file
    std::ifstream snp_file;
    snp_file.open(c_input.c_str());
    if(!snp_file.is_open()){
    		std::string error_message = "Cannot open base file: "+c_input;
    		throw std::runtime_error(error_message);
    }
    // Some QC counts
    size_t num_duplicated = 0;
    size_t num_stat_not_convertible = 0;
	size_t num_p_not_convertible = 0;
	size_t num_indel = 0;
    std::string line;
    int max_index = index.back();
    bool se_error = false;
	if(!c_commander.index()) std::getline(snp_file, line);
	bool read_error = false;
	// Actual reading the file, will do a bunch of QC
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
	    					if(pvalue < 0.0 || pvalue > 1.0){
	    						read_error =true;
	    						fprintf(stderr, "ERROR: %s's p-value is %f\n", rs_id.c_str(), pvalue);
	    					}
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
		   		}
		   		double se = 0.0;
		    		if(index[6] >= 0){
		    			try{
		    				se = misc::convert<double>(token[index[6]]);
		    			}
		    			catch(const std::runtime_error &error){
		    				se_error = true; // This does nothing
		    			}
		    		}
		    		size_t loc = 0;
		    		if(index[5]>=0){
		    			int temp = atoi(token[index[5]].c_str());
		    			if(temp <0){
		    				read_error=true;
		    				fprintf(stderr, "ERROR: %s has negative loci\n", rs_id.c_str());
		    			}
		    			else loc = temp;
		    		}
		    		if(ref_allele.compare("-")==0 || ref_allele.compare("I") == 0 || ref_allele.compare("D")==0 ||
		    				ref_allele.size()>1){
		    			num_indel++;
		    		}
		    		else if(!alt_allele.empty() &&
		    				(alt_allele.compare("-")==0 || alt_allele.compare("I") == 0 ||
		    						alt_allele.compare("D")==0 || alt_allele.size()>1)){
		    			num_indel++;
		    		}
		    		else{
		    			snp_list.push_back(new SNP(rs_id, chr, loc, ref_allele, alt_allele, stat, se, pvalue, new SNP::long_type[1], region.size()));
		    		}
		   	}
		}
		if(read_error) throw std::runtime_error("Please check if you have the correct input");
	}
	snp_file.close();
	// This may seems unnecessary and time consuming
	// but hopefully this should speed up the region discovery step
	// The concept is, by performing sorting, we can ensure that the
	// region and SNP is sorted in the same way. Thus as long as they
	// are from the same genome build (e.g. not hg19 vs b37), then
	// we should always be able to assume that the SNPs are in
	// the correct order
	snp_list.sort(SNP::sort_snp);
	// now write in the index
	for(size_t i_snp = 0; i_snp < snp_list.size(); ++i_snp){
		snp_index[snp_list[i_snp].get_rs_id()]=i_snp;
		snp_list[i_snp].set_flag(region.check(snp_list[i_snp].get_chr(), snp_list[i_snp].get_loc()));
	}

	// Now output the statistics. Might want to improve the outputs
	fprintf(stderr, "Number of SNPs from base  : %zu\n", snp_list.size());
	if(num_indel!=0) fprintf(stderr, "Number of Indels          : %zu\n", num_indel);
	if(num_duplicated!=0) fprintf(stderr, "Number of duplicated SNPs : %zu\n", num_duplicated);
	if(num_stat_not_convertible!=0) fprintf(stderr, "Failed to convert %zu OR/beta\n", num_stat_not_convertible);
	if(num_p_not_convertible!=0) fprintf(stderr, "Failed to convert %zu p-value\n", num_p_not_convertible);
}

void PRSice::calculate_score(const Commander &c_commander, bool target_binary,
		const size_t c_i_target, const std::map<std::string, size_t> &inclusion,
		const boost::ptr_vector<SNP> &snp_list, const Region &c_region)
{
	// Might want to add additional parameter for the region output
	// keep it for now
	// First, get the phenotype and covariate matrix
	// fam_index is used for storing the index of each individual on the matrix
	// it is used because there might be missing sample which we will exclude from
	// the matrix
	std::map<std::string, size_t> fam_index;
	// This is a vector of name, so that we can quickly copy it to the result
	// vectors and declare the correct size
    std::vector<std::pair<std::string, double> > prs_fam;
    // check whether if regression is required
    bool no_regress =c_commander.no_regression();
    // below is only required if regression is performed
    Eigen::VectorXd phenotype;
    std::string target = c_commander.get_target(c_i_target);
    std::string pheno_file = c_commander.get_pheno(c_i_target);
    Eigen::MatrixXd independent_variables;
    if(!no_regress){
        // This should generate the phenotype matrix by reading from the fam/pheno file
    		phenotype = gen_pheno_vec(	target, pheno_file,
    													target_binary, fam_index, prs_fam);
    		// This should generate the covariate matrix
    		independent_variables = gen_cov_matrix(	target, c_commander.get_cov_file(),
    													c_commander.get_cov_header(), fam_index);
    }else{
//  		We will need to initialize prs_fam
    		std::string fam_name = target+".fam";
    		std::ifstream fam;
    		fam.open(fam_name.c_str());
    		if(!fam.is_open()){
    			std::string error_message = "ERROR: Cannot open fam file: " + fam_name;
    			throw std::runtime_error(error_message);
    		}
    		std::string line;
    		while(std::getline(fam, line)){
    			misc::trim(line);
    			if(!line.empty()){
    				std::vector<std::string> token = misc::split(line);
    				if(token.size() < 6) throw std::runtime_error("Malformed fam file, should contain at least 6 columns");
    				prs_fam.push_back(std::pair<std::string, double>(token[1], 0.0));
    			}
    		}
    		fam.close();
    }

    // Declare the variables, might need to implement fastscore here
    bool fastscore = c_commander.fastscore();
	double bound_start =  (fastscore)? 0.001: c_commander.get_lower();
	double bound_end = (fastscore)? 0.5:c_commander.get_upper();
	double bound_inter = c_commander.get_inter();
	// we will use a "special" vector to speed up the bim reading
	// the concept is such that we can skip lines efficiently
    std::string bim_name = target+".bim";
    // This should contain some important information to speed up the bed file reading
    // should contain string, size_t, int, size_t
    // string = rsid
    // size_t = line number of bed file
    // int = p-value threshold category, -1 occurs when bound_start != 0
    // size_t = index of the SNP on snp_list, avoid repeat finding from map_inclusion
    std::vector<p_partition > quick_ref;
    std::ifstream bim;
    bim.open(bim_name.c_str());
    if(!bim.is_open()){
    		std::string error_message = "Cannot open bim file: " +bim_name;
    		throw std::runtime_error(error_message);
    }

    std::string line;
    size_t cur_line = 0;
    bool pre_run = false;
    while(getline(bim, line)){
    		misc::trim(line);
    		if(!line.empty()){
    			std::vector<std::string> token = misc::split(line);
    			if(token.size() < 6) throw std::runtime_error("Malformed bim file, should contain at least 6 columns");
    			if(inclusion.find(token[1])!=inclusion.end()){
    				double p = snp_list.at(inclusion.at(token[1])).get_p_value();
    				if(p< bound_start){
    					pre_run=true;
    					quick_ref.push_back(p_partition(token[1], cur_line, -1,inclusion.at(token[1])));
    				}
    				else if(p<bound_end){
    					int category = -1;
    					if(fastscore){
    						category = c_commander.get_category(p);
    						if(category ==-2){
    							throw std::runtime_error("Undefined category!");
    						}
    					}
    					else category = (int)((p-bound_start)/bound_inter);
    					if(category <0) pre_run=true;
    					quick_ref.push_back(p_partition(token[1], cur_line, (category<0)?-1:category ,inclusion.at(token[1])));
    				}
    			}
    			cur_line++;
    		}
    }
    bim.close();
    if(quick_ref.size() ==0){
    		fprintf(stderr, "None of the SNPs met the threshold\n");
    		return;
    }
    // now we sort quick_ref such that the smaller p_value categories are always in the front
    std::sort(begin(quick_ref), end(quick_ref),
        [](PRSice::p_partition const &t1, PRSice::p_partition const &t2) {
            if(std::get<2>(t1)==std::get<2>(t2)) return std::get<1>(t1)<std::get<1>(t2);
            else return std::get<2>(t1)<std::get<2>(t2);
        }
    );
    // as proxy only affect clumping, we don't need to handle it here
    std::vector<std::vector<prs_score> > region_prs_score;
    std::vector<std::vector<prs_score> > region_best_prs_score;
    for(size_t i_region=0; i_region < c_region.size(); ++i_region){
    		region_prs_score.push_back(prs_fam);
    }
    // cur_start_index indiciates how much of quick_ref has been processed
	size_t cur_start_index = 0;
	// num_snp_included should contain the number of SNPs included for each individual
	// regions
	std::vector<size_t> num_snp_included(c_region.size());
    // first read everything that are smaller than 0
	if(pre_run) get_prs_score(quick_ref, snp_list, target, region_prs_score, num_snp_included, cur_start_index);
//	This will update the best score to equal to the base. When there is no valid threshold
//	this should be used to provide the PRS score output
	region_best_prs_score = region_prs_score;
	// if people require only the PRS score, then we will open this no_regress_out file
    	std::string output_name = c_commander.get_out()+"."+target+".all.score";
    	std::ofstream no_regress_out;
    	if(no_regress){
    		no_regress_out.open(output_name.c_str());
    		if(!no_regress_out.is_open()){
    			std::string error_message = "Cannot open file "+output_name+" for write";
    			throw std::runtime_error(error_message);
    		}
    		no_regress_out << "Threshold\tRegion";
    		for(size_t i = 0; i < prs_fam.size(); ++i) no_regress_out << "\t" << std::get<0>(prs_fam[i]);
//    	just in case we have already got something, we will print the no_regress here too
    		if(pre_run){
				for(size_t i_region=0; i_region < region_prs_score.size(); ++i_region){
					for(size_t i_reg_score=0; i_reg_score < region_prs_score[i_region].size(); ++i_reg_score){
						no_regress_out << "\t" << std::get<1>(region_prs_score[i_region][i_reg_score])/(double)num_snp_included[i_region];
					}
					no_regress_out << std::endl;
				}
    		}
    	}
    	// Now we prepare for the PRS analysis
    	double current_upper =0.0;
    	// need to initialize this to use multithreading with EIGEN library
    	Eigen::initParallel();
    	// This is the storage for therad
    	std::vector<std::thread> thread_store;
    	// getting the number of thread
    	size_t n_thread = c_commander.get_thread();
    	// some initialization. TBH, don't think this should slow down too much even if we put it in the loop
    	double p_value=0.0, r2=0.0, r2_adjust = 0.0;
    	// this should hold the PRS results
    	// The contents are threshold, r2, r2_adjust, p-vaule, num_snps
    	std::vector<std::vector<PRSice_result> > prs_results;
    	// this contain the information of the best cutoff
    	//r2, threshold at best r2 and size at best r2
    	std::vector<PRSice_best> region_best_threshold(c_region.size());
    	// try to initialize them w.r.t regions such that we don't need to worry about the
    	// index in the later analysi
    	for(size_t i = 0; i < region_best_threshold.size(); ++i){
    		region_best_threshold[i] = PRSice_best(0,0,0);
    		prs_results.push_back(std::vector<PRSice_result>(0));
    	}
    	// now start going through all the thresholds
    	// when cur_start_index == quick_ref.size(), all SNPs should have processed
    	if(cur_start_index == quick_ref.size()){
    		fprintf(stderr, "There are no valid threshold to test\n");
    		fprintf(stderr, "All SNPs have p-value lower than --lower and higher than --upper\n");
    		fprintf(stderr, "We will output the PRS score for you just to be nice\n");
    	}
//    	Calculate the null here if required. This should speed things up a lot when there are a large amount
//    	of covariates
    	double null_r2 = 0.0, null_r2_adjust=0.0, null_p=0.0;
    	if(independent_variables.cols()>2){
    		Eigen::MatrixXd covariates_only;
    		covariates_only = independent_variables;
    		covariates_only.block(0,1,covariates_only.rows(),covariates_only.cols()-2) = covariates_only.topRightCorner(covariates_only.rows(),covariates_only.cols()-2);
    		covariates_only.conservativeResize(covariates_only.rows(),covariates_only.cols()-1);
    		if(target_binary){
    			Regression::glm(phenotype, covariates_only, null_p, null_r2, 25, n_thread, true);
    		}
    		else{
    			Regression::linear_regression(phenotype, covariates_only, null_p, null_r2, null_r2_adjust, n_thread, true);
    		}

    	}
    	while(cur_start_index != quick_ref.size()){
    		// getting the current cutoff
    		current_upper = std::min((std::get<2>(quick_ref[cur_start_index])+1)*bound_inter+bound_start, bound_end);
    		fprintf(stderr, "\rProcessing %f", current_upper);

    		// now calculate the PRS for each region
    		bool reg = get_prs_score(quick_ref, snp_list, target, region_prs_score,
    				num_snp_included, cur_start_index);
    		// if regression is not required, we will simply output the score
    		if(no_regress){
    			for(size_t i_region=0; i_region < region_prs_score.size(); ++i_region){
    				no_regress_out << current_upper << "\t" << c_region.get_name(i_region);
    				for(size_t i_reg_score=0; i_reg_score < region_prs_score[i_region].size(); ++i_reg_score){
    					no_regress_out << "\t" << std::get<1>(region_prs_score[i_region][i_reg_score])/(double)num_snp_included[i_region];
    				}
    				no_regress_out << std::endl;
    			}
    		}
    		// update the boolean, basically, if we have added new SNPs AND require regression
    		// we will perform the regression analysis
    		reg=reg&&!no_regress;
    		if(reg){
    			// here we do multithreading
    			if(n_thread == 1 || region_prs_score.size()==1){
    				thread_score(independent_variables, phenotype, region_prs_score, num_snp_included, fam_index,
    						region_best_threshold, region_best_prs_score, prs_results, 0, region_prs_score.size(),
							target_binary, current_upper,n_thread);
    			}
    			else{
    				// perform multi threading
    				if(c_region.size() < n_thread){
    					for(size_t i_region = 0; i_region < c_region.size(); ++i_region){
    						thread_store.push_back(std::thread(&PRSice::thread_score, this, std::ref(independent_variables),
    								std::cref(phenotype), std::cref(region_prs_score), std::cref(num_snp_included),
									std::cref(fam_index), std::ref(region_best_threshold), std::ref(region_best_prs_score),
									std::ref(prs_results), i_region, i_region+1, target_binary, current_upper,1));
    					}
    				}else{
    					int job_size = c_region.size()/n_thread;
    					int remain = c_region.size()%n_thread;
    					size_t start =0;
    					for(size_t i_thread = 0; i_thread < n_thread; ++i_thread){
    						size_t ending = start+job_size+(remain>0);
    						ending = (ending>c_region.size())? c_region.size(): ending;
    						 thread_store.push_back(std::thread(&PRSice::thread_score, this, std::ref(independent_variables),
    								 std::cref(phenotype), std::cref(region_prs_score), std::cref(num_snp_included),
									 std::cref(fam_index), std::ref(region_best_threshold), std::ref(region_best_prs_score),
									 std::ref(prs_results), start, ending, target_binary, current_upper,1));
    						start=ending;
    						remain--;
    					}
    				}
    				// joining the threads
    				for(size_t i_thread = 0; i_thread < thread_store.size(); ++i_thread){
    					thread_store[i_thread].join();
    				}
    			}
    		}
    	}
    	fprintf(stderr, "\n");
    	no_regress_out.close();
    	if(no_regress) return;
    	// now perform the output for all the best scores and PRSice results
    	// unfortunately this is not as simple, as the base and target might contain some path information
    	// instead, we need to tokenize them before we do this

	std::string output_prefix = c_commander.get_out()+"."+m_current_base+"."+target;
    	for(size_t i_region = 0; i_region< prs_results.size(); ++i_region){
    		output_name = output_prefix+"."+c_region.get_name(i_region);
    		std::string out_best = output_name+".best";
    		std::string out_prsice = output_name+".prsice";
    		std::string out_pheno = output_name+".pheno";
    		std::ofstream best_out, prsice_out, pheno_out;
    		best_out.open(out_best.c_str());
    		prsice_out.open(out_prsice.c_str());
    		pheno_out.open(out_pheno.c_str());
    		if(!best_out.is_open()){
    			std::string error_message = "ERROR: Cannot open file: " +out_best+" to write";
    			throw std::runtime_error(error_message);
    		}
    		if(!prsice_out.is_open()){
    			std::string error_message = "ERROR: Cannot open file: " +out_prsice+" to write";
    			throw std::runtime_error(error_message);
    		}
    		if(!pheno_out.is_open()){
    			// we won't need this once we have incorporate the plotting into PRSice
    			std::string error_message = "ERROR: Cannot open file: "+out_pheno+" to write";
    			throw std::runtime_error(error_message);
    		}
    		best_out << "IID\tprs_"<<std::get<1>(region_best_threshold[i_region]) << std::endl;
    		prsice_out << "Threshold\tR2\tP\tNum_SNP" << std::endl;
    		pheno_out << "IID\tPheno";
    		// We want to skip the intercept for now
    		for(size_t i_col=2; i_col < independent_variables.cols(); ++i_col) pheno_out << "\tCov"<<std::to_string(i_col);
    		pheno_out << std::endl;
    		for(size_t i_prsice=0; i_prsice< prs_results[i_region].size();++i_prsice){
    			prsice_out << std::get<0>(prs_results[i_region][i_prsice]) << "\t" <<
    					std::get<1>(prs_results[i_region][i_prsice])-null_r2 << "\t" <<
						std::get<3>(prs_results[i_region][i_prsice])<< "\t" <<
						std::get<4>(prs_results[i_region][i_prsice]) << std::endl;
    		}
    		for(size_t i_score=0; i_score < region_best_prs_score[i_region].size(); ++i_score){
    			pheno_out << std::get<0>(region_best_prs_score[i_region][i_score]) << "\t" << phenotype(i_score);
    			for(size_t i_col=2; i_col < independent_variables.cols(); ++i_col){
    				pheno_out << "\t" << independent_variables(i_score,i_col);
    			}
    			pheno_out << std::endl;
    			best_out << std::get<0>(region_best_prs_score[i_region][i_score]) << "\t" <<
    					std::get<1>(region_best_prs_score[i_region][i_score])/std::get<2>(region_best_threshold[i_region])<< std::endl;

    		}
    		prsice_out.close();
    		best_out.close();
    		pheno_out.close();
    	}
}

void PRSice::thread_score( Eigen::MatrixXd &independent_variables, const Eigen::VectorXd &c_pheno,
        		const std::vector<std::vector<PRSice::prs_score > > &c_region_prs_score,
        		const std::vector<size_t> &c_num_snp_included, const std::map<std::string, size_t> &c_fam_index,
        		std::vector<PRSice::PRSice_best > &region_best_threshold,
        		std::vector<std::vector<PRSice::prs_score> > & region_best_prs_score,
			std::vector<std::vector<PRSice::PRSice_result> > &region_result,
        		size_t region_start, size_t region_end, bool target_binary, double threshold,
			size_t thread){

	Eigen::MatrixXd X;
	bool thread_safe=false;
	// so we will only copy the matrix when it is not thread safe to do so
	if(region_start==0 && region_end == c_region_prs_score.size()) thread_safe = true;
	else X = independent_variables;
	// c_prs_region_score = prs score of the region at the current threshold
	// c_num_snp_include = num of SNP in region at current threshold
	// c_fam_index = index for each individual on the matrix (mainly to deal with missing data)
	// prs_best_info = threshold info with best score. p here is the threshold, not the p-value for PRSice
	// prs_best_score = store the best score
	// region_result = PRSice result, containing all the required information for output
	//std::vector<PRSice::PRSice_result> temp_region_result;
	double r2 = 0.0, r2_adjust=0.0, p_value = 0.0;
	for(size_t iter = region_start; iter < region_end; ++iter){
		std::vector<PRSice::prs_score>::const_iterator prs_iter = c_region_prs_score.at(iter).begin();
		std::vector<PRSice::prs_score>::const_iterator prs_end = c_region_prs_score.at(iter).end();
		for(;prs_iter != prs_end; ++prs_iter){
			std::string sample = std::get<0>(*prs_iter);
			if(c_fam_index.find(sample)!=c_fam_index.end()){
				if(thread_safe) independent_variables(c_fam_index.at(sample), 1) = std::get<1>(*prs_iter)/(double)c_num_snp_included[iter];
				else X(c_fam_index.at(sample), 1) = std::get<1>(*prs_iter)/(double)c_num_snp_included[iter];
			}
		}
		if(target_binary){
			try{
				if(thread_safe) Regression::glm(c_pheno, independent_variables, p_value, r2, 25, thread, true);
				else Regression::glm(c_pheno, X, p_value, r2, 25, thread, true);
			}
			catch(const std::runtime_error &error){
				// This should only happen when the glm doesn't converge.
				// Let's hope that won't happen...
				fprintf(stderr, "ERROR: GLM model did not converge!\n");
				fprintf(stderr, "       Please send me the DEBUG files\n");
				std::ofstream debug;
				debug.open("DEBUG");
				if(thread_safe) debug << independent_variables << std::endl;
				else 	debug << X<< std::endl;
				debug.close();
				debug.open("DEBUG.y");
				debug << c_pheno << std::endl;
				debug.close();
				std::cerr << "ERROR: " << error.what() << std::endl;
				exit(-1);
			}
		}
		else{
			if(thread_safe) Regression::linear_regression(c_pheno, independent_variables, p_value, r2, r2_adjust, thread, true);
			else Regression::linear_regression(c_pheno, X, p_value, r2, r2_adjust, thread, true);
		}

		// This should be thread safe as each thread will only mind their own region
        // now add the PRS result to the vectors (hopefully won't be out off scope
		region_result.at(iter).push_back(PRSice::PRSice_result(threshold, r2, r2_adjust, p_value, c_num_snp_included.at(iter)));
		// It this is the best r2, then we will add it
		if(std::get<0>(region_best_threshold[iter]) < r2){
			region_best_threshold[iter] = PRSice::PRSice_best(r2, threshold, c_num_snp_included.at(iter));
			region_best_prs_score[iter] = c_region_prs_score.at(iter);
		}

	}

}


Eigen::VectorXd PRSice::gen_pheno_vec(const std::string &c_target,
		const std::string c_pheno, bool target_binary,
		std::map<std::string, size_t> &fam_index,
		std::vector<std::pair<std::string, double> > &prs_score)
{
	std::vector<double> phenotype_store;
	std::ifstream pheno_file;
	std::string fam_name = c_target+".fam";
	std::string line;
	size_t cur_index = 0;
	// add this information just because I have made the mistake regarding the case
	// control label before
	size_t num_case = 0, num_control = 0;
	if(c_pheno.empty()){
		// Use fam file
		pheno_file.open(fam_name.c_str());
		while(std::getline(pheno_file, line)){
			misc::trim(line);
			if(!line.empty()){
				std::vector<std::string> token = misc::split(line);
				if(token.size() < 6) throw std::runtime_error("Malformed fam file, should contain at least 6 columns");
				prs_score.push_back(std::pair<std::string, double>(token[1], 0.0));
				if(token[5].compare("NA")!=0){
					try{
						if(target_binary){
							double temp = misc::convert<int>(token[5]);
							if(temp-1>=0 && temp-1<2){
								fam_index[token[1]]=cur_index;
								phenotype_store.push_back(temp-1);
								cur_index++;
								if(temp==2) num_case++;
								if(temp==1) num_control++;
							}
							// anything other than 1 or 2 will be treated as missing
						}
						else{
							double temp = misc::convert<double>(token[5]);
							fam_index[token[1]]=cur_index;
							phenotype_store.push_back(temp);
							cur_index++;
						}
					}
					catch(const std::runtime_error &error){
						// if we can't handle it, it is missing
					}
				}
			}
		}
		pheno_file.close();
	}
	else{
		// Use pheno file
		std::ifstream fam;
		fam.open(fam_name.c_str());
		pheno_file.open(c_pheno.c_str());
		// First form the map using the fam
		std::map<std::string, std::string> pheno_info;
		while(std::getline(pheno_file, line)){
			misc::trim(line);
			if(!line.empty()){
				std::vector<std::string> token = misc::split(line);
				if(token.size() < 2) std::runtime_error("Malformed pheno file, should contain at least 2 columns");
				pheno_info[token[0]] = token[1];
			}
		}
		pheno_file.close();
		while(std::getline(fam, line)){
			misc::trim(line);
			if(!line.empty()){
				std::vector<std::string> token = misc::split(line);
				if(token.size() < 6) std::runtime_error("Malformed fam file, should contain at least 6 columns");
				if(pheno_info.find(token[1])!= pheno_info.end()){
					std::string p = pheno_info[token[1]];
					if(p.compare("NA")!=0){
						try{
							if(target_binary){
								double temp = misc::convert<int>(p);
								if(temp-1>=0 && temp-1<=2){
									fam_index[token[1]]=cur_index;
									phenotype_store.push_back(temp-1);
									cur_index++;
									if(temp==2) num_case++;
									if(temp==1) num_control++;
								}
								// anything other than 1 or 2 will be treated as missing
							}
							else{
								double temp = misc::convert<double>(p);
								fam_index[token[1]]=cur_index;
								phenotype_store.push_back(temp);
								cur_index++;
							}
						}
						catch(const std::runtime_error &error){ }
					}
				}
			}
		}
		fam.close();
	}
	if(phenotype_store.size()==0) throw std::runtime_error("No phenotypes present");
	Eigen::Map<Eigen::VectorXd> res(phenotype_store.data(), phenotype_store.size());
	if(target_binary){
		if(num_control==0) throw std::runtime_error("There are no control samples");
		if(num_case==0) throw std::runtime_error("There are no cases");
	    fprintf(stderr,"Number of controls : %zu\n", num_control);
	    fprintf(stderr,"Number of cases : %zu\n", num_case);
	}
	else{
	    fprintf(stderr,"Number of sample(s) with phenotype  : %zu\n", res.rows());
	}
	return res;
}

Eigen::MatrixXd PRSice::gen_cov_matrix(const std::string &target, const std::string &c_cov_file, const std::vector<std::string> &c_cov_header, std::map<std::string, size_t> &fam_index){
	size_t num_sample = fam_index.size();
	if(c_cov_file.empty()){
		// changed this to 2 to avoid needing to add the intercept for every regression
		return Eigen::MatrixXd::Ones(num_sample,2);
	}
	else{
		// First read in the header
		std::ifstream cov;
		cov.open(c_cov_file.c_str());
		if(!cov.is_open()){
			std::string error_message = "ERROR: Cannot open covariate file: "+c_cov_file;
			throw std::runtime_error(error_message);
		}
		std::string line;
		std::vector<size_t> cov_index;
		int max_index = 0;
		std::getline(cov, line);
		if(!line.empty()){
			std::vector<std::string> token = misc::split(line);
			if(c_cov_header.size() == 0){
				// if no header is provided, we will use all the covariates included
				for(size_t i = 1; i < token.size(); ++i) cov_index.push_back(i);
				max_index = cov_index.size()-1;
			}
			else{
				std::map<std::string, bool> include;
				// if specific headers are provided, we should only include them
				for(size_t i = 0; i < c_cov_header.size(); ++i) include[c_cov_header[i]]= true;
				for(size_t i = 1; i < token.size(); ++i){
					if(include.find(token[i])!=include.end()){
						cov_index.push_back(i);
						if(i > max_index) max_index=i;
					}
				}
			}
		}
		else throw std::runtime_error("First line of covariate file is empty!");
		// now we know how much we are working with
		// need to include the space for the polygenic risk score
		// +2 because the first row is for intercept and the second row is for the PRS
		Eigen::MatrixXd result = Eigen::MatrixXd::Ones(num_sample, cov_index.size()+2);
		while(std::getline(cov, line)){
			misc::trim(line);
			if(!line.empty()){
				std::vector<std::string> token = misc::split(line);
				if(token.size() <= max_index){
					std::string error_message = "ERROR: Malformed covariate file, should contain at least "+std::to_string(max_index+1)+" column!";
					throw std::runtime_error(error_message);
				}
				if(fam_index.find(token[1])!= fam_index.end()){
					int index = fam_index[token[1]];
					for(size_t i = 0; i < cov_index.size(); ++i){
						try{
							double temp = misc::convert<double>(token[cov_index[i]]);
							result(index, i+2) = temp;
						}
						catch(const std::runtime_error &error){
							// All missing values are treated as 0
							result(index, i+2) = 0;
						}
					}
				}
			}
		}
		return result;
	}
}

// basically update the score vector to contain the new polygenic score
bool PRSice::get_prs_score(const std::vector<PRSice::p_partition> &quick_ref,
		const boost::ptr_vector<SNP> &snp_list, const std::string &target,
	std::vector< std::vector<std::pair<std::string, double> > > &prs_score,
	std::vector<size_t> &num_snp_included, size_t &cur_index)
{
	// Here is the actual calculation of the PRS
	if(quick_ref.size()==0) return false; // nothing to do
	int prev_index =prev_index = std::get<2>(quick_ref[cur_index]);
	int end_index = 0;
	bool ended =false;
	for(size_t i = cur_index; i < quick_ref.size(); ++i){
		if(std::get<2>(quick_ref[i]) != prev_index && std::get<2>(quick_ref[i])>=0 ){
			end_index = i;
			ended=true;
			break;
		}
		else if(std::get<2>(quick_ref[i])!=prev_index) prev_index=std::get<2>(quick_ref[i]); // only when the category is still negative
		// Use as part of the output
		for(size_t i_region=0; i_region< num_snp_included.size(); ++i_region){
			if(snp_list.at(std::get<3>(quick_ref[i])).in(i_region)) num_snp_included[i_region]++;
		}
	}
	if(!ended) end_index = quick_ref.size();
	PLINK prs(target);
	prs.initialize();
	prs.get_score(quick_ref, snp_list, prs_score, cur_index, end_index);

	cur_index = end_index;
	return true;
}



void PRSice::update_inclusion(std::map<std::string, size_t> &inclusion, const std::string &c_target_bim_name,
		boost::ptr_vector<SNP> &snp_list, const std::map<std::string, size_t> &c_snp_index){
    std::ifstream target_file;
    target_file.open(c_target_bim_name.c_str());
    if(!target_file.is_open()){
        std::string error_message = "Cannot open target bim file: "+c_target_bim_name;
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
                    std::string error_message = "Negative coordinate of SNP in "+c_target_bim_name;
                    throw std::runtime_error(error_message);
                }
                loc = temp;
            }
            catch(std::runtime_error &error){
                std::string error_message = "Non-numeric coordinate of SNP in "+c_target_bim_name;
                throw std::runtime_error(error_message);
            }
            std::string ref_allele = token[4];
            std::string alt_allele = token[5];
            if(inclusion.find(rsid) != inclusion.end() && c_snp_index.find(rsid)!=c_snp_index.end()){
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
                    // This suggest there are mismatching between the LD
                    // and the target file
                    if(inclusion.find(rsid)!=inclusion.end()) inclusion.erase(rsid);
                }
                else{
                    // not ambiguous, now do soft checking
                    size_t index = c_snp_index.at(rsid);
                    bool same = snp_list[index].check_loc(chr, loc, ref_allele, alt_allele);
                    if(!same){
                        fprintf(stderr, "WARNING: %s differ between target and base file\n", rsid.c_str());
                        fprintf(stderr, "         It is advised that you check the files are \n");
                        fprintf(stderr, "         From the same genome build\n");
                    }
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




PRSice::PRSice()
{
    //ctor
}

PRSice::~PRSice()
{
    //dtor
}

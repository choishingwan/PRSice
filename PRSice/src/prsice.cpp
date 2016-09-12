#include "prsice.hpp"

void PRSice::get_snp(boost::ptr_vector<SNP> &snp_list,
		std::map<std::string, size_t> &snp_index, const std::string &c_input, bool beta,
		const Commander &c_commander, Region &region){
	std::vector<int> index = SNP::get_index(c_commander, c_input);
    std::ifstream snp_file;
    snp_file.open(c_input.c_str());
    if(!snp_file.is_open()){
    		std::string error_message = "Cannot open base file: "+c_input;
    		throw std::runtime_error(error_message);
    }
    size_t num_duplicated = 0;
    size_t num_stat_not_convertible = 0;
	size_t num_p_not_convertible = 0;
    std::string line;
    int max_index = index.back();
    bool se_error = false;
	if(!c_commander.index()) std::getline(snp_file, line);
	bool read_error = false;
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
		    		snp_index[rs_id] =snp_list.size();
		      	snp_list.push_back(new SNP(rs_id, chr, loc, ref_allele, alt_allele, stat, se, pvalue, region.check(chr, loc), region.size()));
		   	}
		}
		if(read_error) throw std::runtime_error("Please check if you have the correct input");
	}
	snp_file.close();
	fprintf(stderr, "Number of duplicated SNPs : %zu\n", num_duplicated);
	fprintf(stderr, "Number of SNPs from base  : %zu\n", snp_list.size());
	if(num_stat_not_convertible!=0) fprintf(stderr, "Failed to convert %zu OR/beta\n", num_stat_not_convertible);
	if(num_p_not_convertible!=0) fprintf(stderr, "Failed to convert %zu p-value\n", num_p_not_convertible);
}

void PRSice::calculate_score(const Commander &c_commander, bool target_binary,
		const std::string c_target, const std::map<std::string, size_t> &inclusion,
		const boost::ptr_vector<SNP> &snp_list)
{
	// First, get the phenotype and covariate matrix
	std::map<std::string, size_t> fam_index;
    std::vector<std::pair<std::string, double> > prs_score, prs_best_score;
	Eigen::VectorXd phenotype = gen_pheno_vec(	c_target, c_commander.get_pheno(),
												target_binary, fam_index, prs_score);
	Eigen::MatrixXd covariates = gen_cov_matrix(	c_target, c_commander.get_cov_file(),
												c_commander.get_cov_header(), fam_index);
	std::string output_name = c_commander.get_out()+"."+c_target+".prsice";
    std::ofstream prs_out, prs_best;
    prs_out.open(output_name.c_str());
    if(!prs_out.is_open()){
    		std::string error_message= "Cannot open file "+output_name+" for write!";
    		throw std::runtime_error(error_message);
    }
    prs_out << "#Threshold\tR2\tR2_Adjusted\tP-value\tNum_Snp"<< std::endl;
	double bound_start = c_commander.get_lower();
	double bound_end = c_commander.get_upper();
	double bound_inter = c_commander.get_inter();
    double current_lower = 0.0, p_value=0.0, r2=0.0, r2_adjust = 0.0, best_r2 =0.0;
    size_t num_snp_included = 0;
    // First, construct a SNP vector,
    std::string bim_name = c_target+".bim";
    // Use for speeding up the bed file reading
    std::vector<p_partition > quick_ref;
    std::ifstream bim;
    bim.open(bim_name.c_str());
    if(!bim.is_open()){
    		std::string error_message = "Cannot open bim file: " +bim_name;
    		throw std::runtime_error(error_message);
    }
    std::string line;
    size_t cur_line = 0;
    while(getline(bim, line)){
    		misc::trim(line);
    		if(!line.empty()){
    			std::vector<std::string> token = misc::split(line);
    			if(token.size() < 6) throw std::runtime_error("Malformed bim file, should contain at least 6 columns");
    			if(inclusion.find(token[1])!=inclusion.end()){
    				double p = snp_list.at(inclusion.at(token[1])).get_p_value();
    				if(p<bound_end && p>bound_start){
    					quick_ref.push_back(p_partition(token[1], cur_line, (int)((p-bound_start)/bound_inter),inclusion.at(token[1])));
    				}
    			}
    			cur_line++;
    		}
    }
    bim.close();
    std::sort(begin(quick_ref), end(quick_ref),
        [](PRSice::p_partition const &t1, PRSice::p_partition const &t2) {
            if(std::get<2>(t1)==std::get<2>(t2)) return std::get<1>(t1)<std::get<1>(t2);
            else return std::get<2>(t1)<std::get<2>(t2);
        }
    );
    // Now the quick_ref can be use for fast SNP reading
    // indicate the start of the quick_ref;
    size_t cur_start_index = 0;
    double current_upper=0.0;
    while(cur_start_index!=quick_ref.size()){
		current_upper = std::min((std::get<2>(quick_ref[cur_start_index])+1)*bound_inter, bound_end);
		fprintf(stderr, "\rProcessing %f", current_upper);
    		bool reg = get_prs_score(quick_ref, snp_list, c_target, prs_score,
    				num_snp_included, cur_start_index);
    		if(reg){
    			for(size_t i = 0; i < prs_score.size(); ++i){
    				std::string sample = std::get<0>(prs_score[i]);
    				if(fam_index.find(sample)!=fam_index.end()){
    					covariates(fam_index[sample], 0) = std::get<1>(prs_score[i]);
    				}
    			}
    		}
    		if(reg && target_binary){
    			try{
    			Regression::glm(phenotype, covariates, p_value, r2, 25, c_commander.get_thread(), true);
    			}
    			catch(const std::runtime_error &error){
    				std::ofstream debug;
    				debug.open("DEBUG");
    				debug << covariates<< std::endl;
    				debug.close();
    				debug.open("DEBUG.y");
    				debug << phenotype << std::endl;
    				debug.close();
					std::cerr << "ERROR: " << error.what() << std::endl;
    				exit(-1);
    			}
    			double null_p, null_r2;
    			if(covariates.cols() > 1){
    				Regression::glm(phenotype, covariates.block(0,1,covariates.rows(),
    						covariates.cols()-1), null_p, null_r2, 25, c_commander.get_thread(), true);
    				r2-=null_r2;
    			}
    			prs_out << current_upper << "\t" << r2 << "\tNA\t" << p_value  <<"\t" << num_snp_included << std::endl;
    		}
    		else if(reg){
    			Regression::linear_regression(phenotype, covariates, p_value, r2, r2_adjust,
    					c_commander.get_thread(), true);
    			prs_out << current_upper << "\t" << r2 << "\t" << r2_adjust << "\t" << p_value <<"\t" << num_snp_included << std::endl;
    		}
    		if(r2 > best_r2){
    			best_r2 = r2;
    			prs_best_score = prs_score;
    		}
    }
    prs_out.close();
    output_name = c_commander.get_out()+"."+c_target+".best.prsice";
    prs_best.open(output_name.c_str());
    if(!prs_best.is_open()){
    		std::string error_message = "ERROR: Cannot open file "+output_name+" for write";
    		throw std::runtime_error(error_message);
    }
    prs_best << "#SampleID\tPRS" << std::endl;
    for(size_t i = 0; i < prs_best_score.size(); ++i){
    		std::string sample = std::get<0>(prs_best_score[i]);
    		if(fam_index.find(sample)!=fam_index.end()){
    			prs_best << sample << "\t" << std::get<1>(prs_best_score[i]) << std::endl;
    		}
    		else prs_best << sample << "\tNA" << std::endl;
    }
    prs_best.close();
    fprintf(stderr, "\n");
}

void PRSice::process(const std::string &c_input, bool beta, const Commander &c_commander, Region &region){
	//We will now use the ptr_vector
	boost::ptr_vector<SNP> snp_list;
    std::map<std::string, size_t> snp_index;
    get_snp(snp_list, snp_index, c_input, beta, c_commander, region);
    std::vector<std::string> target = c_commander.get_target();
	fprintf(stderr,"\nClumping Parameters: \n");
    fprintf(stderr,"==============================\n");
    fprintf(stderr,"P-Threshold  : %f\n", c_commander.get_clump_p());
    fprintf(stderr,"R2-Threshold : %f\n", c_commander.get_clump_r2());
    fprintf(stderr,"Window Size  : %zu\n", c_commander.get_clump_kb());

    for(size_t i_target = 0; i_target < target.size(); ++i_target){
    		fprintf(stderr,"\nStart processing: %s\n", target[i_target].c_str());
    	    fprintf(stderr,"==============================\n");
    	    // Key = rsid, Value = index w.r.t. snp_list
        std::map<std::string, size_t> inclusion;
        bool has_ld = !c_commander.ld_prefix().empty();
        std::string ld_file = (has_ld)? c_commander.ld_prefix(): target[i_target];
        PLINK clump(ld_file, c_commander.get_thread());
        if(has_ld) fprintf(stderr,"\nIn LD Reference %s\n", ld_file.c_str());
        else fprintf(stderr,"\nStart performing clumping\n");
        // after initialize, we have a list of SNPs to include
        // also, the bim file should be closed
        clump.initialize(inclusion, snp_list, snp_index);
        if(has_ld){
        		std::string target_bim_name = target[i_target]+".bim";
        		// perform additional filtering
            fprintf(stderr,"\nIn target %s\n", ld_file.c_str());
        		update_inclusion(inclusion, target_bim_name, snp_list, snp_index);
        }
        clump.start_clumping(inclusion, snp_list, snp_index, c_commander.get_clump_p(),
        		c_commander.get_clump_r2(), c_commander.get_clump_kb());
        // inclusion now represent the index SNPs
        calculate_score(c_commander, c_commander.get_target_binary(i_target),
        		target[i_target], inclusion, snp_list);
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
		return Eigen::MatrixXd::Zero(num_sample,1);
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
				for(size_t i = 1; i < token.size(); ++i) cov_index.push_back(i);
				max_index = cov_index.size()-1;
			}
			else{
				std::map<std::string, bool> include;
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
		Eigen::MatrixXd result = Eigen::MatrixXd::Zero(num_sample, cov_index.size()+1);
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
							result(index, i+1) = temp;
						}
						catch(const std::runtime_error &error){
							// Here, we handle all missing covariates as 0
							result(index, i+1) = 0;
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
	std::vector<std::pair<std::string, double> > &prs_score, size_t &num_snp_included, size_t &cur_index)
{
	// Here we will make a different inclusion for the inclusion
	if(quick_ref.size()==0) return false; // again, nothing to do
	size_t prev_index =prev_index = std::get<2>(quick_ref[cur_index]);;
	size_t end_index = 0;
	bool ended =false;
	for(size_t i = cur_index; i < quick_ref.size(); ++i){
		if(std::get<2>(quick_ref[i]) != prev_index){
			end_index = i;
			ended=true;
			break;
		}
		// Use as part of the output
		num_snp_included++;
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


void PRSice::run(const Commander &c_commander, Region &region){
    std::vector<std::string> base = c_commander.get_base();
    int num_base = base.size();
    if(num_base == 0) throw std::runtime_error("There is no base case to run");
    if(num_base < 0) throw std::runtime_error("Negative number of base");
    else{
        for(size_t i = 0; i < num_base; ++i){
        		region.reset();
            process(base[i], c_commander.get_base_binary(i), c_commander, region);
        }
    }
    // completed when we reach here
}

PRSice::PRSice()
{
    //ctor
}

PRSice::~PRSice()
{
    //dtor
}

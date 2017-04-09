#include "prsice.hpp"

std::mutex PRSice::score_mutex;


void PRSice::pheno_check(const Commander &c_commander) {
	std::vector < std::string > pheno_header = c_commander.get_pheno_col();
	std::string pheno_file = c_commander.get_pheno();
	if (pheno_header.size() != 0 && pheno_file.empty()) {
		throw std::runtime_error( "You must provide a phenotype file for multiple phenotype analysis");
	}
	if (!pheno_file.empty()) {
		pheno_info.use_pheno = false;
		pheno_info.binary.push_back(c_commander.is_binary(0));
	} else {
		std::ifstream pheno;
		pheno.open(pheno_file.c_str());
		if (!pheno.is_open()) {
			std::string error_message = "Cannot open phenotype file: " + pheno_file;
			throw std::runtime_error(error_message);
		}
		std::string line;
		std::getline(pheno, line);
		if (line.empty()) {
			throw std::runtime_error( "Cannot have empty header line for phenotype file!");
		}
		pheno.close();
		misc::trim(line);
		std::vector < std::string > col = misc::split(line);
		bool found = false;
		std::unordered_map<std::string, bool> dup_col;
		if (pheno_header.size() == 0)
		{
			// use the second column from the pheno file
			pheno_info.use_pheno = true;
			pheno_info.col.push_back(1+!m_ignore_fid);
			pheno_info.name.push_back("");
			pheno_info.order.push_back(0);
			pheno_info.binary.push_back(c_commander.is_binary(0));
		}
		else
		{
			for (size_t i_pheno = 0; i_pheno < pheno_header.size(); ++i_pheno) {
				if (dup_col.find(pheno_header[i_pheno]) == dup_col.end()) {
					found = false;
					dup_col[pheno_header[i_pheno]] = true;
					// start from 1+!m_ignore_fid to skip the iid and fid part
					for (size_t i_column = 1+!m_ignore_fid; i_column < col.size(); ++i_column) {
						if (col[i_column].compare(pheno_header[i_pheno]) == 0) {
							found = true;
							pheno_info.use_pheno=true;
							pheno_info.col.push_back(i_column);
							pheno_info.name.push_back(pheno_header[i_pheno]);
							pheno_info.order.push_back(i_pheno);
							pheno_info.binary.push_back(c_commander.is_binary(i_pheno));
							break;
						}
					}
					if (!found) {
						fprintf(stderr, "Phenotype: %s cannot be found in phenotype file\n",
								pheno_header[i_pheno].c_str());
					}
				}
			}
		}

	}
	size_t num_pheno= (pheno_info.use_pheno)? pheno_info.col.size() : 1;
	fprintf(stderr, "There are a total of %zu phenotype to process\n",num_pheno);
}

void PRSice::init_matrix(const commander &c_commander, const size_t pheno_index, Genotype &target,
		const bool prslice)
{
	m_null_r2 = 0.0;
	m_phenotype = Eigen::VectorXd::Zero(0);
	m_independent_variables.resize(0,0);
	bool no_regress = c_commander.no_regress();
	bool all = c_commander.all();
	std::string pheno_file = c_commander.pheno_file();
	std::string output_name = c_commander.out();

	std::ofstream all_out;
	bool multi = pheno_info.col.size()>1;
	if(all && !prslice)
	{
		std::string all_out_name = output_name;
		if(multi)
		{
			all_out_name.append("."+pheno_info.name[pheno_index]);
		}
		all_out_name.append(".all.score");
		all_out.open(all_out_name.c_str());
		if(!all_out.is_open())
		{
			std::string error_message = "Cannot open file: "+all_out_name+" for write";
			throw std::runtime_error(error_message);
		}
	}
	m_sample_names = target.sample_names();
	gen_pheno_vec(pheno_file, pheno_index, !no_regress);
	if (!no_regress)
	{
		gen_cov_matrix(c_commander.get_cov_file(), c_commander.get_cov_header());
	}


}

void PRSice::gen_pheno_vec(const std::string &pheno_file_name, const int pheno_index, bool regress)
{
	std::vector<double> pheno_store(m_sample_names.size());
	bool binary = pheno_info.binary[pheno_index];
	int max_num = 0;
	int num_case =0;
	int num_control =0;
	size_t num_not_found = 0;
	if(pheno_info.use_pheno) // use phenotype file
	{
		std::ifstream pheno_file;
		pheno_file.open(pheno_file_name.c_str());
		if(!pheno_file.is_open())
		{
			std::string error_message = "Cannot open phenotype file: " + pheno_file_name;
			throw std::runtime_error(error_message);
		}

		std::unordered_map<std::string, std::string> phenotype_info;
		while (std::getline(pheno_file, line))
		{
			misc::trim(line);
			if (line.empty()) continue;
			std::vector < std::string > token = misc::split(line);
			if (token.size() < pheno_index + 1)
			{
				std::string error_message = "Malformed pheno file, should contain at least "
						+ std::to_string(pheno_index + 1) + " columns";
				throw std::runtime_error(error_message);
			}
			std::string id =(m_ignore_fid)? token[0]:token[0]+"_"+token[1];
			phenotype_info[id] = token[pheno_index];
		}
		pheno_file.close();
		// now go through the sample information
		size_t cur_index = 0;
		for(auto &&sample : m_sample_names)
		{
			std::string id = (m_ignore_fid)? sample.IID : sample.FID+"_"+sample.IID;
			if(phenotype_info.find(id)!=phenotype_info.end())
			{
				try{
					if(binary)
					{
						int temp = misc::convert<int>(phenotype_info[id]);
						if(temp >=0 && temp <= 2)
						{
							m_sample_with_phenotypes[id]=cur_index;
							phenotype_store[cur_index++] = temp;
							max_num = (temp>max_num)?temp:max_num;
							num_case+=(temp==1);
							num_control+=(temp==0);
						}
						else
						{
							sample.included = false;
						}
					}
					else
					{
						m_sample_with_phenotypes[id]=cur_index;
						phenotype_store[cur_index++]=misc::convert<double>(phenotype_info[id]);
					}
				}catch(const std::runtime_error &error){
					sample.included=false;
				}
			}
			else
			{
				sample.included =false;
				num_not_found++;
			}
		}
	}
	else
	{
		// directly extract it from the sample_name stuff
		size_t cur_index =0;
		for(auto &&sample: m_sample_names)
		{
			if(sample.phenotype.compare("NA")==0){
				sample.included=false;
				continue;
			}
			try{
				if(binary)
				{
					int temp = misc::convert<int>(sample.phenotype);
					if(temp >=0 && temp <= 2)
					{
						m_sample_with_phenotypes[m_ignore_fid?sample.IID:sample.FID+"_"+sample.IID]=cur_index;
						phenotype_store[cur_index++] = temp;
						max_num = (temp>max_num)?temp:max_num;
						num_case+=(temp==1);
						num_control+=(temp==0);
					}
					else
					{
						sample.included = false;
					}
				}
				else
				{
					m_sample_with_phenotypes[m_ignore_fid?sample.IID:sample.FID+"_"+sample.IID]=cur_index;
					phenotype_store[cur_index++]=misc::convert<double>(sample.phenotype);
				}
			}catch(const std::runtime_error &error){
				sample.included = false;
			}
		}
	}
	if (num_not_found != 0)
	{
		fprintf(stderr, "Number of missing samples: %zu\n", n_not_found);
	}
	bool error = false;
	if(max_num > 1 && binary)
	{
		num_case = 0;
		num_control = 0;
		for(auto &&pheno : phenotype_store.size())
		{
			pheno--;
			if(pheno < 0) error = true;
			else (pheno==1)? num_case++: num_control++;
		}
	}
	if(error)
	{
		throw std::runtime_error("Mixed encoding! Both 0/1 and 1/2 encoding found!");
	}
	if (phenotype_store.size() == 0) throw std::runtime_error("No phenotype presented");
	m_phenotype = Eigen::Map<Eigen::VectorXd>(phenotype_store.data(), phenotype_store.size());
	if (binary) {
		fprintf(stderr, "Number of controls : %zu\n", num_control);
		fprintf(stderr, "Number of cases : %zu\n", num_case);
		if (regress) {
			if (num_control == 0) throw std::runtime_error("There are no control samples");
			if (num_case == 0) throw std::runtime_error("There are no cases");
		}
	} else {
		fprintf(stderr, "Number of sample(s) with phenotype  : %zu\n", m_phenotype.rows());
	}
}


std::vector<size_t> PRSice::get_cov_index(const std::string &c_cov_file,
		const std::vector<std::string> &c_cov_header)
{
	std::vector<size_t> cov_index;
	std::ifstream cov;
	cov.open(c_cov_file.c_str());
	if (!cov.is_open())
	{
		std::string error_message = "ERROR: Cannot open covariate file: " + c_cov_file;
		throw std::runtime_error(error_message);
	}
	std::string line;
	size_t num_valid = 0;
	std::getline(cov, line);
	// obtain the header information of the covariate file
	if(line.empty()) throw std::runtime_error("First line of covariate file is empty!");
	std::vector < std::string > token = misc::split(line);
	if (c_cov_header.size() == 0)
	{
		// if no header is provided, we will use all the covariates included
		for (size_t i = 1+!m_ignore_fid; i < token.size(); ++i) cov_index.push_back(i); //FID, therefore+1
	}
	else
	{
		std::unordered_set<std::string> included;
		for (auto &&cov: c_cov_header)
		{
			if(included.find(cov) == included.end())// to avoid duplicated covariance headers
			{
				included.insert(cov);
			}
		}
		// same, +1 when fid is include
		for (size_t i_header = 1+!m_ignore_fid; i_header < token.size(); ++i_header)
		{
			if (included.find(token[i_header]) != included.end())
			{
				cov_index.push_back(i_header);
			}
		}
	}
	std::sort(cov_index.begin(), cov_index.end());
	return cov_index;
}

void PRSice::gen_cov_matrix(const std::string &c_cov_file,
		const std::vector<std::string> &c_cov_header)
{
	size_t num_sample = m_sample_with_phenotypes.size();
	std::vector<size_t> cov_index = get_cov_index(c_cov_file, c_cov_header);
	fprintf(stderr, "\nStart processing the covariates\n");
	fprintf(stderr, "==============================\n");
	std::vector < std::pair<std::string, size_t> > valid_sample_index;
	m_independent_variables = Eigen::MatrixXd::Ones(num_sample, cov_index.size() + 2);
	bool valid=true;
	std::ifstream cov;
	cov.open(c_cov_file.c_str());
	if (!cov.is_open())
	{
		std::string error_message = "ERROR: Cannot open covariate file: " + c_cov_file;
		throw std::runtime_error(error_message);
	}
	std::string line;
	std::getline(cov, line); // remove header
	int max_index = cov_index.back()+1;
	while(std::getline(cov, line))
	{
		misc::trim(line);
		if(line.empty()) continue;
		valid = true;
		std::vector<std::string> token = misc::split(line);
		if(token.size() < max_index)
		{
			std::string error_message = "ERROR: Malformed covariate file, should contain at least "
					+ std::to_string(max_index) + " column!";
			throw std::runtime_error(error_message);
		}
		std::string id =(m_ignore_fid)? token[0]:token[0]+"_"+token[1];
		if (m_sample_with_phenotypes.find(id) != m_sample_with_phenotypes.end())
		{ // sample is found in the phenotype vector
			int index = m_sample_with_phenotypes[id];
			for (size_t i_cov = 0; i_cov < cov_index.size(); ++i_cov)
			{
				try {
					double temp = misc::convert<double>( token[cov_index[i_cov]]);
					m_independent_variables(index, i_cov + 2) = temp; // + 2 because first line = intercept, second line = PRS
				} catch (const std::runtime_error &error) {
					valid = false;
					m_independent_variables(index, i_cov + 2) = 0;
				}
			}
			if (valid)
			{
				valid_sample_index.push_back( std::pair<std::string, size_t>(id, index));
				num_valid++;
			}
		}
	}
	// now we need to handle the situation where there are a different number of samples

	if (valid_sample_index.size() != num_sample && num_sample != 0) {
		int removed = num_sample - valid_sample_index.size();
		fprintf(stderr, "Number of samples with invalid covariate: %d\n", removed);
		double portion = (double) removed / (double) num_sample;
		if(valid_sample_index.size() == 0)
		{
			throw std::runtime_error("All samples removed due to missingness in covariate file!");
		}
		if (portion > 0.05) {
			fprintf(stderr, "WARNING! More than %03.2f%% of the samples were removed!\n", portion * 100);
			fprintf(stderr, "         Do check if your covariate file is correct\n");
		}
		std::sort(begin(valid_sample_index), end(valid_sample_index),
				[](std::pair<std::string, size_t> const &t1, std::pair<std::string, size_t> const &t2)
				{
			if(std::get<1>(t1)==std::get<1>(t2)) return std::get<0>(t1).compare(std::get<0>(t2)) < 0;
			else return std::get<1>(t1) < std::get<1>(t2);
				});

		// update the m_phenotype and m_independent
		m_sample_with_phenotypes.clear();
		for (size_t cur_index = 0; cur_index < valid_sample_index.size(); ++cur_index) {
			std::string name = std::get < 0 > (valid_sample_index[cur_index]);
			m_sample_with_phenotypes[name]  = cur_index;
			size_t update_index = std::get < 1 > (valid_sample_index[cur_index]);
			if (update_index != cur_index) {
				m_phenotype(cur_index, 0) = m_phenotype(update_index, 0);
				for (size_t i_cov = 0; i_cov < cov_index.size(); ++i_cov) {
					m_independent_variables(cur_index, i_cov + 2) = m_independent_variables(update_index, i_cov + 2);
				}
			}
		}
		m_independent_variables.conservativeResize(valid_sample_index.size(),
				m_independent_variables.cols());
		m_phenotype.conservativeResize(valid_sample_index.size(), 1);

		fprintf(stderr, "\nFinal number of samples: %zu\n\n", valid_sample_index.size());
	}
}

void PRSice::init_matrix(const Commander &c_commander, const size_t c_pheno_index, const bool prslice) {
	m_null_r2 = 0.0;

	// Clean up the matrix
	m_phenotype = Eigen::VectorXd::Zero(0);
	m_independent_variables.resize(0, 0);
	bool no_regress = c_commander.no_regression();
	bool all = c_commander.all();

	std::string pheno_file = c_commander.get_pheno();
	std::string output_name = c_commander.get_out();
	std::ofstream all_out;
	bool multi = m_pheno_names.size() > 1;
	if (all && !prslice) // we don't want this output for PRSlice
	{
		std::string all_out_name = output_name + "." + m_base_name;
		if (multi)
		{
			all_out_name.append("." + std::get < pheno_store::NAME > (m_pheno_names.at(c_pheno_index)));
		}
		all_out_name.append(".all.score");
		all_out.open(all_out_name.c_str());
		if (!all_out.is_open())
		{
			std::string error_message = "Cannot open file " + all_out_name + " for write";
			throw std::runtime_error(error_message);
		}
	}
	gen_pheno_vec(
			std::get < pheno_store::FILE_NAME > (m_pheno_names[c_pheno_index]),
			std::get < pheno_store::INDEX > (m_pheno_names[c_pheno_index]),
			std::get < pheno_store::ORDER > (m_pheno_names[c_pheno_index]),
			!no_regress);
	if (!no_regress)
	{
		gen_cov_matrix(c_commander.get_cov_file(), c_commander.get_cov_header());
	}


	if (all && !prslice) {
		all_out << "Threshold\tRegion";
		for (auto &&sample : m_sample_names)
			all_out << "\t" << std::get < +PRS::IID > (sample);
		all_out << std::endl;
		all_out.close();
	}
	double null_r2_adjust = 0.0, null_p = 0.0, null_coeff = 0.0;
	// calculate the null r2
	int n_thread = c_commander.get_thread();
	if (m_independent_variables.cols() > 2) {
		Eigen::MatrixXd covariates_only;
		covariates_only = m_independent_variables;
		covariates_only.block(0, 1, covariates_only.rows(),
				covariates_only.cols() - 2) = covariates_only.topRightCorner(
				covariates_only.rows(), covariates_only.cols() - 2);
		covariates_only.conservativeResize(covariates_only.rows(),
				covariates_only.cols() - 1);
		if (m_target_binary[c_pheno_index]) {
			Regression::glm(m_phenotype, covariates_only, null_p, m_null_r2,
					null_coeff, 25, n_thread, true);
		} else {
			Regression::linear_regression(m_phenotype, covariates_only, null_p,
					m_null_r2, null_r2_adjust, null_coeff, n_thread, true);
		}
	}
}



void PRSice::categorize(const Commander &c_commander) {
	m_partition.clear();
	bool fastscore = c_commander.fastscore();
	double bound_start = c_commander.get_lower();
	double bound_end = c_commander.get_upper();
	double bound_inter = c_commander.get_inter();
	bool full_model = c_commander.full();
	std::vector < std::string > file_names;
	if (m_target.find("#") != std::string::npos)
	{
		for (auto &&chr : m_chr_list)
		{
			std::string name = m_target;
			misc::replace_substring(name, "#", chr);
			file_names.push_back(name);
		}
	}
	else if(file_names.size()==0) // in the case where the file name really do have the # sign
	{
		file_names.push_back(m_target);
	}
	// WARNING: In someway, this should be the same as the plink read sequence because both uses the same
	// chromosome list. But then, we should always be away that it is possible for the two to be out of
	// sync.

	for (auto &&name : file_names)
	{
		std::ifstream bim;
		std::string bim_name = name + ".bim";
		bim.open(bim_name.c_str());
		if (!bim.is_open())
		{
			std::string error_message = "Cannot open bim file: " + bim_name;
			throw std::runtime_error(error_message);
		}
		std::string line;
		size_t cur_line = 0;
		while (std::getline(bim, line))
		{
			misc::trim(line);
			if (!line.empty()) {
				std::vector < std::string > token = misc::split(line);
				if (token.size() < 6)
					throw std::runtime_error( "Malformed bim file, should contain at least 6 columns");
				if (m_include_snp.find(token[+BIM::RS]) != m_include_snp.end()) {
					size_t cur_snp_index = m_include_snp[token[+BIM::RS]]; //because found, we knwo it is ok
					double p = m_snp_list[cur_snp_index].get_p_value();
					p_partition part;
					std::get < +PRS::RS > (part) = token[+BIM::RS];
					std::get < +PRS::LINE > (part) = cur_line;
					std::get < +PRS::INDEX > (part) = cur_snp_index;
					std::get < +PRS::FILENAME > (part) = name;
					int category = -1;
					if (fastscore)
					{
						category = c_commander.get_category(p);
						std::get < +PRS::CATEGORY > (part) = category;
						std::get < +PRS::P_THRES > (part) = c_commander.get_threshold(category);
						m_partition.push_back(part);
					}
					else
					{
						// calculate the threshold instead
						if (p > bound_end && full_model)
						{
							std::get < +PRS::CATEGORY > (part) = std::ceil((bound_end + 0.1 - bound_start) / bound_inter);
							std::get < +PRS::P_THRES > (part) = 1.0;
							m_partition.push_back(part);
						}
						else
						{
							category = std::ceil((p - bound_start) / bound_inter);
							category = (category < 0) ? 0 : category;
							std::get < +PRS::CATEGORY > (part) = category;
							std::get < +PRS::P_THRES > (part) = category * bound_inter + bound_start;
							m_partition.push_back(part);
						}
					}
				}
				cur_line++;
			}
		}
		bim.close();
	}
	if (m_partition.size() == 0) {
		throw std::runtime_error("None of the SNPs met the threshold\n");
	}
	std::sort(begin(m_partition), end(m_partition),
			[](p_partition const &t1, p_partition const &t2)
			{
				if(std::get<+PRS::CATEGORY>(t1)==std::get<+PRS::CATEGORY>(t2))
				{
					if(std::get<+PRS::FILENAME>(t1).compare(std::get<+PRS::FILENAME>(t2))==0)
					{
						return std::get<+PRS::LINE>(t1)<std::get<+PRS::LINE>(t2);
					}
					else return std::get<+PRS::FILENAME>(t1).compare(std::get<+PRS::FILENAME>(t2))<0;
				}
				else return std::get<+PRS::CATEGORY>(t1)<std::get<+PRS::CATEGORY>(t2);
			});
}


void PRSice::prsice(const Commander &c_commander, const Region &c_region,
		const size_t c_pheno_index, bool prslice) {
	if (m_partition.size() == 0)
	{
		throw std::runtime_error("None of the SNPs fall into the threshold\n");
	}
	// not allowed for prslice
	bool no_regress = c_commander.no_regression() && !prslice;
	bool require_all = c_commander.all() && !prslice;
	bool multi = m_pheno_names.size() > 1;
	std::ofstream all_out;
	if (require_all)
	{
		std::string all_out_name = c_commander.get_out() + "." + m_base_name;
		if (multi)
		{
			all_out_name.append("." + std::get < pheno_store::NAME > (m_pheno_names.at(c_pheno_index)));
		}
		all_out_name.append(".all.score");
		all_out.open(all_out_name.c_str(), std::ofstream::app);
		if (!all_out.is_open())
		{
			std::string error_message = "Cannot open file " + all_out_name + " for write";
			throw std::runtime_error(error_message);
		}
	}

	Eigen::initParallel();
	std::vector < std::thread > thread_store;
	size_t n_thread = c_commander.get_thread();
	m_best_threshold.clear();
	m_current_prs.clear();
	m_prs_results.clear();
	// below is also initialization
	// no need to worry about PRSlice as that should always follow
	m_num_snp_included = std::vector < size_t > (m_region_size);
	for (size_t i_region = 0; i_region < m_region_size; ++i_region)
	{
		m_current_prs.push_back(m_sample_names);
		PRSice_best cur_best;
		std::get<+PRS::THRESHOLD>(cur_best)=0;
		std::get<+PRS::R2>(cur_best)=0;
		std::get<+PRS::NSNP>(cur_best)=0;
		std::get<+PRS::COEFF>(cur_best)=0;
		std::get<+PRS::P>(cur_best)=0;
		std::get<+PRS::EMPIRICAL_P>(cur_best)=-1;
		m_best_threshold.push_back(cur_best);
		m_prs_results.push_back(std::vector < PRSice_result > (0));
	}

	m_best_score = m_current_prs;

	// avoid 100% before complete
	int max_category = std::get < +PRS::CATEGORY > (m_partition.back()) + 1;

	// With the new PLINK class, we should start the plink here instead of initializing it
	// in the get_prs_score
	// rewind should work like magic?
	PLINK score_plink(m_target, false, c_commander.get_thread());

	size_t partition_size = m_partition.size();
	double cur_threshold = 0.0;
	size_t cur_partition_index = 0;
	bool reg = false;
	while (cur_partition_index < partition_size)
	{
		int cur_category = std::get < +PRS::CATEGORY > (m_partition[cur_partition_index]);
		cur_threshold = std::get< +PRS::P_THRES > (m_partition[cur_partition_index]);
		if (!prslice)
			fprintf(stderr, "\rProcessing %03.2f%%", (double) cur_category / (double) (max_category) * 100.0);

		reg = get_prs_score(cur_partition_index, score_plink);
		if (require_all && all_out.is_open()) {
			for (size_t i_region = 0; i_region < m_region_size; ++i_region)
			{
				all_out << cur_threshold << "\t" << c_region.get_name(i_region);
				for (auto &&prs : m_current_prs[i_region])
				{
					all_out << "\t" << std::get < +PRS::PRS > (prs) / (double) std::get < +PRS::NNMISS > (prs);
				}
				all_out << std::endl;
			}
		}
		reg = reg && !no_regress;
		if (reg)
		{
			if (n_thread == 1 || m_region_size == 1)
			{
				thread_score(0, m_region_size, cur_threshold, n_thread, c_pheno_index);
			}
			else
			{
				if (m_region_size < n_thread)
				{
					for (size_t i_region = 0; i_region < m_region_size; ++i_region)
					{
						thread_store.push_back( std::thread(&PRSice::thread_score, this,
								i_region, i_region + 1, cur_threshold,
								1, c_pheno_index));
					}
				}
				else
				{
					int job_size = m_region_size / n_thread;
					int remain = m_region_size % n_thread;
					size_t start = 0;
					for (size_t i_thread = 0; i_thread < n_thread; ++i_thread)
					{
						size_t ending = start + job_size + (remain > 0);
						ending = (ending > m_region_size) ? m_region_size : ending;
						thread_store.push_back( std::thread(&PRSice::thread_score, this, start,
								ending, cur_threshold, 1,
								c_pheno_index));
						start = ending;
						remain--;
					}
				}
				// joining the threads
				for (auto &&thread : thread_store) thread.join();
				thread_store.clear();
			}
		}
		score_plink.clear();
	}
	if (all_out.is_open()) all_out.close();
	if (!prslice) fprintf(stderr, "\rProcessing %03.2f%%\n", 100.0);
}

bool PRSice::get_prs_score(size_t &cur_index, PLINK &score_plink)
{
	if (m_partition.size() == 0) return false; // nothing to do
	int prev_index = std::get < +PRS::CATEGORY > (m_partition[cur_index]);
	int end_index = 0;
	bool ended = false;
	for (size_t i = cur_index; i < m_partition.size(); ++i)
	{
		if (std::get < +PRS::CATEGORY > (m_partition[i]) != prev_index
				&& std::get < +PRS::CATEGORY > (m_partition[i]) >= 0)
		{
			end_index = i;
			ended = true;
			break;
		}
		else if (std::get < +PRS::CATEGORY > (m_partition[i]) != prev_index)
		{
			prev_index = std::get < +PRS::CATEGORY > (m_partition[i]); // only when the category is still negative
		}
//		// Use as part of the output
		for (size_t i_region = 0; i_region < m_region_size; ++i_region)
		{
			if (m_snp_list[std::get < +PRS::INDEX > (m_partition[i])].in( i_region)) m_num_snp_included[i_region]++;
		}
	}
	if (!ended) end_index = m_partition.size();
	score_plink.get_score(m_partition, m_snp_list, m_current_prs, cur_index, end_index, m_region_size, m_score);

	cur_index = end_index;
	return true;
}

void PRSice::thread_score(size_t region_start, size_t region_end,
		double threshold, size_t thread, const size_t c_pheno_index)
{

	Eigen::MatrixXd X;
	bool thread_safe = false;
	if (region_start == 0 && region_end == m_current_prs.size())
		thread_safe = true;
	else
		X = m_independent_variables;
	double r2 = 0.0, r2_adjust = 0.0, p_value = 0.0, coefficient = 0.0;
	for (size_t iter = region_start; iter < region_end; ++iter)
	{
		// The m_prs size check is just so that the back will be valid
		// m_prs will only be empty for the first run
		if (m_num_snp_included[iter] == 0 ||
				(m_prs_results[iter].size() != 0 &&
						m_num_snp_included[iter] == std::get < +PRS::NSNP > (m_prs_results[iter].back())
				)
		) continue; // don't bother when there is no additional SNPs added
		for (auto &&prs : m_current_prs[iter])
		{
			std::string sample = (m_ignore_fid)? std::get <+PRS::IID > (prs):
					std::get<+PRS::FID>(prs)+"_"+std::get<+PRS::IID>(prs);
			// The reason why we need to udpate the m_sample_with_phenotypes matrix
			if (m_sample_with_phenotypes.find(sample) != m_sample_with_phenotypes.end()) {

				if(std::get < +PRS::NNMISS >(prs)!= 0)
				{
					if (thread_safe)
						m_independent_variables(m_sample_with_phenotypes.at(sample), 1) = std::get < +PRS::PRS > (prs) / (double) std::get < +PRS::NNMISS >(prs) ;
					else
						X(m_sample_with_phenotypes.at(sample), 1) = std::get < +PRS::PRS > (prs) / (double) std::get < +PRS::NNMISS >(prs);
				}
				else{
					if (thread_safe)
						m_independent_variables(m_sample_with_phenotypes.at(sample), 1) = 0 ;
					else
						X(m_sample_with_phenotypes.at(sample), 1) =0;
				}
			}
		}
		int num_better = -1;
		if (m_target_binary[c_pheno_index])
		{
			try {
				if (thread_safe)
					Regression::glm(m_phenotype, m_independent_variables, p_value, r2, coefficient, 25, thread, true);
				else
					Regression::glm(m_phenotype, X, p_value, r2, coefficient, 25, thread, true);
			} catch (const std::runtime_error &error) {
				// This should only happen when the glm doesn't converge.
				// Let's hope that won't happen...
				fprintf(stderr, "ERROR: GLM model did not converge!\n");
				fprintf(stderr, "       Please send me the DEBUG files\n");
				std::ofstream debug;
				debug.open("DEBUG");
				if (thread_safe)
					debug << m_independent_variables << std::endl;
				else
					debug << X << std::endl;
				debug.close();
				debug.open("DEBUG.y");
				debug << m_phenotype << std::endl;
				debug.close();
				fprintf(stderr, "ERROR: %s\n", error.what());
				exit(-1);
			}
		} else {
			if (thread_safe)
				Regression::linear_regression(m_phenotype, m_independent_variables, p_value, r2, r2_adjust,
						coefficient, thread, true);
			else
				Regression::linear_regression(m_phenotype, X, p_value, r2, r2_adjust, coefficient, thread, true);
		}


		// It this is the best r2, then we will add it
		if (std::get < +PRS::R2 > (m_best_threshold[iter]) < r2) {
			num_better = 0;
			if(m_target_binary[c_pheno_index]){
				Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm( m_phenotype.rows());
				for (size_t i_perm = 0; i_perm < m_perm; ++i_perm)
				{
					perm.setIdentity();
					std::random_shuffle(perm.indices().data(), perm.indices().data() + perm.indices().size());
					Eigen::MatrixXd A_perm = perm * m_phenotype; // permute columns
					double perm_p, perm_r2, perm_coefficient;
					if (thread_safe)
						Regression::glm(A_perm, m_independent_variables, perm_p, perm_r2, perm_coefficient, 25, thread, true);
					else
						Regression::glm(A_perm, X, perm_p, perm_r2, perm_coefficient, 25, thread, true);
					if (perm_p < p_value)num_better++;
				}
			}
			else{
				Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm( m_phenotype.rows());
				for (size_t i_perm = 0; i_perm < m_perm; ++i_perm) {
					perm.setIdentity();
					std::random_shuffle(perm.indices().data(), perm.indices().data() + perm.indices().size());
					Eigen::MatrixXd A_perm = perm * m_phenotype; // permute columns
					double perm_p, perm_r2, perm_coefficient, perm_r2_adj;
					if (thread_safe)
						Regression::linear_regression(A_perm, m_independent_variables, perm_p, perm_r2,
								perm_r2_adj, perm_coefficient, thread, true);
					else
						Regression::linear_regression(A_perm, X, perm_p, perm_r2, perm_r2_adj, perm_coefficient, thread, true);
					if (perm_p < p_value) num_better++;
				}
			}



			PRSice_best best;
			std::get < +PRS::THRESHOLD > (best) = threshold;
			std::get < +PRS::R2 > (best) = r2;
			std::get < +PRS::NSNP > (best) = m_num_snp_included[iter];
			std::get < +PRS::COEFF > (best) = coefficient;
			std::get < +PRS::P > (best) = p_value;
			std::get < +PRS::EMPIRICAL_P > (best) = num_better;
			m_best_threshold[iter] = best;
			m_best_score[iter] = m_current_prs.at(iter);
		}
		// This should be thread safe as each thread will only mind their own region
		// now add the PRS result to the vectors (hopefully won't be out off scope
		PRSice_result res;
		std::get < +PRS::THRESHOLD > (res) = threshold;
		std::get < +PRS::R2 > (res) = r2;
		std::get < +PRS::NSNP > (res) = m_num_snp_included[iter];
		std::get < +PRS::R2ADJ > (res) = r2_adjust;
		std::get < +PRS::P > (res) = p_value;
		std::get < +PRS::COEFF > (res) = coefficient;
		std::get < +PRS::EMPIRICAL_P > (res) = num_better;
		m_prs_results[iter].push_back(res);
	}
}

void PRSice::output(const Commander &c_commander, const Region &c_region,
		size_t pheno_index) const {
	// this is ugly, need to make it better
	std::string pheno_name = std::get < pheno_store::NAME > (m_pheno_names[pheno_index]);
	std::string output_prefix = c_commander.get_out() + "." + m_base_name;
	if (!pheno_name.empty()) output_prefix.append("." + pheno_name);
	size_t total_perm = c_commander.get_perm();
	bool perm = total_perm > 0;
	std::string output_name = output_prefix;
	for(size_t i_region = 0; i_region < m_region_size; ++i_region)
	{
		if(m_region_size > 1) output_name = output_prefix+"."+c_region.get_name(i_region);
		std::string out_best = output_name + ".best";
		std::string out_prsice = output_name + ".prsice";
		std::string out_snp = output_name +".snps";
		std::ofstream best_out, prsice_out, snp_out;
		prsice_out.open(out_prsice.c_str());
		if (!prsice_out.is_open())
		{
			std::string error_message = "ERROR: Cannot open file: " + out_prsice + " to write";
			throw std::runtime_error(error_message);
		}
		prsice_out << "Threshold\tR2\tP\tCoefficient\tNum_SNP";
		if (perm) prsice_out << "\tEmpirical_P";
		prsice_out << std::endl;
		for (auto &&prs : m_prs_results[i_region]) {
			prsice_out << std::get < +PRS::THRESHOLD > (prs) << "\t"
					<< std::get < +PRS::R2 > (prs) - m_null_r2 << "\t"
					<< std::get < +PRS::P > (prs) << "\t"
					<< std::get < +PRS::COEFF > (prs) << "\t"
					<< std::get < +PRS::NSNP > (prs);
			if (perm) prsice_out << "\t" << (double) (std::get < +PRS::EMPIRICAL_P > (prs) + 1.0) / (double) (total_perm + 1.0);
			prsice_out << std::endl;
		}
		prsice_out.close();

		best_out.open(out_best.c_str());
		if (!best_out.is_open())
		{
			std::string error_message = "ERROR: Cannot open file: " + out_best + " to write";
			throw std::runtime_error(error_message);
		}
		best_out << "FID\tIID\tIncluded\tprs_" << std::get < +PRS::THRESHOLD > (m_best_threshold[i_region]) << std::endl;
		int best_snp_size = std::get < +PRS::NSNP > (m_best_threshold[i_region]);
		if (best_snp_size == 0) {
			fprintf(stderr, "ERROR: Best R2 obtained when no SNPs were included\n");
			fprintf(stderr, "       Cannot output the best PRS score\n");
		} else {
			for (auto &&prs : m_best_score[i_region]) {
				std::string id = (m_ignore_fid)? std::get < +PRS::IID > (prs) :
						std::get<+PRS::FID>(prs)+"_"+std::get < +PRS::IID > (prs);
				char in = (m_sample_with_phenotypes.find(id)==m_sample_with_phenotypes.end())? 'N':'Y';
				best_out << std::get<+PRS::FID>(prs) << "\t"
						<< std::get < +PRS::IID > (prs) << "\t"
						<< in << "\t"
						<< std::get< +PRS::PRS> (prs) / (double) best_snp_size
						<< std::endl;
			}
		}
		best_out.close();
		if(c_commander.print_snp())
		{
			snp_out.open(out_snp);
			if (!snp_out.is_open()) {
				std::string error_message = "ERROR: Cannot open file: " + out_snp + " to write";
				throw std::runtime_error(error_message);
			}
			size_t num_snp = std::get < +PRS::NSNP> (m_best_threshold[i_region]);
			for(auto snp : m_partition)
			{
				if(m_snp_list.at(m_include_snp.at(std::get< +PRS::RS >(snp))).in(i_region))
				{
						snp_out  << std::get< +PRS::RS >(snp) << std::endl;;
				}
			}
			snp_out.close();
		}
		if(!c_commander.print_all()) break;
	}

	if(m_region_size > 1)
	{
		// now print the group information
		std::string out_region = output_prefix + ".prset";
		std::ofstream region_out;
		region_out.open(out_region.c_str());
		region_out << "Region\tThreshold\tR2\tCoefficient\tP\tNum_SNP";
		if (perm) region_out << "\tEmpirical_P";
		region_out	<< std::endl;
		size_t i_region = 0;
		for (auto best_region : m_best_threshold) {
			region_out << c_region.get_name(i_region) << "\t" <<
					std::get< +PRS::THRESHOLD > (best_region) << "\t"
					<< std::get< +PRS::R2 > (best_region) - m_null_r2 << "\t"
					<< std::get < +PRS::COEFF > (best_region) << "\t"
					<< std::get < +PRS::P> (best_region) << "\t"
					<< std::get < +PRS::NSNP> (best_region);
			if(perm) region_out << "\t" << (double) (std::get < +PRS::EMPIRICAL_P > (best_region) + 1.0) / (double) (total_perm + 1.0);
			region_out << std::endl;
			i_region++;
		}
		region_out.close();
	}
}

PRSice::~PRSice() {
	//dtor
}

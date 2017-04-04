#include "prsice.hpp"

std::mutex PRSice::score_mutex;

std::vector<SNP> read_base(unordered_map<std::string, int> &snp_index,
		const Commander &c_commander, const std::vector<int32_t> &xymt_codes,
		const uint32_t max_code)
{
	// do region later on
	const std::string input = c_commander.base_name();
	const bool beta = c_commander.beta();
	const bool fastscore = c_commander.fastscore();
	const bool full = c_commander.full();
	std::vector<int> index = c_commander.index();
	std::ifstream snp_file;
	snp_file.open(input.c_str());
	if(!snp_file.is_open())
	{
		std::string error_message = "ERROR: Cannot open base file: " +input;
		throw std::runtime_error(error_message);
	}
	int max_index = index[+BASE_INDEX::MAX]; // a struct will actually be better, but nvm
	std::string line;
	if (!c_commander.has_index()) std::getline(snp_file, line);

	// category related stuff
	double threshold = (c_commander.fastscore())? c_commander.bar_upper() : c_commander.upper();
	double bound_start = c_commander.lower();
	double bound_end = c_commander.upper();
	double bound_inter = c_commander.inter();
	threshold = (full)? 1.0 : threshold;
	std::vector < std::string > token;

	bool exclude = false;
	// Some QC countss
	size_t num_duplicated = 0;
	size_t num_excluded = 0;
	size_t num_stat_not_converted = 0; // this is for NA
	size_t num_p_not_converted = 0; // this is for NA
	size_t num_negative_stat = 0;
	size_t num_ambiguous = 0;
	std::unordered_set<std::string> dup_index;
	std::vector<SNP> base_snps;
	while (std::getline(snp_file, line))
	{
		misc::trim(line);
		if (line.empty()) continue;

		exclude = false;
		token = misc::split(line);
		if (token.size() <= max_index)
			throw std::runtime_error("More index than column in data");
		std::string rs_id = token[index[+BASE_INDEX::RS]];
		if(dup_index.find(rs_id)==dup_index.end())
		{
			dup_index.insert(rs_id);
			int32_t chr_code = -1;
			if (index[+BASE_INDEX::CHR] >= 0)
			{
				chr_code = get_chrom_code_raw(token[index[+BASE_INDEX::CHR]].c_str());
				if (((const uint32_t)chr_code) > max_code) {
					if (chr_code != -1) {
						if (chr_code >= MAX_POSSIBLE_CHROM) {
							chr_code= xymt_codes[chr_code - MAX_POSSIBLE_CHROM];
						}
						else
						{
							std::string error_message ="ERROR: Cannot parse chromosome code: "
									+ token[index[+BASE_INDEX::CHR]];
							throw std::runtime_error(error_message);
						}
					}
				}
			}
			std::string ref_allele = (index[+BASE_INDEX::REF] >= 0) ?
					token[index[+BASE_INDEX::REF]] : "";
			std::string alt_allele = (index[+BASE_INDEX::ALT] >= 0) ?
					token[index[+BASE_INDEX::ALT]] : "";
			int loc = -1;
			if (index[+BASE_INDEX::BP] >= 0)
			{
				// obtain the SNP coordinate
				try {
					loc = misc::convert<int>( token[index[+BASE_INDEX::BP]].c_str());
					if (loc < 0)
					{
						std::string error_message = "ERROR: "+rs_id+" has negative loci!\n";
						throw std::runtime_error(error_message);
					}
				} catch (const std::runtime_error &error) {
					std::string error_message = "ERROR: Non-numeric loci for "+rs_id+"!\n";
					throw std::runtime_error(error_message);
				}
			}
			double pvalue = 2.0;
			try{
				misc::convert<double>( token[index[+BASE_INDEX::P]]);
				if (pvalue < 0.0 || pvalue > 1.0)
				{
					std::string error_message = "ERROR: Invalid p-value for "+rs_id+"!\n";
					throw std::runtime_error(error_message);
				}
				else if (pvalue > threshold)
				{
					exclude = true;
					num_excluded++;
				}
			}catch (const std::runtime_error& error) {
				exclude = true;
				num_p_not_converted++;
			}
			double stat = 0.0;
			try {
				stat = misc::convert<double>( token[index[+BASE_INDEX::STAT]]);
				if(stat <0 && !beta)
				{
					num_negative_stat++;
					exclude = true;
				}
				else if (!beta) stat = log(stat);
			} catch (const std::runtime_error& error) {
				num_stat_not_converted++;
				exclude = true;
			}
			if(!alt_allele.empty() && SNP::ambiguous(ref_allele, alt_allele)){
				num_ambiguous++;
				exclude= true;
			}
			if(!exclude)
			{
				int category = -1;
				double pthres = 0.0;
				if (fastscore)
				{
					category = c_commander.get_category(pvalue);
					pthres = c_commander.get_threshold(category);
				}
				else
				{
					// calculate the threshold instead
					if (pvalue > bound_end && full)
					{
						category = std::ceil((bound_end + 0.1 - bound_start) / bound_inter);
						pthres = 1.0;
					}
					else
					{
						category = std::ceil((pvalue - bound_start) / bound_inter);
						category = (category < 0) ? 0 : category;
						pthres = category * bound_inter + bound_start;
					}
				}
				// ignore the SE as it currently serves no purpose
				snp_index[rs_id] = base_snps.size();
				base_snps.push_back(SNP(rs_id, chr_code, loc, ref_allele, alt_allele, stat, 0,
						pvalue, category, pthres));
			}
		}
		else
		{
			num_duplicated++;
		}
	}
	snp_file.close();
	if (num_excluded != 0) fprintf(stderr, "Number of SNPs excluded: %zu\n", num_excluded);
	if (num_ambiguous != 0) fprintf(stderr, "%zu ambiguous SNPs\n", num_ambiguous);
	if (num_duplicated != 0) fprintf(stderr, "Number of duplicated SNPs : %d\n", num_duplicated);
	if (num_stat_not_converted != 0) fprintf(stderr, "Failed to convert %zu OR/beta\n", num_stat_not_converted);
	if (num_p_not_converted != 0) fprintf(stderr, "Failed to convert %zu p-value\n", num_p_not_converted);
	if (num_negative_stat != 0) fprintf(stderr, "%zu negative OR\n", num_negative_stat);
		fprintf(stderr, "Number of SNPs from base  : %zu\n", base_snps.size());
	return base_snps;
}

void PRSice::get_snp(const Commander &c_commander, Region &region) {
	const std::string input = c_commander.get_base(m_base_index);
	const bool beta = c_commander.get_base_binary(m_base_index);
	// just issue the warning. would not terminate
	if (beta && c_commander.statistic().compare("OR") == 0) fprintf(stderr, "WARNING: OR detected but user suggest the input is beta!\n");
	// First, we need to obtain the index of different columns from the base files
	// NOTE: -1 means missing and index is hard coded such that each index will represent
	//       specific header
	std::vector<int> index = SNP::get_index(c_commander, input);
	if ((c_commander.get_target().find("#") != std::string::npos && index[+SNP_Index::CHR] < 0)
			|| (c_commander.ld_prefix().find("#") != std::string::npos && index[+SNP_Index::CHR] < 0))
	{
		std::string error_message = "To use chromosome separated PLINK input, you must provide"
				" the CHR header as we use the CHR information form the base file to substitute #";
		throw std::runtime_error(error_message);
	}
	else if (region.size() > 1 && (index[+SNP_Index::CHR] < 0 || index[+SNP_Index::BP] < 0))
	{
		std::string error_message = "To perform PRSet, you must provide the CHR and LOC header such"
				" that we can determine the set membership";
		throw std::runtime_error(error_message);
	}
	else if (c_commander.prslice() > 0.0 && (index[+SNP_Index::CHR] < 0 || index[+SNP_Index::BP] < 0))
	{
		std::string error_message = "To perform PRSlice, you must provide the CHR and LOC header such"
				" that we can perform the slicing";
		throw std::runtime_error(error_message);
	}
	// Open the file
	std::ifstream snp_file;
	snp_file.open(input.c_str());
	if (!snp_file.is_open())
	{
		std::string error_message = "Cannot open base file: " + input;
		throw std::runtime_error(error_message);
	}
	// Some QC countss
	int num_duplicated = 0;
	size_t num_stat_not_convertible = 0;
	size_t num_p_not_convertible = 0;
	size_t num_indel = 0;
	size_t num_se_not_convertible = 0;
	size_t num_exclude = 0;
	int max_index = index[+SNP_Index::MAX];
	std::string line;
	// remove header if index is not provided
	if (!c_commander.index()) std::getline(snp_file, line);
	bool read_error = false;
	bool not_converted = false;
	bool exclude = false;
	bool chr_error =false;
	// category related stuff
	double threshold = (c_commander.fastscore())? c_commander.get_bar_upper() : c_commander.get_upper();
	threshold = (c_commander.full())? 1.0 : threshold;
	std::vector < std::string > token;
	// Actual reading the file, will do a bunch of QC
	while (std::getline(snp_file, line))
	{
		misc::trim(line);
		if (!line.empty())
		{
			not_converted = false;
			exclude = false;
			token = misc::split(line);
			if (token.size() <= max_index)
				throw std::runtime_error("More index than column in data");
			else
			{
				std::string rs_id = token[index[+SNP_Index::RS]];
				std::string chr = (index[+SNP_Index::CHR] >= 0)? token[index[+SNP_Index::CHR]] : "";
				int32_t chr_code = get_chrom_code_raw(chr.c_str());
				if (((const uint32_t)chr_code) > PLINK::g_max_code) {
					if (chr_code != -1) {
						if (chr_code >= MAX_POSSIBLE_CHROM) {
							chr_code= PLINK::g_xymt_codes[chr_code - MAX_POSSIBLE_CHROM];
						}
						else
						{
							std::string error_message ="ERROR: Cannot parse chromosome code: " + chr;
							throw std::runtime_error(error_message);
						}
					}
				}
				if(chr_code == (uint32_t)PLINK::g_xymt_codes[X_OFFSET] ||
						chr_code == (uint32_t)PLINK::g_xymt_codes[Y_OFFSET])
				{
					exclude= true;
					if(!chr_error)
					{
						fprintf(stderr, "WARNING: Sex chromosome detected. They will be ignored\n");
						chr_error = true;
					}
				}
				std::string ref_allele = (index[+SNP_Index::REF] >= 0) ? token[index[+SNP_Index::REF]] : "";
				std::string alt_allele = (index[+SNP_Index::ALT] >= 0) ? token[index[+SNP_Index::ALT]] : "";
				double pvalue = 0.0;
				if (index[+SNP_Index::P] >= 0)
				{
					try {
						// obtain the p-value and calculate the corresponding threshold & category
						pvalue = misc::convert<double>(
								token[index[+SNP_Index::P]]);
						if (pvalue < 0.0 || pvalue > 1.0)
						{
							read_error = true;
							fprintf(stderr, "ERROR: %s's p-value is %f\n", rs_id.c_str(), pvalue);
						}
						else if (pvalue > threshold)
						{
							exclude = true;
							num_exclude++;
						}
					} catch (const std::runtime_error &error) {
						num_p_not_convertible++;
						not_converted = true;
					}
				}
				double stat = 0.0;
				if (index[+SNP_Index::STAT] >= 0)
				{
					//Obtain the test statistic
					try {
						stat = misc::convert<double>( token[index[+SNP_Index::STAT]]);
						if (!beta) stat = log(stat);
					} catch (const std::runtime_error& error) {
						num_stat_not_convertible++;
						not_converted = true;
					}
				}
				double se = 0.0;
				if (index[+SNP_Index::SE] >= 0)
				{
					// obtain the standard error (though it is currently useless)
					try {
						se = misc::convert<double>( token[index[+SNP_Index::SE]]);
					} catch (const std::runtime_error &error) {
						num_se_not_convertible++;
					}
				}
				int loc = -1;
				if (index[+SNP_Index::BP] >= 0)
				{
					// obtain the SNP coordinate
					try {
						int temp = misc::convert<int>( token[index[+SNP_Index::BP]].c_str());
						if (temp < 0)
						{
							read_error = true;
							fprintf(stderr, "ERROR: %s has negative loci\n", rs_id.c_str());
						}
						else loc = temp;
					} catch (const std::runtime_error &error) { }
				}
				if (!SNP::valid_snp(ref_allele))
				{
					num_indel++;
				}
				else if (!alt_allele.empty() && !SNP::valid_snp(alt_allele)) num_indel++;
				else if(!alt_allele.empty() && SNP::ambiguous(ref_allele, alt_allele)){
					num_exclude++;
				}
				else if (!not_converted && !exclude)
				{
					m_snp_list.push_back( new SNP(rs_id, chr, loc, ref_allele, alt_allele, stat, se, pvalue));
				} else {
//		    			We skip any SNPs with non-convertible stat and p-value as we don't know how to
//		    			handle them. Most likely those will be NA, which should be ignored anyway
				}
			}
		}
		if (read_error) throw std::runtime_error( "Please check if you have the correct input");
	}
	snp_file.close();

	m_snp_list.sort();

	// The problem of this filtering is that if the same SNP is duplicated but with different
	// p-value etc, they will not be removed
	size_t before = m_snp_list.size();
	m_snp_list.erase(std::unique(m_snp_list.begin(), m_snp_list.end()), m_snp_list.end());
	size_t after = m_snp_list.size();

	std::unordered_set<std::string> unique_chr;
	std::unordered_set<std::string> unique_rsid;

	for (size_t i_snp = 0; i_snp < m_snp_list.size(); ++i_snp)
	{
		m_snp_index[m_snp_list[i_snp].get_rs_id()] = i_snp;
		std::string cur_chr = m_snp_list[i_snp].get_chr();
		if (unique_chr.find(cur_chr) == unique_chr.end())
		{
			unique_chr.insert(cur_chr);
			m_chr_list.push_back(cur_chr);
		}
		std::string rs = m_snp_list[i_snp].get_rs_id();
		if(unique_rsid.find(rs)==unique_rsid.end())
		{
			unique_rsid.insert(rs);
		}
		else
		{
			std::string error_message = "WARNING: Duplicated SNP ID: " + rs;
			throw std::runtime_error(error_message);
		}
		if (index[+SNP_Index::CHR] >= 0 && index[+SNP_Index::BP] >= 0) // if we have chr and bp information
		{
			m_snp_list[i_snp].set_flag( region.check(cur_chr, m_snp_list[i_snp].get_loc()));
		}
	}

	num_duplicated = (int) before - (int) after;
	if (num_indel != 0) fprintf(stderr, "Number of invalid SNPs    : %zu\n", num_indel);
	if (num_exclude != 0) fprintf(stderr, "Number of SNPs excluded: %zu\n", num_exclude);
	if (num_duplicated != 0) fprintf(stderr, "Number of duplicated SNPs : %d\n", num_duplicated);
	if (num_stat_not_convertible != 0) fprintf(stderr, "Failed to convert %zu OR/beta\n", num_stat_not_convertible);
	if (num_p_not_convertible != 0) fprintf(stderr, "Failed to convert %zu p-value\n", num_p_not_convertible);
	if (num_se_not_convertible != 0) fprintf(stderr, "Failed to convert %zu SE\n", num_se_not_convertible);
	fprintf(stderr, "Final Number of SNPs from base  : %zu\n", m_snp_list.size());
	PLINK::set_chromosome(m_chr_list); // update the chromosome information for PLINK
}



void PRSice::perform_clump(const Commander &c_commander) {
	// Main reason why we don't read this earlier is because of allele matching
	// otherwise, in theory, inclusion only need to be read once
	// if we don't want clumping, then the ld file is useless
	bool has_ld = !c_commander.ld_prefix().empty() && c_commander.no_clump();
	std::string ld_file = (has_ld) ? c_commander.ld_prefix() : m_target;

	size_t num_ambig = 0, not_found = 0, num_duplicate = 0;
	std::string target_bim_name = m_target + ".bim";
	if (target_bim_name.find("#") != std::string::npos) {
		for (auto &&chr : m_chr_list) {
			std::string target_chr_bim_name = target_bim_name;
			misc::replace_substring(target_chr_bim_name, "#", chr);
			check_inclusion(target_chr_bim_name, num_ambig, not_found, num_duplicate);
		}
	} else check_inclusion(target_bim_name, num_ambig, not_found, num_duplicate);


	fprintf(stderr, "\nIn Target File\n");
	fprintf(stderr, "==============================\n");

	if (has_ld) {
		fprintf(stderr, "\nIn LD Reference %s\n", ld_file.c_str());
		fprintf(stderr, "==============================\n");
		if (ld_file.find("#") != std::string::npos) {
			for (auto &&chr : m_chr_list) {
				std::string ld_chr_bim_name = ld_file;
				misc::replace_substring(ld_chr_bim_name, "#", chr);
				check_inclusion(ld_chr_bim_name, num_ambig, not_found, num_duplicate);
			}
		} else check_inclusion(target_bim_name, num_ambig, not_found, num_duplicate);
	}

	if (num_ambig != 0) fprintf(stderr, "Number of ambiguous SNPs      : %zu\n", num_ambig);
	if (num_duplicate != 0) fprintf(stderr, "Number of duplicated SNPs     : %zu\n", num_duplicate);
	fprintf(stderr, "Final Number of SNPs included : %zu\n", m_include_snp.size());


	// Thing is, this might be the LD file, which only need to be read once
	PLINK clump(ld_file, true, c_commander.get_thread(), m_include_snp);
	if (!c_commander.no_clump()) {
		fprintf(stderr, "\nStart performing clumping\n");

		clump.start_clumping(m_include_snp, m_snp_list,
				c_commander.get_clump_p(), c_commander.get_clump_r2(),
				c_commander.get_clump_kb(), c_commander.get_proxy());

	}
}

// Currently does not check indel. Might want to also check it?
void PRSice::check_inclusion(const std::string &c_target_bim_name,
		size_t &num_ambig, size_t &not_found, size_t &num_duplicate) {
	std::ifstream target_file;
	target_file.open(c_target_bim_name.c_str());
	if (!target_file.is_open()) {
		std::string error_message = "Cannot open target bim file: " + c_target_bim_name;
		throw std::runtime_error(error_message);
	}

	std::string line;
	size_t num_line = 0;
	size_t num_diff = 0;
	while (std::getline(target_file, line)) {
		misc::trim(line);
		if (!line.empty()) {
			std::vector < std::string > token = misc::split(line);
			if (token.size() < 6)
				throw std::runtime_error( "Malformed bim file. Should contain at least 6 column");
			std::string chr = token[+BIM::CHR];
			std::string rsid = token[+BIM::RS];
			int loc = -1;
			int temp = 0;
			try {
				temp = misc::convert<int>(token[+BIM::BP]);
				if (temp < 0) {
					std::string error_message = "Negative coordinate of SNP in " + c_target_bim_name;
					throw std::runtime_error(error_message);
				}
				loc = temp;
			} catch (std::runtime_error &error) {
				std::string error_message = "Non-numeric coordinate of SNP in " + c_target_bim_name;
				throw std::runtime_error(error_message);
			}
			std::string ref_allele = token[+BIM::A1];
			std::string alt_allele = token[+BIM::A2];

			if (m_include_snp.find(rsid) == m_include_snp.end() && m_snp_index.find(rsid) != m_snp_index.end()) {
				// will do some soft checking, will issue warning if there are any problem
				// first check if ambiguous
				if (SNP::ambiguous(ref_allele, alt_allele) || SNP::ambiguous(alt_allele, ref_allele)) {
					num_ambig++;
				} else {
					// not ambiguous, now do soft checking
					size_t index = m_snp_index.at(rsid);
					if (loc != -1 && m_snp_list[index].get_loc() != -1) {
						bool same = m_snp_list[index].check_loc(chr, loc, ref_allele, alt_allele);
						if (!same) {
							num_diff++;
						}
					}
					m_include_snp[rsid] = index;
				}
			}
			else if (m_include_snp.find(rsid) != m_include_snp.end())  num_duplicate++;
			else not_found++;
			num_line++;
		}
	}
	target_file.close();
	double portion = (double)num_diff/(double)num_line;
	if(num_diff > 0) fprintf(stderr, "WARNING: %zu snp has different information between target and base file\n", num_diff);
	if(portion > 0.05){
		fprintf(stderr, "         This account for %03.2f%% of all SNPs\n", portion);
		fprintf(stderr, "         It is strongly advised that you check the files are \n");
		fprintf(stderr, "         from the same genome build\n");
	}
}

void PRSice::pheno_check(const Commander &c_commander) {
	std::vector < std::string > pheno_header = c_commander.get_pheno_col();
	std::string pheno_file = c_commander.get_pheno();
	if (pheno_header.size() != 0 && pheno_file.empty()) {
		throw std::runtime_error( "You must provide a phenotype file for multiple phenotype analysis");
	}
	if (pheno_file.empty()) {
		std::string fam = m_target + ".fam";
		if(fam.find("#")!= std::string::npos)
		{
			misc::replace_substring(fam, "#", m_chr_list.front());
		}
		pheno_storage temp;
		std::get < pheno_store::FILE_NAME > (temp) = fam;
		std::get < pheno_store::INDEX > (temp) = +FAM::PHENOTYPE;
		std::get < pheno_store::NAME > (temp) = "";
		std::get < pheno_store::ORDER > (temp) = 0;
		m_pheno_names.push_back(temp);
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
			pheno_storage temp;
			std::get < pheno_store::FILE_NAME > (temp) = pheno_file;
			std::get < pheno_store::INDEX > (temp) = 1+!m_ignore_fid; // should +1 if fid is here
			std::get < pheno_store::NAME > (temp) = "";
			std::get < pheno_store::ORDER > (temp) = 0;
			m_pheno_names.push_back(temp);
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
							pheno_storage temp;
							std::get < pheno_store::FILE_NAME > (temp) = pheno_file;
							std::get < pheno_store::INDEX > (temp) = i_column;
							std::get < pheno_store::NAME > (temp) = pheno_header[i_pheno];
							std::get < pheno_store::ORDER > (temp) = i_pheno;
							m_pheno_names.push_back(temp);
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
	fprintf(stderr, "There are a total of %zu phenotype to process\n",
			m_pheno_names.size());
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


void PRSice::init_matrix(const Commander &c_commander, const size_t c_pheno_index, const bool prslice) {
	m_null_r2 = 0.0;
	m_sample_names.clear();
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



void PRSice::gen_pheno_vec(const std::string c_pheno, const int pheno_index,
		const int col_index, bool regress) {
	std::vector<double> phenotype_store;
	std::ifstream pheno_file;
	std::string fam_name = m_target + ".fam";
	bool binary = m_target_binary.at(col_index);
	if (fam_name.find("#") != std::string::npos && m_chr_list.size() > 0) // prevent undefine behaviour when chr_list ==0
	{
		misc::replace_substring(fam_name, "#", m_chr_list.front());
	}
	std::string line;
	size_t cur_index = 0;
	size_t num_case = 0, num_control = 0;
	size_t n_not_found = 0;
	size_t max_num = 0;

	if (c_pheno.empty() || fam_name == c_pheno) // the phenotype file = fam file
	{
		pheno_file.open(fam_name.c_str());
		if (!pheno_file.is_open())
		{
			std::string error_message = "Cannot open phenotype file: " + fam_name;
			throw std::runtime_error(error_message);
		}
		while (std::getline(pheno_file, line))
		{
			misc::trim(line);
			if (!line.empty()) {
				std::vector < std::string > token = misc::split(line);
				if (token.size() < 6)
					throw std::runtime_error( "Malformed fam file, should contain at least 6 columns");
				prs_score cur_score;
				std::get<+PRS::FID>(cur_score) = token[+FAM::FID];
				std::get<+PRS::IID>(cur_score) = token[+FAM::IID];
				std::get<+PRS::PRS>(cur_score) = 0.0;
				std::get<+PRS::NNMISS>(cur_score) = 0;
				m_sample_names.push_back(cur_score);
				std::string id =(m_ignore_fid)? token[+FAM::IID]:token[+FAM::FID]+"_"+token[+FAM::IID];
				if (token[+FAM::PHENOTYPE] != "NA")
				{
					try {
						if (binary)
						{ //current coding must be 1 2
							int temp = misc::convert<int>( token[+FAM::PHENOTYPE]);
							if (temp >= 0 && temp <= 2)
							{
								max_num = (temp>max_num)? temp:max_num;
								m_sample_with_phenotypes[id] = cur_index;
								phenotype_store.push_back(temp);
								cur_index++;
								if(temp == 1) num_case++;
								if(temp == 0) num_control++;
							}
						}
						else
						{
							double temp = misc::convert<double>( token[+FAM::PHENOTYPE]);
							m_sample_with_phenotypes[id] = cur_index;
							phenotype_store.push_back(temp);
							cur_index++;
						}
					} catch (const std::runtime_error &error) { }
				}
			}
		}
		pheno_file.close();
	}
	else
	{
		std::ifstream fam;
		fam.open(fam_name.c_str());
		pheno_file.open(c_pheno.c_str());
		if (!fam.is_open())
		{
			std::string error_message = "Cannot open fam file: " + fam_name;
			throw std::runtime_error(error_message);
		}
		if (!pheno_file.is_open())
		{
			std::string error_message = "Cannot open phenotype file: " + fam_name;
			throw std::runtime_error(error_message);
		}
		// Main problem: we want the order following the fam file
		std::unordered_map < std::string, std::string > phenotype_info;
		cur_index = 0;
		while (std::getline(pheno_file, line))
		{
			misc::trim(line);
			if (!line.empty())
			{
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
		}
		pheno_file.close();
		while (std::getline(fam, line))
		{
			misc::trim(line);
			if (!line.empty())
			{
				std::vector < std::string > token = misc::split(line);
				if (token.size() < 6)
					std::runtime_error( "Malformed fam file, should contain at least 6 columns");
				prs_score cur_score;
				std::get<+PRS::FID>(cur_score) = token[+FAM::FID];
				std::get<+PRS::IID>(cur_score) = token[+FAM::IID];
				std::get<+PRS::PRS>(cur_score) = 0.0;
				std::get<+PRS::NNMISS>(cur_score) = 0;
				m_sample_names.push_back(cur_score);
				std::string id =(m_ignore_fid)? token[+FAM::IID]:token[+FAM::FID]+"_"+token[+FAM::IID];
				if (phenotype_info.find(id) != phenotype_info.end())
				{
					std::string p = phenotype_info[id];
					if (p.compare("NA") != 0)
					{
						try {
							if (binary)
							{
								int temp = misc::convert<int>(p);
								if (temp  >= 0 && temp <= 2)
								{
									max_num = (temp>max_num)? temp:max_num;
									m_sample_with_phenotypes[id] = cur_index;
									phenotype_store.push_back(temp);
									cur_index++;
									if(temp == 1) num_case++;
									if(temp == 0) num_control++;
								}
							}
							else
							{
								double temp = misc::convert<double>(p);
								m_sample_with_phenotypes[id] = cur_index;
								phenotype_store.push_back(temp);
								cur_index++;
							}
						} catch (const std::runtime_error &error) {}
					}
				}
				else
				{
					n_not_found++;
				}
			}
		}
		fam.close();
	}
	if (n_not_found != 0)
	{
		fprintf(stderr, "Number of missing samples: %zu\n", n_not_found);
	}
	bool error = false;
	if(max_num > 1 && binary)
	{
		num_case = 0;
		num_control = 0;
		for(size_t i = 0; i < phenotype_store.size(); ++i)
		{
			phenotype_store[i]--;
			if(phenotype_store[i] < 0) error = true;
			else (phenotype_store[i]==1)? num_case++: num_control++;
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

void PRSice::gen_cov_matrix(const std::string &c_cov_file,
		const std::vector<std::string> &c_cov_header) {
	size_t num_sample = m_sample_with_phenotypes.size();
	if (c_cov_file.empty())
	{
		m_independent_variables = Eigen::MatrixXd::Ones(num_sample, 2);
		return;
	}
	std::ifstream cov;
	cov.open(c_cov_file.c_str());
	if (!cov.is_open())
	{
		std::string error_message = "ERROR: Cannot open covariate file: " + c_cov_file;
		throw std::runtime_error(error_message);
	}
	std::string line;
	std::vector < size_t > cov_index;
	int max_index = 0;
	size_t num_valid = 0;
	std::getline(cov, line);
	// obtain the header information of the covariate file
	if (!line.empty())
	{
		std::vector < std::string > token = misc::split(line);
		if (c_cov_header.size() == 0)
		{
			// if no header is provided, we will use all the covariates included
			for (size_t i = 1+!m_ignore_fid; i < token.size(); ++i) cov_index.push_back(i); //FID, therefore+1
			max_index = cov_index.size() - 1;
		}
		else
		{
			std::unordered_set<std::string> included;
			// if specific headers are provided, we should only include them
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
					if (i_header > max_index) max_index = i_header;
				}
			}
		}
	}
	else throw std::runtime_error("First line of covariate file is empty!");
	fprintf(stderr, "\nStart processing the covariates\n");
	fprintf(stderr, "==============================\n");
	std::vector < std::pair<std::string, size_t> > valid_samples;
	m_independent_variables = Eigen::MatrixXd::Ones(num_sample, cov_index.size() + 2);
	bool valid = true;
	while (std::getline(cov, line))
	{
		misc::trim(line);
		if (!line.empty())
		{
			valid = true;
			std::vector < std::string > token = misc::split(line);
			if (token.size() <= max_index)
			{
				std::string error_message = "ERROR: Malformed covariate file, should contain at least "
						+ std::to_string(max_index + 1) + " column!";
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
					valid_samples.push_back( std::pair<std::string, size_t>(id, index));
					num_valid++;
				}
			}
		}
	}
	// now we need to handle the situation where there are a different number of samples

	if (valid_samples.size() != num_sample && num_sample != 0) {
		int removed = num_sample - valid_samples.size();
		fprintf(stderr, "Number of samples with invalid covariate: %d\n", removed);
		double portion = (double) removed / (double) num_sample;
		if(valid_samples.size() == 0)
		{
			throw std::runtime_error("All samples removed due to missingness in covariate file!");
		}
		if (portion > 0.05) {
			fprintf(stderr, "WARNING! More than %03.2f%% of the samples were removed!\n", portion * 100);
			fprintf(stderr, "         Do check if your covariate file is correct\n");
		}
		std::sort(begin(valid_samples), end(valid_samples),
				[](std::pair<std::string, size_t> const &t1, std::pair<std::string, size_t> const &t2)
				{
					if(std::get<1>(t1)==std::get<1>(t2)) return std::get<0>(t1).compare(std::get<0>(t2)) < 0;
					else return std::get<1>(t1) < std::get<1>(t2);
				});

		// update the m_phenotype and m_independent
		m_sample_with_phenotypes.clear();
		for (size_t cur_index = 0; cur_index < valid_samples.size(); ++cur_index) {
			std::string name = std::get < 0 > (valid_samples[cur_index]);
			m_sample_with_phenotypes[name]  = cur_index;
			size_t update_index = std::get < 1 > (valid_samples[cur_index]);
			if (update_index != cur_index) {
				m_phenotype(cur_index, 0) = m_phenotype(update_index, 0);
				for (size_t i_cov = 0; i_cov < cov_index.size(); ++i_cov) {
					m_independent_variables(cur_index, i_cov + 2) = m_independent_variables(update_index, i_cov + 2);
				}
			}
		}
		m_independent_variables.conservativeResize(valid_samples.size(), m_independent_variables.cols());
		m_phenotype.conservativeResize(valid_samples.size(), 1);

		fprintf(stderr, "\nFinal number of samples: %zu\n\n", valid_samples.size());
	}
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

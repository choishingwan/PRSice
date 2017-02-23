/*
 * plink.cpp
 *
 *  Created on: 19 Feb 2017
 *      Author: shingwan
 */

#include "plink.hpp"
std::vector<std::string> PLINK::g_chr_list;
std::mutex PLINK::clump_mtx;
PLINK::PLINK(){}

PLINK::PLINK(std::string prefix, bool verbose, const size_t thread, const catelog &inclusion):m_thread(thread){
	// TODO Auto-generated constructor stub
	// All three file is CLOSED after this
	if(prefix.find("#")!=std::string::npos)
	{
		for(auto &&chr: g_chr_list)
		{
			std::string name = prefix;
			misc::replace_substring(name, "#", chr);
			m_prefix.push_back(name);
		}
	}
	if(m_prefix.size()==0) m_prefix.push_back(prefix);

	int32_t retval = 0;
	retval = load_bim(inclusion);
	retval = load_fam();

	// marker_ct is specific to individual bim file
	m_marker_ct = m_unfiltered_marker_ct - m_marker_exclude_ct; // seems reasonable
	m_unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
	m_unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
	uint32_t uii = BITCT_TO_WORDCT(m_unfiltered_marker_ct);
	//m_marker_reverse = new uintptr_t[uii];
	//std::memset(m_marker_reverse, 0x0, uii*sizeof(uintptr_t));

	retval = load_bed();
	if(verbose)
	{
		fprintf(stderr, "%zu people (%zu males, %zu females) loaded from .fam\n", m_unfiltered_sample_ct, m_num_male, m_num_female);
		fprintf(stderr, "%zu variants included\n", m_marker_ct);
	}
}

PLINK::~PLINK() {
	// TODO Auto-generated destructor stub
	// Unfortunately, as I only have partial understanding of the plink code
	// and are far worst a programmer when compared to Chris, I don't know
	// how Chris free the memory in plink. So we will have to live with
	// the fact that the plink code will have memory leak
	if(m_marker_exclude != nullptr) delete [] m_marker_exclude;
	if(m_sex_male != nullptr) delete [] m_sex_male;
	if(m_founder_info != nullptr) delete [] m_founder_info;
	if(m_sample_exclude != nullptr) delete [] m_sample_exclude;
	//delete [] m_marker_reverse;
}

int32_t PLINK::load_bed(const std::string &bedname)
{
	uint32_t uii = 0;
	int64_t llxx = 0;
	int64_t llyy = 0;
	int64_t llzz = 0;
	m_bedfile = fopen(bedname.c_str(), FOPEN_RB);
	if (fseeko(m_bedfile, 0, SEEK_END)) {
		std::string error_message = "Cannot read bed file: "+bedname;
		throw std::runtime_error(error_message);
	}
	llxx = ftello(m_bedfile);
	if (!llxx) {
		throw std::runtime_error("Error: Empty .bed file.");
	}
	rewind(m_bedfile);
	// will let the g_textbuf stay for now
	char version_check[3];
	uii = fread(version_check, 1, 3, m_bedfile);
	llyy = ((uint64_t)m_unfiltered_sample_ct4) * m_unfiltered_marker_ct;
	llzz = ((uint64_t)m_unfiltered_sample_ct) * ((m_unfiltered_marker_ct + 3) / 4);
	bool sample_major = false;
	// compare only the first 3 bytes
	if ((uii == 3) && (!memcmp(version_check, "l\x1b\x01", 3)))
	{
		llyy += 3;
	}
	else if ((uii == 3) && (!memcmp(version_check, "l\x1b", 3)))
	{
		// v1.00 sample-major
		sample_major=true;
		llyy = llzz + 3;
		m_bed_offset = 2;
	}
	else if (uii && (*version_check == '\x01'))
	{
		// v0.99 SNP-major
		llyy += 1;
		m_bed_offset = 1;
	}
	else if (uii && (!(*version_check)))
	{
		// v0.99 sample-major
		sample_major=true;
		llyy = llzz + 1;
		m_bed_offset = 2;
	}
	else
	{
		// pre-v0.99, sample-major, no header bytes
		sample_major=true;
		if (llxx != llzz)
		{
			// probably not PLINK-format at all, so give this error instead of
			// "invalid file size"
			throw std::runtime_error("Error: Invalid header bytes in .bed file.");
		}
		llyy = llzz;
		m_bed_offset = 2;
	}
	if (llxx != llyy)
	{
		if ((*version_check == '#') || ((uii == 3) && (!memcmp(version_check, "chr", 3))))
		{
			throw std::runtime_error("Error: Invalid header bytes in PLINK 1 .bed file.  (Is this a UCSC Genome\nBrowser BED file instead?)");
		}
		else
		{
			throw std::runtime_error("Error: Invalid .bed file size.");
		}
	}
	if(sample_major)
	{
		throw std::runtime_error("Error: Currently do not support sample major format");
	}
	fclose(m_bedfile);
	m_bedfile = nullptr;
    return 0;
}



int32_t PLINK::load_bed()
{
	uint32_t uii = 0;
	int64_t llxx = 0;
	int64_t llyy = 0;
	int64_t llzz = 0;
	for(auto &&prefix : m_prefix)
	{
		std::string bedname = prefix+".bed";
		m_bedfile = fopen(bedname.c_str(), FOPEN_RB);
		if (fseeko(m_bedfile, 0, SEEK_END)) {
			std::string error_message = "Cannot read bed file: "+bedname;
			throw std::runtime_error(error_message);
		}
		llxx = ftello(m_bedfile);
		if (!llxx) {
		  throw std::runtime_error("Error: Empty .bed file.");
		}
		rewind(m_bedfile);
		// will let the g_textbuf stay for now
		char version_check[3];
		uii = fread(version_check, 1, 3, m_bedfile);
		llyy = ((uint64_t)m_unfiltered_sample_ct4) * m_unfiltered_marker_ct;
		llzz = ((uint64_t)m_unfiltered_sample_ct) * ((m_unfiltered_marker_ct + 3) / 4);
		bool sample_major = false;
		// compare only the first 3 bytes
		if ((uii == 3) && (!memcmp(version_check, "l\x1b\x01", 3)))
		{
			llyy += 3;
		}
		else if ((uii == 3) && (!memcmp(version_check, "l\x1b", 3)))
		{
			// v1.00 sample-major
			sample_major=true;
			llyy = llzz + 3;
			m_bed_offset = 2;
		}
		else if (uii && (*version_check == '\x01'))
		{
			// v0.99 SNP-major
			llyy += 1;
			m_bed_offset = 1;
		}
		else if (uii && (!(*version_check)))
		{
			// v0.99 sample-major
			sample_major=true;
			llyy = llzz + 1;
			m_bed_offset = 2;
		}
		else
		{
			// pre-v0.99, sample-major, no header bytes
			sample_major=true;
			if (llxx != llzz)
			{
				// probably not PLINK-format at all, so give this error instead of
				// "invalid file size"
				throw std::runtime_error("Error: Invalid header bytes in .bed file.");
			}
			llyy = llzz;
			m_bed_offset = 2;
		}
		if (llxx != llyy)
		{
			if ((*version_check == '#') || ((uii == 3) && (!memcmp(version_check, "chr", 3))))
			{
				throw std::runtime_error("Error: Invalid header bytes in PLINK 1 .bed file.  (Is this a UCSC Genome\nBrowser BED file instead?)");
			}
			else
			{
				throw std::runtime_error("Error: Invalid .bed file size.");
			}
		}
		if(sample_major)
		{
			throw std::runtime_error("Error: Currently do not support sample major format");
		}
		fclose(m_bedfile);
		m_bedfile = nullptr;
	}
    return 0;
}

int32_t PLINK::load_bim(const catelog &inclusion)
{
	m_unfiltered_marker_ct = 0;
	uintptr_t unfiltered_marker_ctl;
	std::ifstream bimfile;
	size_t num_snp = 0;
	std::vector<size_t> filter;
	for(auto &&prefix : m_prefix)
	{
		std::string bimname = prefix+".bim";
		bimfile.open(bimname.c_str());
		if(!bimfile.is_open())
		{
			std::string error_message = "Error: Cannot open bim file: "+bimname;
			throw std::runtime_error(error_message);
		}
		std::string line;
		int num_line = 0;
		while(std::getline(bimfile, line))
		{

			misc::trim(line);
			if(!line.empty())
			{
				std::vector<std::string> token = misc::split(line);
				if(token.size() < 6)
				{
					fprintf(stderr, "Error: Malformed bim file. Less than 6 column on line: %i\n",num_line);
					throw std::runtime_error("");
				}
				//	for filtering snps
				//	SAM: with my way of memory control, this will likely cause problem
				//	SET_BIT(marker_uidx, marker_exclude);
				if(!inclusion.empty() && //we don't want this when inclusion isn't provided
						inclusion.find(token[+BIM::RS])==inclusion.end()) // this avoid reading the file twice
				{
					filter.push_back(num_snp);
					m_marker_exclude_ct++;
				}
				else if(inclusion.find(token[+BIM::RS])!=inclusion.end())
				{
					snp_link info;
					std::get < +FILE_INFO::FILE >(info) =prefix;
					std::get < +FILE_INFO::INDEX >(info) =inclusion.at(token[+BIM::RS]);
					std::get < +FILE_INFO::LINE >(info) =num_line;;
					m_snp_link.push_back(info);
				}
			}
			num_snp++;
			m_unfiltered_marker_ct++; //add in the checking later on
			num_line++;
		}
		bimfile.close();
	}
	unfiltered_marker_ctl = BITCT_TO_WORDCT(m_unfiltered_marker_ct);
	m_marker_exclude = new uintptr_t[unfiltered_marker_ctl];
	std::memset(m_marker_exclude, 0x0, unfiltered_marker_ctl*sizeof(uintptr_t));

	if(m_unfiltered_marker_ct > 2147483645)
	{
		throw std::runtime_error("Error: PLINK does not suport more than 2^31 -3 variants. "
			"As we are using PLINK for some of our functions, we might encounter problem too. "
			"Sorry.");
	}
	if(m_unfiltered_marker_ct==m_marker_exclude_ct)
	{
		throw std::runtime_error("Error: All variants excluded.");
	}
	// only bother doing this if we have the correct number of markers to work with
	for(auto &index : filter)
	{
		SET_BIT(index, m_marker_exclude);
	}
	return 0;
}


int32_t PLINK::load_fam(){
	std::string famname = m_prefix.front()+".fam";
	m_unfiltered_sample_ct = 0;

	uintptr_t sample_uidx;
	uintptr_t unfiltered_sample_ctl;

	// this should give us the correct information
	std::ifstream famfile;
	famfile.open(famname.c_str());
	if(!famfile.is_open())
	{
		std::string error_message = "Error: Cannot open fam file: "+famname;
		throw std::runtime_error(error_message);
	}
	std::string line;
	sample_uidx=0;
	//first pass to get the number of samples
	while(std::getline(famfile, line))
	{
		misc::trim(line);
		if(!line.empty())
		{
			std::vector<std::string> token = misc::split(line);
			if(token.size() < 6)
			{
				fprintf(stderr, "Error: Malformed fam file. Less than 6 column on line: %zu\n",sample_uidx+1);
				throw std::runtime_error("");
			}
			m_unfiltered_sample_ct++;
			sample_uidx++;
		}
	}
	famfile.clear();
	famfile.seekg(0);
	sample_uidx=0;
	unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
	m_sex_male = new uintptr_t[unfiltered_sample_ctl];
	std::memset(m_sex_male, 0x0, unfiltered_sample_ctl*sizeof(uintptr_t));
	m_founder_info = new uintptr_t[unfiltered_sample_ctl];
	std::memset(m_founder_info, 0x0, unfiltered_sample_ctl*sizeof(uintptr_t));
	m_sample_exclude = new uintptr_t[unfiltered_sample_ctl];
	std::memset(m_sample_exclude, 0x0, unfiltered_sample_ctl*sizeof(uintptr_t));
	m_num_male = 0, m_num_female = 0, m_num_ambig_sex=0;
	while(std::getline(famfile, line))
	{
		misc::trim(line);
		if(!line.empty())
		{
			std::vector<std::string> token = misc::split(line);
			if(token.size() < 6)
			{
				fprintf(stderr, "Error: Malformed fam file. Less than 6 column on line: %zu\n",sample_uidx+1);
				throw std::runtime_error("");
			}
			if(token[+FAM::FATHER].compare("0")==0 && token[+FAM::MOTHER].compare("0")==0)
			{
				m_founder_ct++;
				SET_BIT(sample_uidx, m_founder_info); // if individual is founder e.g. 0 0, then set bit
			}
			if(token[+FAM::SEX].compare("1")==0)
			{
				m_num_male++;
				SET_BIT(sample_uidx, m_sex_male); // if that individual is male, need to set bit
			}
			else if(token[+FAM::SEX].compare("2")==0)
			{
				 m_num_female++;
			}
			else
			{
				 m_num_ambig_sex++;
			}
			sample_uidx++;
		}

	}
	famfile.close();
	return 0;
}

void PLINK::get_score(const std::vector<p_partition> &partition,
                      const boost::ptr_vector<SNP> &snp_list, std::vector< std::vector<prs_score> > &prs_score,
                      size_t start_index, size_t end_bound, size_t num_region, SCORING scoring)
{
    size_t prev =0;
    uint32_t uii;
    uint32_t ujj;
    uint32_t ukk;
    uintptr_t ulii = 0;
    uintptr_t sample_idx;
    int32_t delta1 = 1;
	int32_t delta2 = 2;
	int32_t deltam = 0;
    uintptr_t final_mask = get_final_mask(m_founder_ct);
	// a lot of code in PLINK is to handle the sex chromosome
	// which suggest that PRS can be done on sex chromosome
	// that should be something later
    std::string prev_name = "";
    uint32_t max_reverse = BITCT_TO_WORDCT(m_unfiltered_marker_ct);
	for(size_t i_region=0; i_region < num_region; ++i_region)
	{
		if(prs_score[i_region].size() < m_unfiltered_sample_ct)
		{
			throw std::runtime_error("Size of vector doesn't match number of samples!!");
		}
	}

	std::vector<bool> in_region;
	// index is w.r.t. partition, which contain all the information
    uintptr_t* genotype = new uintptr_t[m_unfiltered_sample_ctl*2];
    uintptr_t* tmp_genotype = new uintptr_t[m_unfiltered_sample_ctl*2];
    for(size_t i_snp = start_index; i_snp < end_bound; ++i_snp)
    { // for each SNP
        if(prev_name.empty() || prev_name.compare(std::get<+PRS::FILENAME>(partition[i_snp]))!=0)
        {
        	if(m_bedfile!=nullptr) fclose(m_bedfile);
            prev_name= std::get<+PRS::FILENAME>(partition[i_snp]);
            std::string bedname = prev_name+".bed";
            load_bed(bedname);
        	m_bedfile = fopen(bedname.c_str(), FOPEN_RB);
            prev=0;
        }

        int snp_index = std::get<+PRS::INDEX>(partition[i_snp]);

    	for(size_t i_region = 0; i_region <num_region; ++i_region)
    	{
    		in_region.push_back(snp_list[snp_index].in(i_region));
    	}
        size_t cur_line = std::get<+PRS::LINE>(partition[i_snp]);
        if (fseeko(m_bedfile, m_bed_offset + (cur_line* ((uint64_t)m_unfiltered_sample_ct4))
        		, SEEK_SET)) {
        	throw std::runtime_error("ERROR: Cannot read the bed file!");
        }
        //loadbuf_raw is the temporary
        //loadbuff is where the genotype will be located
        std::memset(genotype, 0x0, m_unfiltered_sample_ctl*2*sizeof(uintptr_t));
        std::memset(tmp_genotype, 0x0, m_unfiltered_sample_ctl*2*sizeof(uintptr_t));
        if(load_and_collapse_incl(m_unfiltered_sample_ct, m_founder_ct, m_founder_info, final_mask,
        		false, m_bedfile, tmp_genotype, genotype))
        {
        	throw std::runtime_error("ERROR: Cannot read the bed file!");
        }
        uintptr_t* lbptr = genotype;
        uii = 0;
        // normally, sample_ct is the number of samples remain "after" excluding them
        // but as PRSice does not allow sample exclusion, we can just use m_unfiltered_sample_ct
        // Also, it is tempting to use PLINK's version of PRS counting. But that'd be too much
        // for me


        std::vector<size_t> missing_samples;
        double stat = snp_list[snp_index].get_stat();
        std::vector<double> genotypes(m_unfiltered_sample_ct);
        int total_num = 0;
        do {
        	ulii = ~(*lbptr++);
        	if (uii + BITCT2 > m_unfiltered_sample_ct) {
        		ulii &= (ONELU << ((m_unfiltered_sample_ct & (BITCT2 - 1)) * 2)) - ONELU;
        	}
        	while (ulii) {
        		ujj = CTZLU(ulii) & (BITCT - 2);
        		ukk = (ulii >> ujj) & 3;
        		sample_idx = uii + (ujj / 2);
        		if(ukk!=1) // Because 01 is coded as missing
        		{
        			int flipped_geno = snp_list[snp_index].geno(ukk!=1);
        			total_num+=flipped_geno;
        			genotypes[sample_idx] = flipped_geno;
        		}
        		else
        		{
        			missing_samples.push_back(sample_idx);
        		}
        		ulii &= ~((3 * ONELU) << ujj);
        	}
        	uii += BITCT2;
        } while (uii < m_unfiltered_sample_ct);
        size_t i_missing = 0;
        double center_score = stat*((double)total_num/((double)m_unfiltered_sample_ct*2.0));
        size_t num_miss = missing_samples.size();
        for(size_t i_sample=0; i_sample < m_unfiltered_sample_ct; ++i_sample)
        {
        	if(i_missing < num_miss && i_sample == missing_samples[i_missing])
        	{

        		for(size_t i_region; i_region < num_region; ++i_region)
        		{

        			if(in_region[i_region])
        			{
        				if(scoring == SCORING::MEAN_IMPUTE) std::get<+PRS::PRS>(prs_score[i_region][i_sample]) += center_score;
        				if(scoring != SCORING::SET_ZERO) std::get<+PRS::NNMISS>(prs_score[i_region][i_sample])++;
        			}
        		}
        		i_missing++;
        	}
        	else
        	{ // not missing sample
        		for(size_t i_region=0; i_region < num_region; ++i_region)
        		{
        			if(in_region[i_region])
        			{
        				if(scoring == SCORING::CENTER)
        				{
        					// if centering, we want to keep missing at 0
        					std::get<+PRS::PRS>(prs_score[i_region][i_sample]) -= center_score;
        				}
        				std::get<+PRS::PRS>(prs_score[i_region][i_sample]) += genotypes[i_sample]*stat*0.5;
        				std::get<+PRS::NNMISS>(prs_score[i_region][i_sample]) ++;
        			}
        		}
        	}
        }
    }
    delete [] genotype;
    delete [] tmp_genotype;
}
/*
        char *genotype_list = new char[m_num_bytes];
        m_bed.read((char*)genotype_list, m_num_bytes);
        if (!m_bed) throw std::runtime_error("Problem with the BED file...has the FAM/BIM file been changed?");
        prev++;
        size_t sample_index = 0;
        int snp_index = std::get<+PRS::INDEX>(partition[i_snp]);
        if(snp_index >= snp_list.size()) throw std::runtime_error("Out of bound! In PRS score calculation");
        std::vector<bool> in_region;
        int num_region = 0;
        for(size_t i_region = 0; i_region < prs_score.size(); ++i_region)
        {
        		num_region++;
        		in_region.push_back(snp_list[snp_index].in(i_region));
        }
        double stat = snp_list[snp_index].get_stat();
        std::vector<size_t> missing_samples;
        std::vector<double> genotypes(m_num_sample);
        int total_num = 0;
        for(size_t i_byte = 0; i_byte < m_num_bytes; ++i_byte)
        {

        		size_t geno_bit = 0;
    			int geno_batch = static_cast<int>(genotype_list[i_byte]);
        		while(geno_bit < 7 && sample_index < m_num_sample)
        		{
        			int geno = geno_batch>>geno_bit & 3; // This will access the corresponding genotype
        			if(geno!=1) // Because 01 is coded as missing
        			{
        				int flipped_geno = snp_list[snp_index].geno(geno);
        				total_num+=flipped_geno;
        				genotypes[sample_index] = flipped_geno;
        			}
        			else
        			{
        				missing_samples.push_back(sample_index);
        			}
        			sample_index++;
        			geno_bit+=2;
        		}
        }
        delete[] genotype_list;


		size_t i_missing = 0;
		double center_score = stat*((double)total_num/((double)m_num_sample*2.0));
		size_t num_miss = missing_samples.size();
		for(size_t i_sample=0; i_sample < m_num_sample; ++i_sample)
		{
			if(i_missing < num_miss && i_sample == missing_samples[i_missing])
			{
				for(size_t i_region; i_region < num_region; ++i_region)
				{
					if(in_region[i_region])
					{
						if(m_scoring == SCORING::MEAN_IMPUTE) std::get<+PRS::PRS>(prs_score[i_region][i_sample]) += center_score;
						if(m_scoring != SCORING::SET_ZERO) std::get<+PRS::NNMISS>(prs_score[i_region][i_sample])++;
					}
				}
				i_missing++;
			}
			else
			{ // not missing sample
				for(size_t i_region=0; i_region < num_region; ++i_region)
				{
					if(in_region[i_region])
					{
						if(m_scoring == SCORING::CENTER){
							std::get<+PRS::PRS>(prs_score[i_region][i_sample]) -= center_score;
						}
						std::get<+PRS::PRS>(prs_score[i_region][i_sample]) += genotypes[i_sample]*stat*0.5;
    						std::get<+PRS::NNMISS>(prs_score[i_region][i_sample]) ++;
					}
				}
			}
		}
    }
    */

void PLINK::start_clumping(catelog& inclusion, boost::ptr_vector<SNP> &snp_list, double p_threshold, double r2_threshold,
		size_t kb_threshold, double proxy_threshold){
	assert(m_snp_link.size()!=0);
	// The m_snp_link vector should contain the bim file name and the line number for
	// the SNP at certain index in the m_snp_list vector
	// This is cryptic but then, hopefully that should work
	uintptr_t unfiltered_sample_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(m_unfiltered_sample_ct);
	std::deque<size_t> clump_snp_index; // Index for SNP within the clumping region
	size_t snp_id_in_list = std::get<+FILE_INFO::INDEX>(m_snp_link.front());
	std::string prev_chr= snp_list[snp_id_in_list].get_chr();

	std::string prev_file=std::get<+FILE_INFO::FILE>(m_snp_link.front());;
	std::string bedname = prev_file+".bed";
	m_bedfile = fopen(bedname.c_str(), FOPEN_RB);
	size_t bp_of_core =snp_list[snp_id_in_list].get_loc();
    size_t core_genotype_index=0; //index of the core SNP on our genotype deque
    bool require_clump=false; // Whether if the current interval contain the core snp
    uintptr_t final_mask = get_final_mask(m_founder_ct);
    size_t num_snp = m_snp_link.size();
    for(size_t i_info = 0; i_info < m_snp_link.size(); ++i_info)
	{
    	size_t cur_snp_index = std::get<+FILE_INFO::INDEX>(m_snp_link[i_info]);
    	size_t cur_line_num = std::get<+FILE_INFO::LINE>(m_snp_link[i_info]);
    	std::string cur_chr = snp_list[cur_snp_index].get_chr();
    	size_t cur_loc = snp_list[cur_snp_index].get_loc();
    	if(prev_chr.compare(cur_chr)!=0)
		{
    		perform_clump(clump_snp_index, snp_list, core_genotype_index, require_clump, p_threshold,
					r2_threshold, kb_threshold, cur_chr, cur_loc);
			if(prev_file.empty() || prev_file.compare(std::get<+FILE_INFO::FILE>(m_snp_link[i_info]))!=0)
			{
				prev_file = std::get<+FILE_INFO::FILE>(m_snp_link[i_info]);
				if (m_bedfile != nullptr)
				{
					fclose(m_bedfile);
					m_bedfile=nullptr;
				}
				std::string bed_name = prev_file+".bed";
				load_bed(bed_name);
				m_bedfile = fopen(bedname.c_str(), FOPEN_RB);
			}
			prev_chr = cur_chr;
		}
        else if((cur_loc-bp_of_core) > kb_threshold){
			perform_clump(clump_snp_index, snp_list, core_genotype_index, require_clump, p_threshold,
					r2_threshold, kb_threshold, cur_chr, cur_loc);
        }
        if (fseeko(m_bedfile, m_bed_offset + (cur_line_num* ((uint64_t)m_unfiltered_sample_ct4))
        		, SEEK_SET)) {
        	throw std::runtime_error("ERROR: Cannot read the bed file!");
        }
        //loadbuf_raw is the temporary
        //loadbuff is where the genotype will be located
        uintptr_t* genotype = new uintptr_t[m_unfiltered_sample_ctl*2];
        std::memset(genotype, 0x0, m_unfiltered_sample_ctl*2*sizeof(uintptr_t));
        uintptr_t* tmp_genotype = new uintptr_t[m_unfiltered_sample_ctl*2];
        std::memset(tmp_genotype, 0x0, m_unfiltered_sample_ctl*2*sizeof(uintptr_t));
        if(load_and_collapse_incl(m_unfiltered_sample_ct, m_founder_ct, m_founder_info, final_mask,
        		false, m_bedfile, tmp_genotype, genotype))
        {
        	throw std::runtime_error("ERROR: Cannot read the bed file!");
        }
        m_genotype.push_back(genotype);
        delete [] tmp_genotype;// don't need the temporary now

        clump_snp_index.push_back(cur_snp_index);
        if(!require_clump && snp_list[cur_snp_index].get_p_value() < p_threshold)
        {
        	bp_of_core =snp_list[cur_snp_index].get_loc();
        	core_genotype_index=m_genotype.size()-1; // Should store the index on genotype
        	require_clump= true;
        }
    	fprintf(stderr, "\rClumping Progress: %03.2f%%", (double) i_info / (double) (num_snp) * 100.0);
	}
    if(clump_snp_index.size()!=0)
    {
    	// this make sure this will be the last
    	perform_clump(clump_snp_index, snp_list, core_genotype_index, require_clump, p_threshold,
    						r2_threshold, kb_threshold, prev_chr+"_", bp_of_core+2*kb_threshold);
    }

	fprintf(stderr, "\rClumping Progress: %03.2f%%\n\n", 100.0);


	std::unordered_map<std::string, size_t> inclusion_backup = inclusion;
	inclusion.clear();
	std::vector<size_t> p_sort_order = SNP::sort_by_p(snp_list);
	bool proxy = proxy_threshold > 0.0;
	for(auto &&i_snp : p_sort_order){
		if(inclusion_backup.find(snp_list.at(i_snp).get_rs_id()) != inclusion_backup.end() &&
				snp_list[i_snp].get_p_value() < p_threshold)
		{
			if(proxy && !snp_list[i_snp].clumped() )
			{
				snp_list[i_snp].proxy_clump(snp_list, proxy_threshold);
				inclusion[snp_list[i_snp].get_rs_id()]=i_snp;
			}
			else if(!snp_list[i_snp].clumped())
			{
				snp_list[i_snp].clump(snp_list);
				inclusion[snp_list[i_snp].get_rs_id()]=i_snp;
			}
		}
		else if(snp_list[i_snp].get_p_value() >= p_threshold) break;

	}
	fprintf(stderr, "Number of SNPs after clumping : %zu\n\n", inclusion.size());
}

void PLINK::lerase(int num){
	if(num <0)
	{
		std::string error_message = "Number of removed SNPs cannot be less than 1: "+std::to_string(num);
		throw std::runtime_error(error_message);
	}
	if(num > m_genotype.size())
	{
		std::string error_message = "Number of removed SNPs exceed number of SNPs available "+std::to_string(num)+" "+std::to_string(m_genotype.size());
		throw std::runtime_error(error_message);
	}
	for(size_t i = 0; i < num; ++i)
	{
		delete [] m_genotype[i];
	}
	if(num==m_genotype.size())
	{
		m_genotype.clear();
	}
	else
	{
		m_genotype.erase(m_genotype.begin(), m_genotype.begin()+num);
	}

}


void PLINK::compute_clump( size_t core_snp_index, size_t i_start, size_t i_end, boost::ptr_vector<SNP> &snp_list,
		const std::deque<size_t> &clump_snp_index, const double r2_threshold, uintptr_t* geno1,
		bool nm_fixed, uint32_t* tot1)
{
	uintptr_t* loadbuf;
	uintptr_t founder_ctl = BITCT_TO_WORDCT(m_founder_ct);
	uint32_t founder_ctv3 = BITCT_TO_ALIGNED_WORDCT(m_founder_ct);
	uint32_t founder_ctsplit = 3 * founder_ctv3; // Required
	uint32_t marker_idx2_maxw =  m_marker_ct - 1;
	uint32_t counts[18];
	double freq11;
	double freq11_expected;
	double freq1x;
	double freq2x;
	double freqx1;
	double freqx2;
	double dxx;
	double r2 =0.0;
	bool zmiss2=false;
	uintptr_t ulii = founder_ctsplit * sizeof(intptr_t) + 2 * sizeof(int32_t) + marker_idx2_maxw * 2 * sizeof(double);
	size_t max_size = clump_snp_index.size();
    size_t ref_index = clump_snp_index.at(core_snp_index);
    double ref_p_value = snp_list.at(ref_index).get_p_value();
    std::vector<double> r2_store;
    std::vector<size_t> target_index_store; // index we want to push into the current index
	uintptr_t* ulptr =  new uintptr_t[3*founder_ctsplit +founder_ctv3];
    uintptr_t* dummy_nm = new uintptr_t[founder_ctl];
    for(size_t i_snp = i_start; i_snp < i_end && i_snp < max_size; ++i_snp)
    {
        zmiss2 = false;
        size_t target_index = clump_snp_index[i_snp];
        if(i_snp != core_snp_index && snp_list[target_index].get_p_value() > ref_p_value)
        {
            // only calculate r2 if more significant
        	std::memset(ulptr, 0x0, (3*founder_ctsplit +founder_ctv3)*sizeof(uintptr_t));
            std::memset(dummy_nm, ~0, founder_ctl*sizeof(uintptr_t)); // set all bits to 1
        	load_and_split3(m_genotype[i_snp], m_founder_ct, ulptr, dummy_nm, dummy_nm,
        			founder_ctv3, 0, 0, 1, &ulii);
        	uintptr_t uiptr[3];
        	uiptr[0] = popcount_longs(ulptr, founder_ctv3);
        	uiptr[1] = popcount_longs(&(ulptr[founder_ctv3]), founder_ctv3);
        	uiptr[2] = popcount_longs(&(ulptr[2 * founder_ctv3]), founder_ctv3);
        	if (ulii == 3) {
        		zmiss2 = true;
        	}

    		if (nm_fixed) {
    			two_locus_count_table_zmiss1(geno1, ulptr, counts, founder_ctv3, zmiss2);
    			if (zmiss2) {
    				counts[2] = tot1[0] - counts[0] - counts[1];
    				counts[5] = tot1[1] - counts[3] - counts[4];
    			}
    			counts[6] = uiptr[0] - counts[0] - counts[3];
    			counts[7] = uiptr[1] - counts[1] - counts[4];
    			counts[8] = uiptr[2] - counts[2] - counts[5];
    		} else {
    			two_locus_count_table(geno1, ulptr, counts, founder_ctv3, zmiss2);
    			if (zmiss2) {
    				counts[2] = tot1[0] - counts[0] - counts[1];
    				counts[5] = tot1[1] - counts[3] - counts[4];
    				counts[8] = tot1[2] - counts[6] - counts[7];
    			}
    		}
    		if(em_phase_hethet_nobase(counts, false, false, &freq1x, &freq2x, &freqx1, &freqx2, &freq11))
    		{
    			r2 = -1;
    		}
    		else
    		{
    			freq11_expected = freqx1 * freq1x;
    			dxx = freq11 - freq11_expected;
    			if (fabs(dxx) < SMALL_EPSILON)
    			{
    				r2 = 0.0;
    			}
    			else
    			{
    				r2 = fabs(dxx) * dxx / (freq11_expected * freq2x * freqx2);
    				//dxx/= MINV(freqx1 * freq2x, freqx2 * freq1x);
    				//dprime = dxx;

    			}
    		}
    		if(r2 >= r2_threshold)
            {
            	target_index_store.push_back(target_index);
   	            r2_store.push_back(r2);
            }
        }
    }
    delete [] ulptr;
	delete [] dummy_nm;
    PLINK::clump_mtx.lock();
    snp_list[ref_index].add_clump(target_index_store);
    snp_list[ref_index].add_clump_r2(r2_store);
    PLINK::clump_mtx.unlock();
}

void PLINK::clump_thread(const size_t c_core_index, const std::deque<size_t> &c_clump_snp_index,
		boost::ptr_vector<SNP> &snp_list, const double c_r2_threshold)
{
	// do this without the clumping first
	size_t wind_size = c_clump_snp_index.size();
	if(wind_size <=1 ) return; // nothing to do

	uintptr_t founder_ctl = BITCT_TO_WORDCT(m_founder_ct);
	uint32_t founder_ctv3 = BITCT_TO_ALIGNED_WORDCT(m_founder_ct);
	uint32_t founder_ctsplit = 3 * founder_ctv3; // Required
	uint32_t marker_idx2_maxw =  m_marker_ct - 1;
	bool nm_fixed=false;
	uint32_t tot1[6];
	uintptr_t ulii = founder_ctsplit * sizeof(intptr_t) + 2 * sizeof(int32_t) + marker_idx2_maxw * 2 * sizeof(double);
	uintptr_t* geno1 = new uintptr_t[3*founder_ctsplit +founder_ctv3];
	std::memset(geno1, 0x0, (3*founder_ctsplit +founder_ctv3)*sizeof(uintptr_t));
	uintptr_t* dummy_nm = new uintptr_t[founder_ctl];
	std::memset(dummy_nm, ~0, founder_ctl*sizeof(uintptr_t)); // set all bits to 1
	load_and_split3(m_genotype[c_core_index], m_founder_ct,
			geno1, dummy_nm, dummy_nm, founder_ctv3, 0, 0, 1, &ulii);
	delete [] dummy_nm;

	tot1[0] = popcount_longs(geno1, founder_ctv3);
	tot1[1] = popcount_longs(&(geno1[founder_ctv3]), founder_ctv3);
	tot1[2] = popcount_longs(&(geno1[2 * founder_ctv3]), founder_ctv3);

	if (ulii == 3)
	{
		nm_fixed = true;
	}

	std::vector<std::thread> thread_store;
	if((wind_size-1) < m_thread)
	{
		for(size_t i_snp = 0; i_snp < wind_size; ++i_snp)
		{
			if(c_clump_snp_index[i_snp]!=c_core_index)
			{
				thread_store.push_back(std::thread(&PLINK::compute_clump, this,
						c_core_index,i_snp, i_snp+1, std::ref(snp_list),
						std::cref(c_clump_snp_index), c_r2_threshold,
						std::ref(geno1), nm_fixed, std::ref(tot1)));
			}

		}
	}
	else
	{
		int num_snp_per_thread =(int)(wind_size) / (int)m_thread;  //round down
		int remain = (int)(wind_size) % (int)m_thread;
		int cur_start = 0;
		int cur_end = num_snp_per_thread;
		for(size_t i_thread = 0; i_thread < m_thread; ++i_thread)
		{
			thread_store.push_back(std::thread(&PLINK::compute_clump, this, c_core_index, cur_start,
					cur_end+(remain>0), std::ref(snp_list), std::cref(c_clump_snp_index),c_r2_threshold,
					std::ref(geno1), nm_fixed, std::ref(tot1)));
			cur_start = cur_end+(remain>0);
			cur_end+=num_snp_per_thread+(remain>0);
			if(cur_end>wind_size) cur_end =wind_size;
			remain--;
		}
	}
	for(auto &&thread_runner : thread_store) thread_runner.join();
	thread_store.clear();

	delete [] geno1;
}

void PLINK::perform_clump(std::deque<size_t> &clump_snp_index, boost::ptr_vector<SNP> &snp_list,
		size_t &core_snp_index, bool &require_clump, double p_threshold, double r2_threshold,
		size_t kb_threshold, std::string next_chr, size_t next_loc){
	// The next_chr and next_loc basically = currently start_clump is pointing to this snp
	if(clump_snp_index.size()==0) return; // got nothing to do
	std::string core_chr = snp_list[core_snp_index].get_chr();
	size_t core_loc = snp_list[core_snp_index].get_loc();
	size_t infinite_guard = 0;
	size_t max_possible = clump_snp_index.size();
	while(require_clump && (core_chr.compare(next_chr)!=0 || (next_loc - core_loc) > kb_threshold))
	{ // as long as we still need to perform clumping
		clump_thread(core_snp_index, clump_snp_index, snp_list, r2_threshold);
		require_clump = false;
		for(size_t core_finder = core_snp_index+1; core_finder < clump_snp_index.size(); ++core_finder)
		{
			if(snp_list[clump_snp_index[core_finder]].get_p_value() < p_threshold)
			{
				core_snp_index=core_finder;
				core_chr = snp_list[core_snp_index].get_chr();
				core_loc = snp_list[core_snp_index].get_loc();
				require_clump= true;
				break;
			}
		}
		// New core found, need to clean things up a bit
		if(require_clump)
		{
			size_t num_remove = 0;
			for(size_t remover = 0; remover < core_snp_index; ++remover)
			{
				if(core_loc-snp_list[clump_snp_index[remover]].get_loc() > kb_threshold) num_remove++;
				else break;
			}
			if(num_remove!=0)
			{
				lerase(num_remove);
				clump_snp_index.erase(clump_snp_index.begin(), clump_snp_index.begin()+num_remove);
				core_snp_index-=num_remove;
			}
		}
		infinite_guard++;
		if(infinite_guard>max_possible) throw std::logic_error("ERROR: While loop run longer than expected");
	}
	// for this to be true, the require clump should be false or the core_snp is now within reach of the
	// new snp
	if(core_chr.compare(next_chr)!=0)
	{ 	// new chromosome
		// just remove everything
		lerase(m_genotype.size());
		clump_snp_index.clear();
	}
	else if(!require_clump){
		//remove anything that is too far ahead
		size_t num_remove = 0;
		for(auto &&remover : clump_snp_index)
		{
			if(next_loc-snp_list[remover].get_loc()>kb_threshold) num_remove++;
			else break;
		}
		if(num_remove!=0)
		{
			lerase(num_remove);
			clump_snp_index.erase(clump_snp_index.begin(), clump_snp_index.begin()+num_remove);
		}
	}
}

/*
void PLINK::initialize(){
	intptr_t malloc_size_mb = 0; //supposed to be an option in plink for people to select amount of memory use
	intptr_t default_alloc_mb = 0;
	int32_t mib[2];
	int64_t llxx = 0;
	size_t sztmp;
	unsigned char* bigstack_ua = nullptr; // ua = unaligned
	unsigned char* bigstack_initial_base;
	char* bubble = nullptr; // seems like another crazy trick by Chris
#ifdef _WIN32
	SYSTEM_INFO sysinfo;
	MEMORYSTATUSEX memstatus;
	DWORD windows_dw;
#endif
	bubble = (char*)malloc(NON_BIGSTACK_MIN * sizeof(char));
	if(!bubble)
	{
		throw std::runtime_error("Error: Not enough memory for any process!");
	}
	// 	see e.g. http://nadeausoftware.com/articles/2012/09/c_c_tip_how_get_physical_memory_size_system .
#ifdef __APPLE__
	mib[0] = CTL_HW;
	mib[1] = HW_MEMSIZE;
	llxx = 0;
	sztmp = sizeof(int64_t);
	sysctl(mib, 2, &llxx, &sztmp, nullptr, 0);
	llxx /= 1048576;
#else
#ifdef _WIN32
	memstatus.dwLength = sizeof(memstatus);
	GlobalMemoryStatusEx(&memstatus);
	llxx = memstatus.ullTotalPhys / 1048576;
#else
	llxx = ((uint64_t)sysconf(_SC_PHYS_PAGES)) * ((size_t)sysconf(_SC_PAGESIZE)) / 1048576;
#endif
#endif
	if (!llxx) {
		default_alloc_mb = BIGSTACK_DEFAULT_MB;
	} else if (llxx < (BIGSTACK_MIN_MB * 2)){
		default_alloc_mb = BIGSTACK_MIN_MB;
	} else {
		default_alloc_mb = llxx / 2;
	}
	malloc_size_mb = default_alloc_mb;
#ifndef __LP64__
	if (malloc_size_mb > 2047){
		malloc_size_mb = 2047;
	}
#endif
	bigstack_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
	while (!bigstack_ua) {
		malloc_size_mb = (malloc_size_mb * 3) / 4;
		if (malloc_size_mb < BIGSTACK_MIN_MB) {
			malloc_size_mb = BIGSTACK_MIN_MB;
		}
		bigstack_ua = (unsigned char*)malloc(malloc_size_mb * 1048576 * sizeof(char));
		if (malloc_size_mb == BIGSTACK_MIN_MB) {
			throw std::runtime_error("Error: Not enough memory for data structures!");
		}
	}
	bigstack_initial_base = (unsigned char*)round_up_pow2((uintptr_t)bigstack_ua, CACHELINE);
	// as the g suggest... global variable, declared in plink_common
	g_bigstack_base = bigstack_initial_base;
	g_bigstack_end = &(bigstack_initial_base[(malloc_size_mb * 1048576 - (uintptr_t)(bigstack_initial_base - bigstack_ua)) & (~(CACHELINE - ONELU))]);
	free(bubble);
}

*/

double PLINK::calc_lnlike(double known11, double known12, double known21, double known22, double center_ct_d,
		double freq11, double freq12, double freq21, double freq22, double half_hethet_share, double freq11_incr) {
	double lnlike;
	freq11 += freq11_incr;
	freq22 += freq11_incr;
	freq12 += half_hethet_share - freq11_incr;
	freq21 += half_hethet_share - freq11_incr;
	lnlike = center_ct_d * log(freq11 * freq22 + freq12 * freq21);
	if (known11 != 0.0) {
		lnlike += known11 * log(freq11);
	}
	if (known12 != 0.0) {
		lnlike += known12 * log(freq12);
	}
	if (known21 != 0.0) {
		lnlike += known21 * log(freq21);
	}
	if (known22 != 0.0) {
		lnlike += known22 * log(freq22);
	}
	return lnlike;
}

// This is where the magic happens
uint32_t PLINK::em_phase_hethet(double known11, double known12, double known21, double known22, uint32_t center_ct,
		double* freq1x_ptr, double* freq2x_ptr, double* freqx1_ptr, double* freqx2_ptr, double* freq11_ptr,
		uint32_t* onside_sol_ct_ptr) {
	// Returns 1 if at least one SNP is monomorphic over all valid observations;
	// returns 0 otherwise, and fills all frequencies using the maximum
	// likelihood solution to the cubic equation.
	// (We're discontinuing most use of EM phasing since better algorithms have
	// been developed, but the two marker case is mathematically clean and fast
	// enough that it'll probably remain useful as an input for some of those
	// better algorithms...)
	double center_ct_d = (int32_t)center_ct;
	double twice_tot = known11 + known12 + known21 + known22 + 2 * center_ct_d;
	uint32_t sol_start_idx = 0;
	uint32_t sol_end_idx = 1;
	double solutions[3];
	double twice_tot_recip;
	double half_hethet_share;
	double freq11;
	double freq12;
	double freq21;
	double freq22;
	double prod_1122;
	double prod_1221;
	double incr_1122;
	double best_sol;
	double best_lnlike;
	double cur_lnlike;
	double freq1x;
	double freq2x;
	double freqx1;
	double freqx2;
	double lbound;
	double dxx;
	uint32_t cur_sol_idx;
	// shouldn't have to worry about subtractive cancellation problems here
	if (twice_tot == 0.0) {
		return 1;
	}
	twice_tot_recip = 1.0 / twice_tot;
	freq11 = known11 * twice_tot_recip;
	freq12 = known12 * twice_tot_recip;
	freq21 = known21 * twice_tot_recip;
	freq22 = known22 * twice_tot_recip;
	prod_1122 = freq11 * freq22;
	prod_1221 = freq12 * freq21;
	half_hethet_share = center_ct_d * twice_tot_recip;
	// the following four values should all be guaranteed nonzero except in the
	// NAN case
	freq1x = freq11 + freq12 + half_hethet_share;
	freq2x = 1.0 - freq1x;
	freqx1 = freq11 + freq21 + half_hethet_share;
	freqx2 = 1.0 - freqx1;
	if (center_ct) {
		if ((prod_1122 != 0.0) || (prod_1221 != 0.0)) {
			sol_end_idx = cubic_real_roots(0.5 * (freq11 + freq22 - freq12 - freq21 - 3 * half_hethet_share),
					0.5 * (prod_1122 + prod_1221 + half_hethet_share * (freq12 + freq21 - freq11 - freq22 + half_hethet_share)),
					-0.5 * half_hethet_share * prod_1122, solutions);
			while (sol_end_idx && (solutions[sol_end_idx - 1] > half_hethet_share + SMALLISH_EPSILON))
			{
				sol_end_idx--;
				assert(sol_end_idx && sol_end_idx-1 >= 0);
			}
			while ((sol_start_idx < sol_end_idx) && (solutions[sol_start_idx] < -SMALLISH_EPSILON))
			{
				sol_start_idx++;
				assert((sol_start_idx < sol_end_idx) &&sol_start_idx < 3);
			}
			if (sol_start_idx == sol_end_idx)
			{
				// Lost a planet Master Obi-Wan has.  How embarrassing...
				// lost root must be a double root at one of the boundary points, just
				// check their likelihoods
				sol_start_idx = 0;
				sol_end_idx = 2;
				solutions[0] = 0;
				solutions[1] = half_hethet_share;
			}
			else
			{
				if (solutions[sol_start_idx] < 0)
				{
					solutions[sol_start_idx] = 0;
				}
				// checking here
				if (solutions[sol_end_idx-1] > half_hethet_share)
				{
					solutions[sol_end_idx-1] = half_hethet_share;
				}
			}
		} else {
			solutions[0] = 0;
			if ((freq22 + SMALLISH_EPSILON < half_hethet_share + freq21) && (freq21 + SMALLISH_EPSILON < half_hethet_share + freq22)) {
				sol_end_idx = 3;
				solutions[1] = (half_hethet_share + freq21 - freq22) * 0.5;
				solutions[2] = half_hethet_share;
			} else {
				sol_end_idx = 2;
				solutions[1] = half_hethet_share;
			}
		}
		best_sol = solutions[sol_start_idx];
		if (sol_end_idx > sol_start_idx + 1) {
			// select largest log likelihood
			best_lnlike = calc_lnlike(known11, known12, known21, known22, center_ct_d, freq11, freq12, freq21, freq22, half_hethet_share, best_sol);
			cur_sol_idx = sol_start_idx + 1;
			do {
				incr_1122 = solutions[cur_sol_idx];
				cur_lnlike = calc_lnlike(known11, known12, known21, known22, center_ct_d, freq11, freq12, freq21, freq22, half_hethet_share, incr_1122);
				if (cur_lnlike > best_lnlike) {
					cur_lnlike = best_lnlike;
					best_sol = incr_1122;
				}
			} while (++cur_sol_idx < sol_end_idx);
		}
		if (onside_sol_ct_ptr && (sol_end_idx > sol_start_idx + 1)) {
			if (freqx1 * freq1x >= freq11) {
				dxx = freq1x * freqx1 - freq11;
				if (dxx > half_hethet_share) {
					dxx = half_hethet_share;
				}
			} else {
				dxx = 0.0;
			}
			// okay to NOT count suboptimal boundary points because they don't permit
			// direction changes within the main interval
			// this should exactly match haploview_blocks_classify()'s D sign check
			if ((freq11 + best_sol) - freqx1 * freq1x >= 0.0) {
				if (best_sol > dxx + SMALLISH_EPSILON) {
					lbound = dxx + SMALLISH_EPSILON;
				} else {
					lbound = dxx;
				}
				if (best_sol < half_hethet_share - SMALLISH_EPSILON) {
					half_hethet_share -= SMALLISH_EPSILON;
				}
			} else {
				if (best_sol > SMALLISH_EPSILON) {
					lbound = SMALLISH_EPSILON;
				} else {
					lbound = 0.0;
				}
				if (best_sol < dxx - SMALLISH_EPSILON) {
					half_hethet_share = dxx - SMALLISH_EPSILON;
				} else {
					half_hethet_share = dxx;
				}
			}
			for (cur_sol_idx = sol_start_idx; cur_sol_idx < sol_end_idx; cur_sol_idx++) {
				if (solutions[cur_sol_idx] < lbound) {
					sol_start_idx++;
				}
				if (solutions[cur_sol_idx] > half_hethet_share) {
					break;
				}
			}
			if (cur_sol_idx >= sol_start_idx + 2) {
				*onside_sol_ct_ptr = cur_sol_idx - sol_start_idx;
			}
		}
		freq11 += best_sol;
	} else if ((prod_1122 == 0.0) && (prod_1221 == 0.0)) {
		return 1;
	}
	*freq1x_ptr = freq1x;
	*freq2x_ptr = freq2x;
	*freqx1_ptr = freqx1;
	*freqx2_ptr = freqx2;
	*freq11_ptr = freq11;
	return 0;
}

uint32_t PLINK::em_phase_hethet_nobase(uint32_t* counts, uint32_t is_x1, uint32_t is_x2, double* freq1x_ptr,
		double* freq2x_ptr, double* freqx1_ptr, double* freqx2_ptr, double* freq11_ptr) {
	// if is_x1 and/or is_x2 is set, counts[9]..[17] are male-only counts.
	double known11 = (double)(2 * counts[0] + counts[1] + counts[3]);
	double known12 = (double)(2 * counts[2] + counts[1] + counts[5]);
	double known21 = (double)(2 * counts[6] + counts[3] + counts[7]);
	double known22 = (double)(2 * counts[8] + counts[5] + counts[7]);
	if (is_x1 || is_x2) {
		if (is_x1 && is_x2) {
			known11 -= (double)((int32_t)counts[9]);
			known12 -= (double)((int32_t)counts[11]);
      	  known21 -= (double)((int32_t)counts[15]);
      	  known22 -= (double)((int32_t)counts[17]);
		} else if (is_x1) {
			known11 -= ((double)(2 * counts[9] + counts[10])) * (1.0 - SQRT_HALF);
			known12 -= ((double)(2 * counts[11] + counts[10])) * (1.0 - SQRT_HALF);
			known21 -= ((double)(2 * counts[15] + counts[16])) * (1.0 - SQRT_HALF);
			known22 -= ((double)(2 * counts[17] + counts[16])) * (1.0 - SQRT_HALF);
		} else {
			known11 -= ((double)(2 * counts[9] + counts[12])) * (1.0 - SQRT_HALF);
			known12 -= ((double)(2 * counts[11] + counts[12])) * (1.0 - SQRT_HALF);
			known21 -= ((double)(2 * counts[15] + counts[14])) * (1.0 - SQRT_HALF);
			known22 -= ((double)(2 * counts[17] + counts[14])) * (1.0 - SQRT_HALF);
		}
	}
	return em_phase_hethet(known11, known12, known21, known22, counts[4], freq1x_ptr, freq2x_ptr, freqx1_ptr, freqx2_ptr, freq11_ptr, nullptr);
}


uint32_t PLINK::load_and_split3(uintptr_t* rawbuf, uint32_t unfiltered_sample_ct, uintptr_t* casebuf, uintptr_t* pheno_nm, uintptr_t* pheno_c, uint32_t case_ctv, uint32_t ctrl_ctv, uint32_t do_reverse, uint32_t is_case_only, uintptr_t* nm_info_ptr)
{
	uintptr_t* rawbuf_end = &(rawbuf[unfiltered_sample_ct / BITCT2]);
	uintptr_t* ctrlbuf = &(casebuf[3 * case_ctv]);
	uintptr_t case_words[4];
	uintptr_t ctrl_words[4];
	uint32_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
	uint32_t case_rem = 0;
	uint32_t ctrl_rem = 0;
	uint32_t read_shift_max = BITCT2;
	uint32_t sample_uidx = 0;
	uint32_t offset0_case = do_reverse * 2 * case_ctv;
	uint32_t offset2_case = (1 - do_reverse) * 2 * case_ctv;
	uint32_t offset0_ctrl = do_reverse * 2 * ctrl_ctv;
	uint32_t offset2_ctrl = (1 - do_reverse) * 2 * ctrl_ctv;
	uint32_t read_shift;
	uintptr_t read_word;
	uintptr_t ulii;

	case_words[0] = 0;
	case_words[1] = 0;
	case_words[2] = 0;
	case_words[3] = 0;
	ctrl_words[0] = 0;
	ctrl_words[1] = 0;
	ctrl_words[2] = 0;
	ctrl_words[3] = 0;
	while (1) {
		while (rawbuf < rawbuf_end) {
			read_word = *rawbuf++;
			for (read_shift = 0; read_shift < read_shift_max; sample_uidx++, read_shift++) {
				if (is_set(pheno_nm, sample_uidx)) {
					ulii = read_word & 3;
					if (is_set(pheno_c, sample_uidx)) { // Both is_set is always true, because dummy_nm is set
						case_words[ulii] |= ONELU << case_rem;
						if (++case_rem == BITCT) {
							casebuf[offset0_case] = case_words[0];
							casebuf[case_ctv] = case_words[2];
							casebuf[offset2_case] = case_words[3];
							casebuf++;
							case_words[0] = 0;
							case_words[2] = 0;
							case_words[3] = 0;
							case_rem = 0;
						}
					} else if (!is_case_only) {
						ctrl_words[ulii] |= ONELU << ctrl_rem;
						if (++ctrl_rem == BITCT) {
							ctrlbuf[offset0_ctrl] = ctrl_words[0];
							ctrlbuf[ctrl_ctv] = ctrl_words[2];
							ctrlbuf[offset2_ctrl] = ctrl_words[3];
							ctrlbuf++;
							ctrl_words[0] = 0;
							ctrl_words[2] = 0;
							ctrl_words[3] = 0;
							ctrl_rem = 0;
						}
					}
				}
				read_word >>= 2;
			}
		}
		if (sample_uidx == unfiltered_sample_ct) {
			if (case_rem) {
				casebuf[offset0_case] = case_words[0];
				casebuf[case_ctv] = case_words[2];
				casebuf[offset2_case] = case_words[3];
			}
			if (ctrl_rem) {
				ctrlbuf[offset0_ctrl] = ctrl_words[0];
				ctrlbuf[ctrl_ctv] = ctrl_words[2];
				ctrlbuf[offset2_ctrl] = ctrl_words[3];
			}
			ulii = 3;
			if (case_words[1]) {
				ulii -= 1;
			}
			if (ctrl_words[1]) {
				ulii -= 2;
			}
			*nm_info_ptr = ulii;
			return 0;
		}
		rawbuf_end++;
		read_shift_max = unfiltered_sample_ct % BITCT2;
  }
}



void PLINK::two_locus_count_table(uintptr_t* lptr1, uintptr_t* lptr2, uint32_t* counts_3x3, uint32_t sample_ctv3,
		uint32_t is_zmiss2) {
#ifdef __LP64__
  uint32_t uii;
  fill_uint_zero(9, counts_3x3);
  if (!is_zmiss2) {
    two_locus_3x3_tablev((__m128i*)lptr1, (__m128i*)lptr2, counts_3x3, sample_ctv3 / 2, 3);
  } else {
    two_locus_3x3_tablev((__m128i*)lptr2, (__m128i*)lptr1, counts_3x3, sample_ctv3 / 2, 2);
    uii = counts_3x3[1];
    counts_3x3[1] = counts_3x3[3];
    counts_3x3[3] = uii;
    counts_3x3[6] = counts_3x3[2];
    counts_3x3[7] = counts_3x3[5];
  }
#else
  counts_3x3[0] = popcount_longs_intersect(lptr2, lptr1, sample_ctv3);
  counts_3x3[3] = popcount_longs_intersect(lptr2, &(lptr1[sample_ctv3]), sample_ctv3);
  counts_3x3[6] = popcount_longs_intersect(lptr2, &(lptr1[2 * sample_ctv3]), sample_ctv3);
  lptr2 = &(lptr2[sample_ctv3]);
  counts_3x3[1] = popcount_longs_intersect(lptr2, lptr1, sample_ctv3);
  counts_3x3[4] = popcount_longs_intersect(lptr2, &(lptr1[sample_ctv3]), sample_ctv3);
  counts_3x3[7] = popcount_longs_intersect(lptr2, &(lptr1[2 * sample_ctv3]), sample_ctv3);
  if (!is_zmiss2) {
    lptr2 = &(lptr2[sample_ctv3]);
    counts_3x3[2] = popcount_longs_intersect(lptr2, lptr1, sample_ctv3);
    counts_3x3[5] = popcount_longs_intersect(lptr2, &(lptr1[sample_ctv3]), sample_ctv3);
    counts_3x3[8] = popcount_longs_intersect(lptr2, &(lptr1[2 * sample_ctv3]), sample_ctv3);
  }
#endif
}

void PLINK::two_locus_count_table_zmiss1(uintptr_t* lptr1, uintptr_t* lptr2, uint32_t* counts_3x3,
		uint32_t sample_ctv3, uint32_t is_zmiss2) {

#ifdef __LP64__
  fill_uint_zero(6, counts_3x3);
  if (is_zmiss2) {
    two_locus_3x3_zmiss_tablev((__m128i*)lptr1, (__m128i*)lptr2, counts_3x3, sample_ctv3 / 2);
  } else {
    two_locus_3x3_tablev((__m128i*)lptr1, (__m128i*)lptr2, counts_3x3, sample_ctv3 / 2, 2);
  }
#else
  counts_3x3[0] = popcount_longs_intersect(lptr1, lptr2, sample_ctv3);
  counts_3x3[1] = popcount_longs_intersect(lptr1, &(lptr2[sample_ctv3]), sample_ctv3);
  if (!is_zmiss2) {
    counts_3x3[2] = popcount_longs_intersect(lptr1, &(lptr2[2 * sample_ctv3]), sample_ctv3);
    counts_3x3[5] = popcount_longs_intersect(&(lptr1[sample_ctv3]), &(lptr2[2 * sample_ctv3]), sample_ctv3);
  }
  lptr1 = &(lptr1[sample_ctv3]);
  counts_3x3[3] = popcount_longs_intersect(lptr1, lptr2, sample_ctv3);
  counts_3x3[4] = popcount_longs_intersect(lptr1, &(lptr2[sample_ctv3]), sample_ctv3);
#endif
}


#ifdef __LP64__
void PLINK::two_locus_3x3_tablev(__m128i* vec1, __m128i* vec2, uint32_t* counts_3x3, uint32_t sample_ctv6,
		uint32_t iter_ct) {
  const __m128i m1 = {FIVEMASK, FIVEMASK};
  const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
  const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
  __m128i* vec20;
  __m128i* vec21;
  __m128i* vec22;
  __m128i* vend1;
  __m128i loader1;
  __m128i loader20;
  __m128i loader21;
  __m128i loader22;
  __m128i count10;
  __m128i count11;
  __m128i count12;
  __m128i count20;
  __m128i count21;
  __m128i count22;
  __univec acc0;
  __univec acc1;
  __univec acc2;
  uint32_t ct;
  uint32_t ct2;
  while (iter_ct--) {
    ct = sample_ctv6;
    vec20 = vec2;
    vec21 = &(vec20[sample_ctv6]);
    vec22 = &(vec20[2 * sample_ctv6]);
    while (ct >= 30) {
      ct -= 30;
      vend1 = &(vec1[30]);
      acc0.vi = _mm_setzero_si128();
      acc1.vi = _mm_setzero_si128();
      acc2.vi = _mm_setzero_si128();
      do {
      two_locus_3x3_tablev_outer:
	loader1 = *vec1++;
	loader20 = *vec20++;
	loader21 = *vec21++;
	loader22 = *vec22++;
	count10 = _mm_and_si128(loader1, loader20);
	count11 = _mm_and_si128(loader1, loader21);
	count12 = _mm_and_si128(loader1, loader22);
	count10 = _mm_sub_epi64(count10, _mm_and_si128(_mm_srli_epi64(count10, 1), m1));
	count11 = _mm_sub_epi64(count11, _mm_and_si128(_mm_srli_epi64(count11, 1), m1));
	count12 = _mm_sub_epi64(count12, _mm_and_si128(_mm_srli_epi64(count12, 1), m1));
      two_locus_3x3_tablev_two_left:
        // unlike the zmiss variant, this apparently does not suffer from
	// enough register spill to justify shrinking the inner loop
	loader1 = *vec1++;
	loader20 = *vec20++;
	loader21 = *vec21++;
	loader22 = *vec22++;
	count20 = _mm_and_si128(loader1, loader20);
	count21 = _mm_and_si128(loader1, loader21);
	count22 = _mm_and_si128(loader1, loader22);
	count20 = _mm_sub_epi64(count20, _mm_and_si128(_mm_srli_epi64(count20, 1), m1));
	count21 = _mm_sub_epi64(count21, _mm_and_si128(_mm_srli_epi64(count21, 1), m1));
	count22 = _mm_sub_epi64(count22, _mm_and_si128(_mm_srli_epi64(count22, 1), m1));
      two_locus_3x3_tablev_one_left:
	loader1 = *vec1++;
	loader20 = *vec20++;
	loader21 = _mm_and_si128(loader1, loader20); // half1
	loader22 = _mm_and_si128(_mm_srli_epi64(loader21, 1), m1); // half2
	count10 = _mm_add_epi64(count10, _mm_and_si128(loader21, m1));
	count20 = _mm_add_epi64(count20, loader22);
	loader20 = *vec21++;
	loader21 = _mm_and_si128(loader1, loader20);
	loader22 = _mm_and_si128(_mm_srli_epi64(loader21, 1), m1);
	count11 = _mm_add_epi64(count11, _mm_and_si128(loader21, m1));
	count21 = _mm_add_epi64(count21, loader22);
	loader20 = *vec22++;
	loader21 = _mm_and_si128(loader1, loader20);
	loader22 = _mm_and_si128(_mm_srli_epi64(loader21, 1), m1);
	count12 = _mm_add_epi64(count12, _mm_and_si128(loader21, m1));
	count22 = _mm_add_epi64(count22, loader22);

	count10 = _mm_add_epi64(_mm_and_si128(count10, m2), _mm_and_si128(_mm_srli_epi64(count10, 2), m2));
	count11 = _mm_add_epi64(_mm_and_si128(count11, m2), _mm_and_si128(_mm_srli_epi64(count11, 2), m2));
	count12 = _mm_add_epi64(_mm_and_si128(count12, m2), _mm_and_si128(_mm_srli_epi64(count12, 2), m2));
	count10 = _mm_add_epi64(count10, _mm_add_epi64(_mm_and_si128(count20, m2), _mm_and_si128(_mm_srli_epi64(count20, 2), m2)));
	count11 = _mm_add_epi64(count11, _mm_add_epi64(_mm_and_si128(count21, m2), _mm_and_si128(_mm_srli_epi64(count21, 2), m2)));
	count12 = _mm_add_epi64(count12, _mm_add_epi64(_mm_and_si128(count22, m2), _mm_and_si128(_mm_srli_epi64(count22, 2), m2)));
	acc0.vi = _mm_add_epi64(acc0.vi, _mm_add_epi64(_mm_and_si128(count10, m4), _mm_and_si128(_mm_srli_epi64(count10, 4), m4)));
	acc1.vi = _mm_add_epi64(acc1.vi, _mm_add_epi64(_mm_and_si128(count11, m4), _mm_and_si128(_mm_srli_epi64(count11, 4), m4)));
	acc2.vi = _mm_add_epi64(acc2.vi, _mm_add_epi64(_mm_and_si128(count12, m4), _mm_and_si128(_mm_srli_epi64(count12, 4), m4)));
      } while (vec1 < vend1);
      const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
      acc0.vi = _mm_add_epi64(_mm_and_si128(acc0.vi, m8), _mm_and_si128(_mm_srli_epi64(acc0.vi, 8), m8));
      acc1.vi = _mm_add_epi64(_mm_and_si128(acc1.vi, m8), _mm_and_si128(_mm_srli_epi64(acc1.vi, 8), m8));
      acc2.vi = _mm_add_epi64(_mm_and_si128(acc2.vi, m8), _mm_and_si128(_mm_srli_epi64(acc2.vi, 8), m8));
      counts_3x3[0] += ((acc0.u8[0] + acc0.u8[1]) * 0x1000100010001LLU) >> 48;
      counts_3x3[1] += ((acc1.u8[0] + acc1.u8[1]) * 0x1000100010001LLU) >> 48;
      counts_3x3[2] += ((acc2.u8[0] + acc2.u8[1]) * 0x1000100010001LLU) >> 48;
    }
    if (ct) {
      vend1 = &(vec1[ct]);
      ct2 = ct % 3;
      acc0.vi = _mm_setzero_si128();
      acc1.vi = _mm_setzero_si128();
      acc2.vi = _mm_setzero_si128();
      ct = 0;
      if (ct2) {
	count10 = _mm_setzero_si128();
	count11 = _mm_setzero_si128();
	count12 = _mm_setzero_si128();
	if (ct2 == 2) {
	  goto two_locus_3x3_tablev_two_left;
	}
	count20 = _mm_setzero_si128();
	count21 = _mm_setzero_si128();
	count22 = _mm_setzero_si128();
	goto two_locus_3x3_tablev_one_left;
      }
      goto two_locus_3x3_tablev_outer;
    }
    counts_3x3 = &(counts_3x3[3]);
  }
}

#endif

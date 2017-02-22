/*
 * plink.cpp
 *
 *  Created on: 19 Feb 2017
 *      Author: shingwan
 */

#include "plink.hpp"
std::vector<std::string> PLINK::g_chr_list;


PLINK::PLINK(std::string prefix, const size_t thread, const catelog &inclusion):m_thread(thread){
	// TODO Auto-generated constructor stub

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

	m_marker_ct = m_unfiltered_marker_ct - m_marker_exclude_ct; // seems reasonable
	m_unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
	m_unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
	retval = load_bed();
	fprintf(stderr, "%zu people (%zu males, %zu females) loaded from .fam\n", m_unfiltered_sample_ct, m_num_male, m_num_female);
	fprintf(stderr, "%zu variants included\n", m_marker_ct);
}

PLINK::~PLINK() {
	// TODO Auto-generated destructor stub
	// Unfortunately, as I only have partial understanding of the plink code
	// and are far worst a programmer when compared to Chris, I don't know
	// how Chris free the memory in plink. So we will have to live with
	// the fact that the plink code will have memory leak
	delete [] m_marker_exclude;
	delete [] m_sex_male;
	delete [] m_founder_info;
	delete [] m_sample_exclude;
	delete [] m_marker_reverse;
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
	}
	uii = BITCT_TO_WORDCT(m_unfiltered_marker_ct);
	m_marker_reverse = new uintptr_t[uii];
	std::memset(m_marker_reverse, 0x0, uii*sizeof(uintptr_t));
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
					std::cerr << token.size() << std::endl;
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
	return 0;
}


void PLINK::start_clumping(boost::ptr_vector<SNP> &snp_list, double p_threshold, double r2_threshold,
		size_t kb_threshold, double proxyy){
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
	size_t bp_of_core =0;
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
				if (m_bedfile != NULL) fclose(m_bedfile);
				std::string bed_name = prev_file+".bed";
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
        uintptr_t* genotype = new uintptr_t[unfiltered_sample_ctv2];
        std::memset(genotype, 0x0, unfiltered_sample_ctv2*sizeof(uintptr_t));
        uintptr_t* tmp_genotype = new uintptr_t[unfiltered_sample_ctv2];
        std::memset(tmp_genotype, 0x0, unfiltered_sample_ctv2*sizeof(uintptr_t));
        if(load_and_collapse_incl(m_unfiltered_sample_ct, m_founder_ct, m_founder_info, final_mask,
        		IS_SET(m_marker_reverse, cur_line_num), m_bedfile, tmp_genotype, genotype))
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

	fprintf(stderr, "\rClumping Progress: %03.2f%%\n\n", 100.0);
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


void PLINK::clump_thread(const size_t c_core_index, const std::deque<size_t> &c_clump_snp_index,
		boost::ptr_vector<SNP> &snp_list, const double c_r2_threshold)
{

	// do this without the clumping first
	uintptr_t founder_ctl = BITCT_TO_WORDCT(m_founder_ct);

	uint32_t founder_ctv3 = BITCT_TO_ALIGNED_WORDCT(m_founder_ct);
	uint32_t founder_ctsplit = 3 * founder_ctv3; // Required

	uintptr_t ulii = founder_ctsplit * sizeof(intptr_t) + 2 * sizeof(int32_t) + marker_idx2_maxw * 2 * sizeof(double);


	uintptr_t* geno1 = new uintptr_t[founder_ctsplit];
	std::memset(geno1, 0x0, founder_ctsplit*sizeof(uintptr_t))

	uintptr_t* dummy_nm = new uintptr_t[founder_ctl];
	std::memset(dummy_nm, ~0, founder_ctl*sizeof(uintptr_t)); // set all bits to 1

	load_and_split3(nullptr, m_genotype[c_core_index], m_founder_ct,
			geno_1, dummy_nm, dummy_nm, founder_ctv3, 0, 0, 1, &ulii);

	if (ulii == 3) {
		SET_BIT(block_idx1, g_epi_zmiss1); // some missingness observed
	}
	/*
	size_t snp_in_region = c_clump_snp_index.size();
    if(snp_in_region <=1 ) return; // nothing to do
    std::vector<std::thread> thread_store;
    if((snp_in_region-1) < m_thread)
    {
        for(size_t i_snp = 0; i_snp < c_snp_index.size(); ++i_snp)
        {
            if(c_snp_index[i_snp]!=c_core_index)
            	{
            		thread_store.push_back(std::thread(&PLINK::compute_clump, this,
            				c_core_index,i_snp, i_snp+1, std::ref(snp_list), std::cref(c_clump_snp_index),
							c_r2_threshold));
            	}

        }
    }
    else
    {
        int num_snp_per_thread =(int)(snp_in_region) / (int)m_thread;  //round down
        int remain = (int)(snp_in_region) % (int)m_thread;
        int cur_start = 0;
        int cur_end = num_snp_per_thread;
        for(size_t i_thread = 0; i_thread < m_thread; ++i_thread)
        {
            thread_store.push_back(std::thread(&PLINK::compute_clump, this, c_core_index, cur_start,
            		cur_end+(remain>0), std::ref(snp_list), std::cref(c_clump_snp_index),c_r2_threshold ));
            cur_start = cur_end+(remain>0);
            cur_end+=num_snp_per_thread+(remain>0);
            if(cur_end>snp_in_region) cur_end =snp_in_region;
            remain--;
        }
    }
    for(auto &&thread_runner : thread_store) thread_runner.join();
    thread_store.clear();
    */
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
		clump_thread(core_snp_index, snp_index, snp_list, r2_threshold);
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
			sol_end_idx = cubic_real_roots(0.5 * (freq11 + freq22 - freq12 - freq21 - 3 * half_hethet_share), 0.5 * (prod_1122 + prod_1221 + half_hethet_share * (freq12 + freq21 - freq11 - freq22 + half_hethet_share)), -0.5 * half_hethet_share * prod_1122, solutions);
			while (sol_end_idx && (solutions[sol_end_idx - 1] > half_hethet_share + SMALLISH_EPSILON)) {
				sol_end_idx--;
			}
			while ((sol_start_idx < sol_end_idx) && (solutions[sol_start_idx] < -SMALLISH_EPSILON)) {
				sol_start_idx++;
			}
			if (sol_start_idx == sol_end_idx) {
				// Lost a planet Master Obi-Wan has.  How embarrassing...
				// lost root must be a double root at one of the boundary points, just
				// check their likelihoods
				sol_start_idx = 0;
				sol_end_idx = 2;
				solutions[0] = 0;
				solutions[1] = half_hethet_share;
			} else {
				if (solutions[sol_start_idx] < 0) {
					solutions[sol_start_idx] = 0;
				}
				if (solutions[sol_end_idx] > half_hethet_share) {
					solutions[sol_end_idx] = half_hethet_share;
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

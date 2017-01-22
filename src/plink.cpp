//
//  plink.cpp
//  plink
//
//  Created by Shing Wan Choi on 18/08/2016.
//  Copyright © 2016 Shing Wan Choi. All rights reserved.
//

#include "plink.hpp"

#define MULTIPLEX_LD 1920
#define BITCT 64
#define BITCT2 (BITCT / 2)

std::mutex PLINK::clump_mtx;

void PLINK::initialize()
{
    std::string fam_name = m_prefix+".fam";
    // Start processing the fam file
    if(m_prefix.find("#")!=std::string::npos)
    {
        if(m_chr_list.size()==0)
        {
            std::string error_message = "# is reserved for chromosome number. Chromosome information must be provided in order to use the chromosome separated PLINK file!";
            throw std::runtime_error(error_message);
        }
        for(auto &&chr: m_chr_list)
        {
            std::string name = m_prefix;
            misc::replace_substring(name, "#", chr);;
            m_names.push_back(name);
        }
        misc::replace_substring(fam_name, "#", m_chr_list.front());
    }
    else
    {
        m_names.push_back(m_prefix);
    }
    std::cerr << "name get?" << std::endl;
//  This should be ok for sample size as that will always be the same for all file
    std::ifstream fam;
    fam.open(fam_name.c_str());
    if(!fam.is_open())
    {
        std::string error_message = "Cannot open fam file: "+fam_name;
        throw std::runtime_error(error_message);
    }
    std::cerr << "Fam open" << std::endl;
    std::string line;
    while(std::getline(fam, line))
        if(!misc::trimmed(line).empty()) m_num_sample++;
    fam.close();
    std::cerr << "Fam counted" << std::endl;
    // Check whether if the bed file is correct
    for(auto &&bed:m_names)
    {
        std::string bed_name = bed+".bed";
        bool snp_major = openPlinkBinaryFile(bed_name, m_bed);
        if(!snp_major)
        {
            std::string error_message = "Currently does not support sample major format\n";
            error_message.append(bed+" is in sample major format");
            throw std::runtime_error(error_message);
        }
        m_bed.close();
    }
    std::cerr << "bed check" << std::endl;
    std::string bed_name = m_names.front()+".bed";
    openPlinkBinaryFile(bed_name, m_bed);
    // Check whether if the bim file is correct
    std::cerr << "File open" << std::endl;
    for(auto &&bim: m_names)
    {
        std::string bim_name = bim+".bim";
        m_num_snp.push_back(0);
        m_bim.open(bim_name.c_str());
        if(!m_bim.is_open())
        {
            std::string error_message = "Cannot open bim file: "+bim;
            throw std::runtime_error(error_message);
        }
        std::cerr << "start getting the snp id" << std::endl;
        std::cerr << bim_name.c_str() << std::endl;
        while(std::getline(m_bim, line))
        {
        	std::cerr << "Adding stuff?" << std::endl;
            if(!misc::trimmed(line).empty()) m_num_snp.back()++;
            std::cerr << "Check size: " << misc::split(line).size() << std::endl;
            std::cerr << "Size ok" << std::endl;
            std::vector<std::string> token = misc::split(line);
            if(token.size() > 2)
            {
            	std::cerr << "Adding" << std::endl;
            	 m_snp_id.push_back(token[1]);
            }
            std::cerr << "Finished push back" << std::endl;
        }
        m_bim.close();
        std::cerr << "bim completed" << std::endl;
    }
    std::string bim_name = m_names.front()+".bim";
    m_bim.open(bim_name.c_str());
    m_name_index=0;
    std::cerr << "Cal bytes" << std::endl;
    m_num_bytes=ceil((double)m_num_sample/4.0);
    std::cerr << "???" << std::endl;
    m_required_bit = m_num_sample*2;
    m_snp_iter=0;
    m_init = true;
}

void PLINK::clump_initialize(const std::unordered_map<std::string, size_t> &inclusion)
{
    std::string fam_name = m_prefix+".fam";
    if(m_prefix.find("#")!=std::string::npos)
    {
        if(m_chr_list.size()==0)
        {
            std::string error_message = "# is reserved for chromosome number. Chromosome information must be provided in order to use the chromosome separated PLINK file!";
            throw std::runtime_error(error_message);
        }
        for(auto &&chr: m_chr_list)
        {
            std::string name = m_prefix;
            misc::replace_substring(name, "#", chr);
            m_names.push_back(name);
        }
        misc::replace_substring(fam_name, "#", m_chr_list.front());
    }
    else m_names.push_back(m_prefix);
    // Start processing the fam file

    std::ifstream fam;
    fam.open(fam_name.c_str());
    if(!fam.is_open())
    {
        std::string error_message = "Cannot open fam file: "+fam_name;
        throw std::runtime_error(error_message);
    }
    std::string line;
    while(std::getline(fam, line)) if(!misc::trimmed(line).empty()) m_num_sample++;
    fam.close();

    // Check whether if the bed file is correct
    for(auto &&bed:m_names)
    {
        std::string bed_name = bed+".bed";
        bool snp_major = openPlinkBinaryFile(bed_name, m_bed);
        if(!snp_major)
        {
            std::string error_message = "Currently does not support sample major format\n";
            error_message.append(bed+" is in sample major format");
            throw std::runtime_error(error_message);
        }
        m_bed.close();
    }
    std::string bed_name = m_names.front()+".bed";
    openPlinkBinaryFile(bed_name, m_bed);
    // Check whether if the bim file is correct
    size_t cur_num_line=0;
    for(auto &&bim:m_names)
    {
        cur_num_line=0;
        std::string bim_name = bim+".bim";
        m_bim.open(bim_name.c_str());
        if(!m_bim.is_open())
        {
            std::string error_message = "Cannot open bim file: "+bim_name;
            throw std::runtime_error(error_message);
        }
        m_num_snp.push_back(0);
        while(std::getline(m_bim, line))
        {
            cur_num_line++;
            if(!misc::trimmed(line).empty())
            {
                std::vector<std::string> token = misc::split(line);
                if(token.size() < 6) throw std::runtime_error("Malformed bim file. Should contain at least 6 column");
                m_num_snp.back()++;
                std::string rsid = token[+BIM::RS];
                m_snp_id.push_back(rsid);
                if(inclusion.find(rsid)!=inclusion.end())
                {
                    m_clump_ref.push_back(file_info(bim, cur_num_line, inclusion.at(rsid)));
                }
            }
        }
        m_bim.close();
    }
    std::string bim_name = m_names.front()+".bim";
    m_bim.open(bim_name.c_str());
    m_name_index=0;
    m_num_bytes=ceil((double)m_num_sample/4.0);
    m_required_bit = m_num_sample*2;
    m_snp_iter=0;
    m_init = true;
}

void PLINK::clump_initialize(std::unordered_map<std::string, size_t> &inclusion, boost::ptr_vector<SNP> &snp_list, const std::unordered_map<std::string, size_t> &c_snp_index)
{
    std::string fam_name = m_prefix+".fam";
    if(m_prefix.find("#")!=std::string::npos)
    {
        if(m_chr_list.size()==0)
        {
            std::string error_message = "# is reserved for chromosome number. Chromosome information must be provided in order to use the chromosome separated PLINK file!";
            throw std::runtime_error(error_message);
        }
        for(auto &&chr: m_chr_list)
        {
            std::string name = m_prefix;
            misc::replace_substring(name, "#", chr);
            m_names.push_back(name);
        }
        misc::replace_substring(fam_name, "#", m_chr_list.front());
    }
    else m_names.push_back(m_prefix);
    // Start processing the fam file

    std::ifstream fam;
    fam.open(fam_name.c_str());
    if(!fam.is_open())
    {
        std::string error_message = "Cannot open fam file: "+fam_name;
        throw std::runtime_error(error_message);
    }
    std::string line;
    while(std::getline(fam, line)) if(!misc::trimmed(line).empty()) m_num_sample++;
    fam.close();

    // Check whether if the bed file is correct
    for(auto &&bed:m_names)
    {
        std::string bed_name = bed+".bed";
        bool snp_major = openPlinkBinaryFile(bed_name, m_bed);
        if(!snp_major)
        {
            std::string error_message = "Currently does not support sample major format\n";
            error_message.append(bed+" is in sample major format");
            throw std::runtime_error(error_message);
        }
        m_bed.close();
    }
    std::string bed_name = m_names.front()+".bed";
    openPlinkBinaryFile(bed_name, m_bed);
    // Check whether if the bim file is correct
    size_t num_ambig=0, not_found=0;
    size_t num_line=0, cur_num_line=0;
    std::unordered_map<std::string, bool> dup_check;
    for(auto &&bim:m_names)
    {
        cur_num_line=0;
        std::string bim_name = bim+".bim";
        m_bim.open(bim_name.c_str());
        if(!m_bim.is_open())
        {
            std::string error_message = "Cannot open bim file: "+bim_name;
            throw std::runtime_error(error_message);
        }
        m_num_snp.push_back(0);
        while(std::getline(m_bim, line))
        {
            num_line++;
            cur_num_line++;
            if(!misc::trimmed(line).empty())
            {
                std::vector<std::string> token = misc::split(line);
                if(token.size() < 6) throw std::runtime_error("Malformed bim file. Should contain at least 6 column");
                m_num_snp.back()++;
                std::string chr = token[+BIM::CHR];
                std::string rsid = token[+BIM::RS];
                m_snp_id.push_back(rsid);
                int loc = -1;
                int temp = 0;
                try
                {
                    temp =misc::convert<int>(token[+BIM::BP]);
                    if(temp < 0)
                    {
                        std::string error_message = "Negative coordinate of SNP in "+bim;
                        throw std::runtime_error(error_message);
                    }
                    loc = temp;
                }
                catch(std::runtime_error &error)
                {
                    std::string error_message = "Non-numeric coordinate of SNP in "+bim;
                    throw std::runtime_error(error_message);
                }
                std::string ref_allele = token[+BIM::A1];
                std::string alt_allele = token[+BIM::A2];
                if(c_snp_index.find(rsid)!=c_snp_index.end() && inclusion.find(rsid)!=inclusion.end()
                        && (dup_check.find(rsid)==dup_check.end() || !dup_check[rsid]) )
                {
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
                        dup_check[rsid] = false;
                    }
                    else
                    {
                        // not ambiguous, now do soft checking
                        size_t index = c_snp_index.at(rsid);
                        bool same = snp_list[index].check_loc(chr, loc, ref_allele, alt_allele);
                        if(snp_list[index].get_loc() == -1 && loc!=-1)
                        {
                            snp_list[index].set_loc(loc);
                        }
                        if(!same && (snp_list[index].get_loc()!=-1 || loc!=-1))
                        {
                            fprintf(stderr, "WARNING: %s differ between LD and base file\n", rsid.c_str());
                            fprintf(stderr, "         It is advised that you check the files are \n");
                            fprintf(stderr, "         From the same genome build\n");
                        }
                        //inclusion[rsid] = index; // Doesn't need to add the index as it should already be there
                        if(dup_check.find(rsid)!=dup_check.end() && !dup_check[rsid]) num_ambig--;
                        m_clump_ref.push_back(file_info(bim, cur_num_line, inclusion[rsid]));
                        dup_check[rsid] = true;
                    }
                }
                else
                {
                    not_found++;
                    if(inclusion.find(rsid)!=inclusion.end()) inclusion.erase(rsid);
                }
            }
        }
        m_bim.close();
    }
    // don't calculate the number of duplication as that might be time consuming
    int num_dup = num_line-dup_check.size();
    if(num_dup!=0) fprintf(stderr, "Number of duplicated SNPs : %d\n", num_dup);
    if(num_ambig != 0)	fprintf(stderr, "Number of ambiguous SNPs  : %zu\n", num_ambig);
    if(not_found != 0)	fprintf(stderr, "Number of SNPs not found  : %zu\n", not_found);
    fprintf(stderr, "Number of SNPs included   : %zu\n", inclusion.size());

    std::string bim_name = m_names.front()+".bim";
    m_bim.open(bim_name.c_str());
    m_name_index=0;
    m_num_bytes=ceil((double)m_num_sample/4.0);
    m_required_bit = m_num_sample*2;
    m_snp_iter=0;
    m_init = true;
}


void PLINK::start_clumping(std::unordered_map<std::string, size_t> &inclusion,
                           boost::ptr_vector<SNP> &snp_list, double p_threshold,
                           double r2_threshold, size_t kb_threshold, double proxy_threshold)
{
    std::deque<size_t> snp_index; // SNPs within the region
    std::string prev_chr="";
    std::string prev_file="";
    size_t read_snps=0; //<- past tense =P
    size_t bp_of_core =0;
    size_t core_genotype_index=0; //index of the core SNP on our genotype deque
    bool require_clump=false; // Whether if the current interval contain the core snp
    for(size_t i_info = 0; i_info < m_clump_ref.size(); ++i_info)
    {
        size_t cur_snp_index = std::get<+FILE_INFO::INDEX>(m_clump_ref[i_info]);
        size_t cur_line_num = std::get<+FILE_INFO::LINE>(m_clump_ref[i_info]);
        std::string cur_chr = snp_list[cur_snp_index].get_chr();
        size_t cur_loc = snp_list[cur_snp_index].get_loc();
        if(prev_chr.empty() || prev_chr.compare(cur_chr)!=0)
        {
			perform_clump(snp_index, snp_list, core_genotype_index, require_clump, p_threshold,
							 r2_threshold, kb_threshold, cur_chr, cur_loc);
			if(prev_file.empty()||prev_file.compare(std::get<+FILE_INFO::FILE>(m_clump_ref[i_info]))!=0)
			{
                prev_file = std::get<+FILE_INFO::FILE>(m_clump_ref[i_info]);
				read_snps = 0;
                if(m_bed.is_open()) m_bed.close();
                if(m_bim.is_open()) m_bim.close();
                std::string bed_name = prev_file+".bed";
                openPlinkBinaryFile(bed_name, m_bed);
			}
			prev_chr = cur_chr;
        }
        else if((cur_loc-bp_of_core) > kb_threshold){
			perform_clump(snp_index, snp_list, core_genotype_index, require_clump, p_threshold,
							 r2_threshold, kb_threshold, cur_chr, cur_loc);
        }
        if((cur_line_num-read_snps)!=0)
        {
            m_bed.seekg((cur_line_num-read_snps)*m_num_bytes, m_bed.cur);
            read_snps=cur_line_num;
        }
		read_snp(1, true);
        read_snps++;
		snp_index.push_back(cur_snp_index);
        if(!require_clump && snp_list[cur_snp_index].get_p_value() < p_threshold)
        {
        		bp_of_core =snp_list[cur_snp_index].get_loc();
        		core_genotype_index=m_genotype.size()-1; // Should store the index on genotype
            require_clump= true;
        }
    }

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
    fprintf(stderr, "Number of SNPs after clumping : %zu\n", inclusion.size());
}

void PLINK::perform_clump(std::deque<size_t> &snp_index, boost::ptr_vector<SNP> &snp_list, size_t &core_snp_index,
						bool &require_clump, double p_threshold, double r2_threshold, size_t kb_threshold,
						std::string next_chr, size_t next_loc)
{
	if(snp_index.size()==0) return; // got nothing to do
    std::string core_chr = snp_list[core_snp_index].get_chr();
    size_t core_loc = snp_list[core_snp_index].get_loc();
    while(require_clump && (core_chr.compare(next_chr)!=0 || (next_loc - core_loc) > kb_threshold))
	{
		clump_thread(core_snp_index, snp_index, snp_list, r2_threshold);
		require_clump = false;
        for(size_t core_finder = core_snp_index+1; core_finder < snp_index.size(); ++core_finder)
		{
			if(snp_list[snp_index[core_finder]].get_p_value() < p_threshold)
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
				if(core_loc-snp_list[snp_index[remover]].get_loc() > kb_threshold) num_remove++;
                else break;
			}
			if(num_remove!=0)
			{
				lerase(num_remove);
				snp_index.erase(snp_index.begin(), snp_index.begin()+num_remove);
				core_snp_index-=num_remove;
			}
		}
    }
    // for this to be true, require_clump must be false, otherwise it will still be within the loop
    if(core_chr.compare(next_chr)!=0)
	{
		// just remove everything
		lerase(m_genotype.size());
		snp_index.clear();
    }
    else if(!require_clump){
    		//remove anything that is too far ahead
    		size_t num_remove = 0;
    		for(auto &&remover : snp_index)
    		{
    			if(next_loc-snp_list[remover].get_loc()>kb_threshold) num_remove++;
    			else break;
    		}
    		if(num_remove!=0)
    		{
    			lerase(num_remove);
    			snp_index.erase(snp_index.begin(), snp_index.begin()+num_remove);
    		}
    }
}



void PLINK::clump_thread(const size_t c_core_index, const std::deque<size_t> &c_snp_index, boost::ptr_vector<SNP> &snp_list, const double c_r2_threshold)
{
	size_t snp_in_region = c_snp_index.size();
    if(snp_in_region <=1 ) return; // nothing to do
    std::vector<std::thread> thread_store;
    if((snp_in_region-1) < m_thread)
    {
        for(size_t i_snp = 0; i_snp < c_snp_index.size(); ++i_snp)
        {
            if(c_snp_index[i_snp]!=c_core_index)
            	{
            		thread_store.push_back(std::thread(&PLINK::compute_clump, this,
            				c_core_index,i_snp, i_snp+1, std::ref(snp_list), std::cref(c_snp_index), c_r2_threshold));
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
            		cur_end+(remain>0), std::ref(snp_list), std::cref(c_snp_index),c_r2_threshold ));
            cur_start = cur_end+(remain>0);
            cur_end+=num_snp_per_thread+(remain>0);
            if(cur_end>snp_in_region) cur_end =snp_in_region;
            remain--;
        }
    }
    for(auto &&thread_runner : thread_store) thread_runner.join();
    thread_store.clear();
}

void PLINK::compute_clump( size_t core_snp_index, size_t i_start, size_t i_end, boost::ptr_vector<SNP> &snp_list,
		const std::deque<size_t> &snp_index_list, const double r2_threshold)
{
    size_t ref_index = snp_index_list[core_snp_index];
    double ref_p_value = snp_list[ref_index].get_p_value();
    std::vector<double> r2_store;
    std::vector<size_t> target_index_store; // index we want to push into the current index
    for(size_t i_snp = i_start; i_snp < i_end && i_snp < snp_index_list.size(); ++i_snp)
    {
        size_t target_index = snp_index_list[i_snp];
        if(i_snp != core_snp_index && snp_list[target_index].get_p_value() > ref_p_value)
        {
            // only calculate r2 if more significant
            double r2 = get_r2(i_snp, core_snp_index);
            if(r2 >= r2_threshold)
            {
            		target_index_store.push_back(target_index);
                r2_store.push_back(r2);
            }
        }
    }
    PLINK::clump_mtx.lock();
    snp_list[ref_index].add_clump(target_index_store);
    snp_list[ref_index].add_clump_r2(r2_store);
    PLINK::clump_mtx.unlock();
}

PLINK::~PLINK()
{
    for(size_t i = 0; i < m_genotype.size(); ++i)
    {
        delete [] m_genotype[i];
        delete [] m_missing[i];
    }
}

#ifdef __LP64__
//This is obtained from plink
double PLINK::get_r2(const size_t i, const size_t j, bool adjust)
{
    uintptr_t founder_ct_mld = (m_num_sample + MULTIPLEX_LD - 1) / MULTIPLEX_LD;
    uint32_t founder_ct_mld_m1 = ((uint32_t)founder_ct_mld) - 1;
    uint32_t founder_ct_mld_rem = (MULTIPLEX_LD / 192) - (founder_ct_mld * MULTIPLEX_LD - m_num_sample) / 192;
    uint32_t fixed_missing_ct;
    uint32_t fixed_non_missing_ct;
    uint32_t non_missing_ct;
    uintptr_t founder_ctwd = m_num_sample / BITCT2;
    uintptr_t founder_ctwd12 = founder_ctwd / 12;
    uintptr_t founder_ctwd12_rem = founder_ctwd - (12 * founder_ctwd12);
    uintptr_t lshift_last = 2 * ((0x7fffffc0 - m_num_sample) % BITCT2);
    long_type* vec1 = m_genotype[j];
    long_type* vec2 = m_genotype[i];
    long_type* mask1 = m_missing[j];
    long_type* mask2 = m_missing[i];
    int32_t dp_result[5];
    fixed_missing_ct = m_num_missing[i];
    fixed_non_missing_ct = m_num_sample - fixed_missing_ct;
    non_missing_ct = fixed_non_missing_ct - m_num_missing[j];
    if (fixed_missing_ct && m_num_missing[j])
    {
        non_missing_ct += ld_missing_ct_intersect(mask1, mask2, founder_ctwd12, founder_ctwd12_rem, lshift_last);
    }
    dp_result[0] = m_num_sample;
    dp_result[1] = -fixed_non_missing_ct;
    dp_result[2] = (int)m_num_missing[j] - (int)m_num_sample;
    dp_result[3] = dp_result[1];
    dp_result[4] = dp_result[2];
    while (founder_ct_mld_m1--)
    {
        ld_dot_prod_batch((__m128i*)vec1, (__m128i*)vec2, (__m128i*)mask1, (__m128i*)mask2, dp_result, MULTIPLEX_LD / 192);
        vec1 = &(vec1[MULTIPLEX_LD / BITCT2]);
        vec2 = &(vec2[MULTIPLEX_LD / BITCT2]);
        mask1 = &(mask1[MULTIPLEX_LD / BITCT2]);
        mask2 = &(mask2[MULTIPLEX_LD / BITCT2]);
    }
    ld_dot_prod_batch((__m128i*)vec1, (__m128i*)vec2, (__m128i*)mask1, (__m128i*)mask2, dp_result, founder_ct_mld_rem);
    double non_missing_ctd = (double)((int32_t)non_missing_ct);
    double dxx = dp_result[1];
    double dyy = dp_result[2];
    double cov12 = dp_result[0] * non_missing_ctd - dxx * dyy;
    dxx = (dp_result[3] * non_missing_ctd + dxx * dxx) * (dp_result[4] * non_missing_ctd + dyy * dyy);
    dxx = (cov12 * cov12) / dxx;
//	if(adjust){
//		std::cerr << m_num_sample << "\t" << dp_result[0] << "\t"<< dp_result[1] << "\t"<< dp_result[2] << "\t"<< dp_result[3] << "\t"<< dp_result[4] << std::endl;
//		exit(-1);
//		size_t n = 0;
//		for(size_t check = 0; check < (m_required_bit /(m_bit_size))+1; ++check){
//			n += __builtin_popcountll(m_missing[i][check]&m_missing[j][check])/2;
//		}
//		dxx = dxx-(1.0-dxx)/(double)((n>2)?n-2:n);
//	}
    return dxx;
}


#else
// This should work for uint32_t but not uint64_t because of the special popcount he used
double PLINK::get_r2(const size_t i, const size_t j)
{
    double r2 =0.0;
    if(i >= m_genotype.size() || j >=m_genotype.size()) throw std::runtime_error("Out of bound error! In R2 calculation.");
    // Crazy stuff of plink here
    size_t range = (m_required_bit /(m_bit_size))+1;
    long_type loader1, loader2, sum1, sum2, sum11, sum12, sum22;
    long_type final_sum1 = 0;
    long_type final_sum2 = 0;
    long_type final_sum11 = 0;
    long_type final_sum22 = 0;
    long_type final_sum12 = 0;
    double return_vals[5];
    return_vals[0] = (double) m_num_sample;
    return_vals[1] = -(double) m_num_missing[j];
    return_vals[2] = -(double) m_num_missing[i];
    return_vals[3] = return_vals[1];
    return_vals[4] = return_vals[2];
    size_t N =0;
    for(size_t i_geno = 0; i_geno < range;)
    {
        loader1 = m_genotype[i][i_geno];
        loader2 = m_genotype[j][i_geno];
        sum1 = m_missing[j][i_geno];
        sum2 = m_missing[i][i_geno];
        i_geno++;
        N+= __builtin_popcountll(sum1&sum2)/2;
        sum12 = (loader1 | loader2) & FIVEMASK;
        sum1 = sum1 & loader1;
        sum2 = sum2 & loader2;
        loader1 = (loader1 ^ loader2) & (AAAAMASK - sum12);
        sum12 = sum12 | loader1;
        sum11 = sum1 & FIVEMASK;
        sum22 = sum2 & FIVEMASK;
        sum1 = (sum1 & THREEMASK) + ((sum1 >> 2) & THREEMASK);
        sum2 = (sum2 & THREEMASK) + ((sum2 >> 2) & THREEMASK);
        sum12 = (sum12 & THREEMASK) + ((sum12 >> 2) & THREEMASK);
        long_type tmp_sum1=0 , tmp_sum2=0;
        if(i_geno < range)
        {
            loader1 = m_genotype[i][i_geno];
            loader2 = m_genotype[j][i_geno];
            tmp_sum1 = m_missing[j][i_geno];
            tmp_sum2 = m_missing[i][i_geno];
            N+= __builtin_popcountll(tmp_sum1&tmp_sum2)/2;
        }
        else
        {
            loader1 = 0;
            loader2 = 0;
        }
        i_geno++;
        long_type tmp_sum12 = (loader1 | loader2) & FIVEMASK;
        tmp_sum1 = tmp_sum1 & loader1;
        tmp_sum2 = tmp_sum2 & loader2;
        loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
        tmp_sum12 = tmp_sum12 | loader1;
        sum11 += tmp_sum1 & FIVEMASK;
        sum22 += tmp_sum2 & FIVEMASK;
        sum1 += (tmp_sum1 & THREEMASK) + ((tmp_sum1 >> 2) & THREEMASK);
        sum2 += (tmp_sum2 & THREEMASK) + ((tmp_sum2 >> 2) & THREEMASK);
        sum12 += (tmp_sum12 & THREEMASK) + ((tmp_sum12 >> 2) & THREEMASK);
        if(i_geno < range)
        {
            loader1 = m_genotype[i][i_geno];
            loader2 = m_genotype[j][i_geno];
            tmp_sum1 = m_missing[j][i_geno];
            tmp_sum2 = m_missing[i][i_geno];
            N+= __builtin_popcountll(tmp_sum1&tmp_sum2)/2;
        }
        else
        {
            loader1=0;
            loader2=0;
            tmp_sum1=0;
            tmp_sum2=0;
        }
        i_geno++;
        tmp_sum12 = (loader1 | loader2) & FIVEMASK;
        tmp_sum1 = tmp_sum1 & loader1;
        tmp_sum2 = tmp_sum2 & loader2;
        loader1 = (loader1 ^ loader2) & (AAAAMASK - tmp_sum12);
        tmp_sum12 = tmp_sum12 | loader1;
        sum11 += tmp_sum1 & FIVEMASK;
        sum22 += tmp_sum2 & FIVEMASK;
        sum1 += (tmp_sum1 & THREEMASK) + ((tmp_sum1 >> 2) & THREEMASK);
        sum2 += (tmp_sum2 & THREEMASK) + ((tmp_sum2 >> 2) & THREEMASK);
        sum11 = (sum11 & THREEMASK) + ((sum11 >> 2) & THREEMASK);
        sum22 = (sum22 & THREEMASK) + ((sum22 >> 2) & THREEMASK);
        sum12 += (tmp_sum12 & THREEMASK) + ((tmp_sum12 >> 2) & THREEMASK);
        sum1 = (sum1 & OFMASK) + ((sum1 >> 4) & OFMASK);
        sum2 = (sum2 & OFMASK) + ((sum2 >> 4) & OFMASK);
        sum11 = (sum11 & OFMASK) + ((sum11 >> 4) & OFMASK);
        sum22 = (sum22 & OFMASK) + ((sum22 >> 4) & OFMASK);
        sum12 = (sum12 & OFMASK) + ((sum12 >> 4) & OFMASK);
        final_sum1 += (sum1 * ONEZEROMASK) >> 24;
        final_sum2 += (sum2 * ONEZEROMASK) >> 24;
        final_sum11 += (sum11 * ONEZEROMASK) >> 24;
        final_sum22 += (sum22 * ONEZEROMASK) >> 24;
        final_sum12 += (sum12 * ONEZEROMASK) >> 24;
    }

    return_vals[0] -= final_sum12;
    return_vals[1] += final_sum1;
    return_vals[2] += final_sum2;
    return_vals[3] += final_sum11;
    return_vals[4] += final_sum22;

    double dxx = return_vals[1];
    double dyy = return_vals[2];
    double n = N;
    double cov12 = return_vals[0] * n - dxx * dyy;
    dxx = (return_vals[3] * n + dxx * dxx) * (return_vals[4] * n + dyy * dyy);
    if(dxx !=0.0) r2 =(cov12 * cov12) / dxx;
    return r2;
}
#endif

void PLINK::lerase(int num)
{
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
        delete [] m_missing[i];
    }
    if(num==m_genotype.size())
    {
        m_genotype.clear();
        m_missing.clear();
//    		m_chr_list.clear();
//    		m_cm_list.clear();
        m_bp_list.clear();
//    		m_ref_allele.clear();
//    		m_alt_allele.clear();
        m_maf.clear();
        m_num_missing.clear();
    }
    else
    {
        m_maf.erase(m_maf.begin(), m_maf.begin()+num);
        m_num_missing.erase(m_num_missing.begin(), m_num_missing.begin()+num);
        m_genotype.erase(m_genotype.begin(), m_genotype.begin()+num);
        m_missing.erase(m_missing.begin(), m_missing.begin()+num);
    }
}

//The return value should be the number of remaining SNPs
int PLINK::read_snp(int num_snp, bool ld)
{
    if(!m_init) throw std::runtime_error("Class uninitialize! Must initialize before use!");
    if(num_snp <= 0)
    {
        std::string error_message = "Number of required SNPs cannot be less than 1: "+std::to_string(num_snp);
        throw std::runtime_error(error_message);
    }
    std::string line;
    //First get the information of the SNPs
    size_t cur_iter = 0;
    while(cur_iter < num_snp)
    {
        for(; m_snp_iter < m_num_snp[m_name_index] && cur_iter<num_snp; ++m_snp_iter)
        {
            cur_iter++;
            char genotype_list[m_num_bytes];
            m_bed.read(genotype_list, m_num_bytes);
            size_t i_genotype = 0;
            size_t total_allele = 0;
            size_t num_missing = 0;
            uintptr_t founder_ct_mld = (m_num_sample + MULTIPLEX_LD - 1) / MULTIPLEX_LD;
            uint32_t founder_ct_mld_m1 = ((uint32_t)founder_ct_mld) - 1;
            uint32_t founder_ct_mld_rem = (MULTIPLEX_LD / 192) - (founder_ct_mld * MULTIPLEX_LD - m_num_sample) / 192;
            size_t range = (founder_ct_mld_m1*(MULTIPLEX_LD / 192)*6+founder_ct_mld_rem*6)*2/(sizeof(long_type)/4);
            long_type *genotype = new long_type[range];
            long_type *missing = new long_type[range];
            std::memset(genotype, 0x0,(range)*sizeof(long_type));
            std::memset(missing, 0x0,(range)*sizeof(long_type));
            for(size_t byte_runner= 0; byte_runner < m_num_bytes;)
            {
#ifdef __LP64__
                long_type current_genotypes = 0ULL;
#else
                long_type current_genotypes=0UL;
#endif
                for(int byte_set = 0; byte_set < sizeof(long_type)/sizeof(char) && byte_runner < m_num_bytes; ++byte_set)
                {
                		long_type current_byte = static_cast<long_type>(genotype_list[byte_runner]) << ((sizeof(long_type)-1)*CHAR_BIT) >> (((sizeof(long_type)-1)-byte_set)*CHAR_BIT);
                    current_genotypes |= current_byte;
                    byte_runner++;
                }
                long_type five_masked_geno = current_genotypes & FIVEMASK;
                long_type inter = (five_masked_geno & (current_genotypes>>1)) ^ five_masked_geno;
                long_type current_missing = inter | (inter << 1);
                missing[i_genotype] = ~current_missing;
                if(!ld){
                		genotype[i_genotype] = current_genotypes;
                }
                else
                {
                    genotype[i_genotype] =(current_genotypes &(five_masked_geno <<1));
                    genotype[i_genotype] |= (five_masked_geno^((current_genotypes &(FIVEMASK*2))>>1));
                    total_allele += __builtin_popcountll(current_genotypes & (~current_missing));
                    num_missing +=__builtin_popcountll(inter); // because inter only contain one bit for each missing
                    i_genotype++;
                }
            }
            m_genotype.push_back(genotype);
            m_missing.push_back(missing);
            double maf= 0.0;
            if(m_required_bit - num_missing != 0) maf = (double)total_allele/((double)m_required_bit-(double)num_missing);
            maf = (maf > 0.5)? 1.0-maf: maf;
            m_maf.push_back(maf);
            m_num_missing.push_back(num_missing);
        }
        if(m_snp_iter >= m_num_snp[m_name_index])
        {
            m_snp_iter = 0;
            m_name_index++;
            m_bed.close();
            m_bim.close();
            if(m_name_index< m_names.size())
            {
                std::string bed_name = m_names[m_name_index]+".bed";
                std::string bim_name = m_names[m_name_index]+".bim";
                openPlinkBinaryFile(bed_name, m_bed);
                m_bim.open(bim_name.c_str());
            }
            else if(cur_iter < num_snp)
            {
                throw std::runtime_error("All SNP read");
            }
            else
            {
                return 0; // everything is read
            }
        }
    }

    return m_num_snp[m_name_index]-m_snp_iter;
}

// Because of the per chromosome method, we will no longer be using the bool and the whole structure of this class
// will change
//This initialization will also perform the filtering and flipping
void PLINK::get_score(const std::vector<p_partition> &partition,
                      const boost::ptr_vector<SNP> &snp_list, std::vector< std::vector<std::pair<std::string, double> > > &prs_score,
                      size_t start_index, size_t end_bound)
{

    size_t prev =0;
//	This allow for consistence at least in this specific use case
    if(m_bed.is_open()) m_bed.close();
    if(m_bim.is_open()) m_bim.close();
    std::string prev_name = "";
    for(size_t i_snp = start_index; i_snp < end_bound; ++i_snp)
    {
        if(prev_name.empty() || prev_name.compare(std::get<+PRS::FILENAME>(partition[i_snp]))!=0)
        {
            m_bed.close();
            prev_name= std::get<+PRS::FILENAME>(partition[i_snp]);
            std::string bed_name = prev_name+".bed";
            openPlinkBinaryFile(bed_name, m_bed);
            prev=0;
        }
        size_t cur_line = std::get<+PRS::LINE>(partition[i_snp]);
        if((cur_line-prev)!=0)
        {
            // Skip snps
            m_bed.seekg((std::get<+PRS::LINE>(partition[i_snp])-prev)*m_num_bytes, m_bed.cur);
            prev=std::get<+PRS::LINE>(partition[i_snp]);
        }
        //read_snp(1, false);
        char genotype_list[m_num_bytes];
        m_bed.read(genotype_list, m_num_bytes);
        if (!m_bed) throw std::runtime_error("Problem with the BED file...has the FAM/BIM file been changed?");
        prev++;
        size_t sample_index = 0;
        int snp_index = std::get<+PRS::INDEX>(partition[i_snp]);
        if(snp_index >= snp_list.size()) throw std::runtime_error("Out of bound! In PRS score calculation");
        for(auto &&byte : genotype_list)
        {
        		size_t geno_bit = 0;
    			int geno_batch = static_cast<int>(byte);
        		while(geno_bit < 7 && sample_index < m_num_sample)
        		{
        			int geno = geno_batch>>geno_bit & 3; // This will access the corresponding genotype
        			if(geno!=1) // Because 01 is coded as missing
        			{
        				for(size_t i_region = 0; i_region < prs_score.size(); ++i_region)
        				{
        					if(snp_list[snp_index].in(i_region))
        					{
        						prs_score[i_region][sample_index].second += snp_list[snp_index].score(geno);
        					}
        				}
        			}
        			sample_index++;
        			geno_bit+=2;
        		}
        }
    }
}

void PLINK::get_score(const std::vector<p_partition> &partition,
                      const boost::ptr_vector<SNP> &snp_list, std::vector< std::vector<std::pair<std::string, double> > > &prs_score,
                      size_t start_index, size_t end_bound, size_t i_region)
{

    size_t prev =0;
//	This allow for consistence at least in this specific use case
    if(m_bed.is_open()) m_bed.close();
    if(m_bim.is_open()) m_bim.close();
    std::string prev_name = "";
    for(size_t i_snp = start_index; i_snp < end_bound; ++i_snp)
    {
        if(prev_name.empty() || prev_name.compare(std::get<+PRS::FILENAME>(partition[i_snp]))!=0)
        {
            m_bed.close();
            prev_name= std::get<+PRS::FILENAME>(partition[i_snp]);
            std::string bed_name = prev_name+".bed";
            openPlinkBinaryFile(bed_name, m_bed);
            prev=0;
        }
        size_t cur_line = std::get<+PRS::LINE>(partition[i_snp]);
        if((cur_line-prev)!=0)
        {
            // Skip snps
            m_bed.seekg((std::get<+PRS::LINE>(partition[i_snp])-prev)*m_num_bytes, m_bed.cur);
            prev=std::get<+PRS::LINE>(partition[i_snp]);
        }
        //read_snp(1, false);
        char genotype_list[m_num_bytes];
        m_bed.read(genotype_list, m_num_bytes);
        if (!m_bed) throw std::runtime_error("Problem with the BED file...has the FAM/BIM file been changed?");
        prev++;
        size_t sample_index = 0;
        int snp_index = std::get<+PRS::INDEX>(partition[i_snp]);
        if(snp_index >= snp_list.size()) throw std::runtime_error("Out of bound! In PRS score calculation");
        for(auto &&byte : genotype_list)
        {
        		size_t geno_bit = 0;
    			int geno_batch = static_cast<int>(byte);
        		while(geno_bit < 7 && sample_index < m_num_sample)
        		{
        			int geno = geno_batch>>geno_bit & 3; // This will access the corresponding genotype
        			if(geno!=1) // Because 01 is coded as missing
        			{
        				if(snp_list[snp_index].in(i_region))
        				{
        					prs_score[i_region][sample_index].second += snp_list[snp_index].score(geno);
        				}
        			}
        			sample_index++;
        			geno_bit+=2;
        		}
        }
    }
}

bool PLINK::openPlinkBinaryFile(const std::string s, std::ifstream & BIT)
{
    BIT.open(s.c_str(), std::ios::in | std::ios::binary);
    if(!BIT.is_open())
    {
        std::string error_message= "Cannot open the bed file: "+s;
        throw std::runtime_error(error_message);
    }
    // 1) Check for magic number
    // 2) else check for 0.99 SNP/Ind coding
    // 3) else print warning that file is too old
    char ch[3];
    BIT.read(ch,3);
    bool bfile_SNP_major = false;
    bool v1_bfile = true;
    // If v1.00 file format
    // Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
    // check magic number
    if(static_cast<uint32_t>(ch[0])== 108 && static_cast<uint32_t>(ch[1])==27)
        bfile_SNP_major=(static_cast<uint32_t>(ch[2])==1);
    else
        v1_bfile = false;
    // Reset file if < v1
    if ( ! v1_bfile )
    {
        std::cerr << "Warning, old BED file <v1.00 : will try to recover..." << std::endl;
        std::cerr << "  but you should --make-bed from PED )" << std::endl;
        BIT.close();
        BIT.clear();
        BIT.open(s.c_str(), std::ios::in | std::ios::binary);
        BIT.read(ch,1);
        uint32_t file_info = static_cast<uint32_t>(ch[0]);
        if(file_info != 1 && file_info!=0)
        {
            std::cerr << std::endl <<   " *** Possible problem: guessing that BED is < v0.99       *** " << std::endl;
            std::cerr <<                " *** High chance of data corruption, spurious results     *** " << std::endl;
            std::cerr <<                " *** Unless you are _sure_ this really is an old BED file *** " << std::endl;
            std::cerr <<                " *** you should recreate PED -> BED                       *** " << std::endl << std::endl;
            bfile_SNP_major = false;
            BIT.close();
            BIT.clear();
            BIT.open(s.c_str(), std::ios::in | std::ios::binary);
        }
        else
        {
            std::cerr << "Binary PED file is v0.99" << std::endl;
            bfile_SNP_major = (file_info==1);
        }
    }
    return bfile_SNP_major;
}

void PLINK::ld_dot_prod_batch(__m128i* vec1, __m128i* vec2, __m128i* mask1, __m128i* mask2, int32_t* return_vals, uint32_t iters)
{
    const __m128i m1 = {FIVEMASK,FIVEMASK};
    const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
    const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
    __m128i loader1;
    __m128i loader2;
    __m128i sum1;
    __m128i sum2;
    __m128i sum11;
    __m128i sum22;
    __m128i sum12;
    __m128i tmp_sum1;
    __m128i tmp_sum2;
    __m128i tmp_sum12;
    __univec acc;
    __univec acc1;
    __univec acc2;
    __univec acc11;
    __univec acc22;
    acc.vi = _mm_setzero_si128();
    acc1.vi = _mm_setzero_si128();
    acc2.vi = _mm_setzero_si128();
    acc11.vi = _mm_setzero_si128();
    acc22.vi = _mm_setzero_si128();
    do
    {
        loader1 = *vec1++;
        loader2 = *vec2++;
        sum1 = *mask2++;
        sum2 = *mask1++;
        sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
        sum1 = _mm_and_si128(sum1, loader1);
        sum2 = _mm_and_si128(sum2, loader2);
        sum11 = _mm_and_si128(sum1, m1);
        sum22 = _mm_and_si128(sum2, m1);
        loader1 = _mm_andnot_si128(_mm_add_epi64(m1, sum12), _mm_xor_si128(loader1, loader2));
        sum12 = _mm_or_si128(sum12, loader1);

        sum1 = _mm_add_epi64(_mm_and_si128(sum1, m2), _mm_and_si128(_mm_srli_epi64(sum1, 2), m2));
        sum2 = _mm_add_epi64(_mm_and_si128(sum2, m2), _mm_and_si128(_mm_srli_epi64(sum2, 2), m2));
        sum12 = _mm_add_epi64(_mm_and_si128(sum12, m2), _mm_and_si128(_mm_srli_epi64(sum12, 2), m2));

        loader1 = *vec1++;
        loader2 = *vec2++;
        tmp_sum1 = *mask2++;
        tmp_sum2 = *mask1++;

        tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
        tmp_sum1 = _mm_and_si128(tmp_sum1, loader1);
        tmp_sum2 = _mm_and_si128(tmp_sum2, loader2);
        sum11 = _mm_add_epi64(sum11, _mm_and_si128(tmp_sum1, m1));
        sum22 = _mm_add_epi64(sum22, _mm_and_si128(tmp_sum2, m1));
        loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
        tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);

        sum1 = _mm_add_epi64(sum1, _mm_add_epi64(_mm_and_si128(tmp_sum1, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum1, 2), m2)));
        sum2 = _mm_add_epi64(sum2, _mm_add_epi64(_mm_and_si128(tmp_sum2, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum2, 2), m2)));
        sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));

        loader1 = *vec1++;
        loader2 = *vec2++;

        tmp_sum1 = *mask2++;
        tmp_sum2 = *mask1++;

        tmp_sum12 = _mm_and_si128(_mm_or_si128(loader1, loader2), m1);
        tmp_sum1 = _mm_and_si128(tmp_sum1, loader1);
        tmp_sum2 = _mm_and_si128(tmp_sum2, loader2);
        sum11 = _mm_add_epi64(sum11, _mm_and_si128(tmp_sum1, m1));
        sum22 = _mm_add_epi64(sum22, _mm_and_si128(tmp_sum2, m1));
        loader1 = _mm_andnot_si128(_mm_add_epi64(m1, tmp_sum12), _mm_xor_si128(loader1, loader2));
        tmp_sum12 = _mm_or_si128(loader1, tmp_sum12);
//	    	if(bug) std::cerr << "Loader1: " << std::hex << loader1[0] << "\t" << loader1[1] << std::endl;
//	    	if(bug) std::cerr << "Loader2: " << std::hex << loader2[0] << "\t" << loader2[1] << std::endl;
//	    	if(bug) std::cerr << "Sum1b: " << std::hex << sum12[0] << "\t" << sum12[1] << std::endl;
//	    	if(bug) std::cerr << "Sum1t: " << std::hex << tmp_sum12[0] << "\t" << tmp_sum12[1] << std::endl;
        sum1 = _mm_add_epi64(sum1, _mm_add_epi64(_mm_and_si128(tmp_sum1, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum1, 2), m2)));
        sum2 = _mm_add_epi64(sum2, _mm_add_epi64(_mm_and_si128(tmp_sum2, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum2, 2), m2)));
        sum11 = _mm_add_epi64(_mm_and_si128(sum11, m2), _mm_and_si128(_mm_srli_epi64(sum11, 2), m2));
        sum22 = _mm_add_epi64(_mm_and_si128(sum22, m2), _mm_and_si128(_mm_srli_epi64(sum22, 2), m2));
        sum12 = _mm_add_epi64(sum12, _mm_add_epi64(_mm_and_si128(tmp_sum12, m2), _mm_and_si128(_mm_srli_epi64(tmp_sum12, 2), m2)));
//	    	if(bug) std::cerr << "Before: " << std::hex << acc.vi[0] << "\t" << acc.vi[1] << std::endl;
//	    	if(bug) std::cerr << "Sum12: " << std::hex << sum12[0] << "\t" << sum12[1] << std::endl;
        acc1.vi = _mm_add_epi64(acc1.vi, _mm_add_epi64(_mm_and_si128(sum1, m4), _mm_and_si128(_mm_srli_epi64(sum1, 4), m4)));
        acc2.vi = _mm_add_epi64(acc2.vi, _mm_add_epi64(_mm_and_si128(sum2, m4), _mm_and_si128(_mm_srli_epi64(sum2, 4), m4)));
        acc11.vi = _mm_add_epi64(acc11.vi, _mm_add_epi64(_mm_and_si128(sum11, m4), _mm_and_si128(_mm_srli_epi64(sum11, 4), m4)));
        acc22.vi = _mm_add_epi64(acc22.vi, _mm_add_epi64(_mm_and_si128(sum22, m4), _mm_and_si128(_mm_srli_epi64(sum22, 4), m4)));
        acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(sum12, m4), _mm_and_si128(_mm_srli_epi64(sum12, 4), m4)));
//	    	if(bug) std::cerr << "After: " << std::hex << acc.vi[0] << "\t" << acc.vi[1] << std::dec << std::endl;

    }
    while (--iters);
    // moved down because we've almost certainly run out of xmm registers
    const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
#if MULTIPLEX_LD > 960
    acc1.vi = _mm_add_epi64(_mm_and_si128(acc1.vi, m8), _mm_and_si128(_mm_srli_epi64(acc1.vi, 8), m8));
    acc2.vi = _mm_add_epi64(_mm_and_si128(acc2.vi, m8), _mm_and_si128(_mm_srli_epi64(acc2.vi, 8), m8));
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
#else
    acc1.vi = _mm_and_si128(_mm_add_epi64(acc1.vi, _mm_srli_epi64(acc1.vi, 8)), m8);
    acc2.vi = _mm_and_si128(_mm_add_epi64(acc2.vi, _mm_srli_epi64(acc2.vi, 8)), m8);
    acc.vi = _mm_and_si128(_mm_add_epi64(acc.vi, _mm_srli_epi64(acc.vi, 8)), m8);
#endif
    acc11.vi = _mm_and_si128(_mm_add_epi64(acc11.vi, _mm_srli_epi64(acc11.vi, 8)), m8);
    acc22.vi = _mm_and_si128(_mm_add_epi64(acc22.vi, _mm_srli_epi64(acc22.vi, 8)), m8);
    return_vals[0] -= ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
    return_vals[1] += ((acc1.u8[0] + acc1.u8[1]) * 0x1000100010001LLU) >> 48;
    return_vals[2] += ((acc2.u8[0] + acc2.u8[1]) * 0x1000100010001LLU) >> 48;
    return_vals[3] += ((acc11.u8[0] + acc11.u8[1]) * 0x1000100010001LLU) >> 48;
    return_vals[4] += ((acc22.u8[0] + acc22.u8[1]) * 0x1000100010001LLU) >> 48;
}

uint32_t PLINK::ld_missing_ct_intersect(long_type* lptr1, long_type* lptr2, uintptr_t word12_ct, uintptr_t word12_rem, uintptr_t lshift_last)
{
    // variant of popcount_longs_intersect()
    uintptr_t tot = 0;
    long_type* lptr1_end2;
#ifdef __LP64__
    const __m128i m1 = {FIVEMASK, FIVEMASK};
    const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
    const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
    const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
    __m128i* vptr1 = (__m128i*)lptr1;
    __m128i* vptr2 = (__m128i*)lptr2;
    __m128i* vend1;
    __m128i loader1;
    __m128i loader2;
    __univec acc;

    while (word12_ct >= 10)
    {
        word12_ct -= 10;
        vend1 = &(vptr1[60]);
ld_missing_ct_intersect_main_loop:
        acc.vi = _mm_setzero_si128();
        do
        {
            loader1 = _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1);
            loader2 = _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1);
            loader1 = _mm_add_epi64(loader1, _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
            loader2 = _mm_add_epi64(loader2, _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
            loader1 = _mm_add_epi64(loader1, _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
            loader2 = _mm_add_epi64(loader2, _mm_andnot_si128(_mm_or_si128(*vptr2++, *vptr1++), m1));
            loader1 = _mm_add_epi64(_mm_and_si128(loader1, m2), _mm_and_si128(_mm_srli_epi64(loader1, 2), m2));
            loader1 = _mm_add_epi64(loader1, _mm_add_epi64(_mm_and_si128(loader2, m2), _mm_and_si128(_mm_srli_epi64(loader2, 2), m2)));
            acc.vi = _mm_add_epi64(acc.vi, _mm_add_epi64(_mm_and_si128(loader1, m4), _mm_and_si128(_mm_srli_epi64(loader1, 4), m4)));
        }
        while (vptr1 < vend1);
        acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8), _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
        tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
    }
    if (word12_ct)
    {
        vend1 = &(vptr1[word12_ct * 6]);
        word12_ct = 0;
        goto ld_missing_ct_intersect_main_loop;
    }
    lptr1 = (long_type*)vptr1;
    lptr2 = (long_type*)vptr2;
#else
    uintptr_t* lptr1_end = &(lptr1[word12_ct * 12]);
    uintptr_t tmp_stor;
    uintptr_t loader1;
    uintptr_t loader2;
    while (lptr1 < lptr1_end)
    {
        loader1 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
        loader2 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
        loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
        loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
        loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
        loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
        loader1 = (loader1 & 0x33333333) + ((loader1 >> 2) & 0x33333333);
        loader1 += (loader2 & 0x33333333) + ((loader2 >> 2) & 0x33333333);
        tmp_stor = (loader1 & 0x0f0f0f0f) + ((loader1 >> 4) & 0x0f0f0f0f);

        loader1 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
        loader2 = (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
        loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
        loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
        loader1 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
        loader2 += (~((*lptr1++) | (*lptr2++))) & FIVEMASK;
        loader1 = (loader1 & 0x33333333) + ((loader1 >> 2) & 0x33333333);
        loader1 += (loader2 & 0x33333333) + ((loader2 >> 2) & 0x33333333);
        tmp_stor += (loader1 & 0x0f0f0f0f) + ((loader1 >> 4) & 0x0f0f0f0f);
        tot += (tmp_stor * 0x01010101) >> 24;
    }
#endif
    lptr1_end2 = &(lptr1[word12_rem]);
    while (lptr1 < lptr1_end2)
    {
        tot += popcount2_long((~((*lptr1++) | (*lptr2++))) & FIVEMASK);
    }
    if (lshift_last)
    {
        tot += popcount2_long(((~((*lptr1) | (*lptr2))) & FIVEMASK) << lshift_last);
    }
    return tot;
}

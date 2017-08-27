// This file is part of PRSice2.0, copyright (C) 2016-2017
// Shing Wan Choi, Jack Euesden, Cathryn M. Lewis, Paul F. O’Reilly
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "binarygen.hpp"

namespace bgenlib=genfile::bgen;

BinaryGen::BinaryGen(std::string prefix, std::string pheno_file,
        bool header, std::string remove_sample, std::string keep_sample,
        std::string extract_snp, std::string exclude_snp,
        bool ignore_fid, int num_auto, bool no_x, bool no_y, bool no_xy,
        bool no_mt, bool keep_ambig, const size_t thread, bool verbose)
{
	m_thread = thread;
	filter.keep_ambig = keep_ambig;
    if(!remove_sample.empty())
    {
    	m_remove_sample = true;
    	m_remove_sample_list = load_ref(remove_sample, ignore_fid);
    }
	else if(!keep_sample.empty())
    {
    	m_keep_sample = true;
    	m_keep_sample_list = load_ref(keep_sample, ignore_fid);
    }
    if(!extract_snp.empty())
    {
    	m_extract_snp = true;
    	m_extract_snp_list = load_snp_list(extract_snp);
    }
    else if(!exclude_snp.empty())
    {
    	m_exclude_snp = true;
    	m_exclude_snp_list = load_snp_list(exclude_snp);
    }

    /** setting the chromosome information **/
    m_xymt_codes.resize(XYMT_OFFSET_CT);
    // we are not using the following script for now as we only support human
    m_haploid_mask = new uintptr_t[CHROM_MASK_WORDS];
    fill_ulong_zero(CHROM_MASK_WORDS, m_haploid_mask);
    m_chrom_mask = new uintptr_t[CHROM_MASK_WORDS];
    fill_ulong_zero(CHROM_MASK_WORDS, m_chrom_mask);
    // now initialize the chromosome
    init_chr(num_auto, no_x, no_y, no_xy, no_mt);

    /** now get the chromosome information we've got by replacing the # in the name **/
    set_genotype_files(prefix);

    m_sample_names = preload_samples(pheno_file, header, ignore_fid);
    // now we update the sample information accordingly
    // unlike binaryPlink, we don't actually want the return value
    // mainly because this is a QC step, checking consistency between
    // the bgen file and the pheno file if the bgen also contain
    // the sample information
    load_samples(ignore_fid);

    /** Load SNP information from the file **/
    m_existed_snps = load_snps();
	m_marker_ct = m_existed_snps.size();

    if(verbose)
    {
    	if(m_num_non_founder!=0) fprintf(stderr, "%zu non-founder sample(s) removed\n", m_num_non_founder);
        fprintf(stderr, "%zu people (%zu males, %zu females) included\n", m_unfiltered_sample_ct, m_num_male, m_num_female);
        if(m_num_ambig!=0 && !keep_ambig) fprintf(stderr, "%u ambiguous variants excluded\n", m_num_ambig);
        else if(m_num_ambig!=0) fprintf(stderr, "%u ambiguous variants kept\n", m_num_ambig);
        fprintf(stderr, "%zu variants included\n", m_marker_ct);
    }

	uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    m_tmp_genotype = new uintptr_t[unfiltered_sample_ctl*2];
    m_remove_sample_list.clear();
    m_keep_sample_list.clear();
    m_extract_snp_list.clear();
    m_exclude_snp_list.clear();
}

BinaryGen::~BinaryGen()
{
    if(m_tmp_genotype!=nullptr) delete [] m_tmp_genotype;
    if(m_bgen_file.is_open()) m_bgen_file.close();
}

std::vector<SNP> BinaryGen::load_snps()
{
    std::vector<SNP> snp_res;
    bool chr_sex_error = false;
    bool chr_error = false;
    size_t chr_index=0;
    size_t expected_total=0;
    m_num_ambig=0;
    m_unfiltered_marker_ct=0;
    for(auto &&info : m_bgen_info)
    {
        expected_total+=info.second.number_of_variants;
    }
    snp_res.resize(expected_total);
    for(auto &&prefix : m_genotype_files)
    {
        std::string bgen_name = prefix+".bgen";
        if(m_bgen_file.is_open()) m_bgen_file.close();
        m_bgen_file.open(bgen_name.c_str(), std::ifstream::binary);
        if(!m_bgen_file.is_open())
        {
            std::string error_message = "ERROR: Cannot open bgen file "+bgen_name;
            throw std::runtime_error(error_message);
        }
        uint32_t offset = m_offset_map[prefix];
        m_bgen_file.seekg(offset+4);
        uint32_t num_snp = m_bgen_info[prefix].number_of_variants;
        uint32_t const layout = m_bgen_info[prefix].flags & genfile::bgen::e_Layout ;
        uint32_t num_sample = m_bgen_info[prefix].number_of_samples;
        std::string prev_chr="";
        int chr_code=0;

        for(size_t i_snp =0; i_snp < num_snp; ++i_snp)
        {
            if(m_unfiltered_marker_ct%1000==0 && m_unfiltered_marker_ct>0)
            {
                fprintf(stderr, "\r%zuK SNPs processed\r", m_unfiltered_marker_ct/1000);
            }
            uint16_t SNPID_size = 0;
            uint16_t RSID_size = 0;
            uint16_t numberOfAlleles = 0 ;
            uint16_t chromosome_size = 0 ;
            uint32_t allele_size = 0;
            std::string allele ;
            std::string SNPID;
            std::string RSID;
            std::string chromosome;
            uint32_t SNP_position;
            if( layout == genfile::bgen::e_Layout1 || layout == genfile::bgen::e_Layout0 ) {
                uint32_t number_of_samples ;
                bgenlib::read_little_endian_integer( m_bgen_file, &number_of_samples ) ;
                if( number_of_samples != num_sample ) {
                    throw std::runtime_error("ERROR: Number of sample doesn't match!");
                }
                bgenlib::read_length_followed_by_data( m_bgen_file, &SNPID_size, &SNPID ) ;
            } else if( layout == genfile::bgen::e_Layout2 ) {
                bgenlib::read_length_followed_by_data( m_bgen_file, &SNPID_size, &SNPID ) ;
            } else {
                assert(0) ;
            }
            bgenlib::read_length_followed_by_data( m_bgen_file, &RSID_size, &RSID ) ;
            bgenlib::read_length_followed_by_data( m_bgen_file, &chromosome_size, &chromosome ) ;
            bgenlib::read_little_endian_integer( m_bgen_file, &SNP_position ) ;
            if( layout == bgenlib::e_Layout2 ) {
                bgenlib::read_little_endian_integer( m_bgen_file, &numberOfAlleles ) ;
            } else {
                numberOfAlleles = 2 ;
            }
            if(numberOfAlleles!=2)
            {
                // after we finish writing up the papers, we can come back and work on this.
                // only use the 2 most common allele and set the remaining to missing
                throw std::runtime_error("ERROR: Currently only support bgen with 2 alleles!");
                // we can, in the future, allow for more than 2 alleles by setting anything
                // but the 2 most common alleles to missing
            }
            std::vector<std::string> final_alleles(numberOfAlleles);
            for( uint16_t i = 0; i < numberOfAlleles; ++i ) {
                bgenlib::read_length_followed_by_data( m_bgen_file, &allele_size, &allele ) ;
                std::transform(allele.begin(), allele.end(),allele.begin(), ::toupper);
                final_alleles[i] = allele;
            }

            if( !m_bgen_file ) {
#if DEBUG_BGEN_FORMAT
                std::cerr << "bgen: layout = " << layout << ", alleles = " << numberOfAlleles << ".\n" << std::flush ;
                std::cerr << *SNPID << ", " << *RSID << ", " << *chromosome << ", " << *SNP_position << ".\n" << std::flush ;
#endif
                throw std::runtime_error("ERROR: Problem reading bgen file!");
            }

            size_t snp_id = (unsigned int)m_bgen_file.tellg();
            std::vector< genfile::byte_t > buffer;
            read_genotype_data_block(m_bgen_file, m_bgen_info[prefix], &buffer); // move read pointer forward
            if(chromosome.compare(prev_chr)!=0)
            {
                prev_chr = chromosome;
                if(m_chr_order.find(chromosome)!= m_chr_order.end())
                {
                    throw std::runtime_error("ERROR: SNPs on the same chromosome must be clustered together!");
                }
                m_chr_order[chromosome] = chr_index++;
                chr_code = get_chrom_code_raw(chromosome.c_str());
                if (((const uint32_t)chr_code) > m_max_code) { // bigger than the maximum code, ignore it
                    if(!chr_error)
                    {
                        fprintf(stderr, "WARNING: SNPs with chromosome number larger than %du\n", m_max_code);
                        fprintf(stderr, "         They will be ignored!\n");
                        chr_error=true;
                        continue;
                    }
                    else if(!chr_sex_error && (is_set(m_haploid_mask, chr_code) ||
                            chr_code==m_xymt_codes[X_OFFSET] ||
                            chr_code==m_xymt_codes[Y_OFFSET]))
                    {
                        fprintf(stderr, "WARNING: Currently not support haploid chromosome and sex chromosomes\n");
                        chr_sex_error=true;
                        continue;
                    }
                }
            }
            if(RSID.compare(".")==0) // when the rs id isn't available, change it to chr:loc coding
            {
            	RSID = std::to_string(chr_code)+":"+std::to_string(SNP_position);
            }
            if((m_extract_snp && m_extract_snp_list.find(RSID) == m_extract_snp_list.end()) ||
            		(m_exclude_snp && m_exclude_snp_list.find(RSID) != m_exclude_snp_list.end())){
            	m_unfiltered_marker_ct++;
            	continue;
            }

            if(m_existed_snps_index.find(RSID)!= m_existed_snps_index.end())
            {
                throw std::runtime_error("ERROR: Duplicated SNP ID detected!\n");
            }
            else if(ambiguous(final_alleles.front(), final_alleles.back()))
            {
                m_num_ambig++;
                if(filter.keep_ambig)
                {
                	m_existed_snps_index[RSID] = m_unfiltered_marker_ct;
                	snp_res[m_unfiltered_marker_ct] = SNP(RSID, chr_code, SNP_position, final_alleles.front(),
                			final_alleles.back(), prefix, snp_id);
                	m_unfiltered_marker_ct++;
                }
            }
            else
            {
                m_existed_snps_index[RSID] = m_unfiltered_marker_ct;
                snp_res[m_unfiltered_marker_ct] = SNP(RSID, chr_code, SNP_position, final_alleles.front(),
                        final_alleles.back(), prefix, snp_id);
                m_unfiltered_marker_ct++;
                // directly ignore all others?
            }
        }
    }
    fprintf(stderr, "\n");
    snp_res.resize(m_unfiltered_marker_ct);// so that it will be more suitable
    if(m_bgen_file.is_open()) m_bgen_file.close();
    return snp_res;
}


void BinaryGen::dosage_score( misc::vec2d<Sample_lite> &current_prs_score,
                        size_t start_index, size_t end_bound)
{
    m_cur_file = "";
    std::vector<bool> in_region(m_region_size);
    std::vector< genfile::byte_t > buffer1, buffer2 ;
    size_t num_included_samples = current_prs_score.cols();
    for(size_t i_snp = start_index; i_snp < end_bound; ++i_snp)
    {
        for(size_t i_region = 0; i_region < m_region_size; ++i_region)
        {
            in_region[i_region]=m_existed_snps[i_snp].in(i_region);
        }
        auto &&snp =  m_existed_snps[i_snp];
        if(m_cur_file.empty() || snp.file_name().compare(m_cur_file)!=0)
        {
            if(m_bgen_file.is_open()) m_bgen_file.close();
            std::string bgen_name = snp.file_name()+".bgen";
            m_bgen_file.open(bgen_name.c_str(), std::ifstream::binary);
            if(!m_bgen_file.is_open())
            {
                std::string error_message = "ERROR: Cannot open bgen file: "+snp.file_name();
                throw std::runtime_error(error_message);
            }
            m_cur_file = snp.file_name();
        }
        m_bgen_file.seekg(snp.snp_id(), std::ios_base::beg);

        Data probability ;
        ProbSetter setter( &probability ) ;
        genfile::bgen::read_and_parse_genotype_data_block< ProbSetter >(
                m_bgen_file,
                m_bgen_info[snp.file_name()],
                setter,
                &buffer1,
                &buffer2
        ) ;
        std::vector<size_t> missing_samples;
        double total = 0.0;
        std::vector<double> score(num_included_samples);
        size_t cur_sample=0;
        for(size_t i_sample=0; i_sample < probability.size(); ++i_sample)
        {
            auto &&prob = probability[i_sample];
            if(prob.size() != 3)
            {
                // this is likely phased
                fprintf(stderr, "ERROR: Currently don't support phased data\n");
                fprintf(stderr, "       (It is because the lack of development time)\n");
                throw std::runtime_error("");
            }
            double expected = 0.0;
            if(IS_SET(m_sample_include,i_sample)) // to ignore non-selected SNPs
            {
                for(int g = 0; g < prob.size(); ++g)
                {
                    if(*max_element(prob.begin(), prob.end()) < filter.hard_threshold)
                    {
                        missing_samples.push_back(i_sample);
                        break;
                    }
                    else expected+=prob[g]*((snp.is_flipped())?abs(g-2) :g);
                }
                score[cur_sample++]=expected;
                total+=expected;
            }
        }
        // now process the missing and clean stuff
        // we divide the mean by 2 so that in situation where there is no dosage,
        // it will behave the same as the genotype data.

        size_t num_miss = missing_samples.size();
        double mean = total/(((double)num_included_samples-(double)num_miss)*2);
        size_t i_missing = 0;
        double stat = snp.stat();
        for(size_t i_sample=0; i_sample < num_included_samples; ++i_sample)
        {
            if(i_missing < num_miss && i_sample == missing_samples[i_missing])
            {
                for(size_t i_region=0; i_region < m_region_size; ++i_region)
                {
                    if(in_region[i_region])
                    {
                        if(m_scoring == SCORING::MEAN_IMPUTE) current_prs_score(i_region,i_sample).prs += stat*mean;
                        if(m_scoring != SCORING::SET_ZERO) current_prs_score(i_region,i_sample).num_snp++;
                    }
                }
                i_missing++;
            }
            else
            { // not missing sample
                for(size_t i_region=0; i_region < m_region_size; ++i_region)
                {
                    if(in_region[i_region])
                    {
                        if(m_scoring == SCORING::CENTER)
                        {
                            // if centering, we want to keep missing at 0
                            current_prs_score(i_region,i_sample).prs -= stat*mean;
                        }
                        // again, so that it will generate the same result as genotype file format
                        // when we are 100% certain of the genotypes
                        current_prs_score(i_region,i_sample).prs += score[i_sample]*stat*0.5;
                        current_prs_score(i_region,i_sample).num_snp++;
                    }
                }
            }
        }
    }
}

void BinaryGen::hard_code_score( misc::vec2d<Sample_lite> &current_prs_score,
                        size_t start_index, size_t end_bound)
{
    m_cur_file = "";
    uint32_t uii;
    uint32_t ujj;
    uint32_t ukk;
    uintptr_t ulii = 0;
	uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
	uintptr_t final_mask = get_final_mask(m_sample_ct);

    size_t num_included_samples = current_prs_score.cols();
    // a lot of code in PLINK is to handle the sex chromosome
    // which suggest that PRS can be done on sex chromosome
    // that should be something later

    std::vector<bool> in_region(m_region_size);
    // index is w.r.t. partition, which contain all the information
    uintptr_t* genotype = new uintptr_t[unfiltered_sample_ctl*2];
    for(size_t i_snp = start_index; i_snp < end_bound; ++i_snp)
    { // for each SNP
       for(size_t i_region = 0; i_region < m_region_size; ++i_region)
        {
            in_region[i_region]=m_existed_snps[i_snp].in(i_region);
        }
        std::fill(genotype, genotype+unfiltered_sample_ctl*2, 0);
        //std::memset(genotype, 0x0, m_unfiltered_sample_ctl*2*sizeof(uintptr_t));
        std::fill(m_tmp_genotype, m_tmp_genotype+unfiltered_sample_ctl*2,0);
        //std::memset(m_tmp_genotype, 0x0, m_unfiltered_sample_ctl*2*sizeof(uintptr_t));
        if(load_and_collapse_incl(m_existed_snps[i_snp].snp_id(), m_existed_snps[i_snp].file_name(),
                m_unfiltered_sample_ct, m_sample_ct, m_sample_include, final_mask,
                false, m_tmp_genotype, genotype))
        {
            throw std::runtime_error("ERROR: Cannot read the bed file!");
        }
        uintptr_t* lbptr = genotype;
        uii = 0;
        std::vector<size_t> missing_samples;
        double stat = m_existed_snps[i_snp].stat();
        bool flipped = m_existed_snps[i_snp].is_flipped();
        std::vector<double> genotypes(num_included_samples);
        int total_num = 0;
        uint32_t sample_idx=0;
        int nmiss = 0;
        // This whole thing was wrong...
        do
        {
            ulii = ~(*lbptr++);
            if (uii + BITCT2 > m_unfiltered_sample_ct)
            {
                ulii &= (ONELU << ((m_unfiltered_sample_ct & (BITCT2 - 1)) * 2)) - ONELU;
            }
            while (ulii)
            {
                ujj = CTZLU(ulii) & (BITCT - 2);
                ukk = (ulii >> ujj) & 3;
                sample_idx = uii + (ujj / 2);
                if(ukk==1 || ukk==3) // Because 01 is coded as missing
                {
                    // 3 is homo alternative
                    //int flipped_geno = snp_list[snp_index].geno(ukk);
                    if(sample_idx < num_included_samples){
                        total_num+=(ukk==3)? 2: ukk;
                        genotypes[sample_idx] = (ukk==3)? 2: ukk;
                    }
                }
                else // this should be 2
                {
                    missing_samples.push_back(sample_idx);
                    nmiss++;
                }
                ulii &= ~((3 * ONELU) << ujj);
            }
            uii += BITCT2;
        } while (uii < num_included_samples);


        size_t i_missing = 0;
        double maf = ((double)total_num/((double)((int)num_included_samples-nmiss)*2.0)); // MAF does not count missing
        if(flipped) maf = 1.0-maf;
        double center_score = stat*maf;
        size_t num_miss = missing_samples.size();
        for(size_t i_sample=0; i_sample < num_included_samples; ++i_sample)
        {
            if(i_missing < num_miss && i_sample == missing_samples[i_missing])
            {
                for(size_t i_region=0; i_region < m_region_size; ++i_region)
                {
                    if(in_region[i_region])
                    {
                        if(m_scoring == SCORING::MEAN_IMPUTE) current_prs_score(i_region,i_sample).prs += center_score;
                        if(m_scoring != SCORING::SET_ZERO) current_prs_score(i_region,i_sample).num_snp++;
                    }
                }
                i_missing++;
            }
            else
            { // not missing sample
                for(size_t i_region=0; i_region < m_region_size; ++i_region)
                {
                    if(in_region[i_region])
                    {
                        if(m_scoring == SCORING::CENTER)
                        {
                            // if centering, we want to keep missing at 0
                            current_prs_score(i_region,i_sample).prs -= center_score;
                        }
                        current_prs_score(i_region,i_sample).prs += ((flipped)?fabs(genotypes[i_sample]-2):genotypes[i_sample])*stat*0.5;
                        current_prs_score(i_region,i_sample).num_snp++;
                    }
                }
            }
        }
    }
    delete [] genotype;
}

void BinaryGen::read_score( misc::vec2d<Sample_lite> &current_prs_score,
               size_t start_index, size_t end_bound)
{
    if(current_prs_score.rows() != m_region_size)
    {
        throw std::runtime_error("Size of Matrix doesn't match number of region!!");
    }

    if(filter.use_hard)
    {
        hard_code_score(current_prs_score, start_index, end_bound);
        return;
    }
    else dosage_score(current_prs_score, start_index, end_bound);
}


Sample BinaryGen::get_sample(std::vector<std::string> &token, bool ignore_fid,
        bool has_sex, int sex_col, std::vector<int> &sex_info)
{
    std::string id = (ignore_fid)? token[0] : token[0]+"_"+token[1];
    // this will pose problem when there are duplicated IID names even if they
    // are from different family. However, we don't know how bgen store the
    // sample information (do they contain the FID?) so we will have to work
    // this way
    if(m_sample_index_check.find((ignore_fid)? token[0] : token[1])==m_sample_index_check.end())
    {
        Sample cur_sample;
        cur_sample.FID=(ignore_fid)? "" : token[0];
        cur_sample.IID=(ignore_fid)? token[0] : token[1];
        cur_sample.pheno = "NA";
        cur_sample.included = false;
        if(m_keep_sample)
        {
            cur_sample.included=(m_keep_sample_list.find(id)!=m_keep_sample_list.end());
        }
        else if(m_remove_sample)
        {
            cur_sample.included = !(m_remove_sample_list.find(id)!=m_remove_sample_list.end());
        }
        else cur_sample.included = true;
        cur_sample.prs = 0;
        cur_sample.num_snp = 0;
        if(has_sex)
        {
            try{
                sex_info.push_back(misc::convert<int>(token[sex_col]));
            }
            catch(const std::runtime_error &er)
            {
                throw std::runtime_error("ERROR: Invalid sex coding!\n");
            }
        }
        return cur_sample;
    }
    else
    {
        std::string error_message = "ERROR: Duplicated sample: "+ id;
        throw std::runtime_error(error_message);
    }
}

std::vector<Sample> BinaryGen::preload_samples(std::string pheno, bool has_header, bool ignore_fid)
{
    std::vector<Sample> sample_res;
    std::ifstream pheno_file;
    pheno_file.open(pheno.c_str());
    if(!pheno_file.is_open())
    {
        std::string error_message = "ERROR: Cannot open phenotype file: "+pheno;
        throw std::runtime_error(error_message);
    }
    std::string line;
    // we will assume we have all
    std::string first_line;
    std::getline(pheno_file, first_line);
    misc::trim(first_line);
    std::vector<std::string> possible_header = misc::split(first_line);
    std::string second_line;
    std::getline(pheno_file, second_line);
    misc::trim(second_line);
    std::vector<std::string> token = misc::split(second_line);
    bool bgen_sample = false;
    if(token.size() > 3) // FID IID and Missing are required, then phenotype
    {
        if(token[0].compare("0")==0 && token[1].compare("0")==0 && token[2].compare("0")==0)
        {
            bgen_sample = true;
            for(size_t i = 3; i < token.size(); ++i)
            {
                if(token[i].compare("D")!= 0&& token[i].compare("C")!= 0&&
                        token[i].compare("P")!= 0 && token[i].compare("B")!= 0)
                {
                    bgen_sample = false;
                    break;
                }
            }
        }
    }
    bool has_sex=false;
    int sex_col = 0;
    std::vector<int> sex_info;
    if(bgen_sample) // then we know the first line is header
    {
        fprintf(stderr, "Detected bgen sample format\n");
        for(size_t i = 3; i < possible_header.size(); ++i)
        {
            if(possible_header[i].compare("Sex")==0 || possible_header[i].compare("sex")==0
                    || possible_header[i].compare("SEX")==0)
            {
                has_sex = true;
                sex_col = i;
                break;
            }
        }
        if(has_sex)
        {
            if(token[sex_col].compare("D")!=0)
            {
                std::string error_message= "ERROR: Sex must be coded as \"D\" in bgen sample file!";
                throw std::runtime_error(error_message);
            }
        }
    }
    else
    {
        // this is just a normal line
        if(!has_header)
        {
            // need to check that for normal pheno file, whether there is
            // a header
            if(possible_header.front().compare("FID")==0 ||
                    (possible_header.size()> 2 && possible_header[1].compare("IID")==0) ||
                    possible_header.front().compare("IID")==0)
            {
                // these are simple test. might not be correct though
            	// if fit, this is header
            }
            else
            {
            	// this isn't a header
                sample_res.push_back(get_sample(possible_header, ignore_fid, has_sex, sex_col, sex_info));
                m_sample_index_check[sample_res.back().IID] = sample_res.size()-1;
            }
        }
        sample_res.push_back(get_sample(token, ignore_fid, has_sex, sex_col, sex_info));
        m_sample_index_check[sample_res.back().IID] = sample_res.size()-1;
    }

    while(std::getline(pheno_file, line))
    {
        misc::trim(line);
        if(line.empty()) continue;
        std::vector<std::string> token = misc::split(line);
        if(token.size() < (has_sex)? (sex_col):(1+!ignore_fid))
        {
            std::string error_message = "ERROR: Header line must contain at least "+
                    std::to_string((has_sex)? (sex_col):(1+!ignore_fid))+" columns!";
        }
        sample_res.push_back(get_sample(token, ignore_fid, has_sex, sex_col, sex_info));
        m_sample_index_check[sample_res.back().IID] = sample_res.size()-1;
    }
    m_unfiltered_sample_ct = sample_res.size();
    uintptr_t unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    //uintptr_t unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;

    //m_sex_male = new uintptr_t[m_unfiltered_sample_ctl];
    //std::memset(m_sex_male, 0x0, m_unfiltered_sample_ctl*sizeof(uintptr_t));

    // assume all are founder
    m_founder_info = new uintptr_t[unfiltered_sample_ctl];
    m_founder_ct = m_unfiltered_sample_ct;
    std::fill(m_founder_info, m_founder_info+unfiltered_sample_ctl,0);
    //std::memset(m_founder_info, 0, m_unfiltered_sample_ctl*sizeof(uintptr_t));
    m_sample_include = new uintptr_t[unfiltered_sample_ctl];
    std::fill(m_sample_include, m_sample_include+unfiltered_sample_ctl,0);
    //std::memset(m_sample_include, 0x0, m_unfiltered_sample_ctl*sizeof(uintptr_t));
	for(size_t i = 0; i < sample_res.size(); ++i)
    {
        if(sample_res[i].included)
        {
            SET_BIT(i, m_founder_info);
        }
    }

    m_num_male = 0, m_num_female = 0, m_num_ambig_sex=(has_sex)? 0 : m_unfiltered_sample_ct;
    for(size_t i = 0; i < sex_info.size() && has_sex; ++i)
    {
        switch(sex_info[i])
        {
            case 1:
                m_num_male++;
                //SET_BIT(i, m_sex_male);
                break;
            case 2:
                m_num_female++;
                break;
            default:
                m_num_ambig_sex++;
        }
    }
    pheno_file.close();
    return sample_res;
}

std::vector<Sample> BinaryGen::load_samples(bool ignore_fid)
{
    std::unordered_set<std::string> dup_check;
    bool first =  true;
    for(auto &&prefix : m_genotype_files)
    {
        if(m_bgen_file.is_open()) m_bgen_file.close();
        std::string bgen_name = prefix+".bgen";
        m_bgen_file.open(bgen_name.c_str(), std::ifstream::binary);
        if(!m_bgen_file.is_open())
        {
            std::string error_message = "ERROR: Cannot open bgen file: "+bgen_name;
            throw std::runtime_error(error_message);
        }
        uint32_t offset = 0;
        bgenlib::read_little_endian_integer(m_bgen_file, &offset);
        uint32_t header_size = 0, number_of_snp_blocks = 0, number_of_samples = 0,
                flags = 0 ;
        char magic[4] ;
        std::size_t fixed_data_size = 20 ;
        std::vector<char> free_data ;
        bgenlib::read_little_endian_integer( m_bgen_file, &header_size ) ;
        assert( header_size >= fixed_data_size ) ;
        bgenlib::read_little_endian_integer( m_bgen_file, &number_of_snp_blocks ) ;
        bgenlib::read_little_endian_integer( m_bgen_file, &number_of_samples ) ;
        m_bgen_file.read( &magic[0], 4 ) ;
        free_data.resize( header_size - fixed_data_size ) ;
        m_bgen_file.read( &free_data[0], free_data.size() ) ;
        bgenlib::read_little_endian_integer( m_bgen_file, &flags ) ;
        if(( magic[0] != 'b' || magic[1] != 'g' || magic[2] != 'e' || magic[3] != 'n' )
                && ( magic[0] != 0 || magic[1] != 0 || magic[2] != 0 || magic[3] != 0 ))
        {
            throw std::runtime_error("ERROR: Incorrect magic string!\nPlease check you have provided a valid bgen file!");
        }
        if( m_bgen_file )
        {
            bgenlib::Context current_context;
            current_context.number_of_samples = number_of_samples ;
            current_context.number_of_variants = number_of_snp_blocks ;
            current_context.magic.assign( &magic[0], &magic[0] + 4 ) ;
            //current_context.free_data.assign( free_data.begin(), free_data.end() ) ;
            current_context.flags = flags ;
            m_bgen_info[prefix] = current_context;
            m_offset_map[prefix] = offset;
        }
        else
        {
            throw std::runtime_error("ERROR: Problem reading bgen file!") ;
        }
        if(m_bgen_info[prefix].number_of_samples!= m_sample_names.size())
        {
            std::string error_message = "ERROR: Number of sample in bgen does not match those in phenotype file! ("
                    +std::to_string(m_bgen_info[prefix] .number_of_samples)+" vs "+
                    std::to_string(m_sample_names.size())+")";
            throw std::runtime_error(error_message);
        }

        uint32_t const compressionType = m_bgen_info[prefix].flags & bgenlib::e_CompressedSNPBlocks;
        if( compressionType == bgenlib::e_ZstdCompression )
        {
            throw std::runtime_error("ERROR: zstd compression currently not supported");
        }
        if((m_bgen_info[prefix] .flags & bgenlib::e_SampleIdentifiers) && first)
        { // only read in the sample information for the first bgen
            first = false;
            uint32_t sample_block_size = 0 ;
            uint32_t actual_number_of_samples = 0 ;
            uint16_t identifier_size ;
            std::string identifier ;
            std::size_t bytes_read = 0 ;
            // the bgen format actually double stored the number of samples and
            // the block size is the sample_block size
            bgenlib::read_little_endian_integer(m_bgen_file, &sample_block_size ) ;
            bgenlib::read_little_endian_integer( m_bgen_file, &actual_number_of_samples ) ;
            bytes_read += 8 ;
            assert( actual_number_of_samples == m_bgen_info[prefix] .number_of_samples ) ;
            for( uint32_t i = 0; i < actual_number_of_samples; ++i ) {
                bgenlib::read_length_followed_by_data( m_bgen_file, &identifier_size, &identifier ) ;
                if( m_bgen_file ) {
                    bytes_read += sizeof( identifier_size ) + identifier_size ;
                    if(m_sample_index_check.find(identifier)== m_sample_index_check.end())
                    {
                        throw std::runtime_error("ERROR: Sample mismatch between bgen and phenotype file!");
                    }
                    else if(m_sample_index_check[identifier] !=i)
                    {
                        throw std::runtime_error("ERROR: Sample sequence differ between bgen and phenotype file!");
                    }
                    if(dup_check.find(identifier)==dup_check.end())
                    {
                        dup_check.insert(identifier);
                        // only for checking
                        if(m_sample_names[i].IID.compare(identifier)!=0)
                        {
                            throw std::runtime_error("ERROR: Sample name mismatch between bgen and phenotype file!");
                        }
                    }
                    else
                    {
                        throw std::runtime_error("ERROR: Duplicated sample within bgen file!");
                    }
                }
                else
                {
                    throw std::runtime_error("ERROR: Problem reading bgen file!") ;
                }
            }
            assert( bytes_read == sample_block_size ) ;
        }
    }
    return std::vector<Sample>();
}
/*
void BinaryGen::get_header()
{

}
*/

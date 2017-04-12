#include "genotype.hpp"


BinaryPlink::BinaryPlink(std::string prefix, int num_auto, bool no_x, bool no_y, bool no_xy, bool no_mt,
		const size_t thread, bool verbose)
{
	m_xymt_codes.resize(XYMT_OFFSET_CT);
	init_chr(num_auto, no_x, no_y, no_xy, no_mt);
	m_thread = thread;
	set_genotype_files(prefix);
	m_sample_names = load_samples();
	m_existed_snps = load_snps();
	if(verbose)
	{
		fprintf(stderr, "%zu people (%zu males, %zu females) included\n", m_unfiltered_sample_ct, m_num_male, m_num_female);
		if(m_num_ambig!=0) fprintf(stderr, "%u ambiguous variants excluded\n", m_num_ambig);
		fprintf(stderr, "%zu variants included\n", m_marker_ct);
	}
	//Genotype(prefix,num_auto, no_x, no_y, no_xy, no_mt, thread, verbose)
	check_bed();
	m_cur_file="";
}

BinaryPlink::~BinaryPlink()
{
	if(m_tmp_genotype != nullptr) delete [] m_tmp_genotype;
}
std::vector<Sample> BinaryPlink::load_samples()
{
	assert(m_genotype_files.size()>0);
	std::string famName = m_genotype_files.front()+".fam";
	std::ifstream famfile;
	famfile.open(famName.c_str());
	if(!famfile.is_open())
	{
		std::string error_message = "ERROR: Cannot open fam file: "+famName;
		throw std::runtime_error(error_message);
	}
	std::string line;
	m_unfiltered_sample_ct = 0;
	uintptr_t sample_uidx=0;
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
	m_unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
	m_unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;
	std::vector<Sample> sample_name(m_unfiltered_sample_ct);

	m_sex_male = new uintptr_t[m_unfiltered_sample_ctl];
	std::memset(m_sex_male, 0x0, m_unfiltered_sample_ctl*sizeof(uintptr_t));

	m_founder_info = new uintptr_t[m_unfiltered_sample_ctl];
	std::memset(m_founder_info, 0x0, m_unfiltered_sample_ctl*sizeof(uintptr_t));

	m_sample_exclude = new uintptr_t[m_unfiltered_sample_ctl];
	std::memset(m_sample_exclude, 0x0, m_unfiltered_sample_ctl*sizeof(uintptr_t));

	m_num_male = 0, m_num_female = 0, m_num_ambig_sex=0;
	std::vector<Sample> sample_names;
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
			Sample cur_sample;
			cur_sample.FID = token[+FAM::FID];
			cur_sample.IID = token[+FAM::IID];
			cur_sample.pheno = token[+FAM::PHENOTYPE];
			cur_sample.included = true;
			cur_sample.prs = 0;
			cur_sample.num_snp = 0;
			sample_name.push_back(cur_sample);
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
				 m_num_ambig_sex++; //currently ignore as we don't do sex chromosome
				 //SET_BIT(sample_uidx, m_sample_exclude); // exclude any samples without sex information
			}
			sample_uidx++;
		}
	}
	famfile.close();
	m_final_mask = get_final_mask(m_founder_ct);
	m_tmp_genotype = new uintptr_t[m_unfiltered_sample_ctl*2];
	return sample_name;
}

std::vector<SNP> BinaryPlink::load_snps()
{
	assert(m_genotype_files.size()>0);
	m_unfiltered_marker_ct = 0;
	m_unfiltered_marker_ctl=0;
	std::ifstream bimfile;
	std::string prev_chr = "";
	int chr_code=0;
	int order = 0;
	bool chr_error = false, chr_sex_error = false;
	m_num_ambig = 0;
	std::vector<SNP> snp_info;
	for(auto &&prefix : m_genotype_files)
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
			if(line.empty()) continue;
			std::vector<std::string> token = misc::split(line);
			if(token.size() < 6)
			{
				fprintf(stderr, "Error: Malformed bim file. Less than 6 column on line: %i\n",num_line);
				throw std::runtime_error("");
			}
			std::string chr = token[+BIM::CHR];
			if(chr.compare(prev_chr)!=0)
			{
				prev_chr = chr;
				if(m_chr_order.find(chr)!= m_chr_order.end())
				{
					throw std::runtime_error("ERROR: SNPs on the same chromosome must be clustered together!");
				}
				m_chr_order[chr] = order++;
				chr_code = get_chrom_code_raw(chr.c_str());
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
			int loc = misc::convert<int>(token[+BIM::BP]);
			if(loc < 0)
			{
				fprintf(stderr, "ERROR: SNP with negative corrdinate: %s:%s\n", token[+BIM::RS].c_str(), token[+BIM::BP].c_str());
				throw std::runtime_error("Please check you have the correct input");
			}
			if(m_existed_snps_index.find(token[+BIM::RS])!= m_existed_snps_index.end())
			{
				throw std::runtime_error("ERROR: Duplicated SNP ID detected!\n");
			}
			else if(ambiguous(token[+BIM::A1], token[+BIM::A2]))
			{
				m_num_ambig++;
			}
			else
			{
				m_existed_snps_index[token[+BIM::RS]] = snp_info.size();
				snp_info.push_back(SNP(token[+BIM::RS], chr_code, loc, token[+BIM::A1],
					token[+BIM::A2], prefix, num_line));
				// directly ignore all others?
			}
			m_unfiltered_marker_ct++; //add in the checking later on
			num_line++;
		}
		bimfile.close();
	}

	m_unfiltered_marker_ctl = BITCT_TO_WORDCT(m_unfiltered_marker_ct);
	m_marker_exclude = new uintptr_t[m_unfiltered_marker_ctl];
	std::memset(m_marker_exclude, 0x0, m_unfiltered_marker_ctl*sizeof(uintptr_t));
	if(m_unfiltered_marker_ct > 2147483645)
	{
		throw std::runtime_error("Error: PLINK does not suport more than 2^31 -3 variants. "
			"As we are using PLINK for some of our functions, we might encounter problem too. "
			"Sorry.");
	}
	m_marker_ct = m_unfiltered_marker_ct - m_marker_exclude_ct;
	return snp_info;
}


void BinaryPlink::check_bed()
{
	uint32_t uii = 0;
	int64_t llxx = 0;
	int64_t llyy = 0;
	int64_t llzz = 0;
	for(auto &&prefix : m_genotype_files)
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
}

void BinaryPlink::read_genotype(uintptr_t* genotype, const uint32_t snp_index, const std::string &file_name)
{
	if(m_cur_file.empty() || m_cur_file.compare(file_name)!=0)
	{
		if (m_bedfile != nullptr)
		{
			fclose(m_bedfile);
			m_bedfile=nullptr;
		}
		std::string bedname = file_name+".bed";
		m_bedfile = fopen(bedname.c_str(), FOPEN_RB); // assume there is no error
	}
	if (fseeko(m_bedfile, m_bed_offset + (snp_index* ((uint64_t)m_unfiltered_sample_ct4))
			, SEEK_SET))
	{
		throw std::runtime_error("ERROR: Cannot read the bed file!");
	}
	genotype = new uintptr_t[m_unfiltered_sample_ctl*2];
	std::memset(genotype, 0x0, m_unfiltered_sample_ctl*2*sizeof(uintptr_t));
	std::memset(m_tmp_genotype, 0x0, m_unfiltered_sample_ctl*2*sizeof(uintptr_t));
	if(load_and_collapse_incl(m_unfiltered_sample_ct, m_founder_ct, m_founder_info, m_final_mask,
			false, m_bedfile, m_tmp_genotype, genotype))
	{
		throw std::runtime_error("ERROR: Cannot read the bed file!");
	}
}

void BinaryPlink::read_score(std::vector< std::vector<Sample_lite> > &current_prs_score, size_t start_index,
		size_t end_bound)
{
    size_t prev =0;
    uint32_t uii;
    uint32_t ujj;
    uint32_t ukk;
    uintptr_t ulii = 0;
    int32_t delta1 = 1;
	int32_t delta2 = 2;
	int32_t deltam = 0;
	//m_final_mask

	// a lot of code in PLINK is to handle the sex chromosome
	// which suggest that PRS can be done on sex chromosome
	// that should be something later
	m_cur_file =""; // just close it
	if(m_bedfile!=nullptr) fclose(m_bedfile);

    uint32_t max_reverse = BITCT_TO_WORDCT(m_unfiltered_marker_ct);
    if(current_prs_score.size() != m_region_size)
    {
    		throw std::runtime_error("Size of Matrix doesn't match number of region!!");
    }
    for(size_t i = 0; i < m_region_size; ++i)
    {
    	 if(current_prs_score[i].size()!=m_unfiltered_sample_ct)
    	 {
    		 throw std::runtime_error("Size of Matrix doesn't match number of samples!!");
    	 }
    }
	std::vector<bool> in_region(m_region_size);
	// index is w.r.t. partition, which contain all the information
	uintptr_t* genotype = new uintptr_t[m_unfiltered_sample_ctl*2];
    for(size_t i_snp = start_index; i_snp < end_bound; ++i_snp)
    { // for each SNP
        if(m_cur_file.empty() || m_cur_file.compare(m_existed_snps[i_snp].file_name())!=0)
        {
        	if(m_bedfile!=nullptr) fclose(m_bedfile);
        	m_cur_file= m_existed_snps[i_snp].file_name();
            std::string bedname = m_cur_file+".bed";
            m_bedfile = fopen(bedname.c_str(), FOPEN_RB); //again, we are assuming that the file is correct
            prev=0;
        }

        for(size_t i_region = 0; i_region < m_region_size; ++i_region)
        {
        	in_region[i_region]=m_existed_snps[i_snp].in(i_region);
        }
        size_t cur_line = m_existed_snps[i_snp].snp_id();
        if (fseeko(m_bedfile, m_bed_offset + (cur_line* ((uint64_t)m_unfiltered_sample_ct4))
        		, SEEK_SET)) {
        	throw std::runtime_error("ERROR: Cannot read the bed file!");
        }
        //loadbuf_raw is the temporary
        //loadbuff is where the genotype will be located

        std::memset(genotype, 0x0, m_unfiltered_sample_ctl*2*sizeof(uintptr_t));
        std::memset(m_tmp_genotype, 0x0, m_unfiltered_sample_ctl*2*sizeof(uintptr_t));
        if(load_and_collapse_incl(m_unfiltered_sample_ct, m_founder_ct, m_sample_exclude, m_final_mask,
        		false, m_bedfile, m_tmp_genotype, genotype))
        {
        	throw std::runtime_error("ERROR: Cannot read the bed file!");
        }
        uintptr_t* lbptr = genotype;
        uii = 0;
        std::vector<size_t> missing_samples;
        double stat = m_existed_snps[i_snp].stat();
        bool flipped = m_existed_snps[i_snp].is_flipped();
        std::vector<double> genotypes(m_unfiltered_sample_ct);
        int total_num = 0;
        uint32_t sample_idx=0;
        int nmiss = 0;
        // This whole thing is wrong...
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
        			total_num+=(ukk==3)? 2: ukk;
        			genotypes[sample_idx] = (ukk==3)? 2: ukk;
        		}
        		else // this should be 2
        		{
        			missing_samples.push_back(sample_idx);
        			nmiss++;
        		}
        		ulii &= ~((3 * ONELU) << ujj);
        	}
        	uii += BITCT2;
        } while (uii < m_unfiltered_sample_ct);


        size_t i_missing = 0;
        double maf = ((double)total_num/((double)((int)m_unfiltered_sample_ct-nmiss)*2.0)); // MAF does not count missing
        if(flipped) maf = 1.0-maf;
        double center_score = stat*maf;
        size_t num_miss = missing_samples.size();
        for(size_t i_sample=0; i_sample < m_unfiltered_sample_ct; ++i_sample)
        {
        	if(i_missing < num_miss && i_sample == missing_samples[i_missing])
        	{
        		for(size_t i_region=0; i_region < m_region_size; ++i_region)
        		{
        			if(in_region[i_region])
        			{
        				if(m_scoring == SCORING::MEAN_IMPUTE) current_prs_score[i_region][i_sample].prs += center_score;
        				if(m_scoring != SCORING::SET_ZERO) current_prs_score[i_region][i_sample].num_snp++;
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
        					current_prs_score[i_region][i_sample].prs -= center_score;
        				}
        				current_prs_score[i_region][i_sample].prs += ((flipped)?fabs(genotypes[i_sample]-2):genotypes[i_sample])*stat*0.5;
        				current_prs_score[i_region][i_sample].num_snp++;
        			}
        		}
        	}
        }
    }
    delete [] genotype;
}


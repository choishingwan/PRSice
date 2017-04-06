#include "genotype.hpp"


BinaryPlink::BinaryPlink(std::string prefix, int num_auto, bool no_x, bool no_y, bool no_xy, bool no_mt,
		const size_t thread, bool verbose):
		Genotype(prefix,num_auto, no_x, no_y, no_xy, no_mt, thread, verbose)
{
	load_bed();
}

void BinaryPlink::load_sample()
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

	m_sex_male = new uintptr_t[m_unfiltered_sample_ctl];
	std::memset(m_sex_male, 0x0, m_unfiltered_sample_ctl*sizeof(uintptr_t));

	m_founder_info = new uintptr_t[m_unfiltered_sample_ctl];
	std::memset(m_founder_info, 0x0, m_unfiltered_sample_ctl*sizeof(uintptr_t));

	m_sample_exclude = new uintptr_t[m_unfiltered_sample_ctl];
	std::memset(m_sample_exclude, 0x0, m_unfiltered_sample_ctl*sizeof(uintptr_t));

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
				 SET_BIT(sample_uidx, m_sample_exclude); // exclude any samples without sex information
			}
			sample_uidx++;
		}
	}
	famfile.close();
}

std::vector<SNP> BinaryPlink::load_snps()
{
	assert(m_genotype_files.size()>0);
	m_unfiltered_marker_ct = 0;
	m_unfiltered_marker_ctl=0;
	std::ifstream bimfile;
	std::vector<SNP> snp_info;
	std::string prev_chr = "";
	int chr_code=0;
	int order = 0;
	bool chr_error = false, chr_sex_error = false;
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
					else if(!chr_sex_error)
					{
						fprintf(stderr, "WARNING: Sex chromosome currently not supported\n");
						chr_sex_error;
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
			// better way is to use struct. But for some reason that doesn't work
			if(m_existed_snps_index.find(token[+BIM::RS])!= m_existed_snps_index.end())
			{
				throw std::runtime_error("ERROR: Duplicated SNP ID detected!\n");
			}
			m_existed_snps_index[token[+BIM::RS]] = m_unfiltered_marker_ct;
			snp_info.push_back(SNP(token[+BIM::RS], chr_code, loc, token[+BIM::A1],
					token[+BIM::A2]));
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

void BinaryPlink::load_bed()
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

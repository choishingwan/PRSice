#include "genotype.hpp"


BinaryPlink::BinaryPlink(std::string prefix, int num_auto, bool x, bool y, bool xy, bool mt,
		const size_t thread, bool verbose)
{
	Genotype(prefix,num_auto, x, y, xy, mt, thread, verbose);
}

void BinaryPlink::load_samples()
{
	assert(m_genotype_files.size()>0);
	std::string famName = m_genotype_files.front()+".fam";
	std::ifstream famfile;
	famfile.open(famName.c_str());
	if(!fam.is_open())
	{
		std::string error_message = "ERROR: Cannot open fam file: "+famName;
		throw std::runtime_error(error_message);
	}
	std::string line;
	m_unfiltered_sample_ct = 0;
	uintptr_t sample_uidx;

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
	m_unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
	m_unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;

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
				 SET_BIT(sample_uidx, m_sample_exclude); // exclude any samples without sex information
			}
			sample_uidx++;
		}
	}
	famfile.close();
	return 0;
}

int32_t BinaryPlink::load_snps()
{
	assert(m_genotype_files.size()>0);
	m_unfiltered_marker_ct = 0;
	m_unfiltered_marker_ctl=0
	std::ifstream bimfile;

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
			if(!line.empty())
			{
				std::vector<std::string> token = misc::split(line);
				if(token.size() < 6)
				{
					fprintf(stderr, "Error: Malformed bim file. Less than 6 column on line: %i\n",num_line);
					throw std::runtime_error("");
				}
			}
			m_unfiltered_marker_ct++; //add in the checking later on
			num_line++;
		}
		bimfile.close();
	}
	m_unfiltered_marker_ctl = BITCT_TO_WORDCT(m_unfiltered_marker_ct);
	m_marker_exclude = new uintptr_t[unfiltered_marker_ctl];
	std::memset(m_marker_exclude, 0x0, unfiltered_marker_ctl*sizeof(uintptr_t));

	if(m_unfiltered_marker_ct > 2147483645)
	{
		throw std::runtime_error("Error: PLINK does not suport more than 2^31 -3 variants. "
			"As we are using PLINK for some of our functions, we might encounter problem too. "
			"Sorry.");
	}
	size_t marker_uidx=0;
	bool chr_error=false;
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
		std::string prev_chr="";
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
				std::string chr = token[+BIM::CHR];
				if(chr.compare(prev_chr)!=0)
				{
					chr_code = get_chrom_code_raw(chr.c_str());
					if (((const uint32_t)chr_code) > m_max_code) { // bigger than the maximum code, ignore it
						if(!chr_error)
						{
							fprintf(stderr, "WARNING: SNPs with chromosome number larger than %zu\n", m_max_code);
							fprintf(stderr, "         They will be ignored!\n");
							chr_error=true;
						}
					}
					else
					{
						int temp = misc::convert<int>(token[+BIM::BP]);
						if(temp < 0)
						{
							fprintf(stderr, "ERROR: SNP with negative corrdinate: %s:%d\n", token[+BIM::RS]. token[+BIM::BP]);
							throw std::runtime_error("Please check you have the correct input");
						}
						size_t loc = temp;
						m_existed_snps[token[+BIM::RS]] = SNP(token[+BIM::CHR], token[+BIM::A1],
								token[+BIM::A2], loc, prefix, num_line, true);
					}
				}
			}
			marker_uidx++;
			num_line++;
		}
		bimfile.close();
	}
	m_marker_ct = m_unfiltered_marker_ct - m_marker_exclude_ct;
	return 0;
}

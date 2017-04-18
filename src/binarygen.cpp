#include "binarygen.hpp"

BinaryGen::BinaryGen(std::string prefix, std::string pheno_file,
        bool header, std::string remove_sample, std::string keep_sample,
        bool ignore_fid, int num_auto, bool no_x, bool no_y, bool no_xy,
        bool no_mt, const size_t thread, bool verbose)
{
    if(remove_sample.empty()) m_remove_sample = false;
    else m_remove_sample_list = load_ref(remove_sample, ignore_fid);
    if(keep_sample.empty()) m_keep_sample = false;
    else m_keep_sample_list = load_ref(keep_sample, ignore_fid);

    m_xymt_codes.resize(XYMT_OFFSET_CT);
    init_chr(num_auto, no_x, no_y, no_xy, no_mt);
    m_thread = thread;
    set_genotype_files(prefix);
    //m_bgen_info = load_bgen_info();
    m_sample_names = preload_samples(pheno_file, header, ignore_fid);
    // here we don't actually want the return value
    load_samples(ignore_fid);

    m_existed_snps = load_snps();
    if(verbose)
    {
        fprintf(stderr, "%zu people (%zu males, %zu females) included\n", m_unfiltered_sample_ct, m_num_male, m_num_female);
        if(m_num_ambig!=0) fprintf(stderr, "%u ambiguous variants excluded\n", m_num_ambig);
        fprintf(stderr, "%zu variants included\n", m_marker_ct);
    }
    //get_header();
    exit(0);

}

std::unordered_map<std::string, BinaryGen::Context> BinaryGen::load_bgen_info()
{
    std::unordered_map<std::string, Context> info;
    for(auto &&prefix : m_genotype_files)
    {
        Context cur_context;
        cur_context.flags =0;
        cur_context.number_of_samples=0;
        cur_context.number_of_variants = 0;
        cur_context.offset = 0;
        info[prefix] = cur_context;
    }
    return info;
}

BinaryGen::~BinaryGen()
{
    std::cerr << "Destroy!" << std::endl;
}

std::vector<SNP> BinaryGen::load_snps()
{
    std::vector<SNP> snp_res;
    bool chr_sex_error = false;
    bool chr_error = false;
    size_t order=0;
    m_num_ambig=0;
    size_t expected_total=0;
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
        uint32_t offset = m_bgen_info[prefix].offset;
        m_bgen_file.seekg(offset+4);
        uint32_t num_snp = m_bgen_info[prefix].number_of_variants;
        uint32_t const layout = m_bgen_info[prefix].flags & e_Layout ;
        uint32_t num_sample = m_bgen_info[prefix].number_of_samples;
        std::string prev_chr="";
        int chr_code=0;

        for(size_t i_snp =0; i_snp < num_snp; ++i_snp)
        {
            if(m_unfiltered_marker_ct%1000==0)
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
            if( layout == e_Layout1 || layout == e_Layout0 ) {
                uint32_t number_of_samples ;
                read_little_endian_integer( &number_of_samples ) ;
                if( number_of_samples != num_sample ) {
                    throw std::runtime_error("ERROR: Number of sample doesn't match!");
                }
                read_length_followed_by_data( &SNPID_size, &SNPID ) ;
            } else if( layout == e_Layout2 ) {
                read_length_followed_by_data( &SNPID_size, &SNPID ) ;
            } else {
                assert(0) ;
            }
            read_length_followed_by_data( &RSID_size, &RSID ) ;
            read_length_followed_by_data( &chromosome_size, &chromosome ) ;
            read_little_endian_integer( &SNP_position ) ;
            if( layout == e_Layout2 ) {
                read_little_endian_integer( &numberOfAlleles ) ;
            } else {
                numberOfAlleles = 2 ;
            }
            if(numberOfAlleles!=2)
            {
                throw std::runtime_error("ERROR: Currently only support bgen with 2 alleles!");
                // we can, in the future, allow for more than 2 alleles by setting anything
                // but the 2 most common alleles to missing
            }
            std::vector<std::string> final_alleles(numberOfAlleles);
            for( uint16_t i = 0; i < numberOfAlleles; ++i ) {
                read_length_followed_by_data( &allele_size, &allele ) ;
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
            std::vector< byte_t > buffer;
            read_genotype_data_block(prefix, &buffer); // move read pointer forward
            if(chromosome.compare(prev_chr)!=0)
            {
                prev_chr = chromosome;
                if(m_chr_order.find(chromosome)!= m_chr_order.end())
                {
                    throw std::runtime_error("ERROR: SNPs on the same chromosome must be clustered together!");
                }
                m_chr_order[chromosome] = order++;
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
            if(m_existed_snps_index.find(RSID)!= m_existed_snps_index.end())
            {
                throw std::runtime_error("ERROR: Duplicated SNP ID detected!\n");
            }
            else if(ambiguous(final_alleles.front(), final_alleles.back()))
            {
                m_num_ambig++;
            }
            else
            {
                if(RSID.compare(".")==0) // when the rs id isn't available, change it to chr:loc coding
                {
                    RSID = std::to_string(chr_code)+":"+std::to_string(SNP_position);
                }
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
    m_marker_ct = snp_res.size();
    return snp_res;
}

void BinaryGen::cleanup()
{
    return;
}

void BinaryGen::read_score( std::vector<std::vector<Sample_lite> > &current_prs_score,
               size_t start_index, size_t end_bound)
{
    std::cerr << "Testing" << std::endl;
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
        if(m_remove_sample)
        {
            if(m_keep_sample && cur_sample.included)
            {
                std::string error_message = "ERROR: Sample ID: "+id+" existed in both remove and keep option!";
                throw std::runtime_error(error_message);
            }
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
            sample_res.push_back(get_sample(possible_header, ignore_fid, has_sex, sex_col, sex_info));
            m_sample_index_check[sample_res.back().IID] = sample_res.size()-1;
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
    m_unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    m_unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;

    m_sex_male = new uintptr_t[m_unfiltered_sample_ctl];
    std::memset(m_sex_male, 0x0, m_unfiltered_sample_ctl*sizeof(uintptr_t));

    // assume all are founder
    m_founder_info = new uintptr_t[m_unfiltered_sample_ctl];
    std::memset(m_founder_info, ~0, m_unfiltered_sample_ctl*sizeof(uintptr_t));

    m_sample_exclude = new uintptr_t[m_unfiltered_sample_ctl];
    std::memset(m_sample_exclude, 0x0, m_unfiltered_sample_ctl*sizeof(uintptr_t));

    m_num_male = 0, m_num_female = 0, m_num_ambig_sex=(has_sex)? 0 : m_unfiltered_sample_ct;
    for(size_t i = 0; i < sex_info.size() && has_sex; ++i)
    {
        switch(sex_info[i])
        {
            case 1:
                m_num_male++;
                SET_BIT(i, m_sex_male);
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
        read_little_endian_integer(&offset);
        uint32_t header_size = 0, number_of_snp_blocks = 0, number_of_samples = 0, flags = 0 ;
        char magic[4] ;
        std::size_t fixed_data_size = 20 ;
        std::vector<char> free_data ;
        read_little_endian_integer( &header_size ) ;
        assert( header_size >= fixed_data_size ) ;
        read_little_endian_integer( &number_of_snp_blocks ) ;
        read_little_endian_integer( &number_of_samples ) ;
        m_bgen_file.read( &magic[0], 4 ) ;
        free_data.resize( header_size - fixed_data_size ) ;
        m_bgen_file.read( &free_data[0], free_data.size() ) ;
        read_little_endian_integer( &flags ) ;
        if(( magic[0] != 'b' || magic[1] != 'g' || magic[2] != 'e' || magic[3] != 'n' )
                && ( magic[0] != 0 || magic[1] != 0 || magic[2] != 0 || magic[3] != 0 ))
        {
            throw std::runtime_error("ERROR: Incorrect magic string!\nPlease check you have provided a valid bgen file!");
        }
        if( m_bgen_file )
        {
            Context current_context;
            current_context.offset = offset;
            current_context.number_of_samples = number_of_samples ;
            current_context.number_of_variants = number_of_snp_blocks ;
            current_context.magic.assign( &magic[0], &magic[0] + 4 ) ;
            //current_context.free_data.assign( free_data.begin(), free_data.end() ) ;
            current_context.flags = flags ;
            m_bgen_info[prefix] = current_context;
        }
        else
        {
            throw std::runtime_error("ERROR: Problem reading bgen file!") ;
        }
        if(m_bgen_info[prefix] .number_of_samples!= m_sample_names.size())
        {
            std::string error_message = "ERROR: Number of sample in bgen does not match those in phenotype file! ("
                    +std::to_string(m_bgen_info[prefix] .number_of_samples)+" vs "+
                    std::to_string(m_sample_names.size())+")";
            throw std::runtime_error(error_message);
        }

        uint32_t const compressionType = m_bgen_info[prefix] .flags & e_CompressedSNPBlocks;
        if( compressionType == e_ZstdCompression )
        {
            throw std::runtime_error("ERROR: zstd compression currently not supported");
        }
        if((m_bgen_info[prefix] .flags & e_SampleIdentifiers) && first)
        { // only read in the sample information for the first bgen
            first = false;
            uint32_t sample_block_size = 0 ;
            uint32_t actual_number_of_samples = 0 ;
            uint16_t identifier_size ;
            std::string identifier ;
            std::size_t bytes_read = 0 ;
            // the bgen format actually double stored the number of samples and
            // the block size is the sample_block size
            read_little_endian_integer( &sample_block_size ) ;
            read_little_endian_integer( &actual_number_of_samples ) ;
            bytes_read += 8 ;
            assert( actual_number_of_samples == m_bgen_info[prefix] .number_of_samples ) ;
            for( uint32_t i = 0; i < actual_number_of_samples; ++i ) {
                read_length_followed_by_data( &identifier_size, &identifier ) ;
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

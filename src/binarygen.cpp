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
    m_sample_names = preload_samples(pheno_file, header, ignore_fid);
    // here we don't actually want the return value
    load_samples(ignore_fid);
    //get_header();
    exit(0);

}


BinaryGen::~BinaryGen()
{
    std::cerr << "Destroy!" << std::endl;
}

std::vector<SNP> BinaryGen::load_snps()
{
    std::vector<SNP> snp_res;
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


void BinaryGen::set_offset()
{
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
        //m_offset_info[bgen_name] = offset;
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
    bool has_sex=false;
    int sex_col = 0;
    if(has_header)
    {
        std::getline(pheno_file, line);
        misc::trim(line);
        if(line.empty()) throw std::runtime_error("ERROR: Empty header line for phenotype file!");
        std::vector<std::string> token = misc::split(line);
        if(token.size() < 1+!ignore_fid)
        {
            std::string error_message = "ERROR: Header line must contain at least "+
                    std::to_string(1+!ignore_fid)+" columns!";
            throw std::runtime_error(error_message);
        }
        for(size_t i = 1+!ignore_fid; i < token.size(); ++i)
        {
            if(token[i].compare("Sex")==0 || token[i].compare("sex")==0
                    || token[i].compare("SEX")==0)
            {
                has_sex = true;
                sex_col = i;
                break;
            }
        }
    }
    size_t index =0;
    bool first_line  = true;
    std::vector<int> sex_info;
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
        if(first_line)
        {
            //this is to check if the input is a bgen sample file
            first_line = false;
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
            if(bgen_sample)
            {
                fprintf(stderr, "Detected bgen sample format");
                if(has_sex)
                {
                    if(token[sex_col].compare("D")!=0)
                    {
                        std::string error_message= "ERROR: Sex must be coded as \"D\" in bgen sample file!";
                        throw std::runtime_error(error_message);
                    }
                }
                continue;
            }
        }
        std::string id = (ignore_fid)? token[0] : token[0]+"_"+token[1];
        if(m_sample_index_check.find(id)==m_sample_index_check.end())
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
            sample_res.push_back(cur_sample);
            m_sample_index_check[cur_sample.IID] = index++;
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
        }
        else
        {
            std::string error_message = "ERROR: Duplicated sample: "+ id;
            throw std::runtime_error(error_message);
        }
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

    m_num_male = 0, m_num_female = 0, m_num_ambig_sex=0;
    for(size_t i = 0; i < sex_info.size(); ++i)
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
        Context current_context;
        current_context.offset = offset;
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
            current_context.number_of_samples = number_of_samples ;
            current_context.number_of_variants = number_of_snp_blocks ;
            current_context.magic.assign( &magic[0], &magic[0] + 4 ) ;
            //current_context.free_data.assign( free_data.begin(), free_data.end() ) ;
            current_context.flags = flags ;
        }
        else
        {
            throw std::runtime_error("ERROR: Problem reading bgen file!") ;
        }
        m_bgen_info[bgen_name] = current_context;
        if(current_context.number_of_samples!= m_sample_names.size())
        {
            std::string error_message = "ERROR: Number of sample in bgen does not match those in phenotype file! ("
                    +std::to_string(current_context.number_of_samples)+" vs "+
                    std::to_string(m_sample_names.size())+")";
            throw std::runtime_error(error_message);
        }
        if(!(current_context.flags & e_SampleIdentifiers))
        {
            throw std::runtime_error("ERROR: BGEN file does not contain sample information!");
        }
        else if(first)
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
            assert( actual_number_of_samples == current_context.number_of_samples ) ;
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

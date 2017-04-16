#include "bgen.hpp"

BGEN::BGEN(std::string prefix,  std::string remove_sample, std::string keep_sample,
        bool ignore_fid, int num_auto, bool no_x, bool no_y, bool no_xy,
        bool no_mt, const size_t thread, bool verbose)
{
    context.number_of_variants = 0;

    if(remove_sample.empty()) m_remove_sample = false;
    else m_remove_sample_list = load_ref(remove_sample, ignore_fid);
    if(keep_sample.empty()) m_keep_sample = false;
    else m_keep_sample_list = load_ref(keep_sample, ignore_fid);

    m_xymt_codes.resize(XYMT_OFFSET_CT);
    init_chr(num_auto, no_x, no_y, no_xy, no_mt);
    m_thread = thread;
    set_genotype_files(prefix);

    read_offset();
    m_sample_names = load_samples(ignore_fid);
    m_existed_snps = load_snps();
    if(verbose)
    {
        fprintf(stderr, "%zu people (%zu males, %zu females) included\n", m_unfiltered_sample_ct, m_num_male, m_num_female);
        if(m_num_ambig!=0) fprintf(stderr, "%u ambiguous variants excluded\n", m_num_ambig);
        fprintf(stderr, "%zu variants included\n", m_marker_ct);
    }

    m_founder_ctl = BITCT_TO_WORDCT(m_founder_ct);
    m_founder_ctv3 = BITCT_TO_ALIGNED_WORDCT(m_founder_ct);
    m_founder_ctsplit = 3 * m_founder_ctv3;
    m_bgen_file.reset();
}

std::vector<SNP> BGEN::load_snps()
{
    std::vector<SNP> snp_result;
    for(auto &&prefix : m_gneotype_files)
    {
        std::string bgen_name = prefix+".bgen";
        m_bgen_file.reset(new std::ifstream(bgen_name.c_str(), std::ifstream::binary ));
        m_bgen_file->seekg( m_offset + 4 ) ;
        std::string SNPID ; // read but ignored in this toy implementation
        std::vector< std::string > alleles; //restrict ourselves to 2 Allele at the moment
        // the main problem is, if there is more than 2 Alleles, then we don't know how to
        // convert it into plink binary
        genfile::bgen::read_snp_identifying_data(
                *m_bgen_file, context,
                &SNPID, rsid, chromosome, position,
                [&alleles]( std::size_t n ) { alleles->resize( n ) ; },
                [&alleles]( std::size_t i, std::string const& allele ) { alleles->at(i) = allele ; }
        )
    }

}

bool BGEN::read_snp_identifying_data(std::istream& aStream,
        std::string* SNPID,
        std::string* RSID,
        std::string* chromosome,
        uint32_t* SNP_position,
        NumberOfAllelesSetter set_number_of_alleles,
        AlleleSetter set_allele)
{

    uint16_t SNPID_size = 0;
    uint16_t RSID_size = 0;
    uint16_t numberOfAlleles = 0 ;
    uint16_t chromosome_size = 0 ;
    uint32_t allele_size = 0;
    std::string allele;
    uint32_t const layout = context.flags & e_Layout ;

    // If we can't read a valid first field we return false; this will indicate EOF.
    // Any other fail to read is an error and an exception will be thrown.
    if( layout == e_Layout1 || layout == e_Layout0 ) {
        uint32_t number_of_samples ;
        try {
            read_little_endian_integer( aStream, &number_of_samples ) ;
        } catch( BGenError const& ) {
            return false ;
        }
        if( number_of_samples != context.number_of_samples ) {
            throw BGenError() ;
        }
        read_length_followed_by_data( aStream, &SNPID_size, SNPID ) ;
    } else if( layout == e_Layout2 ) {
        try {
            read_length_followed_by_data( aStream, &SNPID_size, SNPID ) ;
        } catch( BGenError const& ) {
            return false ;
        }
    } else {
        assert(0) ;
    }

    read_length_followed_by_data( aStream, &RSID_size, RSID ) ;
    read_length_followed_by_data( aStream, &chromosome_size, chromosome ) ;
    read_little_endian_integer( aStream, SNP_position ) ;
    if( layout == e_Layout2 ) {
        read_little_endian_integer( aStream, &numberOfAlleles ) ;
    } else {
        numberOfAlleles = 2 ;
    }
    set_number_of_alleles( numberOfAlleles ) ;
    for( uint16_t i = 0; i < numberOfAlleles; ++i ) {
        read_length_followed_by_data( aStream, &allele_size, &allele ) ;
        set_allele( i, allele ) ;
    }
    if( !aStream ) {
#if DEBUG_BGEN_FORMAT
        std::cerr << "bgen: layout = " << layout << ", alleles = " << numberOfAlleles << ".\n" << std::flush ;
        std::cerr << *SNPID << ", " << *RSID << ", " << *chromosome << ", " << *SNP_position << ".\n" << std::flush ;
#endif
        throw BGenError() ;
    }
    return true ;

}

std::vector<Sample> BGEN::load_samples(bool ignore_fid)
{
    std::string bgen_name = m_genotype_files.front()+".bgen";
    m_bgen_file.reset(new std::ifstream( bgen_name.c_str(), std::ifstream::binary ));
    read_header_block(*m_bgen_file);

    if(!context.flags & e_SampleIdentifiers)
    {
        std::string error_message = "ERROR: Require sample name information!";
        error_message.append("       Information not found in bgen file");
        throw std::runtime_error(error_message);
    }
    uint32_t block_size = 0 ;
    uint32_t number_of_samples = 0 ;
    uint16_t identifier_size ;
    std::string identifier ;
    std::size_t bytes_read = 0 ;

    read_little_endian_integer( *m_bgen_file, &block_size ) ;
    read_little_endian_integer( *m_bgen_file, &number_of_samples ) ;
    bytes_read += 8 ;
    assert( number_of_samples == context.number_of_samples ) ;

    m_unfiltered_sample_ct = context.number_of_samples;
    m_unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct);
    m_unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) / 4;

    m_sex_male = new uintptr_t[m_unfiltered_sample_ctl];
    std::memset(m_sex_male, 0x0, m_unfiltered_sample_ctl*sizeof(uintptr_t));

    // assume all to be founder
    m_founder_info = new uintptr_t[m_unfiltered_sample_ctl];
    std::memset(m_founder_info, ~0, m_unfiltered_sample_ctl*sizeof(uintptr_t));

    m_sample_exclude = new uintptr_t[m_unfiltered_sample_ctl];
    std::memset(m_sample_exclude, 0x0, m_unfiltered_sample_ctl*sizeof(uintptr_t));

    // Sex information was not provided for bgen file
    m_num_male = 0, m_num_female = 0, m_num_ambig_sex=m_unfiltered_sample_ct;
    std::vector<Sample> result;
    for( uint32_t i = 0; i < number_of_samples; ++i )
    {
        read_length_followed_by_data( *m_bgen_file, &identifier_size, &identifier ) ;
        if( *m_bgen_file ) {
            bytes_read += sizeof( identifier_size ) + identifier_size ;
            Sample cur_sample;
            cur_sample.FID="";
            cur_sample.IID=identifier;
            cur_sample.included = false;
            cur_sample.pheno = "NA";
            cur_sample.prs = 0;
            cur_sample.num_snp = 0;
            if(m_keep_sample)
            {
                cur_sample.included=(m_keep_sample_list.find(identifier)!=m_keep_sample_list.end());
            }
            if(m_remove_sample)
            {
                if(m_keep_sample && cur_sample.included)
                {
                    std::string error_message = "ERROR: Sample ID: "+identifier+" existed in both remove and keep option!";
                    throw std::runtime_error(error_message);
                }
                cur_sample.included = !(m_remove_sample_list.find(identifier)!=m_remove_sample_list.end());
            }
            else cur_sample.included = true;
            result.push_back(cur_sample);
        } else {
            throw std::runtime_error("BGEN ERROR!");
        }
    }
    assert( bytes_read == block_size ) ;
    return result;
}


void BGEN::read_offset()
{
    // make sure all file has the same offset
    uint32_t prev_offset=0;
    bool first_offset=true;
    for(auto &&prefix : m_genotype_files)
    {
        std::string bgen_name = prefix+".bgen";
        m_bgen_file.reset(new std::ifstream( bgen_name.c_str(), std::ifstream::binary ));
        if( (m_bgen_file->rdstate() & std::ifstream::failbit ) != 0)
        {
            std::string error_message = "ERROR: Cannot open bgen file: "+bgen_name;
            throw std::runtime_error(error_message);
        }
        read_offset(*m_bgen_file, &m_offset);
        if(first_offset)
        {
            first_offset = false;
            prev_offset = m_offset;
        }
        else if(prev_offset != m_offset)
        {
            std::string error_message = "ERROR: For simplicity, we require all BGEN file to have the same offset!\n";
            error_message.append("      "+bgen_name+" has a different offset!");
            throw std::runtime_error(error_message);
        }
    }
}

std::size_t BGEN::read_header_block(std::istream& aStream)
{
    assert( context != 0 ) ;
    uint32_t
    header_size = 0,
    number_of_snp_blocks = 0,
    number_of_samples = 0,
    flags = 0 ;

    char magic[4] ;
    std::size_t fixed_data_size = 20 ;
    std::vector<char> free_data ;

    read_little_endian_integer( aStream, &header_size ) ;
    assert( header_size >= fixed_data_size ) ;
    read_little_endian_integer( aStream, &number_of_snp_blocks ) ;
    read_little_endian_integer( aStream, &number_of_samples ) ;
    aStream.read( &magic[0], 4 ) ;
    free_data.resize( header_size - fixed_data_size ) ;
    aStream.read( &free_data[0], free_data.size() ) ;
    read_little_endian_integer( aStream, &flags ) ;
    if(
            ( magic[0] != 'b' || magic[1] != 'g' || magic[2] != 'e' || magic[3] != 'n' )
            && ( magic[0] != 0 || magic[1] != 0 || magic[2] != 0 || magic[3] != 0 )
    )
    {
        throw std::runtime_error("BGEN ERROR!");
    }
    if( aStream ) {
        context.number_of_samples = number_of_samples ;
        context.number_of_variants += number_of_snp_blocks ;
        context.magic.assign( &magic[0], &magic[0] + 4 ) ;
        context.free_data.assign( free_data.begin(), free_data.end() ) ;
        context.flags = flags ;

        return( header_size ) ;
    } else {
        throw std::runtime_error("BGEN ERROR!");
    }
}

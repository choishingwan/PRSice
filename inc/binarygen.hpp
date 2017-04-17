
#ifndef BinaryGEN_H
#define BinaryGEN_H

#include "genotype.hpp"

class BinaryGen: public Genotype
{
    public:
        BinaryGen(std::string prefix, std::string pheno_file, bool header,
                std::string remove_sample, std::string keep_sample,
                bool ignore_fid, int num_auto = 22, bool no_x = false,
                bool no_y = false, bool no_xy = false, bool no_mt = false,
                const size_t thread = 1, bool verbose = false);
        ~BinaryGen();
    private:
        enum FlagMask { e_NoFlags = 0, e_CompressedSNPBlocks = 0x3, e_Layout = 0x3C } ;
        enum Layout { e_Layout0 = 0x0, e_Layout1 = 0x4, e_Layout2 = 0x8 } ;
        enum Structure { e_SampleIdentifiers = 0x80000000 } ;
        enum Compression { e_NoCompression = 0, e_ZlibCompression = 1, e_ZstdCompression = 2 } ;
        typedef uint8_t byte_t;
        struct Context {
            uint32_t flags ;
            uint32_t number_of_samples ;
            uint32_t number_of_variants ;
            uint32_t offset;
            std::string magic;
        };
        std::vector<Sample> preload_samples(std::string pheno, bool has_header, bool ignore_fid);
        std::unordered_map<std::string, int> m_sample_index_check;
        std::vector<Sample> load_samples(bool ignore_fid);
        std::vector<SNP> load_snps();
        void cleanup();
        inline void read_genotype(uintptr_t* genotype, const uint32_t snp_index,
                const std::string &file_name)
        {
            std::cerr << "Checking" << std::endl;
        };

        void read_score( std::vector<std::vector<Sample_lite> > &current_prs_score,
                size_t start_index, size_t end_bound);

        void set_offset();
        std::unordered_map<std::string, Context> m_bgen_info;
        std::ifstream m_bgen_file;

/**
 *  DON'T TOUCH (FROM BGEN LIBRARY)
 */

        template< typename IntegerType >
        void read_little_endian_integer(IntegerType* integer_ptr )
        {
            byte_t buffer[ sizeof( IntegerType ) ] ;
            m_bgen_file.read( reinterpret_cast< char* >( buffer ), sizeof( IntegerType )) ;
            if( !m_bgen_file )
            {
                throw std::runtime_error("ERROR: Cannot read bgen file!");
            }
            read_little_endian_integer( buffer, buffer + sizeof( IntegerType ), integer_ptr ) ;
        }

        template< typename IntegerType >
        byte_t const* read_little_endian_integer( byte_t const* buffer, byte_t const* const end,
                IntegerType* integer_ptr )
        {
            assert( end >= buffer + sizeof( IntegerType )) ;
            *integer_ptr = 0 ;
            for( std::size_t byte_i = 0; byte_i < sizeof( IntegerType ); ++byte_i )
            {
                (*integer_ptr) |= IntegerType( *reinterpret_cast< byte_t const* >( buffer++ )) << ( 8 * byte_i ) ;
            }
            return buffer ;
        }
        template< typename IntegerType >
        void read_length_followed_by_data(IntegerType* length_ptr, std::string* string_ptr )
        {
            IntegerType& length = *length_ptr ;
            read_little_endian_integer( length_ptr ) ;
            std::vector< char >buffer ( length ) ;
            m_bgen_file.read( &buffer[0], length ) ;
            if( !m_bgen_file )
            {
                throw std::runtime_error("ERROR: Problem reading bgen file!") ;
            }
            string_ptr->assign( buffer.begin(), buffer.end() ) ;
        }




};

#endif

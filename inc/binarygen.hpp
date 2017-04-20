
#ifndef BinaryGEN_H
#define BinaryGEN_H

#include "genotype.hpp"
#include <zlib.h>

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
        Sample get_sample(std::vector<std::string> &token, bool ignore_fid,
                bool has_sex, int sex_col, std::vector<int> &sex_info);
        std::vector<Sample> preload_samples(std::string pheno, bool has_header, bool ignore_fid);
        std::unordered_map<std::string, int> m_sample_index_check;

        std::vector<Sample> load_samples(bool ignore_fid);


        std::vector<SNP> load_snps();
        void cleanup();
        std::string m_cur_file;
        inline void read_genotype(uintptr_t* genotype, const uint32_t snp_index,
                const std::string &file_name)
        {
            if(m_cur_file.empty() || file_name.compare(m_cur_file)!=0)
            {
                if(m_bgen_file.is_open()) m_bgen_file.close();
                std::string bgen_name = file_name+".bgen";
                m_bgen_file.open(bgen_name.c_str(), std::ifstream::binary);
                if(!m_bgen_file.is_open())
                {
                    std::string error_message = "ERROR: Cannot open bgen file: "+file_name;
                    throw std::runtime_error(error_message);
                }
                m_cur_file = file_name;
            }
            m_bgen_file.seekg(snp_index, std::ios_base::beg);
            std::vector< byte_t > buffer1, buffer2;
            read_genotype_data_block( file_name, &buffer1 ) ;
            uncompress_probability_data( file_name, buffer1, &buffer2 ) ;

            // super difficult here... if we use plink method,
            // we need both the compressed and uncompressed handling
            // but then if we use the bgen library, then it is
            // all mixed
            std::cerr << "Probability check: " << snp_index <<"\t" << buffer2[0] << std::endl;
            exit(0);
        };

        void read_score( std::vector<std::vector<Sample_lite> > &current_prs_score,
                size_t start_index, size_t end_bound);
        std::unordered_map<std::string, Context> load_bgen_info();
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

        void read_genotype_data_block(const std::string &file_name, std::vector< byte_t >* buffer)
        {
            uint32_t payload_size = 0 ;
            uint32_t flags = m_bgen_info[file_name].flags;
            if( (flags & e_Layout) == e_Layout2 || ((flags & e_CompressedSNPBlocks) != e_NoCompression ) )
            {
                read_little_endian_integer( &payload_size ) ;
            }
            else
            {
                payload_size = 6 * m_bgen_info[file_name].number_of_samples; ;
            }
            buffer->resize( payload_size ) ;
            m_bgen_file.read( reinterpret_cast< char* >( &(*buffer)[0] ), payload_size ) ;
        }

        void uncompress_probability_data( const std::string &file_name,
                std::vector< byte_t > const& compressed_data, std::vector< byte_t >* buffer )
        {
            // compressed_data contains the (compressed or uncompressed) probability data
            uint32_t flags = m_bgen_info[file_name].flags;
            uint32_t const compressionType = (flags & e_CompressedSNPBlocks) ;
            if( compressionType != e_NoCompression )
            {
                byte_t const* begin = &compressed_data[0] ;
                byte_t const* const end = &compressed_data[0] + compressed_data.size() ;
                uint32_t uncompressed_data_size = 0 ;
                if( (flags & e_Layout) == e_Layout1 )
                {
                    uncompressed_data_size = 6 * m_bgen_info[file_name].number_of_samples ;
                }
                else
                {
                    begin = read_little_endian_integer( begin, end, &uncompressed_data_size ) ;
                }
                buffer->resize( uncompressed_data_size ) ;
                if( compressionType == e_ZlibCompression )
                {
                    zlib_uncompress( begin, end, buffer ) ;
                } else if( compressionType == e_ZstdCompression ) {
                    throw std::runtime_error("ERROR: Currently don't support BGEN v1.3");
                }
                assert( buffer->size() == uncompressed_data_size ) ;
            }
            else
            {
                // copy the data between buffers.
                buffer->assign( compressed_data.begin(), compressed_data.end() ) ;
            }
        };
        template< typename T >
        void zlib_uncompress(
                byte_t const* begin,
                byte_t const* const end,
                std::vector< T >* dest)
        {
            uLongf const source_size = ( end - begin ) ;
            uLongf dest_size = dest->size() * sizeof( T ) ;
            int result = uncompress(
                    reinterpret_cast< Bytef* >( &dest->operator[]( 0 ) ),
                    &dest_size,
                    reinterpret_cast< Bytef const* >( begin ),
                    source_size
            ) ;
            assert( result == Z_OK ) ;
            assert( dest_size % sizeof( T ) == 0 ) ;
            dest->resize( dest_size / sizeof( T )) ;
        }

};

#endif

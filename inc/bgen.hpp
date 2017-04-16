#ifndef BGEN_H
#define BGEN_H
#include "genotype.hpp"

// ideally, we would like to use the bgen library directly
// but that'd mean we need to include zlib and more boost
// so try to work on it ourselve
class BGEN: public Genotype
{
    public:
        BGEN(std::string prefix, std::string remove_sample,
                std::string keep_sample, bool ignore_fid, int num_auto = 22,
                bool no_x = false, bool no_y = false, bool no_xy = false,
                bool no_mt = false, const size_t thread = 1, bool verbose =
                        false);
        ~BGEN();
    private:
        typedef uint8_t byte_t ;
        std::unique_ptr< std::istream > m_bgen_file; // bgen library use ifstream instead of FILE
        enum Structure { e_SampleIdentifiers = 0x80000000 };
        enum Layout { e_Layout0 = 0x0, e_Layout1 = 0x4, e_Layout2 = 0x8 };
        std::size_t read_header_block(std::istream& aStream);
        struct{
            int32_t number_of_samples ;
            uint32_t number_of_variants ;
            std::string magic ;
            std::string free_data ;
            uint32_t flags ;
        } context;

        void read_offset();
        uint32_t m_offset ;


        std::vector<Sample> load_samples(bool ignore_fid);
        std::vector<SNP> load_snps();
        void cleanup(){};

        inline void read_genotype(uintptr_t* genotype, const uint32_t snp_index,
                const std::string &file_name){}

        void read_score(
                std::vector<std::vector<Sample_lite> > &current_prs_score,
                size_t start_index, size_t end_bound);

        /**
         * DON'T TOUCH AREA (Obtained from bgen library)
         */
        void read_offset( std::istream& iStream, uint32_t* offset )
        {
            read_little_endian_integer( iStream, offset ) ;
        }
        template< typename IntegerType >
        byte_t const* read_little_endian_integer( byte_t const* buffer,
                byte_t const* const end, IntegerType* integer_ptr )
        {
            assert( end >= buffer + sizeof( IntegerType )) ;
            *integer_ptr = 0 ;
            for( std::size_t byte_i = 0; byte_i < sizeof( IntegerType ); ++byte_i )
            {
                (*integer_ptr) |= IntegerType( *reinterpret_cast< byte_t const* >( buffer++ )) << ( 8 * byte_i ) ;
            }
            return buffer ;
        }

        // Read an integer stored in little-endian format into an integer stored in memory.
        // The stream is assumed to have sizeof( Integertype ) readable bytes.
        template< typename IntegerType >
        void read_little_endian_integer( std::istream& in_stream, IntegerType* integer_ptr )
        {
            byte_t buffer[ sizeof( IntegerType ) ] ;
            in_stream.read( reinterpret_cast< char* >( buffer ), sizeof( IntegerType )) ;
            if( !in_stream )
            {
                throw std::runtime_error("BGEN Error!") ;
            }
            read_little_endian_integer( buffer, buffer + sizeof( IntegerType ), integer_ptr ) ;
        }

        template< typename IntegerType >
        void read_length_followed_by_data( std::istream& in_stream, IntegerType* length_ptr, std::string* string_ptr ) {
            IntegerType& length = *length_ptr ;
            read_little_endian_integer( in_stream, length_ptr ) ;
            std::vector< char >buffer ( length ) ;
            in_stream.read( &buffer[0], length ) ;
            if( !in_stream ) {
                throw std::runtime_error("BGEN Error!") ;
            }
            string_ptr->assign( buffer.begin(), buffer.end() ) ;
        }
};
#endif

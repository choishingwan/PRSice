
#ifndef BinaryGEN_H
#define BinaryGEN_H

#include "genotype.hpp"
#include "bgen_lib.hpp"

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

        typedef std::vector< std::vector< double > > Data ;
        std::vector< genfile::byte_t > m_buffer1, m_buffer2 ;
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
            Data probs ;
            ProbSetter setter( &probs ) ;
            genfile::bgen::read_and_parse_genotype_data_block< ProbSetter >(
                    m_bgen_file,
                    m_bgen_info[file_name],
                    setter,
                    &m_buffer1,
                    &m_buffer2
            ) ;

            std::cerr << "Probability check: " << snp_index <<"\t" << probs[0][0] << std::endl;
            exit(0);
        };

        void read_score( std::vector<std::vector<Sample_lite> > &current_prs_score,
                size_t start_index, size_t end_bound);
        std::unordered_map<std::string, genfile::bgen::Context> m_bgen_info;
        std::unordered_map<std::string, uint32_t> m_offset_map;
        std::ifstream m_bgen_file;

        /** DON'T TOUCH      */
        struct ProbSetter {
            ProbSetter( Data* result ):
                m_result( result ),
                m_sample_i(0)
            {}

            // Called once allowing us to set storage.
            void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {
                m_result->clear() ;
                m_result->resize( number_of_samples ) ;
            }

            // If present with this signature, called once after initialise()
            // to set the minimum and maximum ploidy and numbers of probabilities among samples in the data.
            // This enables us to set up storage for the data ahead of time.
            void set_min_max_ploidy( uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries, uint32_t max_entries ) {
                for( std::size_t i = 0; i < m_result->size(); ++i ) {
                    m_result->at( i ).reserve( max_entries ) ;
                }
            }

            // Called once per sample to determine whether we want data for this sample
            bool set_sample( std::size_t i ) {
                m_sample_i = i ;
                // Yes, here we want info for all samples.
                return true ;
            }

            // Called once per sample to set the number of probabilities that are present.
            void set_number_of_entries(
                std::size_t ploidy,
                std::size_t number_of_entries,
                genfile::OrderType order_type,
                genfile::ValueType value_type
            ) {
                assert( value_type == genfile::eProbability ) ;
                m_result->at( m_sample_i ).resize( number_of_entries ) ;
                m_entry_i = 0 ;
            }

            // Called once for each genotype (or haplotype) probability per sample.
            void set_value( uint32_t, double value ) {
                m_result->at( m_sample_i ).at( m_entry_i++ ) = value ;
            }

            // Ditto, but called if data is missing for this sample.
            void set_value( uint32_t, genfile::MissingValue value ) {
                // Here we encode missing probabilities with -1
                m_result->at( m_sample_i ).at( m_entry_i++ ) = -1 ;
            }

            // If present with this signature, called once after all data has been set.
            void finalise() {
                // nothing to do in this implementation.
            }

        private:
            Data* m_result ;
            std::size_t m_sample_i ;
            std::size_t m_entry_i ;
        } ;
        void read_genotype_data_block(
                std::istream& aStream,
                genfile::bgen::Context const& context,
                std::vector< genfile::byte_t >* buffer
        ) {
            uint32_t payload_size = 0 ;
            if( (context.flags & genfile::bgen::e_Layout) == genfile::bgen::e_Layout2 ||
                    ((context.flags & genfile::bgen::e_CompressedSNPBlocks) != genfile::bgen::e_NoCompression ) ) {
                genfile::bgen::read_little_endian_integer( aStream, &payload_size ) ;
            } else {
                payload_size = 6 * context.number_of_samples ;
            }
            buffer->resize( payload_size ) ;
            aStream.read( reinterpret_cast< char* >( &(*buffer)[0] ), payload_size ) ;
        }
};

#endif

#ifndef MOCK_BINARYGEN_HPP
#define MOCK_BINARYGEN_HPP
#include "binarygen.hpp"
#include "reporter.hpp"

class mock_binarygen : public ::BinaryGen
{
public:
    mock_binarygen() {}
    mock_binarygen(GenoFile& geno, Phenotype& pheno, const std::string& delim,
                   Reporter* reporter)
        : BinaryGen(geno, pheno, delim, reporter)
    {
    }
    size_t test_get_sex_col(const std::string& header,
                            const std::string& format_line)
    {
        return get_sex_col(header, format_line);
    }
    void test_handle_pheno_header(std::unique_ptr<std::istream>& sample)
    {
        handle_pheno_header(sample);
    }
    void set_reporter(Reporter* reporter) { m_reporter = reporter; }
    void set_sample_size(uintptr_t sample_size)
    {
        m_unfiltered_sample_ct = sample_size;
    }
    void add_select_sample(const std::string& in)
    {
        m_sample_selection_list.insert(in);
    }
    void change_sample_selection(bool remove) { m_remove_sample = remove; }

    std::vector<uintptr_t> sample_for_ld() const { return m_sample_for_ld; }
    std::vector<uintptr_t> calculate_prs() const { return m_calculate_prs; }
    std::vector<Sample_ID> sample_id() const { return m_sample_id; }
    void gen_bgen_header(const std::string& file_name,
                         uint32_t number_of_snp_blocks,
                         uint32_t number_of_samples, std::string free_data,
                         uint32_t flags)
    {
        std::ofstream dummy(file_name, std::ofstream::binary);
        genfile::bgen::Context context;
        context.number_of_variants = number_of_snp_blocks;
        context.number_of_samples = number_of_samples;
        context.free_data = free_data;
        context.flags = flags;
        genfile::bgen::write_header_block(dummy, context);
        dummy.close();
    }

    void generate_bgen(const std::string& file_name,
                       uint32_t number_of_snp_blocks,
                       uint32_t number_of_samples, std::string free_data,
                       uint32_t flags)
    {
        /*
            void write_header_block(std::ostream & aStream, Context const&
           context); std::size_t write_sample_identifier_block( std::ostream &
           aStream, Context const& context, std::vector<std::string> const&
           sample_ids); byte_t* write_snp_identifying_data( std::vector<byte_t>
           * buffer, Context const& context, std::string SNPID, std::string
           RSID, std::string chromosome, uint32_t position, uint16_t const
           number_of_alleles, AlleleGetter get_allele);


            std::ostringstream outStream;
            a and b are A1 and A2 in string format
            std::vector<genfile::byte_t> buffer;
            std::vector<genfile::byte_t> buffer2;
            genfile::byte_t* end = genfile::bgen::write_snp_identifying_data(
                &buffer, context, SNPID, RSID, chromosome, SNP_position, 2,
                [&a, &b](std::size_t i) {
                    if (i == 0) { return a; }
                    else if (i == 1)
                    {
                        return b;
                    }
                    else
                    {
                        assert(0);
                    }
                });
            outStream.write(reinterpret_cast<char*>(&buffer[0]), end -
           &buffer[0]);
            */
    }
};

#endif // MOCK_BINARYGEN_HPP

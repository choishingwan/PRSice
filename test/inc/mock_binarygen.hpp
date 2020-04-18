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
        genfile::bgen::write_offset(dummy, context.free_data.length()+20);
        genfile::bgen::write_header_block(dummy, context);
        dummy.close();
    }
    size_t test_transverse_bgen_for_snp(
        const std::vector<IITree<size_t, size_t>>& exclusion_regions,
        const std::string mismatch_snp_record_name, const size_t file_idx,
        std::unique_ptr<std::istream> bgen_file,
        std::unordered_set<std::string>& duplicated_snps,
        std::unordered_set<std::string>& processed_snps,
        std::vector<bool>& retain_snp, bool& chr_error, bool& sex_error,
        Genotype* genotype)
    {
        return transverse_bgen_for_snp(
            exclusion_regions, mismatch_snp_record_name, file_idx,
            std::move(bgen_file), duplicated_snps, processed_snps, retain_snp,
            chr_error, sex_error, genotype);
    }
    void set_context(genfile::bgen::Context context, size_t idx)
    {
        if (m_context_map.size() < idx + 1) { m_context_map.resize(idx); }
        m_context_map[idx] = context;
    }
    void load_context(std::istream& in_file, size_t idx=0)
    {
        if(m_context_map.size() <= idx+1){
            m_context_map.resize(idx+1);
        }
        genfile::bgen::Context context;
        uint32_t offset;
        genfile::bgen::read_offset(in_file, &offset);
        genfile::bgen::read_header_block(in_file, &context);
        context.offset = offset;
        m_context_map[idx] = context;
    }
    void test_init_chr(int num_auto = 22, bool no_x = false, bool no_y = false,
                       bool no_xy = false, bool no_mt = false)
    {
        init_chr(num_auto, no_x, no_y, no_xy, no_mt);
    }
    std::string gen_mock_snp(std::vector<SNP>& input,
                             uint32_t number_individual,
                             genfile::OrderType& phased,
                             genfile::bgen::Layout& layout,
                             genfile::bgen::Compression& compress)
    {
        genfile::bgen::Context context;
        context.flags |= compress;
        context.flags |= layout;
        context.number_of_samples = number_individual;
        context.magic = "bgen";
        context.number_of_variants = static_cast<uint32_t>(input.size());
        context.free_data = "Test with automatically generated mock data";
        std::ostringstream mock_file, mock_w_offset;
        genfile::bgen::write_header_block(mock_file, context);
        genfile::bgen::write_little_endian_integer(mock_w_offset, context.header_size());
        for (auto&& snp : input)
        {
            std::string chr= snp.chr() > m_autosome_ct? "chrX": std::to_string(snp.chr());
            gen_ostringstream_for_snp(mock_file, context, snp.rs(), snp.rs(),
                                      chr, snp.loc(),
                                      snp.ref(), snp.alt(), phased);
        }
        mock_w_offset.write(mock_file.str().data(), mock_file.str().size());
        return (mock_w_offset.str());
    }
    void gen_ostringstream_for_snp(std::ostringstream& mock_file,
                                   genfile::bgen::Context& context,
                                   const std::string& snpid,
                                   const std::string& rsid,
                                   const std::string chr, const size_t loc,
                                   const std::string a1, const std::string a2,
                                   const genfile::OrderType& phased)
    {
        std::vector<genfile::byte_t> buffer, buffer2;
        genfile::byte_t* end = genfile::bgen::write_snp_identifying_data(
            &buffer, context, snpid, rsid, chr, loc, 2,
            [&a1, &a2](std::size_t i) {
                if (i == 0) { return a1; }
                return a2;
            });
        mock_file.write(reinterpret_cast<char*>(&buffer[0]), end - &buffer[0]);
        int bits_per_probability = 16;
        genfile::bgen::GenotypeDataBlockWriter writer(
            &buffer, &buffer2, context, bits_per_probability);
        writer.initialise(context.number_of_samples, 2);
        gen_zero_prob_vector(context.number_of_samples, phased, writer);
        mock_file.write(reinterpret_cast<char const*>(writer.repr().first),
                        writer.repr().second - writer.repr().first);

    }
    void gen_zero_prob_vector(const size_t number_of_individuals,
                              const genfile::OrderType& type,
                              genfile::bgen::GenotypeDataBlockWriter& writer)
    {
        for (size_t i = 0; i < number_of_individuals; ++i)
        {

            writer.set_sample(i);
            if (type == genfile::ePerUnorderedGenotype)
            {
                writer.set_number_of_entries(2, 3, type, genfile::eProbability);
                for (uint32_t j = 0; j < 3; ++j) { writer.set_value(j, 0.0); }
            }
            else
            {
                writer.set_number_of_entries( 2, 4, genfile::ePerPhasedHaplotypePerAllele, genfile::eProbability ) ;
               // writer.set_number_of_entries(2, 4, type, genfile::eProbability);
                for (uint32_t j = 0; j < 4; ++j) { writer.set_value(j, 0.0); }
            }
        }
        writer.finalise();
    }
    void manual_load_snp(SNP cur)
    {
        m_existed_snps_index[cur.rs()] = m_existed_snps.size();
        m_existed_snps.emplace_back(cur);
    }
    std::vector<SNP> existed_snps() const { return m_existed_snps; }
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

#ifndef MOCK_REGION_HPP
#define MOCK_REGION_HPP
#include "genotype.hpp"
#include "region.hpp"
#include "reporter.hpp"

class mock_region : public Region
{
public:
    bool test_in_region(const std::string_view& in,
                        const std::vector<std::string>& feature)
    {
        return in_feature(in, feature);
    }
    void set_reporter(Reporter* reporter) { m_reporter = reporter; }
    static bool test_valid_strand(std::string_view input)
    {
        return valid_strand(input);
    }
    static void test_extend_region(std::string_view strand, const size_t wind_5,
                                   const size_t wind_3, size_t& start,
                                   size_t& end)
    {
        extend_region(strand, wind_5, wind_3, start, end);
    }
    static std::tuple<size_t, size_t>
    test_start_end(const std::string& start_str, const std::string& end_str,
                   const bool zero_based)
    {
        return start_end(start_str, end_str, zero_based);
    }
    static std::tuple<size_t, size_t> test_start_end(std::string_view start_str,
                                                     std::string_view end_str,
                                                     const bool zero_based)
    {
        return start_end(start_str, end_str, zero_based);
    }
    std::tuple<std::string, std::string, bool>
    test_get_set_name(const std::string& input)
    {
        return get_set_name(input);
    }
    bool test_duplicated_set(const std::string_view set_name)
    {
        return duplicated_set(set_name);
    }
    bool test_duplicated_set(const std::string& set_name)
    {
        return duplicated_set(set_name);
    }
    static bool
    test_is_bed_header(const std::vector<std::string_view>& bed_line,
                       size_t& column_size)
    {
        return is_bed_header(bed_line, column_size);
    }
    std::tuple<std::string, bool> test_parse_gene_id(std::string_view substr,
                                                     std::string name)
    {
        return parse_gene_id(substr, name);
    }
    bool test_find_gene_info(std::string_view substr, std::string& gene_id,
                             std::string& gene_name, bool& found_id,
                             bool& found_name)
    {
        return find_gene_info(substr, gene_id, gene_name, found_id, found_name);
    }
    void test_parse_attribute(std::string_view attribute_str,
                              std::string& gene_id, std::string& gene_name)
    {
        parse_attribute(attribute_str, gene_id, gene_name);
    }
    void test_load_snp_sets(
        const std::unordered_map<std::string, size_t>& snp_list_idx,
        const std::vector<SNP>& snp_list, const std::string& snp_file,
        size_t& set_idx)
    {
        load_snp_sets(snp_list_idx, snp_list, snp_file, set_idx);
        for (auto&& tree : m_gene_sets) { tree.index(); }
    }
    void test_read_bed(std::unique_ptr<std::istream> bed,
                       std::vector<IITree<size_t, size_t>>& cr,
                       bool& print_bed_strand_warning, const size_t wind_5,
                       const size_t wind_3, const size_t max_chr,
                       const size_t set_idx, const bool ZERO_BASED)
    {
        read_bed(std::move(bed), cr, print_bed_strand_warning, wind_5, wind_3,
                 max_chr, set_idx, ZERO_BASED);
    }
    void test_transverse_snp_file(
        const std::unordered_map<std::string, size_t>& snp_list_idx,
        const std::vector<SNP>& snp_list, const bool is_set_file,
        std::unique_ptr<std::istream> input, size_t& set_idx)
    {
        transverse_snp_file(snp_list_idx, snp_list, is_set_file,
                            std::move(input), set_idx);
        for (auto&& tree : m_gene_sets) { tree.index(); }
    }
    void test_load_msigdb(
        std::unordered_map<std::string, std::vector<size_t>>& msigdb_list,
        std::unique_ptr<std::istream> input, size_t& set_idx)
    {
        load_msigdb(msigdb_list, std::move(input), set_idx);
    }
    bool test_load_bed_regions(const std::string& bed_file,
                               const size_t set_idx, const size_t max_chr)
    {
        return load_bed_regions(bed_file, set_idx, max_chr);
    }
    void test_load_gtf(
        const std::unordered_map<std::string, std::vector<size_t>>& msigdb_list,
        const size_t max_chr)
    {
        return load_gtf(msigdb_list, max_chr);
    }
    std::tuple<size_t, size_t, size_t> test_transverse_gtf(
        const std::unordered_map<std::string, std::vector<size_t>>& msigdb_list,
        const std::streampos file_length, const size_t max_chr,
        const bool gz_input, std::unique_ptr<std::istream> gtf_stream)
    {
        return transverse_gtf(msigdb_list, file_length, max_chr, gz_input,
                              std::move(gtf_stream));
    }
    void test_load_background(
        const std::unordered_map<std::string, size_t>& snp_list_idx,
        const std::vector<SNP>& snp_list, const size_t max_chr,
        std::unordered_map<std::string, std::vector<size_t>>& msigdb_list)
    {
        return load_background(snp_list_idx, snp_list, max_chr, msigdb_list);
    }
    void set_background(const std::string& b) { m_background = b; }
    void set_gtf(const std::string& g) { m_gtf = g; }
    void set_feature(const std::vector<std::string>& in) { m_feature = in; }
    void set_gwas_bk(bool has_bk) { m_genome_wide_background = has_bk; }
    void index()
    {
        for (auto&& tree : m_gene_sets) { tree.index(); }
    }
    void add_duplicated_set(const std::string& in)
    {
        m_processed_sets.insert(in);
    }
    void set_wind(const size_t w5, size_t w3)
    {
        m_window_5 = w5;
        m_window_3 = w3;
    }
};

#endif // MOCK_REGION_HPP

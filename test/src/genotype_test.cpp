#ifndef GENOTYPE_TEST_HPP
#define GENOTYPE_TEST_HPP
#include "genotype.hpp"
#include "gtest/gtest.h"
class GENOTYPE_BASIC : public Genotype, public ::testing::Test
{
public:
    void cleanup()
    {
        m_keep_file = "";
        m_sample_file = "";
        m_remove_file = "";
        m_genotype_file_names.clear();
        m_delim = "";
        m_ignore_fid = false;
    }
};
// This should be the simplest of the three genotype class. Test anything that
// doesn't require binaryplink and binarygen
// set_genotype_files

TEST_F(GENOTYPE_BASIC, CHR_CONVERT)
{
    ASSERT_FALSE(chr_prefix("1"));
    ASSERT_TRUE(chr_prefix("chr1"));
    auto chr = get_chrom_code("1");
    ASSERT_EQ(chr, 1);
    chr = get_chrom_code("chr10");
    ASSERT_EQ(chr, 10);
    chr = get_chrom_code("CHRX");
    ASSERT_EQ(chr, CHROM_X);
    chr = get_chrom_code("M");
    ASSERT_EQ(chr, CHROM_MT);
    chr = get_chrom_code("Y");
    ASSERT_EQ(chr, CHROM_Y);
    chr = get_chrom_code("XY");
    ASSERT_EQ(chr, CHROM_XY);
    chr = get_chrom_code("mt");
    ASSERT_EQ(chr, CHROM_MT);
    chr = get_chrom_code("CHR Happy");
    ASSERT_EQ(chr, -1);
    chr = get_chrom_code("chromosome");
    ASSERT_EQ(chr, -1);
    chr = get_chrom_code("0X");
    ASSERT_EQ(chr, CHROM_X);
    chr = get_chrom_code("0Y");
    ASSERT_EQ(chr, CHROM_Y);
    chr = get_chrom_code("0M");
    ASSERT_EQ(chr, CHROM_MT);
    // TODO: Add more crazy chromosome use cases
}

TEST_F(GENOTYPE_BASIC, INITIALIZE)
{
    GenoFile geno;
    geno.num_autosome = 22;
    Phenotype pheno;
    pheno.ignore_fid = false;
    std::string delim = " ";
    std::string type = "bed";
    Reporter reporter(std::string("LOG"), 60, true);
    // need to check the following
    // m_genotype_file_name is set
    // m_sample_file is set when external file is provided and empty if not
    // m_keep_file, m_remove_file and m_delim should be set
    // there should be no error check
    // Normal input
    geno.file_name = "Genotype";
    initialize(geno, pheno, delim, type, &reporter);
    ASSERT_STREQ(delim.c_str(), m_delim.c_str());
    ASSERT_EQ(m_genotype_file_names.size(), 1);
    ASSERT_STREQ(geno.file_name.c_str(), m_genotype_file_names.front().c_str());
    ASSERT_TRUE(m_keep_file.empty());
    ASSERT_TRUE(m_remove_file.empty());
    ASSERT_TRUE(m_sample_file.empty());
    ASSERT_FALSE(m_ignore_fid);
    cleanup();
    // now with external sample_file
    geno.file_name = "Genotype,External";
    delim = "\t";
    pheno.ignore_fid = true;
    initialize(geno, pheno, delim, type, &reporter);
    ASSERT_STREQ(delim.c_str(), m_delim.c_str());
    ASSERT_EQ(m_genotype_file_names.size(), 1);
    ASSERT_STREQ("Genotype", m_genotype_file_names.front().c_str());
    ASSERT_TRUE(m_keep_file.empty());
    ASSERT_TRUE(m_remove_file.empty());
    ASSERT_TRUE(m_ignore_fid);
    ASSERT_STREQ(m_sample_file.c_str(), "External");
    // now check replacement
    cleanup();
    geno.file_name = "chr#_geno";
    geno.keep = "Hi";
    initialize(geno, pheno, delim, type, &reporter);
    ASSERT_EQ(m_genotype_file_names.size(), geno.num_autosome);
    for (size_t i = 0; i < geno.num_autosome; ++i)
    {
        std::string cur_file = "chr" + std::to_string(i + 1) + "_geno";
        ASSERT_STREQ(m_genotype_file_names[i].c_str(), cur_file.c_str());
    }
    ASSERT_STREQ(m_keep_file.c_str(), "Hi");
    ASSERT_TRUE(m_remove_file.empty());
    ASSERT_TRUE(m_sample_file.empty());
    cleanup();
    // invalid format
    geno.file_name = "chr,x,y";
    try
    {
        initialize(geno, pheno, delim, type, &reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
    cleanup();
    geno.file_name = "chr#_geno,Sample";
    geno.keep = "Hi";
    geno.num_autosome = 30;
    initialize(geno, pheno, delim, type, &reporter);
    ASSERT_EQ(m_genotype_file_names.size(), geno.num_autosome);
    for (size_t i = 0; i < geno.num_autosome; ++i)
    {
        std::string cur_file = "chr" + std::to_string(i + 1) + "_geno";
        ASSERT_STREQ(m_genotype_file_names[i].c_str(), cur_file.c_str());
    }
    ASSERT_STREQ(m_keep_file.c_str(), "Hi");
    ASSERT_TRUE(m_remove_file.empty());
    ASSERT_STREQ(m_sample_file.c_str(), "Sample");
    // finally, with list
    cleanup();
    // we need a dummy file
    std::ofstream dummy("DUMMY");
    std::vector<std::string> expected = {"A", "B", "C D", "E F G", "H#I"};
    for (auto& d : expected) { dummy << d << std::endl; }
    dummy.close();
    geno.file_list = "DUMMY";
    geno.file_name = "";
    initialize(geno, pheno, delim, type, &reporter);
    ASSERT_EQ(m_genotype_file_names.size(), expected.size());
    for (size_t i = 0; i < expected.size(); ++i)
    { ASSERT_STREQ(m_genotype_file_names[i].c_str(), expected[i].c_str()); }
    ASSERT_TRUE(m_sample_file.empty());
    // now try with external sample
    cleanup();
    geno.file_list = "DUMMY,sample";
    initialize(geno, pheno, delim, type, &reporter);
    ASSERT_EQ(m_genotype_file_names.size(), expected.size());
    for (size_t i = 0; i < expected.size(); ++i)
    { ASSERT_STREQ(m_genotype_file_names[i].c_str(), expected[i].c_str()); }
    ASSERT_FALSE(m_sample_file.empty());
    ASSERT_STREQ(m_sample_file.c_str(), "sample");

    // invalid file input
    cleanup();
    geno.file_list = "DUMMY,B,C";
    try
    {
        initialize(geno, pheno, delim, type, &reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
    cleanup();
    std::remove("DUMMY");
    cleanup();
    // should fail here when the file is not found
    try
    {
        initialize(geno, pheno, delim, type, &reporter);
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}

TEST_F(GENOTYPE_BASIC, READ_BASE_PARSE_CHR)
{
    init_chr();
    // we have already check the chr convertion work
    // so we just need to check if the return is correct
    std::string normal = "Chr6 1023 rs1234 A T";
    // wrong code
    std::string wrong = "chromosome 1023 rs1234 A T";
    // too large
    std::string large = "Chr30 1023 rs1234 A T";
    // sex
    std::string sex = "chrX 1234 rs54321 G C";
    std::vector<std::string_view> token = misc::tokenize(normal);
    // if we don't have chr, we return ~size_t(0);
    BaseFile base_file;
    base_file.has_column[+BASE_INDEX::CHR] = false;
    size_t chr;
    ASSERT_EQ(parse_chr(token, base_file, +BASE_INDEX::CHR, chr), 0);
    ASSERT_EQ(chr, ~size_t(0));
    // now we have chr info
    base_file.has_column[+BASE_INDEX::CHR] = true;
    base_file.column_index[+BASE_INDEX::CHR] = 0;
    ASSERT_EQ(parse_chr(token, base_file, +BASE_INDEX::CHR, chr), 0);
    ASSERT_EQ(chr, 6);
    // check wrong code
    token = misc::tokenize(wrong);
    ASSERT_EQ(parse_chr(token, base_file, +BASE_INDEX::CHR, chr), 1);
    // check large code

    token = misc::tokenize(large);
    ASSERT_EQ(parse_chr(token, base_file, +BASE_INDEX::CHR, chr), 2);
    // sex code
    token = misc::tokenize(sex);
    ASSERT_EQ(parse_chr(token, base_file, +BASE_INDEX::CHR, chr), 2);
}

TEST_F(GENOTYPE_BASIC, READ_BASE_ALLELE)
{
    std::string line = "Chr1 1023 rs1234 a c";
    BaseFile base_file;
    std::vector<std::string_view> token = misc::tokenize(line);
    std::string ref_allele;
    parse_allele(token, base_file, +BASE_INDEX::EFFECT, ref_allele);
    ASSERT_TRUE(ref_allele.empty());
    base_file.column_index[+BASE_INDEX::EFFECT] = 3;
    base_file.has_column[+BASE_INDEX::EFFECT] = true;
    parse_allele(token, base_file, +BASE_INDEX::EFFECT, ref_allele);
    ASSERT_STREQ(ref_allele.c_str(), "A");
}

TEST_F(GENOTYPE_BASIC, READ_BASE_LOC)
{
    std::string line = "chr1 1234 rs1234 a c";
    // ok if not provided
    BaseFile base_file;
    std::vector<std::string_view> token = misc::tokenize(line);
    size_t loc;
    ASSERT_TRUE(parse_loc(token, base_file, +BASE_INDEX::BP, loc));
    ASSERT_EQ(loc, ~size_t(0));
    // we are ok as long as it is valid number
    base_file.has_column[+BASE_INDEX::BP] = true;
    base_file.column_index[+BASE_INDEX::BP] = 1;
    ASSERT_TRUE(parse_loc(token, base_file, +BASE_INDEX::BP, loc));
    ASSERT_EQ(loc, 1234);
    // we have already tested the conversion, just need 1 fail case
    line = "chr12 -123 rs12345 a c";
    token = misc::tokenize(line);
    ASSERT_FALSE(parse_loc(token, base_file, +BASE_INDEX::BP, loc));
}

TEST_F(GENOTYPE_BASIC, READ_BASE_FILTER_BY_VALUE)
{
    std::string line = "chr1 1234 rs1234 a c 0.1 0.3 0.9";
    std::vector<std::string_view> token = misc::tokenize(line);
    BaseFile base_file;
    // we will not do filtering when it is not provided
    ASSERT_FALSE(
        base_filter_by_value(token, base_file, 0.05, +BASE_INDEX::MAF));
    base_file.has_column[+BASE_INDEX::MAF] = true;
    base_file.column_index[+BASE_INDEX::MAF] = 6;
    // we will not filter this as 0.05 < 0.3
    ASSERT_FALSE(
        base_filter_by_value(token, base_file, 0.05, +BASE_INDEX::MAF));
    ASSERT_TRUE(base_filter_by_value(token, base_file, 0.4, +BASE_INDEX::MAF));
    // don't filter if it is exact
    ASSERT_FALSE(base_filter_by_value(token, base_file, 0.3, +BASE_INDEX::MAF));
}
TEST_F(GENOTYPE_BASIC, READ_BASE_P)
{
    double pvalue;
    // 0 valid 1 not convert 3 out bound 2 too big
    ASSERT_EQ(parse_pvalue("0.05", 0.67, pvalue), 0);
    ASSERT_DOUBLE_EQ(pvalue, 0.05);
    ASSERT_EQ(parse_pvalue("NA", 0.67, pvalue), 1);
    ASSERT_EQ(parse_pvalue("0.7", 0.67, pvalue), 2);
    ASSERT_EQ(parse_pvalue("2", 0.67, pvalue), 3);
    ASSERT_EQ(parse_pvalue("-1", 0.67, pvalue), 3);
}
TEST_F(GENOTYPE_BASIC, READ_BASE_STAT)
{
    double stat;
    bool OR = true;
    // 0 valid, 1 OR 0 or cannot convert, 2 negative OR
    ASSERT_EQ(parse_stat("0.1236", !OR, stat), 0);
    ASSERT_DOUBLE_EQ(0.1236, stat);
    ASSERT_EQ(parse_stat("0.1236", OR, stat), 0);
    ASSERT_DOUBLE_EQ(log(0.1236), stat);
    ASSERT_EQ(parse_stat("NA", !OR, stat), 1);
    ASSERT_EQ(parse_stat("0", OR, stat), 1);
    ASSERT_EQ(parse_stat("NA", OR, stat), 1);
    ASSERT_EQ(parse_stat("-0.654", !OR, stat), 0);
    ASSERT_DOUBLE_EQ(-0.654, stat);
    ASSERT_EQ(parse_stat("-0.654", OR, stat), 2);
}
TEST_F(GENOTYPE_BASIC, READ_BASE_RS)
{
    std::string line = "Chr1 1023 rs1234 A T";
    std::string line2 = "chr2 1234 rs5432 G T";
    std::vector<std::string_view> token = misc::tokenize(line);
    std::vector<std::string_view> token2 = misc::tokenize(line2);
    BaseFile base_file;
    base_file.has_column[+BASE_INDEX::RS] = false;
    std::unordered_set<std::string> dup_index;
    std::string rs_id;
    // first, check when rs is not provided
    try
    {
        parse_rs_id(token, dup_index, base_file, rs_id);
        FAIL();
    }
    catch (const std::runtime_error&)
    {
        SUCCEED();
    }
    // now we have the correct rs
    base_file.has_column[+BASE_INDEX::RS] = true;
    base_file.column_index[+BASE_INDEX::RS] = 2;
    ASSERT_EQ(parse_rs_id(token, dup_index, base_file, rs_id), 0);
    dup_index.insert(rs_id);
    ASSERT_STREQ(rs_id.c_str(), "rs1234");
    // duplicated SNP
    rs_id.clear();
    ASSERT_EQ(parse_rs_id(token, dup_index, base_file, rs_id), 1);
    // excluding SNP
    m_exclude_snp = true;
    dup_index.clear();
    m_snp_selection_list.insert("rs1234");
    rs_id.clear();
    ASSERT_EQ(parse_rs_id(token, dup_index, base_file, rs_id), 2);
    rs_id.clear();
    ASSERT_EQ(parse_rs_id(token2, dup_index, base_file, rs_id), 0);
    ASSERT_STREQ(rs_id.c_str(), "rs5432");
    // selecting SNP
    rs_id.clear();
    dup_index.clear();
    m_exclude_snp = false;
    ASSERT_EQ(parse_rs_id(token, dup_index, base_file, rs_id), 0);
    ASSERT_STREQ(rs_id.c_str(), "rs1234");
    rs_id.clear();
    ASSERT_EQ(parse_rs_id(token2, dup_index, base_file, rs_id), 2);
}

TEST_F(GENOTYPE_BASIC, SET_FILE_NAME_WITHOUT_HASH)
{
    std::string name = "Test";
    m_genotype_file_names = set_genotype_files(name);
    ASSERT_STREQ(m_genotype_file_names.front().c_str(), name.c_str());
    ASSERT_EQ(m_genotype_file_names.size(), 1);
}

TEST_F(GENOTYPE_BASIC, SET_FILE_NAME_WITH_HASH)
{
    std::string name = "chr#test";
    m_autosome_ct = 22;
    m_genotype_file_names = set_genotype_files(name);
    ASSERT_EQ(m_genotype_file_names.size(), m_autosome_ct);
    for (size_t i = 1; i <= m_autosome_ct; ++i)
    {
        std::string name = "chr" + std::to_string(i) + "test";
        ASSERT_STREQ(m_genotype_file_names[i - 1].c_str(), name.c_str());
    }
}

TEST_F(GENOTYPE_BASIC, SET_FILE_NAME_MULTI_HASH)
{
    // don't think this is a useful case but in theory all # should be
    // substituted
    std::string name = "chr#test#";
    m_autosome_ct = 22;
    m_genotype_file_names = set_genotype_files(name);
    ASSERT_EQ(m_genotype_file_names.size(), m_autosome_ct);
    for (size_t i = 1; i <= m_autosome_ct; ++i)
    {
        std::string name =
            "chr" + std::to_string(i) + "test" + std::to_string(i);
        ASSERT_STREQ(m_genotype_file_names[i - 1].c_str(), name.c_str());
    }
}


TEST_F(GENOTYPE_BASIC, GET_RS_COLUMN)
{
    Reporter reporter(std::string("LOG"), 60, true);
    m_reporter = &reporter;
    // simplest format
    std::string input = "rs1234";
    auto result = get_rs_column(input);
    ASSERT_EQ(result, 0);
    // might want to test the rest, RS.ID RS_ID, RSID, SNP.ID, SNP_ID,
    // Variant_ID, Variant.ID
    input = "A B CHR snp Weight P A1";
    // read from header
    result = get_rs_column(input);
    ASSERT_EQ(result, 3);
    input = "1 rs1234 123 0 A T";
    // this should be bim file
    result = get_rs_column(input);
    ASSERT_EQ(result, 1);
    input = "rs1234 A BC T A D E T";
    // no idea what this is, will take first column
    result = get_rs_column(input);
    ASSERT_EQ(result, 0);
}
TEST_F(GENOTYPE_BASIC, LOAD_SNP_LIST)
{
    Reporter reporter("LOG", 60, true);
    m_reporter = &reporter;
    std::vector<std::string> expected = {"rs1234", "rs76591", "rs139486"};
    std::ofstream dummy("DUMMY");
    for (auto& e : expected) { dummy << e << std::endl; }
    auto res = load_snp_list("DUMMY");
    ASSERT_EQ(res.size(), expected.size());
    for (auto& e : expected) { ASSERT_FALSE(res.find(e) == res.end()); }
    dummy.close();
    std::remove("DUMMY");
    // now different format, just in case the idx is correct (use bim)
    dummy.open("DUMMY");
    for (auto& e : expected) { dummy << "1 " << e << " 3 4 5 6" << std::endl; }
    dummy.close();
    res = load_snp_list("DUMMY");

    ASSERT_EQ(res.size(), expected.size());
    for (auto& e : expected) { ASSERT_FALSE(res.find(e) == res.end()); }
    std::remove("DUMMY");
    // should fail if file not found
    try
    {
        res = load_snp_list("DUMMY");
        FAIL();
    }
    catch (...)
    {
        SUCCEED();
    }
}

// init_chr
// chr_code_check
// load_snp_list
// load_ref

TEST_F(GENOTYPE_BASIC, AMBIGUOUS)
{
    ASSERT_TRUE(ambiguous("A", "A"));
    ASSERT_TRUE(ambiguous("G", "G"));
    ASSERT_TRUE(ambiguous("C", "C"));
    ASSERT_TRUE(ambiguous("T", "T"));
    ASSERT_TRUE(ambiguous("A", "T"));
    ASSERT_TRUE(ambiguous("G", "C"));
    ASSERT_TRUE(ambiguous("C", "G"));
    ASSERT_TRUE(ambiguous("T", "A"));
    ASSERT_FALSE(ambiguous("A", "G"));
    ASSERT_FALSE(ambiguous("G", "A"));
    ASSERT_FALSE(ambiguous("C", "T"));
    ASSERT_FALSE(ambiguous("T", "G"));
    ASSERT_FALSE(ambiguous("A", "C"));
    ASSERT_FALSE(ambiguous("G", "T"));
    ASSERT_FALSE(ambiguous("C", "A"));
    ASSERT_FALSE(ambiguous("T", "C"));
    ASSERT_FALSE(ambiguous("A", ""));
    ASSERT_FALSE(ambiguous("G", ""));
    ASSERT_FALSE(ambiguous("C", ""));
    ASSERT_FALSE(ambiguous("T", ""));
    ASSERT_FALSE(ambiguous("A", "TA"));
    ASSERT_FALSE(ambiguous("G", "CG"));
    ASSERT_FALSE(ambiguous("C", "GC"));
    ASSERT_FALSE(ambiguous("T", "AT"));
}

TEST_F(GENOTYPE_BASIC, CATEGORY)
{
    unsigned long long category;
    double pthres;
    // anything less than lowest threshold is consider 0
    PThresholding mock;
    mock.lower = 0.05;
    mock.inter = 0.01;
    mock.upper = 0.5;
    mock.no_full = true;
    double p_value = 0;
    category = calculate_category(p_value, pthres, mock);
    ASSERT_EQ(category, 0);
    ASSERT_DOUBLE_EQ(pthres, 0.05);
    // anything above the upper threshold is considered as 1.0 and with category
    // higher than the biggest category
    // though in theory, this should never happen as we will always filter SNPs
    // that are bigger than the last required threshold
    // the pthres is meaningless in this situation
    p_value = 0.6;
    category = calculate_category(p_value, pthres, mock);
    ASSERT_GT(category, mock.upper / mock.inter - mock.lower / mock.inter);
    // if we want full model, we will still do the same
    mock.no_full = false;
    category = calculate_category(p_value, pthres, mock);
    ASSERT_GT(category, mock.upper / mock.inter - mock.lower / mock.inter);
    ASSERT_DOUBLE_EQ(pthres, 1);
    mock.no_full = true;
    p_value = 0.05;
    category = calculate_category(p_value, pthres, mock);
    // this should be the first threshold
    ASSERT_EQ(category, 0);
    ASSERT_DOUBLE_EQ(pthres, 0.05);
    p_value = 0.055;
    category = calculate_category(p_value, pthres, mock);
    // this should be the first threshold
    ASSERT_EQ(category, 1);
    ASSERT_DOUBLE_EQ(pthres, 0.06);
}

void category_threshold(unsigned long long category,
                        unsigned long long expected_category, double pthres,
                        double expected_pthres)
{
    ASSERT_EQ(category, expected_category);
    ASSERT_DOUBLE_EQ(pthres, expected_pthres);
}
TEST_F(GENOTYPE_BASIC, BAR_LEVELS)
{
    unsigned long long category;
    double pthres;
    std::vector<double> barlevels = {0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
    // anything less than lowest threshold is consider 0
    category = calculate_category(0.0, barlevels, pthres);
    category_threshold(category, 0, pthres, 0.001);
    category = calculate_category(0.001, barlevels, pthres);
    category_threshold(category, 0, pthres, 0.001);
    category = calculate_category(0.5, barlevels, pthres);
    category_threshold(category, 6, pthres, 0.5);
    category = calculate_category(0.7, barlevels, pthres);
    category_threshold(category, 7, pthres, 1);
}


#endif // GENOTYPE_TEST_HPP

#include "IITree.h"
#include "catch.hpp"
#include "genotype.hpp"
#include "mock_genotype.hpp"
#include "plink_common.hpp"
#include "region.hpp"
#include "reporter.hpp"
#include <algorithm>
#include <limits>
#include <random>
TEST_CASE("base file read")
{
    BaseFile base_file;
    mockGenotype geno;
    geno.test_init_chr();
    std::vector<std::string_view> token;
    size_t chr;
    std::vector<size_t> filter_count;
    SECTION("parse_chr")
    {
        SECTION("did not provide chr")
        {
            base_file.has_column[+BASE_INDEX::CHR] = false;
            REQUIRE(geno.test_parse_chr(token, base_file, filter_count, chr));
            REQUIRE(chr == ~size_t(0));
        }
        SECTION("provided chr")
        {
            std::string tmp = " 19865 rs1234 A T";
            base_file.has_column[+BASE_INDEX::CHR] = true;
            base_file.column_index[+BASE_INDEX::CHR] = 0;
            SECTION("Numeric input")
            {
                auto i = GENERATE(take(50, random(1, 40)));
                std::string line = "chr" + std::to_string(i) + tmp;
                token = misc::tokenize(line);
                // 22 is the default number of autosome
                if (i > 22)
                {
                    REQUIRE_FALSE(geno.test_parse_chr(token, base_file,
                                                      filter_count, chr));
                    REQUIRE(filter_count[+FILTER_COUNT::HAPLOID] == 1);
                }
                else
                {
                    REQUIRE(geno.test_parse_chr(token, base_file, filter_count,
                                                chr));
                    REQUIRE(chr == static_cast<size_t>(i));
                }
            }
            SECTION("non-numeric")
            {
                using record = std::tuple<std::string, size_t>;
                auto in = GENERATE(table<std::string, size_t>(
                    {record {"chrX", +FILTER_COUNT::HAPLOID},
                     record {"0Y", +FILTER_COUNT::HAPLOID},
                     record {"MT", +FILTER_COUNT::HAPLOID},
                     record {"chromosome", +FILTER_COUNT::CHR},
                     record {"wrong", +FILTER_COUNT::CHR}}));
                std::string line = std::get<0>(in) + tmp;
                token = misc::tokenize(line);
                REQUIRE_FALSE(
                    geno.test_parse_chr(token, base_file, filter_count, chr));
                REQUIRE(filter_count[std::get<1>(in)] == 1);
            }
        }
    }
    SECTION("read alleles")
    {
        std::string prefix = "chr1 1234 rs1234 ";
        auto i = GENERATE("A", "C", "G", "T");
        std::string allele;
        SECTION("effective")
        {
            // must be true
            base_file.has_column[+BASE_INDEX::EFFECT] = true;
            base_file.column_index[+BASE_INDEX::EFFECT] = 3;
            prefix.append(i);
            token = misc::tokenize(prefix);
            geno.test_parse_allele(token, base_file, +BASE_INDEX::EFFECT,
                                   allele);
            REQUIRE(allele == i);
        }
        SECTION("non-effective")
        {
            auto j = GENERATE(true, false);
            base_file.has_column[+BASE_INDEX::NONEFFECT] = j;
            base_file.column_index[+BASE_INDEX::NONEFFECT] = 3;
            prefix.append(i);
            token = misc::tokenize(prefix);
            geno.test_parse_allele(token, base_file, +BASE_INDEX::NONEFFECT,
                                   allele);
            if (j) { REQUIRE(allele == i); }
            else
            {
                REQUIRE(allele == "");
            }
        }
    }
    SECTION("parse location")
    {
        std::string pref = "chr1 ";
        std::string suffix = " rs1234 a t";
        // converter tested, doesn't need to test non-numeric
        auto has_column = GENERATE(true, false);
        auto loc = GENERATE(take(100, random(-10000, 10000)));
        base_file.has_column[+BASE_INDEX::BP] = has_column;
        base_file.column_index[+BASE_INDEX::BP] = 1;
        std::string line = pref + std::to_string(loc) + suffix;
        token = misc::tokenize(line);
        size_t bp;
        geno.test_parse_loc(token, base_file, bp);
        if (!has_column) { REQUIRE(bp == ~size_t(0)); }
        else
        {
            if (loc < 0)
            { REQUIRE_FALSE(geno.test_parse_loc(token, base_file, bp)); }
            else
            {
                REQUIRE(geno.test_parse_loc(token, base_file, bp));
                REQUIRE(static_cast<size_t>(loc) == bp);
            }
        }
    }
    SECTION("read p")
    {
        SECTION("numeric")
        {
            auto i = GENERATE(take(100, random(-1.1, 1.1)));
            auto threshold = GENERATE(take(10, random(0.1, 1.0)));
            std::string input = std::to_string(i);
            double expected = misc::Convertor::convert<double>(input);
            double p_value;
            if (expected < 0 || expected > 1)
            {
                REQUIRE_THROWS(geno.test_parse_pvalue(input, threshold,
                                                      filter_count, p_value));
            }
            else if (expected > threshold)
            {
                REQUIRE_FALSE(geno.test_parse_pvalue(input, threshold,
                                                     filter_count, p_value));
                REQUIRE(filter_count[+FILTER_COUNT::P_EXCLUDED] == 1);
            }
            else
            {
                REQUIRE(geno.test_parse_pvalue(input, threshold, filter_count,
                                               p_value));
                REQUIRE(p_value == Approx(expected));
            }
        }
        SECTION("non numeric")
        {
            double p_value;
            REQUIRE_FALSE(
                geno.test_parse_pvalue("NA", 1, filter_count, p_value));
            REQUIRE(filter_count[+FILTER_COUNT::NOT_CONVERT] == 1);
        }
    }
    SECTION("read stat")
    {
        SECTION("numeric")
        {
            auto is_or = GENERATE(true, false);
            auto in = GENERATE(take(10, random(-2.0, 2.0)));
            std::string input = std::to_string(in);
            double expected = misc::convert<double>(input);
            double stat;
            if (is_or && in < 0)
            {
                REQUIRE_FALSE(
                    geno.test_parse_stat(input, is_or, filter_count, stat));
                REQUIRE(filter_count[+FILTER_COUNT::NEGATIVE] == 1);
            }
            else
            {
                REQUIRE(geno.test_parse_stat(input, is_or, filter_count, stat));
                if (is_or) { REQUIRE(stat == Approx(log(expected))); }
                else
                {
                    REQUIRE(stat == Approx(expected));
                }
            }
        }
        SECTION("non-numeric, or zero")
        {
            auto in = GENERATE("0", "NA");
            double stat;
            REQUIRE_FALSE(geno.test_parse_stat(in, true, filter_count, stat));
            REQUIRE(filter_count[+FILTER_COUNT::NOT_CONVERT] == 1);
        }
    }
    SECTION("filter by value")
    {
        using record = std::tuple<size_t, size_t>;
        auto index = GENERATE(table<size_t, size_t>(
            {record {+BASE_INDEX::MAF, +FILTER_COUNT::MAF},
             record {+BASE_INDEX::MAF_CASE, +FILTER_COUNT::MAF},
             record {+BASE_INDEX::INFO, +FILTER_COUNT::INFO}}));
        auto threshold = GENERATE(take(10, random(0.0, 1.0)));
        std::string prefix = "chr1 1234 rs1234 A C ";
        SECTION("valid input")
        {
            auto value = GENERATE(take(10, random(0.0, 1.0)));
            auto has_column = GENERATE(true, false);
            std::string input = prefix + std::to_string(value);
            token = misc::tokenize(input);
            double converted = misc::convert<double>(std::to_string(value));
            base_file.has_column[std::get<0>(index)] = has_column;
            base_file.column_index[std::get<0>(index)] = 5;
            if (has_column && converted < threshold)
            {
                REQUIRE_FALSE(geno.test_base_filter_by_value(
                    token, base_file, threshold, filter_count,
                    std::get<1>(index), std::get<0>(index)));
                REQUIRE(filter_count[std::get<1>(index)] == 1);
            }
            else
            {
                REQUIRE(geno.test_base_filter_by_value(
                    token, base_file, threshold, filter_count,
                    std::get<1>(index), std::get<0>(index)));
            }
        }
        SECTION("invalid input")
        {
            // filter out if invalid
            token = misc::tokenize(prefix + "NA");
            base_file.has_column[std::get<0>(index)] = true;
            REQUIRE_FALSE(geno.test_base_filter_by_value(
                token, base_file, threshold, filter_count, std::get<1>(index),
                std::get<0>(index)));
            REQUIRE(filter_count[std::get<1>(index)] == 1);
        }
    }
    SECTION("parse rs")
    {
        std::unordered_set<std::string> dup_index, processed_rs;
        std::string prefix = "chr1 1023 ";
        std::string suffix = " G T";
        std::string rs_id;

        SECTION("not provided")
        {
            REQUIRE_THROWS(geno.test_parse_rs_id(token, base_file, processed_rs,
                                                 dup_index, filter_count,
                                                 rs_id));
        }
        SECTION("provided input")
        {
            auto rsid = GENERATE("rs1234", "1:134");
            std::string input = prefix + rsid + suffix;
            token = misc::tokenize(input);
            base_file.has_column[+BASE_INDEX::RS] = true;
            base_file.column_index[+BASE_INDEX::RS] = 2;
            SECTION("without selection")
            {
                SECTION("without duplication")
                {
                    REQUIRE(geno.test_parse_rs_id(token, base_file,
                                                  processed_rs, dup_index,
                                                  filter_count, rs_id));
                    REQUIRE(rs_id == rsid);
                }
                SECTION("with duplication")
                {
                    processed_rs.insert("rs1234");
                    processed_rs.insert("1:134");
                    REQUIRE_FALSE(geno.test_parse_rs_id(token, base_file,
                                                        processed_rs, dup_index,
                                                        filter_count, rs_id));
                    REQUIRE(filter_count[+FILTER_COUNT::DUP_SNP] == 1);
                    REQUIRE(dup_index.find(rsid) != dup_index.end());
                }
            }
            SECTION("with selection")
            {
                auto exclude = GENERATE(true, false);
                auto select = GENERATE("rs1234", "1:134");
                geno.add_select_snp(select, exclude);
                if (exclude ^ (select == rsid))
                {
                    // keep this snp
                    REQUIRE(geno.test_parse_rs_id(token, base_file,
                                                  processed_rs, dup_index,
                                                  filter_count, rs_id));
                }
                else
                {
                    // keeping this SNP
                    REQUIRE_FALSE(geno.test_parse_rs_id(token, base_file,
                                                        processed_rs, dup_index,
                                                        filter_count, rs_id));
                    REQUIRE(filter_count[+FILTER_COUNT::SELECT] == 1);
                }
            }
        }
    }
    SECTION("calculate category")
    {
        double pthres;
        SECTION("fastscore")
        {
            std::vector<double> bar_levels = {0.001, 0.05, 0.1, 0.2,
                                              0.3,   0.4,  0.5, 1};
            auto i = GENERATE(take(10, random(0.0, 1.0)));
            unsigned long long expected = 0;
            double exp_pthres = 0.001;
            for (size_t j = 0; j < bar_levels.size(); ++j)
            {
                expected = j;
                exp_pthres = bar_levels[j];
                if (i <= bar_levels[j]) break;
            }
            REQUIRE(geno.test_cal_bar_category(i, bar_levels, pthres)
                    == expected);
            REQUIRE(pthres == Approx(exp_pthres));
        }
        SECTION("high-res")
        {
            PThresholding thres;
            thres.lower = 1e-8;
            thres.upper = 0.5;
            auto no_full = GENERATE(true, false);
            thres.no_full = no_full;
            SECTION("normal range")
            {
                auto pvalue = GENERATE(1e-10, 0.1, 0.7);
                thres.inter = 1e-5;
                if (pvalue == Approx(1e-10))
                {
                    REQUIRE(geno.test_calculate_category(thres, pvalue, pthres)
                            == 0);
                    REQUIRE(pthres == Approx(thres.lower));
                }
                else if (pvalue == Approx(0.7) && !no_full)
                {
                    REQUIRE(geno.test_calculate_category(thres, pvalue, pthres)
                            == 60000);
                    REQUIRE(pthres == Approx(1.0));
                }
                else if (pvalue > thres.upper)
                {
                    // we don't expect this
                }
                else
                {
                    REQUIRE(geno.test_calculate_category(thres, pvalue, pthres)
                            == 10000);
                    REQUIRE(pthres == Approx(0.10000001));
                }
            }
            SECTION("crazy range")
            {
                auto pvalue = GENERATE(0.1, 0.7);
                // When user ask for incredibly fine intervals
                thres.inter =
                    GENERATE(1e-100, std::numeric_limits<double>::denorm_min(),
                             std::numeric_limits<double>::min());
                REQUIRE_THROWS(
                    geno.test_calculate_category(thres, pvalue, pthres));
            }
        }
    }
    SECTION("full parse")
    {
        // just need to ensure the subroutines are working as expected
        geno.add_select_snp("exclude", true);
        base_file.is_or = true;
        QCFiltering base_qc;
        base_qc.info_score = 0.8;
        base_qc.maf = 0.05;
        base_qc.maf_case = 0.04;
        PThresholding threshold_info;
        threshold_info.no_full = true;
        threshold_info.fastscore = GENERATE(true, false);
        threshold_info.bar_levels = {0.5};
        // won't have header as read_base should have dealt with the header
        std::vector<std::string> base = {
            //"CHR BP RS A1 A2 P STAT MAF INFO MAF_CASE",
            "chr1 1234 exclude A C 0.05 1.96 0.1 0.9 0.05",
            "chr1 1234 normal A C 0.05 1.96 0.1 0.9 0.05",
            "chr1 1234 dup A C 0.05 1.96 0.1 0.9 0.05",
            "chr1 1234 dup A C 0.05 1.96 0.1 0.9 0.05",
            "chrX 1234 sex_chr A C 0.07 1.98 0.1 0.9 0.05",
            "chromosome 1234 wrong_chr A C 0.07 1.98 0.1 0.9 0.05",
            "chr1 1234 filter_control A C 0.05 1.96 0.01 0.9 0.05",
            "chr1 1234 filter_case A C 0.05 1.96 0.1 0.9 0.03",
            "chr1 1234 filter_info A C 0.05 1.96 0.1 0.02 0.05",
            "chr1 1234 p_not_convert A C NA 1.96 0.1 0.9 0.05",
            "chr1 1234 p_exclude A C 0.6 1.96 0.1 0.9 0.05",
            "chr1 1234 stat_not_convert A C 0.05 NA 0.1 0.9 0.05",
            "chr1 1234 negative_stat A C 0.05 -1.96 0.1 0.9 0.05",
            "chr6 1234 region_exclude A C 0.05 -1.96 0.1 0.9 0.05",
            "chr1 1234 ambiguous A T 0.05 1.96 0.1 0.9 0.05"};
        std::fill(base_file.has_column.begin(), base_file.has_column.end(),
                  true);
        base_file.column_index[+BASE_INDEX::CHR] = 0;
        base_file.column_index[+BASE_INDEX::BP] = 1;
        base_file.column_index[+BASE_INDEX::RS] = 2;
        base_file.column_index[+BASE_INDEX::EFFECT] = 3;
        base_file.column_index[+BASE_INDEX::NONEFFECT] = 4;
        base_file.column_index[+BASE_INDEX::P] = 5;
        base_file.column_index[+BASE_INDEX::STAT] = 6;
        base_file.column_index[+BASE_INDEX::MAF] = 7;
        base_file.column_index[+BASE_INDEX::INFO] = 8;
        base_file.column_index[+BASE_INDEX::MAF_CASE] = 9;
        base_file.column_index[+BASE_INDEX::MAX] = 9;
        std::vector<size_t> expected(+FILTER_COUNT::MAX, 1);
        expected[+FILTER_COUNT::NUM_LINE] = base.size();
        expected[+FILTER_COUNT::NOT_CONVERT] = 2;
        expected[+FILTER_COUNT::MAF] = 2;
        std::string input_str;
        for (auto&& b : base) { input_str.append(b + "\n"); }
        auto input = std::make_unique<std::istringstream>(input_str);
        std::vector<IITree<size_t, size_t>> exclusion_regions;
        Region::generate_exclusion(exclusion_regions, "chr6:1-2000");
        // set gz input = true so that we ignore the progress output
        auto [filter_count, dup_idx] = geno.test_transverse_base_file(
            base_file, base_qc, threshold_info, exclusion_regions, 10, true,
            std::move(input));

        auto check = geno.existed_snps();
        REQUIRE_THAT(filter_count, Catch::Equals<size_t>(expected));
        REQUIRE_FALSE(dup_idx.find("dup") == dup_idx.end());
        REQUIRE(dup_idx.size() == 1);
        auto idx = geno.existed_snps_idx();
        auto find = idx.find("normal");
        REQUIRE_FALSE(find == idx.end());
        auto snps = geno.existed_snps();
        REQUIRE(snps[find->second].rs() == "normal");
    }
}

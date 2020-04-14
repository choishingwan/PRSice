#include "catch.hpp"
#include "genotype.hpp"
#include "mock_genotype.hpp"
#include "region.hpp"

TEST_CASE("process snp")
{
    mockGenotype geno;
    Reporter reporter("log", 60, true);
    geno.set_reporter(&reporter);
    geno.test_init_chr();
    SECTION("check chromosome")
    {
        bool chr_error, sex_error;
        std::string prev = "chr1";
        size_t chr_num = 1;
        SECTION("valid input")
        {
            SECTION("same chromosome as before")
            {
                size_t prev_num = chr_num;
                std::string new_chr = prev;
                REQUIRE(geno.test_check_chr(new_chr, prev, chr_num, chr_error,
                                            sex_error));
                REQUIRE(prev_num == chr_num);
                REQUIRE(prev == new_chr);
            }
            SECTION("new chromosome")
            {
                std::string new_chr = "chr2";
                REQUIRE(geno.test_check_chr(new_chr, prev, chr_num, chr_error,
                                            sex_error));
                REQUIRE(chr_num == 2);
                REQUIRE(prev == new_chr);
            }
        }
        SECTION("invalid chromosomes")
        {
            auto chr = GENERATE(std::string("chrX"), "chromosome");
            auto chr_error = GENERATE(true, false);
            auto sex_error = GENERATE(true, false);
            bool ori_sex = sex_error;
            bool ori_chr = chr_error;
            REQUIRE_FALSE(
                geno.test_check_chr(chr, prev, chr_num, chr_error, sex_error));
            if (chr == "chrX")
            {
                REQUIRE(sex_error);
                REQUIRE(chr_error == ori_chr);
            }
            else
            {
                REQUIRE(chr_error);
                REQUIRE(sex_error == ori_sex);
            }
        }
    }
    SECTION("check rs")
    {
        std::unordered_set<std::string> processed_snps;
        std::unordered_set<std::string> duplicated_snps;
        std::string snp_id, rs_id, chr_id;
        SECTION("both snp and rs id are . or empty ")
        {
            rs_id = GENERATE(".", "");
            snp_id = GENERATE(".", "");
            REQUIRE_FALSE(geno.test_check_rs(
                snp_id, chr_id, rs_id, processed_snps, duplicated_snps, &geno));
            REQUIRE(geno.base_missed() == 1);
        }
        SECTION("Genotype have not loaded any SNPs")
        {
            // doesn't matter what input it is, as we haven't loaded any SNPs
            rs_id = "rs";
            snp_id = "rs123";
            REQUIRE_FALSE(geno.test_check_rs(
                snp_id, chr_id, rs_id, processed_snps, duplicated_snps, &geno));
            REQUIRE(geno.base_missed() == 1);
        }
        SECTION("Genotype have SNPs loaded")
        {
            geno.load_snp("rs1234");
            SECTION("rs id is found")
            {
                rs_id = "rs1234";
                snp_id = "rs4321";
                REQUIRE(geno.test_check_rs(snp_id, chr_id, rs_id,
                                           processed_snps, duplicated_snps,
                                           &geno));
                REQUIRE(geno.base_missed() == 0);
            }
            SECTION("snp id is found")
            {
                rs_id = "rs4321";
                snp_id = "rs1234";
                REQUIRE(geno.test_check_rs(snp_id, chr_id, rs_id,
                                           processed_snps, duplicated_snps,
                                           &geno));
                // we will update rs_id if snp_id is found instead of rs_id
                REQUIRE(rs_id == snp_id);
                REQUIRE(geno.base_missed() == 0);
            }
            SECTION("duplicated SNPs in target")
            {
                rs_id = "rs1234";
                snp_id = "rs4321";
                processed_snps.insert("rs1234");
                REQUIRE_FALSE(geno.test_check_rs(rs_id, snp_id, chr_id,
                                                 processed_snps,
                                                 duplicated_snps, &geno));
                REQUIRE(duplicated_snps.find(rs_id) != duplicated_snps.end());
                REQUIRE(geno.base_missed() == 0);
            }
        }
    }
    SECTION("check ambig handling")
    {
        auto keep_ambig = GENERATE(true, false);
        geno.keep_ambig(keep_ambig);
        bool flipped = GENERATE(true, false);
        std::string a1 = "A", a2 = "T";
        std::string ref = GENERATE(std::string("T"), "A");
        if (!keep_ambig)
        { REQUIRE_FALSE(geno.test_check_ambig(a1, a2, ref, flipped)); }
        else
        {
            REQUIRE(geno.test_check_ambig(a1, a2, ref, flipped));
            REQUIRE(flipped == (a1 != ref));
        }
        REQUIRE(geno.num_ambig() == 1);
    }
    SECTION("xregion check")
    {
        std::vector<IITree<size_t, size_t>> exclusion_regions;
        Region::generate_exclusion(exclusion_regions, "chr6:1-2000");
        size_t base_chr = GENERATE(~size_t(0), take(1, random(1ul, 22ul)));
        size_t base_loc =
            GENERATE(~size_t(0),
                     take(1, random(1ul, std::numeric_limits<size_t>::max())));
        using record = std::tuple<size_t, size_t, bool>;
        auto extend = GENERATE(table<size_t, size_t, bool>(
            {record {6, 10, true}, record {1, 10, false},
             record(6, 2001, false)}));
        SNP snp("rs", std::get<0>(extend), std::get<1>(extend), "A", "C", 1, 1);
        SNP base("rs", base_chr, base_loc, "A", "C", 1, 1);
        if (base_chr != ~size_t(0) && base_loc != ~size_t(0))
        {
            // always assume already filtered
            REQUIRE(geno.test_not_in_xregion(exclusion_regions, base, snp));
        }
        else
        {
            // we will actually do the check
            REQUIRE(geno.test_not_in_xregion(exclusion_regions, base, snp)
                    != std::get<2>(extend));
            if (std::get<2>(extend)) { REQUIRE(geno.num_xrange() == 1); }
            else
            {
                REQUIRE(geno.num_xrange() == 0);
            }
        }
    }
    SECTION("Full test")
    {
        // this will also test the region exclusion as we don't separate that
        // out
        std::vector<IITree<size_t, size_t>> exclusion_regions;
        Region::generate_exclusion(exclusion_regions, "chr6:1-2000");
        std::string name = "load.mismatch";
        std::string type = "Base";
        // when we go into process, we have already handled chr
        // this allow us to pack the SNP object as a parameter

        std::unordered_set<std::string> processed_snps;
        std::unordered_set<std::string> duplicated_snps;
        std::vector<bool> retain_snp(1, false);
        SECTION("failed at rs")
        {
            // not found
            SNP cur("not found", 1, 1, "A", "C", 1, 1);
            REQUIRE_FALSE(geno.test_process_snp(
                exclusion_regions, name, type, "", cur, processed_snps,
                duplicated_snps, retain_snp, &geno));
            REQUIRE_FALSE(retain_snp[0]);
            REQUIRE(geno.base_missed() == 1);
        }
        SECTION("ok at rs")
        {
            std::string rs = "rs1234";
            size_t chr = 1, loc = 1;
            SECTION("failed match")
            {
                SNP ref(rs, chr, loc, "A", "C", 1, 1);
                geno.load_snp(ref);
                SNP cur(rs, chr + 1, loc + 1, "A", "C", 1, 1);
                REQUIRE_FALSE(geno.test_process_snp(
                    exclusion_regions, name, type, "", cur, processed_snps,
                    duplicated_snps, retain_snp, &geno));
                REQUIRE_FALSE(retain_snp[0]);
            }

            SECTION("ok match")
            {
                SECTION("failed ambig")
                {
                    geno.keep_ambig(false);
                    SNP ref(rs, chr, loc, "A", "T", 1, 1);
                    geno.load_snp(ref);
                    REQUIRE_FALSE(geno.test_process_snp(
                        exclusion_regions, name, type, "", ref, processed_snps,
                        duplicated_snps, retain_snp, &geno));
                    REQUIRE_FALSE(retain_snp[0]);
                    REQUIRE(geno.num_ambig() == 1);
                }
                SECTION("ok ambig")
                {
                    SECTION("filtered by xregion")
                    {
                        SNP ref(rs, ~size_t(0), ~size_t(0), "A", "C", 1, 1);
                        geno.load_snp(ref);
                        SNP cur(rs, 6, 10, "A", "C", 1, 1);
                        REQUIRE_FALSE(geno.test_process_snp(
                            exclusion_regions, name, type, "", cur,
                            processed_snps, duplicated_snps, retain_snp,
                            &geno));
                        REQUIRE_FALSE(retain_snp[0]);
                        REQUIRE(geno.num_xrange() == 1);
                    }
                    SECTION("valid input")
                    {
                        SNP ref(rs, 6, 10, "A", "C", 1, 1);
                        SNP cur(rs, 6, 10, "A", "C", 10, 20);
                        geno.load_snp(ref);
                        SECTION("in target file")
                        {
                            REQUIRE(geno.test_process_snp(
                                exclusion_regions, name, type, "", cur,
                                processed_snps, duplicated_snps, retain_snp,
                                &geno));
                            REQUIRE(retain_snp[0]);
                            REQUIRE(geno.existed_snps().size() == 1);
                            auto res = geno.existed_snps().front();
                            REQUIRE(res.get_file_idx() == cur.get_file_idx());
                            REQUIRE(res.get_byte_pos() == cur.get_byte_pos());
                        }
                        SECTION("in reference")
                        {
                            mockGenotype target;
                            target.set_reporter(&reporter);
                            target.test_init_chr();
                            target.reference();

                            REQUIRE(target.test_process_snp(
                                exclusion_regions, name, type, "", cur,
                                processed_snps, duplicated_snps, retain_snp,
                                &geno));
                            REQUIRE(retain_snp[0]);
                            REQUIRE(geno.existed_snps().size() == 1);
                            REQUIRE(target.existed_snps().size() == 0);
                            auto res = geno.existed_snps().front();
                            REQUIRE(res.get_file_idx(true)
                                    == cur.get_file_idx());
                            REQUIRE(res.get_byte_pos(true)
                                    == cur.get_byte_pos());
                        }
                    }
                }
            }
        }
    }
}

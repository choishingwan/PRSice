#include "catch.hpp"
#include "snp.hpp"
TEST_CASE("Set SNP")
{
    auto rs = "rs1234";
    using geno = std::tuple<std::string, std::string>;
    auto genotypes = geno {"A", "C"};
    auto stats = GENERATE(take(1, random(-1.96, 1.96)));
    auto pvalue = GENERATE(take(1, random(0.0, 1.0)));
    double pthres = 1;
    unsigned long long category = 1;
    auto chr = GENERATE(take(1, random(1ul, 22ul)));
    auto loc = GENERATE(take(1, random(100000ul, 100000000ul)));

    SNP snp(rs, chr, loc, std::get<0>(genotypes), std::get<1>(genotypes), stats,
            pvalue, category, pthres);
    SECTION("Check default")
    {
        SECTION("rs") { REQUIRE(snp.rs() == rs); }
        SECTION("Effective allele")
        {
            REQUIRE(snp.ref() == std::get<0>(genotypes));
        }
        SECTION("Noneffective allele")
        {
            REQUIRE(snp.alt() == std::get<1>(genotypes));
        }
        SECTION("chr") { REQUIRE(snp.chr() == chr); }
        SECTION("bp") { REQUIRE(snp.loc() == loc); }
        SECTION("file_idx")
        {
            for (auto i : {true, false})
            { REQUIRE(snp.get_file_idx(i) == ~size_t(0)); }
        }
        SECTION("byte_pos")
        {
            for (auto i : {true, false}) { REQUIRE(snp.get_byte_pos(i) == 0); }
        }
        SECTION("get_counts")
        {
            size_t homcom, het, homrar, missing;
            for (auto i : {true, false})
            {
                REQUIRE_FALSE(snp.get_counts(homcom, het, homrar, missing, i));
                REQUIRE(homcom == 0);
                REQUIRE(het == 0);
                REQUIRE(homrar == 0);
                REQUIRE(missing == 0);
            }
        }
        SECTION("clumped") { REQUIRE_FALSE(snp.clumped()); }
        SECTION("flipped") { REQUIRE_FALSE(snp.is_flipped()); }
        SECTION("stat") { REQUIRE(snp.stat() == Approx(stats)); }
        SECTION("pvalue") { REQUIRE(snp.p_value() == Approx(pvalue)); }
        SECTION("threshold") { REQUIRE(snp.get_threshold() == Approx(pthres)); }
        SECTION("category") { REQUIRE(snp.category() == category); }
        SECTION("low bound") { REQUIRE(snp.low_bound() == ~size_t(0)); }
        SECTION("up bound") { REQUIRE(snp.up_bound() == ~size_t(0)); }
    }
    SECTION("add information")
    {
        auto is_ref = GENERATE(true, false);
        // default flipping is false, so we just check if it changed to true
        bool flipped = true;
        std::string ref = "G", alt = "T";
        auto file_idx = GENERATE(take(1, random(1ul, 100000ul)));
        auto byte_pos = GENERATE(
            take(1, random(1ll, std::numeric_limits<long long>::max())));
        snp.add_snp_info(file_idx, byte_pos, chr + 1, loc + 1, ref, alt,
                         flipped, is_ref);
        SECTION("updated file_idx")
        {
            REQUIRE(snp.get_file_idx(is_ref) == file_idx);
            if (is_ref) { REQUIRE(snp.get_file_idx(!is_ref) == ~size_t(0)); }
            else
            {
                // for target, we will also update the reference just in case
                REQUIRE(snp.get_file_idx(!is_ref) == file_idx);
            }
        }
        SECTION("updated byte_pos")
        {
            REQUIRE(snp.get_byte_pos(is_ref) == byte_pos);
            if (is_ref) { REQUIRE(snp.get_byte_pos(!is_ref) == 0); }
            else
            {
                // for target, we will also update the reference just in case
                REQUIRE(snp.get_byte_pos(!is_ref) == byte_pos);
            }
        }
        SECTION("updated chr")
        {
            if (!is_ref) { REQUIRE(snp.chr() == chr + 1); }
            else
            {
                REQUIRE(snp.chr() == chr);
            }
        }
        SECTION("updated loc")
        {
            if (!is_ref) { REQUIRE(snp.loc() == loc + 1); }
            else
            {
                REQUIRE(snp.loc() == loc);
            }
        }
        SECTION("updated ref")
        {
            if (!is_ref) { REQUIRE(snp.ref() == ref); }
            else
            {
                REQUIRE(snp.ref() == std::get<0>(genotypes));
            }
        }
        SECTION("updated alt")
        {
            if (!is_ref) { REQUIRE(snp.alt() == alt); }
            else
            {
                REQUIRE(snp.alt() == std::get<1>(genotypes));
            }
        }
        SECTION("flipping")
        {
            if (is_ref)
            {
                REQUIRE_FALSE(snp.is_flipped());
                REQUIRE(snp.is_ref_flipped() == flipped);
            }
            else
            {
                REQUIRE(snp.is_flipped() == flipped);
                REQUIRE_FALSE(snp.is_ref_flipped());
            }
        }
        SECTION("file info update")
        {
            // for intermediate management
            snp.update_file(file_idx + 1, byte_pos + 1, is_ref);
            REQUIRE(snp.get_file_idx(is_ref) == file_idx + 1);
            REQUIRE(snp.get_byte_pos(is_ref) == byte_pos + 1);
            if (is_ref)
            {
                // never touched target, should stay as default
                REQUIRE(snp.get_file_idx(!is_ref) == ~size_t(0));
                REQUIRE(snp.get_byte_pos(!is_ref) == 0);
            }
            else
            {
                // target is now changed, but ref shouldn't
                REQUIRE(snp.get_file_idx(!is_ref) == file_idx);
                REQUIRE(snp.get_byte_pos(!is_ref) == byte_pos);
            }
        }
    }
}

TEST_CASE("SNP Matching")
{
    using geno = std::tuple<std::string, std::string>;
    SECTION("mismatch check")
    {
        auto mref = GENERATE(std::string("A"), "C", "T", "G");
        auto ref = GENERATE(std::string("A"), "C", "T", "G");
        auto malt = GENERATE(std::string("A"), "C", "T", "G");
        auto alt = GENERATE(std::string("A"), "C", "T", "G");
        SNP snp("rs", 1, 1, mref, malt, 0, 0, 0, 0);
        bool flipped;
        if ((ref == mref && malt == alt)
            || (mref == SNP::complement(ref) && (malt == SNP::complement(alt))))
        {
            REQUIRE(snp.matching(1, 1, ref, alt, flipped));
            REQUIRE_FALSE(flipped);
        }
        else if ((alt == mref && ref == malt)
                 || (alt == SNP::complement(mref)
                     && (ref == SNP::complement(malt))))
        {
            REQUIRE(snp.matching(1, 1, ref, alt, flipped));
            REQUIRE(flipped);
        }
        else
        {
            REQUIRE_FALSE(snp.matching(1, 1, ref, alt, flipped));
        }
    }
    SECTION("with missing alternative allele")
    {

        auto chr = GENERATE(~size_t(0), take(1, random(1ul, 22ul)));
        auto loc =
            GENERATE(~size_t(0),
                     take(1, random(1ul, std::numeric_limits<size_t>::max())));
        auto missing_alt = GENERATE(table<std::string, std::string>(
            {geno {"A", ""}, geno {"C", ""}, geno {"T", ""}, geno {"G", ""}}));
        SNP snp("rs", chr, loc, std::get<0>(missing_alt),
                std::get<1>(missing_alt), 0, 0, 0, 0);
        bool flip = false;
        SECTION("same chr and loc")
        {
            auto ref = GENERATE(std::string("A"), "C", "T", "G");
            auto alt = GENERATE(std::string("A"), "C", "T", "G");
            if (ref == std::get<0>(missing_alt)
                || SNP::complement(ref) == std::get<0>(missing_alt)
                || alt == std::get<0>(missing_alt)
                || SNP::complement(alt) == std::get<0>(missing_alt))
            {
                REQUIRE(snp.matching(chr, loc, ref, alt, flip));
                if (ref == std::get<0>(missing_alt)
                    || SNP::complement(ref) == std::get<0>(missing_alt))
                { REQUIRE_FALSE(flip); }
                else
                {
                    REQUIRE(flip);
                }
            }
            else
            {
                REQUIRE_FALSE(snp.matching(chr, loc, ref, alt, flip));
            }
        }
        SECTION("different chr")
        {
            REQUIRE(snp.matching(chr + 1, loc, std::get<0>(missing_alt),
                                 std::get<1>(missing_alt), flip)
                    == (chr == ~size_t(0)));
        }
        SECTION("different loc")
        {
            REQUIRE(snp.matching(chr, loc + 1, std::get<0>(missing_alt),
                                 std::get<1>(missing_alt), flip)
                    == (loc == ~size_t(0)));
        }
    }
}

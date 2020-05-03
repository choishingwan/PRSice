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
            uint32_t homcom, het, homrar, missing;
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
        SNP src("", chr + 1, loc + 1, ref, alt, file_idx, byte_pos);
        snp.add_snp_info(src, flipped, is_ref);
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
    SECTION("with indel")
    {
        // we can't really do complement
        auto mref = GENERATE(std::string("AATT"), "TCG");
        auto malt = GENERATE(std::string("AATT"), "TCG");
        auto ref = GENERATE(std::string("AATT"), "ATCG");
        auto alt = GENERATE(std::string("TCG"), "ATCG");
        SNP snp("rs", 1, 1, mref, malt, 0, 0, 0, 0);
        bool flipped;
        // with indel, we can never do complementary mapping as that just don't
        // work
        if (ref == mref && malt == alt)
        {
            REQUIRE(snp.matching(1, 1, ref, alt, flipped));
            REQUIRE_FALSE(flipped);
        }
        else if (alt == mref && ref == malt)
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

TEST_CASE("set function return")
{
    SNP snp("rs", 1, 1, "a", "c", 1, 1);
    SECTION("change rs")
    {
        snp.rs() = "RS";
        REQUIRE(snp.rs() == "RS");
    }
    SECTION("change ref")
    {
        SECTION("subsitute")
        {
            snp.ref() = "T";
            REQUIRE(snp.ref() == "T");
        }
        SECTION("to upper")
        {
            misc::to_upper(snp.ref());
            REQUIRE(snp.ref() == "A");
        }
    }
    SECTION("change alt")
    {
        SECTION("subsitute")
        {
            snp.alt() = "G";
            REQUIRE(snp.alt() == "G");
        }
        SECTION("to upper")
        {
            misc::to_upper(snp.alt());
            REQUIRE(snp.alt() == "C");
        }
    }
}
TEST_CASE("SNP Clump")
{
    SNP index, target;
    auto num_set = GENERATE(range(127ul, 128ul));
    const auto required_size = BITCT_TO_WORDCT(num_set);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> set_idx(0, num_set - 1);
    auto&& index_flag = index.get_flag();
    auto&& target_flag = target.get_flag();
    index_flag.resize(required_size, 0);
    target_flag.resize(required_size, 0);
    std::unordered_set<size_t> index_set, target_set;
    // randomly assign 80 sets to index
    // randomly assign 80 sets to target
    for (size_t i = 0; i < 80; ++i)
    {
        auto idx = set_idx(gen);
        while (index_set.find(idx) != index_set.end()) { idx = set_idx(gen); }
        index_set.insert(idx);
        SET_BIT(idx, index_flag.data());
        idx = set_idx(gen);
        while (target_set.find(idx) != target_set.end()) { idx = set_idx(gen); }
        target_set.insert(idx);
        SET_BIT(idx, target_flag.data());
    }
    auto ori_index = index_flag;
    auto ori_target = target_flag;
    SECTION("Target was clumped")
    {
        target.set_clumped();
        index.clump(target, 1, false);
        // won't do clumping here
        REQUIRE_THAT(index.get_flag(), Catch::Equals<uintptr_t>(ori_index));
        REQUIRE_THAT(target.get_flag(), Catch::Equals<uintptr_t>(ori_target));
    }
    SECTION("Target not clumped yet")
    {
        using record = std::tuple<bool, double, double, bool>;
        auto instruction = GENERATE(table<bool, double, double, bool>(
            {record {false, 2, 0.5, false}, record {true, 0.4, 0.6, false},
             record {true, 0.4, 0.3, true}}));
        auto want_proxy_clump = std::get<0>(instruction);
        auto r2 = std::get<1>(instruction);
        auto proxy_thres = std::get<2>(instruction);
        auto used_proxy_clump = std::get<3>(instruction);
        index.clump(target, r2, want_proxy_clump, proxy_thres);
        if (!used_proxy_clump)
        {
            // when we didn't use proxy clump, we don't expect index to change
            // at all
            REQUIRE_THAT(index.get_flag(), Catch::Equals<uintptr_t>(ori_index));
            std::vector<uintptr_t> expected(required_size, 0);
            // and we expect target to lost anything found in index
            auto&& cur_target = target.get_flag();
            for (size_t i = 0; i < num_set; ++i)
            {
                if (target_set.find(i) != target_set.end())
                {
                    if (index_set.find(i) == index_set.end())
                    { SET_BIT(i, expected.data()); }
                }
            }
            REQUIRE_THAT(cur_target, Catch::Equals<uintptr_t>(expected));
        }
        else
        {
            // when proxy clumped, we simply set the target to clumped
            REQUIRE(target.clumped());
            // and then index will become the combination of both
            auto&& cur_index = index.get_flag();
            std::vector<uintptr_t> expected(required_size, 0);
            for (size_t i = 0; i < num_set; ++i)
            {
                if (index_set.find(i) != index_set.end()
                    || target_set.find(i) != target_set.end())
                { SET_BIT(i, expected.data()); }
            }
            REQUIRE_THAT(cur_index, Catch::Equals<uintptr_t>(expected));
        }
    }
}

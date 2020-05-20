#include "binarygen.hpp"
#include "catch.hpp"
#include "mock_binarygen.hpp"

TEST_CASE("failed bgen load")
{
    // Require external samples
    GenoFile geno;
    geno.file_name = "pheno_load";
    Phenotype pheno;
    Reporter reporter("log", 60, true);
    REQUIRE_THROWS(mock_binarygen(geno, pheno, " ", &reporter));
}
TEST_CASE("check phenotype header")
{
    // Check if the external sample file contain a header
    // it is technically allow as user can provide one file for external sample
    // and another for phenotype parsing
    mock_binarygen bgen;
    Reporter reporter("log", 60, true);
    bgen.set_reporter(&reporter);
    bgen.set_sample_size(1);
    SECTION("has header")
    {
        std::unique_ptr<std::istream> input =
            std::make_unique<std::istringstream>("FID IID Pheno\n"
                                                 "S1 S1 Case");
        bgen.test_handle_pheno_header(input);
        std::string line;
        std::getline((*input), line);
        REQUIRE(line == "S1 S1 Case");
    }
    SECTION("no header")
    {
        std::unique_ptr<std::istream> input =
            std::make_unique<std::istringstream>("S1 S1 Case");
        bgen.test_handle_pheno_header(input);
        std::string line;
        std::getline((*input), line);
        REQUIRE(line == "S1 S1 Case");
    }
}
TEST_CASE("sample file check")
{
    SECTION("valid formats")
    {
        auto str =
            GENERATE("ID1 ID2 missing\n"
                     "0 0 0\n"
                     "1 1 1",
                     "ID1 ID2 missing Sex\n"
                     "0 0 0 D\n"
                     "1 1 1 F",
                     "ID_1 ID_2 missing sex category binary positive disgrete\n"
                     "0 0 0 D C B P D\n"
                     "1 1 1 F 10.0 1 2 Hi");
        std::unique_ptr<std::istream> input =
            std::make_unique<std::istringstream>(str);
        REQUIRE(mock_binarygen::check_is_sample_format(input));
    }
    SECTION("invalid formats")
    {
        auto str = GENERATE(
            /*simple malform*/
            "ID1 ID2 missing\n"
            "0 0\n"
            "1 1 1",
            /*forgot missing*/
            "ID1 ID2\n"
            "0 0\n"
            "1 1",
            /* unknown format */
            "ID1 ID2 missing sex pheno\n"
            "0 0 0 D N\n"
            "S1 S1 0 F -1"
            /* sex as missing format */
            "ID1 ID2 sex missing pheno\n"
            "0 0 D 0 N\n"
            "S1 S1 0 F -1",
            /*fam file */
            "ID1 ID2 0 0 1 Pheno\n"
            "ID2 ID3 0 0 2 Pheno",
            /*one sample fam*/
            "ID1 ID2 0 0 1 Pheno"
            /*pheno file*/
            "ID1 ID2 Pheno\n"
            "ID3 ID4 Pheno");
        std::unique_ptr<std::istream> input =
            std::make_unique<std::istringstream>(str);
        REQUIRE_FALSE(mock_binarygen::check_is_sample_format(input));
    }
    SECTION("get sex col")
    {
        mock_binarygen bgen;
        Reporter reporter("log", 60, true);
        bgen.set_reporter(&reporter);
        auto ncol = GENERATE(take(2, random(size_t(5), size_t(10))));
        auto sex_col = GENERATE(range(size_t(3), size_t(4)));

        std::string line, format;
        std::vector<std::string> header(ncol, "0");
        auto colname = GENERATE("Sex", "sex", "SEX", "nosex");
        header[static_cast<size_t>(sex_col)] = colname;
        line = header[0];
        for (size_t i = 1; i < ncol; ++i) { line.append("\t" + header[i]); }
        auto code = GENERATE("D", "C");
        header[sex_col] = code;
        format = "0\t0\t0";
        for (size_t i = 3; i < ncol; ++i) { format.append("\t" + header[i]); }

        if (std::string(code) != "D" || std::string(colname) == "nosex")
        { REQUIRE(bgen.test_get_sex_col(line, format) == ~size_t(0)); }
        else
        {
            REQUIRE(bgen.test_get_sex_col(line, format) == sex_col);
        }
    }
}

TEST_CASE("full sample load")
{
    auto ignore_fid = GENERATE(true, false);
    auto is_ref = GENERATE(true, false);
    auto delim = GENERATE("@");
    Reporter reporter("log", 60, true);
    GenoFile geno;
    geno.is_ref = is_ref;
    geno.file_name = "pheno_load,bgen_pheno_load";
    Phenotype pheno;
    pheno.ignore_fid = ignore_fid;
    mock_binarygen bgen(geno, pheno, delim, &reporter);
    bgen.set_reporter(&reporter);
    bgen.gen_bgen_header("pheno_load.bgen", 10, 3, "check phenotype loading",
                         4294967295u);
    if (is_ref) bgen.reference();
    SECTION("with valid input")
    {
        using record = std::tuple<std::string, bool>;
        auto mock_file = GENERATE(
            table<std::string, bool>({record {"ID1 ID2 missing Sex\n"
                                              "0 0 0 D\n"
                                              "S1 s1 1 F\n"
                                              "remove1 remove 1 F\n"
                                              "S2 s2 1 M",
                                              true},
                                      // fam format
                                      record {"S1 s1 0 0 2 Pheno\n"
                                              "remove1 remove 0 0 2 Pheno\n"
                                              "S2 s2 0 0 2 Pheno",
                                              false},
                                      /*pheno file*/
                                      record {"S1 s1 Pheno\n"
                                              "remove1 remove 0 0 2 Pheno\n"
                                              "S2 s2 Pheno\n",
                                              false},
                                      /*pheno file with header*/
                                      record {"FID IID Pheno\n"
                                              "S1 s1 Pheno\n"
                                              "remove1 remove 0 0 2 Pheno\n"
                                              "S2 s2 Pheno",
                                              false}}));
        std::ofstream test("bgen_pheno_load");
        test << std::get<0>(mock_file) << std::endl;
        test.close();
        SECTION("without selection")
        {
            bgen.load_samples(false);
            auto res = bgen.sample_id();
            if (!is_ref)
            {
                std::string expected_iid =
                    std::get<1>(mock_file)
                        ? "remove"
                        : (ignore_fid ? "remove1" : "remove");
                REQUIRE(res.size() == 3);
                REQUIRE(res[0].FID == "S1");
                REQUIRE(
                    res[0].IID
                    == ((std::get<1>(mock_file) || !ignore_fid) ? "s1" : "S1"));
                REQUIRE(res[1].FID == "remove1");
                REQUIRE(res[1].IID == expected_iid);
                REQUIRE(res[2].FID == "S2");
                REQUIRE(
                    res[2].IID
                    == ((std::get<1>(mock_file) || !ignore_fid) ? "s2" : "S2"));
                for (auto&& i : res) { REQUIRE(i.in_regression); }
            }
            else
            {
                REQUIRE(res.empty());
            }
            auto in_prs = bgen.calculate_prs();
            auto sample_ld = bgen.sample_for_ld();
            for (size_t i = 0; i < 3; ++i)
            {
                REQUIRE(IS_SET(in_prs.data(), i));
                REQUIRE(IS_SET(sample_ld.data(), i));
            }
        }
        SECTION("with selection")
        {
            auto remove_sample = GENERATE(true, false);

            std::string expected_iid =
                std::get<1>(mock_file) ? "remove"
                                       : (ignore_fid ? "remove1" : "remove");
            auto dict_sample = !ignore_fid
                                   ? "remove1" + std::string(delim) + "remove"
                                   : expected_iid;
            bgen.add_select_sample(dict_sample);
            bgen.change_sample_selection(remove_sample);
            bgen.load_samples(false);
            auto res = bgen.sample_id();
            if (!is_ref)
            {
                if (remove_sample)
                {
                    REQUIRE(res.size() == 2);
                    REQUIRE(res[0].FID == "S1");
                    REQUIRE(res[0].IID
                            == ((std::get<1>(mock_file) || !ignore_fid)
                                    ? "s1"
                                    : "S1"));
                    REQUIRE(res[1].FID == "S2");
                    REQUIRE(res[1].IID
                            == ((std::get<1>(mock_file) || !ignore_fid)
                                    ? "s2"
                                    : "S2"));
                }
                else
                {
                    REQUIRE(res.size() == 1);
                    REQUIRE(res[0].FID == "remove1");
                    REQUIRE(res[0].IID == expected_iid);
                }
                for (auto&& i : res) { REQUIRE(i.in_regression); }
            }
            else
            {
                REQUIRE(res.empty());
            }
            auto in_prs = bgen.calculate_prs();
            auto sample_ld = bgen.sample_for_ld();
            if (remove_sample)
            {
                REQUIRE_FALSE(IS_SET(in_prs.data(), 1));
                REQUIRE_FALSE(IS_SET(sample_ld.data(), 1));
                for (auto i : {0, 2})
                {
                    REQUIRE(IS_SET(in_prs.data(), i));
                    REQUIRE(IS_SET(sample_ld.data(), i));
                }
            }
            else
            {
                REQUIRE(IS_SET(in_prs.data(), 1));
                REQUIRE(IS_SET(sample_ld.data(), 1));
                for (auto i : {0, 2})
                {
                    REQUIRE_FALSE(IS_SET(in_prs.data(), i));
                    REQUIRE_FALSE(IS_SET(sample_ld.data(), i));
                }
            }
        }
    }
    SECTION("with invalid input")
    {
        std::string mock_file = "ID1 ID2 missing Sex\n"
                                "0 0 0 D\n"
                                "S1 S1 1\n"
                                "S2 S2 1 Pheno\n";
        std::ofstream test("bgen_pheno_load");
        test << mock_file << std::endl;
        test.close();
        if (!is_ref) { REQUIRE_THROWS(bgen.load_samples(false)); }
        else
        {
            // we don't do sample check if we are not using keep / remove for LD
            // samples
            REQUIRE_NOTHROW(bgen.load_samples(false));
        }
    }
    // sample file
    // sample malform
    // not sample file
    // reference
    // reference with selection
    // reference with selection but not sample file
}


TEST_CASE("check sample consistence") {}

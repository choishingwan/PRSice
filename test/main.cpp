#include "bin_plink_test.hpp"
#include "genotype_test.hpp"
#include "order_stat.hpp"
#include "region_test.hpp"
#include "snp_test.hpp"
#include "gtest/gtest.h"

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

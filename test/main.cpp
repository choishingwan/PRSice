#include "global.hpp"
#include "gtest/gtest.h"
std::string path;
#include <iostream>
int main(int argc, char* argv[])
{
    //assert(argc >= 2);
    //path = argv[1];
    // just in case
    path.append("/");
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

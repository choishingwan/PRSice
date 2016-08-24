#include <iostream>
#include <cstdlib>
#include <string>
#include <stdexcept>

#include "../inc/commander.hpp"
#include "../inc/prsice.hpp"

int main(int argc, char *argv[])
{
    Commander commander = Commander();
    try{
        if(!commander.initialize(argc, argv)) return 0; //only require the usage information
    }
    catch (const std::runtime_error& error){
        std::cerr << error.what() << std::endl;
    }
    PRSice prsice = PRSice();
    try{
        prsice.run(commander);
    }
    catch(const std::runtime_error& error){
        std::cerr << error.what() << std::endl;
    }
    
//    boost::ptr_vector<SNP> snp_list;
//    SNP::read_snp(commander, snp_list);
//    PRSice prsice = PRSice();
//    prsice.run();
    return 0;
}

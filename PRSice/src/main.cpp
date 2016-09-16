#include <iostream>
#include <cstdlib>
#include <string>
#include <stdexcept>

#include "commander.hpp"
#include "prsice.hpp"
#include "region.hpp"

int main(int argc, char *argv[])
{
    Commander commander = Commander();
    try{
        if(!commander.initialize(argc, argv)) return 0; //only require the usage information
    }
    catch (const std::runtime_error& error){
        std::cerr << error.what() << std::endl;
		exit(-1);
    }
    Region region = Region();
    try{
    		region.run(commander.get_gtf(), commander.get_msigdb(), commander.get_bed(), commander.get_out(), commander.gen_bed());
    }
    catch(const std::runtime_error &error){
    		std::cerr << error.what() << std::endl;
    		exit(-1);
    }
    PRSice prsice = PRSice();
    try{
        prsice.run(commander, region);
    }
    catch(const std::runtime_error& error){
        std::cerr << error.what() << std::endl;
		exit(-1);
    }

    return 0;
}

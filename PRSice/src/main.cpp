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
    try
    {
        if(!commander.initialize(argc, argv)) return 0; //only require the usage information
    }
    catch (const std::runtime_error& error)
    {
        std::cerr << error.what() << std::endl;
        exit(-1);
    }
    Region region = Region();
    try
    {
        region.run(commander.get_gtf(), commander.get_msigdb(), commander.get_bed(), commander.get_out(), commander.gen_bed());
    }
    catch(const std::runtime_error &error)
    {
        std::cerr << error.what() << std::endl;
        exit(-1);
    }


    // User input should be shown before other stuff to reduce the redundency
    fprintf(stderr,"\nRegion Information\n");
    fprintf(stderr,"==============================\n");
    if(region.size() ==1) fprintf(stderr, "1 region is included\n");
    else if(region.size()>1) fprintf(stderr, "A total of %zu regions are included\n", region.size());
    fprintf(stderr,"\nUser Defined Column Headers\n");
    fprintf(stderr,"==============================\n");
    if(!commander.chr().empty()) fprintf(stderr,"Chr            : %s\n", commander.chr().c_str());
    fprintf(stderr,"SNP            : %s\n", commander.snp().c_str());
    if(!commander.bp().empty()) fprintf(stderr,"BP             : %s\n", commander.bp().c_str());
    fprintf(stderr,"Ref Allele     : %s\n", commander.ref().c_str());
    if(!commander.alt().empty()) fprintf(stderr,"Alt Allele     : %s\n", commander.alt().c_str());
    if(!commander.statistic().empty()) fprintf(stderr,"Statistic      : %s\n", commander.statistic().c_str());
    if(!commander.se().empty()) fprintf(stderr,"Standard Error : %s\n", commander.se().c_str());
    fprintf(stderr,"P-value        : %s\n", commander.p().c_str());
    fprintf(stderr,"\nClumping Parameters: \n");
    fprintf(stderr,"==============================\n");
    fprintf(stderr,"P-Threshold  : %f\n", commander.get_clump_p());
    fprintf(stderr,"R2-Threshold : %f\n", commander.get_clump_r2());
    fprintf(stderr,"Window Size  : %zu\n", commander.get_clump_kb());
    fprintf(stderr,"\nStart processing: %s\n", commander.get_target().c_str());
    fprintf(stderr,"==============================\n");
    bool fastscore = commander.fastscore();
    bool full_model = commander.full();
    double bound_start =  (fastscore)? commander.get_bar_lower(): commander.get_lower();
    double bound_end = (fastscore)? commander.get_bar_upper():commander.get_upper();
    double bound_inter = commander.get_inter();

    std::vector<std::string> base = commander.get_base();
    int num_base = base.size();
    if(num_base == 0) throw std::runtime_error("There is no base case to run");
    if(num_base < 0) throw std::runtime_error("Negative number of base");
    else
    {
        if(num_base > 1) fprintf(stderr, "Multiple base phenotype detected. You might want to run separate instance of PRSice to speed up the process\n");
        for(size_t i_base = 0; i_base < num_base; ++i_base)
        {
            region.reset();
            fprintf(stderr,"\nStart processing: %s\n", base[i_base].c_str());
            fprintf(stderr,"==============================\n");
            //        	Need to handle paths in the name
            std::string base_name=misc::remove_extension<std::string>(misc::base_name<std::string>(base[i_base]));
            try
            {
                PRSice prsice = PRSice(base_name, i_base);
                double threshold = (full_model)? 1.0:bound_end;
                // no point in getting SNPs that we will never use in the model
                prsice.get_snp(commander, region, threshold);
                prsice.clump(commander);
                if(commander.prslice() > 0.0)
                {
                    if(region.size() > 1)
                    {
                        fprintf(stderr, "\nWARNING: Currently, PRSlice cannot be performed together with PRSet. Please run them separately\n");
                    }
                }
                prsice.run_prs(commander, region);
            }
            catch(const std::out_of_range &error)
            {
                std::cerr << error.what() << std::endl;
                exit(-1);
            }
            catch(const std::runtime_error &error)
            {
                std::cerr << error.what() << std::endl;
                exit(-1);
            }
            fprintf(stderr, "\n");
        }
    }
    return 0;
}

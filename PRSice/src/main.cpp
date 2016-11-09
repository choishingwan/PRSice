#include <iostream>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <utility>

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
    // Output the region summary
    std::vector<std::pair<std::string, double> > region_info= region.get_info();
    std::ofstream region_out;
    std::string region_out_name = commander.get_out()+".region";
    region_out.open(region_out_name.c_str());
    if(!region_out.is_open())
    {
    		fprintf(stderr, "Cannot open region information file to write: %s\n", region_out_name.c_str());
    		return -1;
    }
    region_out << "Region\t%Gene info" << std::endl;
    for(auto i : region_info)
    {
    		region_out << std::get<0>(i) << "\t" << std::get<1>(i) << std::endl;
    }
    region_out.close();

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



    bool perform_prslice = commander.prslice() > 0.0;
    bool full_model = commander.full();
    double bound_end = commander.get_upper();
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
            		/**
            		 * Initialize the PRSice object
            		 */
                PRSice prsice = PRSice(base_name, i_base, commander.get_target(), commander.target_is_binary());
                double threshold = (full_model)? 1.0:bound_end;
                /**
                 * Read in SNPs from the base file. We will only include SNPs less than threhsold as
                 * they will be ignored in the whole process anyway
                 */
                prsice.get_snp(commander, region, threshold);
                /**
                 * Perform clumping on the SNPs. This help us to get around the problem of LD
                 */
                prsice.clump(commander);
                /**
                 * Initialize the phenotype information for the target
                 * We can actually perform this outside the loop and set the
                 * phenotype information as static variable. But then this piece
                 * of code should be quick. So might be safer and easier to just
                 * keep it here
                 */
                prsice.init_pheno(commander);
                size_t num_pheno = prsice.num_phenotype();
                if(!perform_prslice)
                {
                		/**
                		 * For PRSet / PRSice, we will perform this set of actions
                		 * Categorize SNPs based on their p-value. This will aid the
                		 * Iterative process of PRSice
                		 * As this only depends on the base file, this should be the
                		 * same for all phenotypes
                		 */
                		prsice.categorize(commander);
                		for(size_t i_pheno=0; i_pheno < num_pheno; ++i_pheno)
                		{
                			/**
                			 * Initialize the matrix. The reason why we do it for each phenotype
                			 * is because of missing data. It is much easier to make it for each
                			 * phenotype than making a full matrix and modifying it
                			 */
                			prsice.init_matrix(commander, i_pheno, perform_prslice);
                			try
                			{
                    			/**
                    			 * Start performing the actual PRSice. The results will all be stored
                    			 * within the class vectors. This help us to reduce the number of
                    			 * parameters required
                    			 */
                				prsice.prsice(commander, region, i_pheno);
                    			fprintf(stderr, "\n");
                    			/**
                    			 * Output the results
                    			 */
                    			prsice.output(commander, region, i_pheno);
                			}
                			catch(const std::runtime_error &error)
                			{
                				fprintf(stderr, "None of the SNPs fall within the threshold\n");
                			}
                		}
                }
                else
                {
                		if(region.size()!=1)
                		{
                			/**
                			 * It doesn't make much sense for PRSet + PRSlice as PRSlice only
                			 * consider the genomic coordinates, which for most of the time,
                			 * it will not fall within the region defined by PRSet anyway.
                			 * It also complicate the algorithm and might take forever to
                			 * run
                			 */
                			std::string error_message= "WARNING: It doesn't make sense to run PRSlice together with "
                					"PRSset. Will only perfrom PRSlice";
                			fprintf(stderr, "%s\n", error_message.c_str());
                		}
                		for(size_t i_pheno=0; i_pheno < num_pheno; ++i_pheno)
                		{
                			/**
                			 * Again, initialize the matrix, which will be used for the whole PRSlice
                			 */
                			prsice.init_matrix(commander, i_pheno, perform_prslice);
                			/**
                			 * Perform PRSice on each window
                			 * region here is only a place holder required by some of
                			 * the functions from PRSice
                			 */
                			prsice.prslice_windows(commander, region);
                			/**
                			 * Now calculate the best window combination
                			 */
                			prsice.prslice(commander, region, i_pheno);
                			/**
                			 * This should produce the output
                			 */
                			prsice.output(commander, i_pheno);
                		}
                }


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

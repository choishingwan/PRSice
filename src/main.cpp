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
		if (!commander.initialize(argc, argv)) return 0; //only require the usage information
	}
	catch (const std::runtime_error& error)
	{
		std::cerr << error.what() << std::endl;
		exit(-1);
	}

	Region region = Region(commander.get_feature());
	try
	{
		region.run(commander.get_gtf(), commander.get_msigdb(),
				commander.get_bed(), commander.get_out());
	}
	catch (const std::runtime_error &error)
	{
		std::cerr << error.what() << std::endl;
		exit(-1);
	}

	std::vector < std::string > base = commander.get_base();
	// Might want to generate a log file?
	region.info();
	commander.user_input();


	bool perform_prslice = commander.prslice() > 0.0;
	bool full_model = commander.full();
	double bound_end = commander.get_upper();
	int num_base = base.size();
	if (num_base == 0) throw std::runtime_error("There is no base case to run");
	if (num_base < 0) throw std::runtime_error("Negative number of base");
	else
	{
		if (num_base > 1)
			fprintf(stderr, "Multiple base phenotype detected. You might want to run separate instance of PRSice to speed up the process\n");
		for (size_t i_base = 0; i_base < num_base; ++i_base)
		{
			region.reset();
			fprintf(stderr, "\nStart processing: %s\n", base[i_base].c_str());
			fprintf(stderr, "==============================\n");
			//        	Need to handle paths in the name
			std::string base_name = misc::remove_extension<std::string>(
					misc::base_name<std::string>(base[i_base]));
			try
			{
				PRSice prsice = PRSice(base_name, i_base,
						commander.get_target(), commander.target_is_binary(),
						commander.get_perm(),
						commander.get_scoring(),
						region.size());
				prsice.get_snp(commander, region);

				std::string region_out_name = commander.get_out() + "." + base_name + ".region";
				region.print_file(region_out_name);
				prsice.perform_clump(commander);
				prsice.pheno_check(commander);
				size_t num_pheno = prsice.num_phenotype();

				if (!perform_prslice) {
					prsice.categorize(commander);
					for (size_t i_pheno = 0; i_pheno < num_pheno; ++i_pheno) {
						prsice.init_matrix(commander, i_pheno, perform_prslice);
						try {
							prsice.prsice(commander, region, i_pheno);

							fprintf(stderr, "\n");
							prsice.output(commander, region, i_pheno);
						} catch (const std::runtime_error &error) {
							fprintf(stderr,
									"None of the SNPs fall within the threshold\n");
						}
					}
				} else {
					if (region.size() != 1) {
						/**
						 * It doesn't make much sense for PRSet + PRSlice as PRSlice only
						 * consider the genomic coordinates, which for most of the time,
						 * it will not fall within the region defined by PRSet anyway.
						 * It also complicate the algorithm and might take forever to
						 * run
						 */
						std::string error_message =
								"WARNING: It doesn't make sense to run PRSlice together with "
										"PRSset. Will only perfrom PRSlice";
						fprintf(stderr, "%s\n", error_message.c_str());
					}
					// clean up the region such that it is easier to handle later on
					region.prslice();
					for (size_t i_pheno = 0; i_pheno < num_pheno; ++i_pheno) {
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
			} catch (const std::out_of_range &error) {
				std::cerr << error.what() << std::endl;
				exit(-1);
			} catch (const std::runtime_error &error) {
				std::cerr << error.what() << std::endl;
				exit(-1);
			}
			fprintf(stderr, "\n");
		}
	}
	return 0;
}

#include <iostream>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <utility>

#include "commander.hpp"
#include "prsice.hpp"
#include "region.hpp"
#include "genotype.hpp"

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

	GenomeFactory factory;
	// change the factory according to the file type
	// to get the file type, we might want to revemp the commander class
	// such that we can have a more elegant handling of the files.
	std::unique_ptr<Genotype> target_file = factory.createGenotype(commander, commander.target_name(), commander.target_type(), true);
	std::unique_ptr<Genotype> ld_file;
	if(!commander.ld_prefix().empty() && commander.ld_prefix().compare(commander.target_name())!=0){
		ld_file =  factory.createGenotype(commander, commander.ld_prefix(), commander.ld_type(), true);
	}

	if(ld_file!=nullptr)
	{
		double matched = target_file->update_existed(*ld_file);
		if(matched==-1)
		{
			fprintf(stderr, "ERROR: None of the SNP matched between the target and LD file\n");
			return -1;
		}
		else if(matched !=1)
		{
			fprintf(stderr, "WARNING: %03.2f%% mismatched SNPs between target and LD reference\n", matched*100);
		}
		ld_file->update_existed(*target_file);
	}


	Region region = Region(commander.feature(), target_file->get_chr_order());
	try
	{
		region.run(commander.gtf(), commander.msigdb(),
				commander.bed(), commander.out());
	}
	catch (const std::runtime_error &error)
	{
		std::cerr << error.what() << std::endl;
		exit(-1);
	}


	std::string base = commander.base_name();
	// Might want to generate a log file?
	region.info();
	commander.user_input();

	bool perform_prslice = commander.perform_prslice();
	fprintf(stderr, "\nStart processing: %s\n", commander.base_name().c_str());
	fprintf(stderr, "==============================\n");
	//        	Need to handle paths in the name
	std::string base_name = misc::remove_extension<std::string>(
			misc::base_name<std::string>(commander.base_name()));
	try
	{
		target_file->read_base(commander, region);
		if(commander.filter_mind())
		{
			target_file->filter_mind(commander.mind());
		}
		//target_file->clump();
		/*
		PRSice prsice = PRSice(base_name, i_base,
				commander.get_target(), commander.target_is_binary(),
				commander.get_perm(),
				commander.get_scoring(),
				region.size(),
				commander.ignore_fid());
		prsice.get_snp(commander, region);
		std::string region_out_name = commander.get_out() + "." + base_name + ".region";
		region.print_file(region_out_name);


		prsice.perform_clump(commander);

		prsice.pheno_check(commander);
		size_t num_pheno = prsice.num_phenotype();

		if (!perform_prslice) {
			prsice.categorize(commander);
			for (size_t i_pheno = 0; i_pheno < num_pheno; ++i_pheno) {
				if(num_pheno==0)
				{
					fprintf(stderr, "\nPRSice Analysis\n");
					fprintf(stderr, "==============================\n");
				}
				prsice.init_matrix(commander, i_pheno, perform_prslice);
				try {
					prsice.prsice(commander, region, i_pheno);
					fprintf(stderr, "\n");
					prsice.output(commander, region, i_pheno);
				} catch (const std::runtime_error &error) {
					std::cerr << "Error is: " << error.what() << std::endl;
					fprintf(stderr,
							"None of the SNPs fall within the threshold\n");
				}
			}

		} else {
			// clean up the region such that it is easier to handle later on
			region.prslice();
			for (size_t i_pheno = 0; i_pheno < num_pheno; ++i_pheno) {
				prsice.init_matrix(commander, i_pheno, perform_prslice);
				prsice.prslice_windows(commander, region);
				prsice.prslice(commander, region, i_pheno);
				prsice.output(commander, i_pheno);
			}
		}
		*/
	} catch (const std::out_of_range &error) {
		std::cerr << error.what() << std::endl;
		exit(-1);
	} catch (const std::runtime_error &error) {
		std::cerr << error.what() << std::endl;
		exit(-1);
	}
	fprintf(stderr, "\n");
	}
	return 0;
}

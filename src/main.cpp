// This file is part of PRSice2.0, copyright (C) 2016-2017
// Shing Wan Choi, Jack Euesden, Cathryn M. Lewis, Paul F. O’Reilly
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include "commander.hpp"
#include "genotype.hpp"
#include "genotypefactory.hpp"
#include "prsice.hpp"
#include "region.hpp"

int main(int argc, char* argv[])
{
    Commander commander = Commander();
    try
    {
        if (!commander.init(argc, argv))
            return 0; // only require the usage information
    }
    catch (const std::runtime_error& error)
    {
        std::cerr << error.what() << std::endl;
        exit(-1);
    }

    bool verbose = true;
    // this allow us to generate the appropriate object (i.e. binaryplink /
    // binarygen)
    GenomeFactory factory;
    Genotype* target_file;
    try
    {
        target_file = factory.createGenotype(commander, commander.target_name(),
                                             commander.target_type(), verbose);
    }
    catch (const std::invalid_argument& ia)
    {
        std::cerr << ia.what() << std::endl;
    }
    // calculate the maf and genotype missingness here? This will give us the
    // hh_exist information required for processing sex chromosomes
    bool used_ld = false;
    Genotype* ld_file = nullptr;
    if (!commander.ld_prefix().empty()
        && commander.ld_prefix().compare(commander.target_name()) != 0)
    {
        used_ld = true;
        ld_file = factory.createGenotype(commander, commander.ld_prefix(),
                                         commander.ld_type(), verbose);
    }

    Region region = Region(commander.feature(), target_file->get_chr_order());
    try
    {
        region.run(commander.gtf(), commander.msigdb(), commander.bed(),
                   commander.out());
    }
    catch (const std::runtime_error& error)
    {
        std::cerr << error.what() << std::endl;
        exit(-1);
    }

    // Might want to generate a log file?
    region.info();

    bool perform_prslice = commander.perform_prslice();

    // Need to handle paths in the name
    std::string base_name = misc::remove_extension<std::string>(
        misc::base_name<std::string>(commander.base_name()));
    fprintf(stderr, "\nStart processing: %s\n", base_name.c_str());
    fprintf(stderr, "==============================\n");
    try
    {
        target_file->read_base(commander, region);
        std::string region_out_name = commander.out() + ".region";
        region.print_file(region_out_name);

        if (!commander.no_clump()) {
            // we will perform clumping on all samples that are included
            // ignoring if they have valid phenotype/covariates or not
            target_file->clump((ld_file == nullptr) ? *target_file : *ld_file);
        }
        PRSice prsice = PRSice(base_name, commander.target_name(),
                               commander.is_binary(), commander.get_scoring(),
                               region.size(), commander.ignore_fid());
        prsice.pheno_check(commander);
        size_t num_pheno = prsice.num_phenotype();
        if (!perform_prslice) {
            fprintf(stderr, "\nPRSice Analysis\n");
            fprintf(stderr, "==============================\n");
            if (!target_file->prepare_prsice()) {
                // check if we can successfully sort the SNP vector by the
                // category as required by PRSice
                return -1;
            }
            for (size_t i_pheno = 0; i_pheno < num_pheno; ++i_pheno) {
                prsice.init_matrix(commander, i_pheno, *target_file, false);
                prsice.prsice(commander, region.names(), i_pheno, *target_file);
                if (!commander.no_regress())
                    prsice.output(commander, region, i_pheno, *target_file);
            }
        }
        /*

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
    }
    catch (const std::out_of_range& error)
    {
        std::cerr << error.what() << std::endl;
        exit(-1);
    }
    catch (const std::runtime_error& error)
    {
        std::cerr << error.what() << std::endl;
        exit(-1);
    }
    fprintf(stderr, "\n");
    delete target_file;
    if (used_ld) delete ld_file;
    return 0;
}

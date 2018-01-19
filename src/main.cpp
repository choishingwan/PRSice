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

#include "plink_common.hpp"
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
#include "reporter.hpp"
int main(int argc, char* argv[])
{
    Reporter reporter;
    Commander commander;
    try
    {
        if (!commander.init(argc, argv, reporter))
            return 0; // only require the usage information
    }
    catch (const std::runtime_error& error)
    {
        return -1; // all error messages should have printed
    }
    bool verbose = true;
    // this allow us to generate the appropriate object (i.e. binaryplink /
    // binarygen)
    GenomeFactory factory;
    Genotype *target_file, *reference_file;
    try
    {
        target_file = factory.createGenotype(
            commander.target_name(), commander.target_type(),
            commander.thread(), commander.ignore_fid(), commander.nonfounders(),
            commander.keep_ambig(), reporter, commander);
        target_file->load_samples(commander.keep_sample_file(),
                                  commander.remove_sample_file(), verbose,
                                  reporter);
        target_file->load_snps(commander.out(), commander.extract_file(),
                               commander.exclude_file(), commander.geno(),
                               commander.maf(), commander.info(),
                               commander.hard_threshold(),
                               commander.hard_coded(), verbose, reporter);
    }
    catch (const std::invalid_argument& ia)
    {
        reporter.report(ia.what());
        return -1;
    }
    catch (const std::runtime_error& error)
    {
        reporter.report(error.what());
        return -1;
    }
    // TODO: Revamp Region to make it suitable for prslice too
    Region region = Region(commander.feature(), target_file->get_chr_order());
    try
    {
        region.run(commander.gtf(), commander.msigdb(), commander.bed(),
                   commander.out());
    }
    catch (const std::runtime_error& error)
    {
        reporter.report(error.what());
        return -1;
    }

    // Might want to generate a log file?
    region.info(reporter);

    bool perform_prslice = commander.perform_prslice();

    // Need to handle paths in the name
    std::string base_name = misc::remove_extension<std::string>(
        misc::base_name<std::string>(commander.base_name()));
    std::string message = "Start processing " + base_name + "\n";
    message.append("==============================\n");
    reporter.report(message);
    try
    {
        target_file->set_info(commander);
        target_file->read_base(commander, region, reporter);

        if (!commander.ref_name().empty()) {
            reporter.report("Loading reference panel\n");
            reference_file = factory.createGenotype(
                commander.ref_name(), commander.ref_type(), commander.thread(),
                commander.ignore_fid(), commander.nonfounders(),
                commander.keep_ambig(), reporter, commander);
            reference_file->load_samples(commander.ld_keep_file(),
                                         commander.ld_remove_file(), true,
                                         reporter);
            // only load SNPs that can be found in the target file index
            reference_file->load_snps(
                commander.out(), target_file->index(), commander.ld_geno(),
                commander.ld_maf(), commander.ld_info(), commander.ld_hard_threshold(),
                commander.hard_coded(), true, reporter);
        }
        // get the sort by p inex vector for target
        // so that we can still find out the relative coordinates of each SNPs
        if (!target_file->sort_by_p())
        {
            std::string error_message = "No SNPs left for PRSice processing";
            reporter.report(error_message);
            return -1;
        }
        // we no longer need the region boundaries
        // as we don't allow multiple base file input
        region.clean();
        std::string region_out_name = commander.out() + ".region";
        // output the number of SNPs observed in each sets
        region.print_file(region_out_name);
        // perform clumping (Main problem with memory here)
        if (!commander.no_clump()) {
            target_file->efficient_clumping(
                (commander.ref_name().empty()) ? *target_file : *reference_file,
                reporter);
            // immediately free the memory if needed
            if (!commander.ref_name().empty()) delete reference_file;
        }
        if (!target_file->prepare_prsice()) {
            std::string error_message = "No SNPs left for PRSice processing";
            reporter.report(error_message);
            return -1;
        }
        // initialize PRSice class
        PRSice prsice = PRSice(base_name, commander, region.size() > 1,
                               target_file->num_sample(), reporter);
        // check the phenotype input columns
        prsice.pheno_check(commander, reporter);
        size_t num_pheno = prsice.num_phenotype();
        if (!perform_prslice) {
            for (size_t i_pheno = 0; i_pheno < num_pheno; ++i_pheno) {
                // initialize the phenotype & independent variable matrix
                prsice.init_matrix(commander, i_pheno, *target_file, reporter);
                // go through each region separately
                // this should reduce the memory usage
                if (region.size() > 1) {
                    fprintf(stderr, "\rProcessing %03.2f%% of sets", 0.0);
                }
                for (size_t i_region = 0; i_region < region.size(); ++i_region)
                {
                    prsice.run_prsice(commander, region.get_name(i_region),
                                      i_pheno, i_region, *target_file);
                    if (region.size() > 1) {
                        fprintf(stderr, "\rProcessing %03.2f%% of sets",
                                (double) i_region / (double) region.size()
                                    * 100.0);
                    }
                    if (!commander.no_regress())
                        prsice.output(commander, region, i_pheno, i_region,
                                      *target_file);
                }
                if (region.size() > 1) {
                    fprintf(stderr, "\rProcessing %03.2f%% of sets\n", 100.0);
                }
            }
            prsice.summarize(commander, reporter);
        }
        else
        {
            std::string error_message =
                "ERROR: We currently have not implemented PRSlice. We will "
                "implement PRSlice once the implementation of PRSice is "
                "stabalized";
            reporter.report(error_message);
            return -1;
        }
    }
    catch (const std::out_of_range& error)
    {
        reporter.report(error.what());
        return -1;
    }
    catch (const std::runtime_error& error)
    {
        reporter.report(error.what());
        return -1;
    }
    delete target_file;
    return 0;
}

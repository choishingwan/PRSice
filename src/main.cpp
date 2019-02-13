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
    // initialize reporter, use to generate log
    Reporter reporter;
    try
    {
        // initialize commander, use to parse command line arguments
        // initialization should help setting the default values
        Commander commander;
        try
        {
            if (!commander.init(argc, argv, reporter)) {
                return 0; // only require the usage information
            }
        }
        catch (...)
        {
            return -1; // all error messages should have printed
        }
        bool verbose = true;
        // parse the exclusion range and put it into the exclusion object
        Region exclusion_region(commander.exclusion_range(), reporter);

        // this allow us to generate the appropriate object (i.e. binaryplink /
        // binarygen)
        double maf, geno, info, hard_threshold;
        bool maf_filter, geno_filter, hard_coded, info_filter;
        // load the filtering parameters for the target file
        maf_filter = commander.target_maf(maf);
        geno_filter = commander.target_geno(geno);
        info_filter = commander.target_info(info);
        commander.target_hard_threshold(hard_threshold);
        hard_coded = commander.hard_coded();

        GenomeFactory factory;
        Genotype *target_file = nullptr, *reference_file = nullptr;
        try
        {
            // initialize the target object using the factory
            target_file = factory.createGenotype(commander, reporter);
            // then we will read in the sample information
            target_file->load_samples(commander.keep_sample_file(),
                                      commander.remove_sample_file(), verbose,
                                      reporter);
            // For bgen or any format that require intermediate generation, we
            // need to know if a reference file is going to be used, therefore
            // decide if we are going to generate the target intermediate or the
            // reference intermediate
            if (commander.use_ref()) target_file->expect_reference();
            // Finally, we can read in the SNP information
            target_file->load_snps(commander.out(), commander.exclude_file(),
                                   commander.extract_file(), maf, geno, info,
                                   hard_threshold, maf_filter, geno_filter,
                                   info_filter, hard_coded, exclusion_region,
                                   verbose, reporter);
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
        // Initialized the region object. This object is responsible for
        // checking if a SNP falls within a genic region
        Region region(commander.feature(), commander.window_5(),
                      commander.window_3(), commander.perform_set_perm(),
                      commander.genome_wide_background());
        try
        {
            // read in all the region inputs and generate the cooresponding
            // boundaries
            region.generate_regions(
                commander.gtf(), commander.msigdb(), commander.bed(),
                commander.single_snp_set(), commander.multi_snp_sets(),
                commander.background(), *target_file, reporter);
        }
        catch (const std::runtime_error& error)
        {
            reporter.report(error.what());
            return -1;
        }

        // Print out the log about the number of region included
        region.print_region_number(reporter);

        // for now, this should always be false
        const bool perform_prslice = commander.perform_prslice();

        // Need to handle paths in the name
        const std::string base_name = misc::remove_extension<std::string>(
            misc::base_name<std::string>(commander.base_name()));
        try
        {
            // initialize PRSice class
            PRSice prsice(commander, region.size() > 1, reporter);
            // check the phenotype input columns
            prsice.pheno_check(commander.pheno_file(), commander.pheno_col(),
                               commander.is_binary(), reporter);
            // set the clumping information to the target file
            target_file->set_info(commander);
            if ((!commander.no_clump() || commander.perform_shrinkage())
                && commander.use_ref())
            {
                // load the reference file if we require it
                reporter.report("Loading reference "
                                "panel\n==============================\n");
                reference_file =
                    factory.createGenotype(commander, reporter, true);

                reference_file->load_samples(commander.ref_keep_file(),
                                             commander.ref_remove_file(),
                                             verbose, reporter);
                // load the filtering for reference panel
                maf_filter = commander.ref_maf(maf);
                geno_filter = commander.ref_geno(geno);
                info_filter = commander.ref_info(info);
                hard_coded = commander.ref_hard_threshold(hard_threshold);
                // only load SNPs that can be found in the target file index
                // we don't need to load the extraction and exclude files as
                // they should already be in effect on the target SNPs
                reference_file->load_snps(
                    commander.out(), "", "", maf, geno, info, hard_threshold,
                    maf_filter, geno_filter, info_filter, hard_coded,
                    exclusion_region, verbose, reporter, target_file);
            }
            std::string message = "Start processing " + base_name + "\n";
            message.append("==============================\n");
            reporter.report(message);
            // now we read in the base file
            target_file->read_base(
                commander.base_name(), commander.out(), commander.index(),
                commander.bar_levels(), commander.lower(), commander.inter(),
                commander.upper(), commander.maf_base_control(),
                commander.maf_base_case(), commander.base_info_score(),
                commander.perform_maf_base_control_filter(),
                commander.perform_maf_base_case_filter(),
                commander.perform_base_info_score_filter(),
                commander.fastscore(), commander.no_full(), commander.beta(),
                commander.is_index(), commander.perform_shrinkage(), region,
                reporter);
            // remove all boundaries from the region object to free up memory
            region.clean();
            // skip clumping if not required
            if (!commander.no_clump()) {
                // get the sort by p index vector for target
                // so that we can still find out the relative coordinates of
                // each SNPs This is only required for clumping
                if (!target_file->sort_by_p()) {
                    std::string error_message =
                        "No SNPs left for PRSice processing";
                    reporter.report(error_message);
                    return -1;
                }
                // now perforrm clumping
                target_file->efficient_clumping(
                    commander.use_ref() ? *reference_file : *target_file,
                    reporter, commander.pearson());
                // immediately free the memory
            }

            if (commander.perform_shrinkage()) {
                // Perform order statistic shrinkage using the reference panel
                double max_p_value =
                    (commander.no_full())
                        ? (commander.fastscore() ? commander.bar_upper()
                                                 : commander.upper())
                        : 1;
                target_file->perform_shrinkage(
                    commander.use_ref() ? *reference_file : *target_file,
                    commander.maf_bin(), commander.base_prevalence(),
                    max_p_value, commander.num_shrinkage_perm(),
                    commander.num_sample(), commander.num_case(),
                    commander.num_control(), commander.base_is_binary(),
                    reporter);
            }

            if (commander.use_ref()) delete reference_file;

            // Prepare the SNP vector in target for PRS calculation
            if (!target_file->prepare_prsice()) {
                std::string error_message =
                    "No SNPs left for PRSice processing";
                reporter.report(error_message);
                return -1;
            }
            // count the number of SNPs in each region so that we can skip
            // regions that does not contain any SNP
            // Also, print out the SNP matrix if required
            // TODO: Maybe consider having an index storage for all set so that
            // we can jump around quicker
            // NOTE: When PRSet isn't performed, this has no effect
            target_file->count_snp_in_region(region, commander.out(),
                                             commander.print_snp());
            // we should never sort the ordering of m_existed_snp from now on or
            // it will distord the background index and cause problem

            // check which region are removed
            std::ofstream removed_regions;
            size_t region_size = region.size();
            // We do not initialize the num_post_clump_snp index if we don't
            // perform PRSet
            for (size_t i = 0; i < region_size && commander.perform_set_perm();
                 ++i)
            {
                if (region.num_post_clump_snp(i) == 0) {
                    // this is not included in the analysis
                    if (!removed_regions.is_open()) {
                        // only generate this file if there are region that are
                        // excluded from the analysis
                        removed_regions.open(
                            std::string(commander.out() + ".excluded_regions")
                                .c_str());
                        if (!removed_regions.is_open()) {
                            fprintf(stderr,
                                    "Error: Cannot open file to write: %s\n",
                                    std::string(commander.out()
                                                + ".excluded_regions")
                                        .c_str());
                            return -1;
                        }
                    }
                    removed_regions << region.get_name(i) << std::endl;
                }
            }
            if (removed_regions.is_open()) removed_regions.close();
            // now we start processing each phenotype
            const intptr_t num_pheno = prsice.num_phenotype();
            if (!perform_prslice) {
                // Initialize the progress bar
                prsice.init_process_count(commander,
                                          static_cast<intptr_t>(region.size()),
                                          target_file->num_threshold());
                // remove background from the region process
                // background is only there if permutation is to be performed
                // (It must be there if permutation is to be performed)
                const size_t num_region_process =
                    region.size()
                    - ((region.size() > 1) && commander.perform_set_perm());
                // go through each phenotype
                for (intptr_t i_pheno = 0; i_pheno < num_pheno; ++i_pheno) {
                    // initialize the phenotype & independent variable matrix
                    fprintf(stderr, "\nProcessing the %zu th phenotype\n",
                            i_pheno + 1);
                    // we now initialize the phenotype and covariance matrix
                    prsice.init_matrix(commander, i_pheno, *target_file,
                                       reporter);
                    // we then prepare the output files. Main complication is
                    // for all score and best score. Others are easier
                    prsice.prep_output(commander.out(), commander.all_scores(),
                                       commander.has_prevalence(), *target_file,
                                       region.names(), i_pheno);
                    // go through each region separately
                    // this should reduce the memory usage
                    for (size_t i_region = 0; i_region < num_region_process;
                         ++i_region)
                    {
                        // we will skip any region without SNPs in it but never
                        // skip the first set (Base  set)
                        if (i_region != 0
                            && region.num_post_clump_snp(i_region) == 0)
                            continue;
                        // now we start running PRSice
                        prsice.run_prsice(commander, i_pheno, i_region,
                                          *target_file);
                        if (!commander.no_regress())
                            // if we performed regression, we'd like to generate
                            // the output file (.prsice)
                            prsice.output(commander, region, i_pheno, i_region);
                    }
                    if (!commander.no_regress() && commander.perform_set_perm())
                    {
                        // only perform permutation if regression is performed
                        // and user request it
                        prsice.run_competitive(*target_file, commander,
                                               i_pheno);
                    }
                }
                // finish the progres bar
                prsice.print_progress(true);
                fprintf(stderr, "\n");
                if (!commander.no_regress())
                    // now generate the summary file
                    prsice.summarize(commander, reporter);
            }
            else
            {
                // TODO: Have some idea for PRSlice, but will be drastically
                // different from our current implementation direction. Test
                // that out first before we start implementing it
                std::string error_message =
                    "Error: We currently have not implemented PRSlice. We will "
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
    }
    catch (const std::exception& ex)
    {
        reporter.report(ex.what());
    }
    catch (...)
    {
        std::string error_message = "Error: Bad Allocation exception detected. "
                                    "This is likely due to insufficient memory "
                                    "for PRSice. You can try re-running PRSice "
                                    "with more memory.";
        reporter.report(error_message);
    }
    return 0;
}

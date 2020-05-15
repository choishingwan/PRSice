// This file is part of PRSice-2, copyright (C) 2016-2019
// Shing Wan Choi, Paul F. O’Reilly
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


#include "IITree.h"
#include "cgranges.h"
#include "commander.hpp"
#include "genotype.hpp"
#include "genotypefactory.hpp"
#include "plink_common.hpp"
#include "prsice.hpp"
#include "region.hpp"
#include "reporter.hpp"
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>


void print_empty_region(
    const std::string& out,
    const std::vector<std::vector<size_t>>& region_membership,
    std::vector<std::string>& region_names);

void initialize_target(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    Genotype* target_file, Commander& commander, Reporter& reporter);
void initialize_reference(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    Genotype* target_file, Genotype* reference_file, Commander& commander,
    Reporter& reporter);

std::tuple<std::vector<std::string>, size_t>
add_gene_set_info(Commander& commander, Genotype* target_file,
                  Reporter& reporter)
{
    Region region(commander.get_set(), &reporter);
    const size_t num_regions = region.generate_regions(
        target_file->included_snps_idx(), target_file->included_snps(),
        target_file->max_chr());
    target_file->add_flags(region.get_gene_sets(), num_regions,
                           commander.get_set().full_as_background);
    return {region.get_names(), num_regions};
}


int main(int argc, char* argv[])
{
    const std::string separator =
        "==================================================";
    // initialize reporter, use to generate log
    Reporter reporter;
    try
    {
        // initialize commander, use to parse command line arguments
        // initialization should help setting the default values
        Commander commander;
        try
        {
            if (!commander.process_command(argc, argv, reporter))
            {
                return 0; // only require the usage information
            }
        }
        catch (...)
        {
            return -1; // all error messages should have printed
        }
        // parse the exclusion range and put it into the exclusion object
        // Generate the exclusion region
        std::vector<IITree<size_t, size_t>> exclusion_regions;
        Region::generate_exclusion(exclusion_regions,
                                   commander.exclusion_range());
        GenomeFactory factory;
        Genotype *target_file = nullptr, *reference_file = nullptr;
        try
        {
            // initialize the target object using the factory
            target_file = factory.createGenotype(commander.get_target(),
                                                 commander.get_pheno(),
                                                 commander.delim(), reporter);
            initialize_target(exclusion_regions, target_file, commander,
                              reporter);
            if (commander.use_ref() && commander.need_ref())
            {
                reference_file = factory.createGenotype(
                    commander.get_reference(), commander.get_pheno(),
                    commander.delim(), reporter);
                initialize_reference(exclusion_regions, target_file,
                                     reference_file, commander, reporter);
            }
            exclusion_regions.clear();
            target_file->calc_freqs_and_intermediate(commander.get_target_qc(),
                                                     commander.out(), true);
            if (reference_file != nullptr)
            {
                reference_file->set_thresholds(commander.get_ref_qc());
                reference_file->calc_freqs_and_intermediate(
                    commander.get_ref_qc(), commander.out(), true, target_file);
            }
            if (target_file->num_snps() == 0)
            {
                reporter.report("No SNPs left for PRSice processing");
                return -1;
            }
            auto [region_names, num_regions] =
                add_gene_set_info(commander, target_file, reporter);

            PRSice prsice(commander.get_prs_instruction(),
                          commander.get_p_threshold(), commander.get_pheno(),
                          commander.get_perm(), commander.out(), &reporter);
            prsice.pheno_check();

            if (!commander.get_clump_info().no_clump)
            {
                target_file->build_clump_windows(
                    commander.get_clump_info().distance);
                target_file->sort_by_p();
                // now perform clumping
                target_file->clumping(commander.get_clump_info(),
                                      commander.use_ref() ? *reference_file
                                                          : *target_file,
                                      commander.get_prs_instruction().thread);
            }
            // immediately free the memory
            if (reference_file != nullptr) { delete reference_file; }
            if (commander.ultra_aggressive())
            {
                // we will do something ultra aggressive here: To load all SNP
                // information into memory (does not work for bgen if hard
                // coding isn't used)
                target_file->load_genotype_to_memory();
            }

            // can do the update structure here
            // Use sparse matrix for space and speed
            // Column = set, row = SNPs (because EIGEN is column major)
            // need to also know the number of threshold included
            std::vector<std::vector<size_t>> region_membership;
            target_file->prepare_prsice();
            // from now on, we are not allow to sort the m_existed_snps
            target_file->build_membership_matrix(region_membership, num_regions,
                                                 commander.out(), region_names,
                                                 commander.print_snp());
            // we can now quickly check if any of the region are empty
            try
            {
                print_empty_region(commander.out(), region_membership,
                                   region_names);
            }
            catch (const std::runtime_error& er)
            {
                reporter.report(er.what());
                return -1;
            }
            // Initialize the progress bar
            prsice.init_progress_count(num_regions,
                                       target_file->get_set_thresholds());
            const size_t num_pheno = prsice.num_phenotype();
            for (size_t i_pheno = 0; i_pheno < num_pheno; ++i_pheno)
            {
                if (prsice.pheno_skip(i_pheno))
                {
                    reporter.simple_report("Skipping the "
                                           + std::to_string(i_pheno + 1)
                                           + " th phenotype");
                    continue;
                }
                reporter.simple_report("Processing the "
                                       + std::to_string(i_pheno + 1)
                                       + " th phenotype");
                prsice.new_phenotype(*target_file);
                if (!commander.get_prs_instruction().no_regress)
                {
                    prsice.init_matrix(i_pheno, commander.delim(),
                                       *target_file);
                }
                fprintf(stderr, "Preparing Output Files\n");
                prsice.prep_output(*target_file, region_names, i_pheno,
                                   commander.all_scores());
                // go through each region
                fprintf(stderr, "\nStart Processing\n");
                for (size_t i_region = 0; i_region < num_regions; ++i_region)
                {
                    // always skip background region
                    if (i_region == 1) continue;
                    if (!prsice.run_prsice(i_pheno, i_region, region_membership,
                                           commander.all_scores(),
                                           *target_file))
                    {
                        // did not run
                        continue;
                    }
                    prsice.output(region_names, i_pheno, i_region);
                }
                if (!commander.get_prs_instruction().no_regress)
                {
                    prsice.print_best(*target_file, region_names, i_pheno);
                    if (commander.get_perm().run_set_perm
                        && region_names.size() > 2)
                    {
                        // only perform permutation if regression is performed
                        // and user request it
                        assert(region_membership.size() >= 2);
                        prsice.run_competitive(
                            *target_file, region_membership[1].begin(),
                            region_membership[1].end(), i_pheno);
                    }
                }
            }
            prsice.print_progress(true);
            fprintf(stderr, "\n");
            if (!commander.get_prs_instruction().no_regress)
                // now generate the summary file
                prsice.summarize();
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
        catch (const std::out_of_range& error)
        {
            reporter.report(error.what());
            return -1;
        }
        catch (const std::exception& ex)
        {
            reporter.report(ex.what());
        }
        catch (...)
        {
            std::string error_message =
                "Error: Bad Allocation exception detected. "
                "This is likely due to insufficient memory "
                "for PRSice. You can try re-running PRSice "
                "with more memory.";
            reporter.report(error_message);
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

void print_empty_region(
    const std::string& out,
    const std::vector<std::vector<size_t>>& region_membership,
    std::vector<std::string>& region_names)
{
    bool has_empty_region = false;
    std::ofstream empty_region;
    std::string empty_region_name = out + ".xregion";
    // region_start_idx size always = num_regions
    // check regions to see if there are any empty regions
    for (size_t region_idx = 2; region_idx < region_membership.size();
         ++region_idx)
    {
        if (region_membership[region_idx].empty())
        {
            if (!has_empty_region)
            {
                empty_region.open(empty_region_name.c_str());
                if (!empty_region.is_open())
                {
                    throw std::runtime_error("Error: Cannot open file: "
                                             + empty_region_name
                                             + " to write!");
                }
                has_empty_region = true;
            }
            empty_region << region_names[region_idx] << std::endl;
            region_names[region_idx] = "";
        }
    }
    if (has_empty_region) empty_region.close();
}

void initialize_reference(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    Genotype* target_file, Genotype* reference_file, Commander& commander,
    Reporter& reporter)
{
    const std::string separator =
        "==================================================";
    reporter.report("Start processing reference\n");
    reference_file =
        &reference_file->reference().intermediate(commander.use_inter());
    reporter.report("Loading Genotype info from reference\n" + separator);
    reference_file->load_samples();
    // load the reference file
    const bool verbose = true;
    reference_file->load_snps(commander.out(), exclusion_regions, verbose,
                              target_file);
    reference_file->set_thresholds(commander.get_ref_qc());
}

void initialize_target(
    const std::vector<IITree<size_t, size_t>>& exclusion_regions,
    Genotype* target_file, Commander& commander, Reporter& reporter)
{

    const std::string separator =
        "==================================================";
    target_file = &target_file->keep_nonfounder(commander.nonfounders())
                       .keep_ambig(commander.keep_ambig())
                       .intermediate(commander.use_inter())
                       .set_prs_instruction(commander.get_prs_instruction())
                       .set_weight();
    const std::string base_name = commander.get_base_name();
    reporter.report("Start processing " + base_name + "\n" + separator);
    target_file->snp_extraction(commander.extract_file(),
                                commander.exclude_file());
    auto [filter_count, dup_rs_id] =
        target_file->read_base(commander.get_base(), commander.get_base_qc(),
                               commander.get_p_threshold(), exclusion_regions);
    target_file->print_base_stat(filter_count, dup_rs_id, commander.out(),
                                 commander.get_base_qc().info_score);
    reporter.report("Loading Genotype info from target\n" + separator);
    if (commander.use_ref()) target_file->expect_reference();
    target_file->load_samples();
    const bool verbose = true;
    target_file->load_snps(commander.out(), exclusion_regions, verbose);
    target_file->set_thresholds(commander.get_target_qc());
}

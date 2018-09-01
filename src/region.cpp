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


#include "region.hpp"
#include "genotype.hpp"

// This is for exclusion region
Region::Region(const std::string& exclusion_range, Reporter& reporter)
{
    // no need to do anything if the input is empty
    if (exclusion_range.empty()) return;

    // First check if it is a simple range in format of
    // chr:start-end,chr:start-end
    std::vector<std::string> region_range = misc::split(exclusion_range, ",");
    bool file_input = false;
    if (region_range.size() == 1) {
        // this can either be: 1 range, or an input bed file
        std::vector<std::string> token = misc::split(region_range[0], ":");
        // a range should contain :, if not, this is likely a bed file
        if (token.size() == 1) {
            file_input = true;
        }
    }
    if (file_input) {
        // we need to use a vector, because process_bed takes a vector of string
        // as input
        std::vector<std::string> bed;
        bed.push_back(exclusion_range);
        process_bed(bed, reporter);
    }
    else
    {
        // we will manually insert the regions
        std::vector<region_bound> current_region;
        for (auto&& range : region_range) {
            // go through each exclusion region, separated by comma
            std::vector<std::string> token = misc::split(range, ":");
            if (token.size() >= 2) {
                // this region is in the form of chr:start[-end]
                int chr = get_chrom_code_raw(token[0].c_str());
                int temp, start, end;
                try
                {
                    std::vector<std::string> coordinates =
                        misc::split(token[1], "-");
                    if (coordinates.size() != 1 && coordinates.size() != 2) {
                        std::string error =
                            "Error: Undefined coordinate format: " + token[1]
                            + ". Format of --x-range must either be chr:start "
                              "or chr:start-end";
                        throw std::runtime_error(error);
                    }
                    // try and see if there is a - in the region, if yes, it is
                    // a start-end format
                    temp = misc::convert<int>(coordinates[0].c_str());
                    if (temp < 0) {
                        std::string error = "Error: Negative Start Coordinate "
                                            "for exclusion range!";
                        throw std::runtime_error(error);
                    }
                    start = temp;
                    // if end is not given, this will generate a end region by
                    // adding 1 to start
                    temp = misc::convert<int>(
                        coordinates[(coordinates.size() > 1)].c_str());
                    if (temp < 0) {
                        std::string error = "Error: Negative End Coordinate "
                                            "for exclusion range!";
                        throw std::runtime_error(error);
                    }
                    end = temp + (!(coordinates.size() > 1));
                }
                catch (const std::runtime_error& error)
                {
                    fprintf(stderr, "Error: Non-numeric coordinate(s)!\n");
                    throw std::runtime_error(error.what());
                }

                if (start > end) {
                    std::string message =
                        "Error: Start coordinate should be smaller than "
                        "end coordinate!\n";
                    message.append("start: " + misc::to_string(start) + "\n");
                    message.append("end: " + misc::to_string(end) + "\n");
                    throw std::runtime_error(message);
                }
                region_bound cur_bound;
                cur_bound.chr = chr;
                cur_bound.start = start;
                cur_bound.end = end;
                current_region.push_back(cur_bound);
            }
            else
            {
                std::string message = "Error: Invalid exclusion range format. "
                                      "Should be chr:start or chr:start-end";
                throw std::runtime_error(message);
            }
        }

        // removing overlapped regions and add the result to the region_list
        m_region_list.push_back(solve_overlap(current_region));
    }
    // initalize the region inforamtion with index start at 0
    m_snp_check_index = std::vector<size_t>(1, 0);
}

// constructor
Region::Region(std::vector<std::string> feature, const int window_5,
               const int window_3, const bool run_perm,
               const bool genome_wide_background)
    : m_gtf_feature(feature)
    , m_5prime(window_5)
    , m_3prime(window_3)
    , m_run_perm(run_perm)
    , m_genome_wide_background(genome_wide_background)
{
    // Always include the base region
    // TODO: Better name for the "base" region? Genome wide?
    m_duplicated_names.insert("Base");
    m_region_name.push_back("Base");
    // initialize the region list, add an empty vector as a place holder,
    // we will handle the base region differently
    m_region_list.push_back(std::vector<region_bound>(1));
    // need to also initialize snp_count vector with one item or it might cause
    // trouble in update_flag if user didn't perform generate_region
    m_region_snp_count.push_back(0);
}

// Function to generate all required regions
void Region::generate_regions(const std::string& gtf, const std::string& msigdb,
                              const std::vector<std::string>& bed,
                              const std::string& snp_set,
                              const std::string& multi_snp_sets,
                              const std::string& background,
                              const Genotype& target, Reporter& reporter)
{

    // ensure we don't add the background flag if we don't want to perform
    // perm set
    if (gtf.empty() && bed.size() == 0 && snp_set.empty()
        && multi_snp_sets.empty())
    {
        // initialize the two vector to only contain 1 set (the base)
        // m_snp_check_index is the storage of vector index for quick SNP
        // region check
        m_snp_check_index = std::vector<size_t>(m_region_name.size(), 0);
        m_region_snp_count = std::vector<int>(m_region_name.size());
        return;
    }
    // Issue a warning to users regarding behaviour of PRSet when we are
    // uncertain of the strand of the region
    // doesn't matter if m_5prime==m_3prime though
    if ((m_5prime > 0 || m_3prime > 0) && (m_5prime != m_3prime)) {
        std::string message = "Warning: We will assume a positive strand for "
                              "any features with unspecific strand information "
                              "e.g. \".\"";
        reporter.report(message);
    }
    // first we process the bed files
    process_bed(bed, reporter);
    // now we want to read in the gtf file. One problem with the gtf file is
    // that it can sometime use ensembl gene ID together with the gene name.
    // This cause problem when matching the msigdb file to the GTF file.
    // Therefore we also need to generate a ensembl ID / gene ID to gene name
    // dictionary, which is the id_to_name
    std::unordered_map<std::string, std::set<std::string>> id_to_name;
    // we want to construct a region boundary for each gene. Thus we store the
    // information wihtin the unordered_map
    std::unordered_multimap<std::string, region_bound> gtf_boundary;
    if (!gtf.empty()) {
        // only process the gtf and msigdb file if we have the gtf file
        reporter.report("Processing the GTF file");
        try
        {
            // we want to avoid spending time reading things outside of what we
            // are going to work on. So we ignore any SNP fall on chromosome
            // that we are not going to work on
            uint32_t max_chr = target.max_chr();
            // Read in the GTF information
            gtf_boundary = process_gtf(gtf, id_to_name, max_chr, reporter);
        }
        catch (const std::runtime_error& error)
        {
            // we have problem reading in the gtf file.
            // Terminate the program
            gtf_boundary.clear();
            throw std::runtime_error(error.what());
        }

        std::string message = "A total of "
                              + misc::to_string(gtf_boundary.size())
                              + " genes found in the GTF file";
        reporter.report(message);
        if (gtf_boundary.size() != 0) {
            // if we have more than one gene read from the GTF, we will start
            // processing the msigdb
            process_msigdb(msigdb, gtf_boundary, id_to_name, reporter);
        }
        else
        {
            throw std::runtime_error("Error: GTF file does not contain "
                                     "information of any gene after feature "
                                     "selection");
        }
    }
    // Now we will read in the snp set(s) and add them to the region
    process_snp_sets(snp_set, multi_snp_sets, target, reporter);
    if (m_run_perm) {
        // only generate the background when we want to run the permutation
        if (background.empty() && !m_genome_wide_background && !gtf.empty()) {
            // We want to generate the background region based on the gtf
            // information
            generate_background(gtf_boundary);
        }
        else if (!background.empty())
        {
            // if we have the background information, we will try to read it
            read_background(background, gtf_boundary, id_to_name, reporter);
        }
        else
        {
            // we must still add the background name as it is used to determine
            // the number of region later on
            m_region_name.push_back("Background");
        }
    }
    // initialize snp_check_index to 0, which is the front of each region
    m_snp_check_index = std::vector<size_t>(m_region_name.size(), 0);
    // we also initialize m_region_snp_count to check if any SNP fall into the
    // region. If none of the SNP fall into the region, we will exclude the
    // region from the analysis
    m_region_snp_count = std::vector<int>(m_region_name.size(), 0);
    m_duplicated_names.clear();
    if (!m_has_background && m_run_perm && !m_genome_wide_background) {
        reporter.report("Warning: Background not provided, will use all SNPs "
                        "as background (only affect --set-perm)");
    }
}
// read in the snp sets
// TODO: Not sure about this, but maybe separate the multi_snp set and single
// snp set into two function for clarity?
void Region::process_snp_sets(const std::string& single_snp_set,
                              const std::string& multi_snp_set,
                              const Genotype& target, Reporter& reporter)
{
    // NOTE: For SNP set, we will not do 3' 5' padding
    if (single_snp_set.empty() && multi_snp_set.empty()) return;
    std::string message = "";
    if (!single_snp_set.empty()) {
        std::ifstream input;
        // read in the file contain one and only one SNP set
        // TODO: Maybe allow user input name?
        input.open(single_snp_set.c_str());
        if (!input.is_open()) {
            message = "Error: " + single_snp_set + " cannot be open!";
            throw std::runtime_error(message);
        }
        if (m_duplicated_names.find(single_snp_set) != m_duplicated_names.end())
        {
            message = "Warning: Set name of " + single_snp_set
                      + " is duplicated, it will be ignored";
            reporter.report(message);
        }
        else
        {
            std::vector<region_bound> current_region;
            std::string line;
            while (std::getline(input, line)) {
                misc::trim(line);
                intptr_t chr, loc;
                // get the SNP coordinate from the target file
                if (target.get_snp_loc(line, chr, loc)) {
                    region_bound cur_bound;
                    cur_bound.chr = chr;
                    cur_bound.start = loc;
                    // end is non-inclusive
                    cur_bound.end = loc + 1;
                    current_region.push_back(cur_bound);
                }
            }
            input.close();
            // then remove duplcate. However, this is unlikely for there to be
            // any duplciate unless multiple SNP at the same location is found
            if (current_region.size() > 0) {
                m_region_list.push_back(solve_overlap(current_region));
                m_region_name.push_back(single_snp_set);
                m_duplicated_names.insert(single_snp_set);
            }
        }
    }
    if (!multi_snp_set.empty()) {
        std::ifstream input;
        input.open(multi_snp_set.c_str());
        // we now process the multi snp set input
        if (!input.is_open()) {
            message = "Error: " + multi_snp_set + " cannot be open!";
            throw std::runtime_error(message);
        }
        std::string line;
        std::vector<std::string> token;
        while (std::getline(input, line)) {
            misc::trim(line);
            if (line.empty()) continue;
            token = misc::split(line);
            if (token.size() <= 1) {
                message = "Error: Multi-set file should contain at least 2 "
                          "columns. Did you want to use --snp-set instead?";
                // can be less stringent
                throw std::runtime_error(message);
            }
            // first column should contain the set name
            if (m_duplicated_names.find(token[0]) != m_duplicated_names.end()) {
                message = "Warning: Set name of " + token[0]
                          + " is duplicated, it will be ignored";
                reporter.report(message);
                continue;
            }
            // now generate ther region
            std::vector<region_bound> current_region;
            for (auto&& snp : token) {
                intptr_t chr, loc;
                if (target.get_snp_loc(snp, chr, loc)) {
                    region_bound cur_bound;
                    cur_bound.chr = chr;
                    cur_bound.start = loc;
                    // end is non-inclusive
                    cur_bound.end = loc + 1;
                    current_region.push_back(cur_bound);
                }
            }
            if (current_region.size() > 0) {
                m_region_list.push_back(solve_overlap(current_region));
                m_region_name.push_back(token[0]);
                m_duplicated_names.insert(token[0]);
            }
        }
        input.close();
    }
}

// read in the bed files
void Region::process_bed(const std::vector<std::string>& bed,
                         Reporter& reporter)
{
    // TODO: Allow user define name by modifying their input (e.g. Bed:Name
    bool print_warning = false;
    std::vector<std::string> token;
    std::string message, name, file, line;
    for (auto& b : bed) {
        //  first, check if the input contain a name
        token = misc::split(b, ":");
        if (token.size() == 2) {
            file = token[0];
            name = token[1];
        }
        else if (token.size() == 1)
        {
            name = file = b;
        }
        else
        {
            std::string error_message =
                "Error: Undefine bed file input format: " + b;
            throw std::runtime_error(error_message);
        }

        message = "Reading: " + b;
        reporter.report(message);
        std::ifstream bed_file;
        bool error = false;
        bed_file.open(file.c_str());
        if (!bed_file.is_open()) {
            // previously we allow user to provide missing bed files. But better
            // way should be to terminate and let user check their input
            message =
                "Error: " + file
                + " cannot be open. Please check you have the correct input";

            throw std::runtime_error(message);
        }
        if (m_duplicated_names.find(name) != m_duplicated_names.end()) {
            message =
                "Error: " + name
                + " is duplicated, please check you have the correct input";
            throw std::runtime_error(message);
        }
        std::vector<region_bound> current_region;
        size_t num_line = 0;

        // now start reading in the bed file
        bool first_read = true;
        bool has_strand = false;
        while (std::getline(bed_file, line)) {
            num_line++;
            misc::trim(line);
            if (line.empty()) continue;
            token = misc::split(line);
            if (token.size() < 3) {

                message = "Error: " + file
                          + " contain less than 3 columns, please check your "
                            "bed files in the correct format";
                throw std::runtime_error(message);
            }
            if (first_read) {
                first_read = false;
                if (token.size() > +BED::STRAND) has_strand = true;
            }
            if (has_strand && token.size() <= +BED::STRAND) {
                message = "Error: line " + misc::to_string(num_line)
                          + " of the bed file: " + file
                          + " contain less than than "
                          + misc::to_string(+BED::STRAND)
                          + " columns. BED file should have the same number of "
                            "column for each row. Please check if you have the "
                            "correct input format!";
                throw std::runtime_error(message);
            }
            if (token.size() <= +BED::STRAND && (m_5prime > 0 || m_3prime > 0)
                && (m_5prime != m_3prime) && !print_warning)
            {
                // as bed file does not necessarily contain the strand
                // information, we will issue a warning and assume positive
                // strand. Though for this situation, it might be better for
                // user to use the same padding for 3' and 5'
                std::string message = "Warning: You bed file does not contain "
                                      "strand information, we will assume all "
                                      "regions are on the positive strand, "
                                      "e.g. start coordinates always on the 5' "
                                      "end";
                reporter.report(message);
                print_warning = true;
            }
            intptr_t start = 0, end = 0;
            message = "";
            try
            {
                start = misc::convert<intptr_t>(token[+BED::START]);
                if (start >= 0)
                    ++start; // That's because bed is 0 based
                else
                {
                    message.append("Error: Negative Start Coordinate at line "
                                   + misc::to_string(num_line) + "!");
                    error = true;
                }
            }
            catch (...)
            {
                message.append("Error: Cannot convert start coordinate! (line: "
                               + misc::to_string(num_line) + ")!");
                error = true;
            }
            try
            {
                end = misc::convert<intptr_t>(token[+BED::END]);
                if (end < 0) {
                    message.append("Error: Negative End Coordinate at line "
                                   + misc::to_string(num_line) + "!");
                    error = true;
                }
                else
                {
                    /* If we have
                     *
                     *  chr1   | T | A | C | C | G |
                     *         | | | | | | | | | | |
                     * 1-base  | 1 | 2 | 3 | 4 | 5 |
                     * 0-base  0   1   2   3   4   5
                     *
                     * For bed file with input of
                     * chr1 1 5
                     * we will want the region_bound variable to be
                     * 1 2 6
                     * Such that any SNP with coordinate of chr1:2-5 on the 1
                     * base coordinate system will be included in the analysis
                     */
                    ++end;
                }
            }
            catch (...)
            {
                message.append("Error: Cannot convert end coordinate! (line: "
                               + misc::to_string(num_line) + ")!");
                error = true;
            }
            if (!error && start > end) {
                // only do this check if there isn't an error before
                error = true;
                message.append("Error: Start coordinate should be smaller than "
                               "end coordinate!\n");
                message.append("start: " + misc::to_string(start) + "\n");
                message.append("end: " + misc::to_string(end) + "\n");
            }
            if (error) {
                message.append("Please check your input is correct");
                throw std::runtime_error(message);
            }
            // only include regions that falls into the chromosome of interest
            if (token.size() > +BED::STRAND) {
                // we haev the strand information, therefore can add the padding
                // accordingly
                if (token[+BED::STRAND] == "-") {
                    // negative strand, so add 3' to start and 5' to end
                    start -= m_3prime;
                    if (start < 1) start = 1;
                    end += m_5prime;
                }
                else if (token[+BED::STRAND] == "+"
                         || token[+BED::STRAND] == ".")
                {
                    // positive or unknown strand, add 5' to start and 3' to end
                    start -= m_5prime;
                    if (start < 1) start = 1;
                    end += m_3prime;
                }
                else
                {
                    std::string error = "Error: Undefined strand "
                                        "information. Possibly a malform "
                                        "BED file: "
                                        + token[+BED::STRAND];
                    throw std::runtime_error(error);
                }
            }
            else
            {
                start -= m_5prime;
                if (start < 1) start = 1;
                end += m_3prime;
            }
            // now add the boundary to the region
            region_bound cur_bound;
            cur_bound.chr = get_chrom_code_raw(token[0].c_str());
            cur_bound.start = start;
            cur_bound.end = end;
            current_region.push_back(cur_bound);
        }

        if (!error) {
            // Then we push_back the non-overlapping region
            m_region_list.emplace_back(solve_overlap(current_region));
            // and then we provide the name information
            m_region_name.push_back(name);
            m_duplicated_names.insert(name);
        }
        else
        {
            throw std::runtime_error(message);
        }
        bed_file.close();
    }
}

// it return an unordered_map where the key is the gene ID, and the value is the
// boundary
std::unordered_multimap<std::string, Region::region_bound> Region::process_gtf(
    const std::string& gtf,
    std::unordered_map<std::string, std::set<std::string>>& id_to_name,
    const uint32_t max_chr, Reporter& reporter)
{

    std::unordered_multimap<std::string, Region::region_bound> result_boundary;
    // do nothing when gtf file is not provided
    if (gtf.empty()) return result_boundary;

    std::string line;
    size_t num_line = 0, exclude_feature = 0;

    bool gz_input = false;
    // we want to allow gz file input (as GTF file can be big)
    GZSTREAM_NAMESPACE::igzstream gz_gtf_file;
    if (gtf.substr(gtf.find_last_of(".") + 1).compare("gz") == 0) {
        gz_gtf_file.open(gtf.c_str());
        if (!gz_gtf_file.good()) {
            std::string error_message =
                "Error: Cannot open GTF (gz) to read!\n";
            throw std::runtime_error(error_message);
        }
        gz_input = true;
    }

    // but user can also provide the depressed file
    std::ifstream gtf_file;
    if (!gz_input) {
        gtf_file.open(gtf.c_str());
        if (!gtf_file.is_open()) {
            std::string error_message = "Cannot open gtf file: " + gtf;
            throw std::runtime_error(error_message);
        }
    }

    std::vector<std::string> token, attribute, extract;
    std::string chr, name, id;
    intptr_t chr_code, start, end;
    // this should ensure we will be reading either from the gz stream or
    // ifstream
    while ((!gz_input && std::getline(gtf_file, line))
           || (gz_input && std::getline(gz_gtf_file, line)))
    {

        num_line++;
        misc::trim(line);
        // skip headers
        if (line.empty() || line[0] == '#') continue;
        token = misc::split(line, "\t");
        // convert chr string into consistent chr_coding
        chr_code = get_chrom_code_raw(token[+GTF::CHR].c_str());
        if (in_feature(token[+GTF::FEATURE]) && chr_code <= max_chr) {
            start = 0;
            end = 0;
            try
            {
                start = misc::convert<intptr_t>(token[+GTF::START]);
                if (start < 0) {
                    // we opt for extreme stringency. Will definitely want the
                    // whole gtf file to contain valid entries
                    std::string error =
                        "Error: Negative Start Coordinate! (line: "
                        + std::to_string(num_line)
                        + "). Will ignore the gtf file\n";
                    id_to_name.clear();
                    throw std::runtime_error(error);
                }
            }
            catch (...)
            {
                std::string error =
                    "Error: Cannot convert the start coordinate! (line: "
                    + misc::to_string(num_line) + "). Will ignore the gtf file";
                id_to_name.clear();
                throw std::runtime_error(error);
            }
            try
            {
                end = misc::convert<intptr_t>(token[+GTF::END]);
                if (end < 0) {
                    std::string error =
                        "Error: Negative End Coordinate! (line: "
                        + misc::to_string(num_line)
                        + "). Will ignore the gtf file\n";
                    result_boundary.clear();
                    id_to_name.clear();
                    throw std::runtime_error(error);
                }
                // end is non-inclusive, so we need to add 1 to it
                ++end;
            }
            catch (...)
            {
                std::string error =
                    "Error: Cannot convert the end coordinate! (line: "
                    + std::to_string(num_line) + "). Will ignore the gtf file";
                id_to_name.clear();
                throw std::runtime_error(error);
            }
            // Now extract the name
            attribute = misc::split(token[+GTF::ATTRIBUTE], ";");
            name = "";
            id = "";
            // It is not required by GTF format to contain Gene ID and Gene
            // Name. In that case, we will just refuse to work on this GTF file
            // as we won't be able to conntect it with the MSigDB file
            for (auto& info : attribute) {
                if (info.find("gene_id") != std::string::npos) {
                    extract = misc::split(info);
                    if (extract.size() > 1) {
                        // WARNING: HARD CODING HERE
                        // we assume this should be of the format gene_id "ID"
                        extract[1].erase(std::remove(extract[1].begin(),
                                                     extract[1].end(), '\"'),
                                         extract[1].end());
                        id = extract[1];
                    }
                }
                else if (info.find("gene_name") != std::string::npos)
                {
                    extract = misc::split(info);
                    if (extract.size() > 1) {
                        // WARNING: HARD CODING HERE
                        // Again, we assume this should be of the format
                        // gene_name "Name"
                        extract[1].erase(std::remove(extract[1].begin(),
                                                     extract[1].end(), '\"'),
                                         extract[1].end());
                        name = extract[1];
                    }
                }
            }
            if (name.empty() && id.empty()) {
                // lack both
                std::string message =
                    "Error: GTF file should contain the "
                    "gene_id field. This GTF file does not contain either the "
                    "gene_id field or gene_name field. Please check if you "
                    "have "
                    "the correct file\n";
                id_to_name.clear();
                throw std::runtime_error(message);
            }
            else if (id.empty())
            {
                // lack ID, but mandate
                std::string message = "Error: GTF file should contain the "
                                      "gene_id field. Please check if you have "
                                      "the correct file\n";
                id_to_name.clear();
                throw std::runtime_error(message);
            }
            // the GTF only mandate the ID field, so GTF can miss out the name
            // field.
            // if we have both the name and ID, we will build the dictionary
            if (!id.empty() && !name.empty()) {
                id_to_name[name].insert(id);
            }
            // it is ok to not check for previous error as they all result in
            // throw
            if (start > end) {
                std::string message = "Error: Start coordinate should be "
                                      "smaller than end coordinate!\n";
                message.append("start: " + misc::to_string(start) + "\n");
                message.append("end: " + misc::to_string(end) + "\n");
                id_to_name.clear();
                throw std::runtime_error(message);
            }
            // add padding
            if (token[+GTF::STRAND].compare("-") == 0) {
                start -= m_3prime;
                if (start < 1) start = 1;
                end += m_5prime;
            }
            else if (token[+GTF::STRAND].compare("+") == 0
                     || token[+GTF::STRAND].compare(".") == 0)
            {
                start -= m_5prime;
                if (start < 1) start = 1;
                end += m_3prime;
            }
            else
            {
                // GTF by definition should contain the strand information. If
                // there isn't, then this is a malformed GTF file
                std::string error = "Error: Undefined strand information. "
                                    "Possibly a malform GTF file: "
                                    + token[+GTF::STRAND];
                id_to_name.clear();
                throw std::runtime_error(error);
            }
            // Now add the information to the map using the id
            region_bound cur_bound;
            cur_bound.chr = chr_code;
            cur_bound.start = start;
            cur_bound.end = end;
            result_boundary.insert(std::make_pair(id, cur_bound));
        }
        else
        {
            exclude_feature++;
        }
    }
    std::string message = "";
    if (exclude_feature == 1) {
        message.append("A total of " + std::to_string(exclude_feature)
                       + " entry removed due to feature selection");
        reporter.report(message);
    }
    if (exclude_feature > 1) {
        message.append("A total of " + std::to_string(exclude_feature)
                       + " entries removed due to feature selection");
        reporter.report(message);
    }
    return result_boundary;
}

void Region::process_msigdb(
    const std::string& msigdb,
    const std::unordered_multimap<std::string, Region::region_bound>& gtf_info,
    const std::unordered_map<std::string, std::set<std::string>>& id_to_name,
    Reporter& reporter)
{
    // Got nothing to do
    if (msigdb.empty() || gtf_info.size() == 0) return;
    // Assume format = Name Gene
    // doesn't matter if the file contain the url or not. As it is unlikely for
    // us to have the url as a gene
    std::ifstream input;
    std::vector<std::string> msigdb_files = misc::split(msigdb, ",");
    bool has_open = false;
    std::string line, name;
    std::vector<std::string> token;
    for (auto&& m : msigdb_files) {
        input.open(m.c_str());
        if (!input.is_open()) {
            std::string error_message = "Warning: Cannot open " + msigdb =
                                            ", will skip this file";
            reporter.report(error_message);
        }
        else
        {
            has_open = true;
            std::vector<region_bound> current_region;
            while (std::getline(input, line)) {
                misc::trim(line);
                if (line.empty()) continue;
                token = misc::split(line);
                if (token.size() < 2) // Will treat the url as gene just in case
                {
                    std::string message =
                        "Error: Each line require at least 2 information\n";
                    message.append(line);
                    reporter.report(message);
                }
                else if (m_duplicated_names.find(token[0])
                         == m_duplicated_names.end())
                {
                    // ignore any gene set with duplicated name
                    name = token[0];
                    // reusing the current_region vector
                    current_region.clear();
                    current_region.shrink_to_fit();
                    current_region.reserve(token.size());
                    for (auto& gene : token) {
                        // for each gene, see if we can find it in the gtf file
                        auto&& gtf_search = gtf_info.equal_range(gene);
                        // for equal range, first==second when not found
                        if (gtf_search.first == gtf_search.second) {
                            // we cannot find the gene name in the gtf
                            // information, check if we can find the id?
                            auto&& gene_name_iter = id_to_name.find(gene);
                            if (gene_name_iter != id_to_name.end()) {
                                // we found a way to convert the gene name
                                // to gene id
                                auto& gene_name = gene_name_iter->second;
                                // when one gene name correspond to
                                // multiple gene id, we will take
                                // all of them
                                for (auto&& translate : gene_name) {
                                    auto&& gene_gtf =
                                        gtf_info.equal_range(translate);
                                    // now read in all the regions
                                    for (auto&& it = gene_gtf.first;
                                         it != gene_gtf.second; ++it)
                                        current_region.push_back(it->second);
                                }
                            }
                        }
                        else
                        {
                            // Read in all the regions. Each gene can have
                            // multiple regions
                            for (auto&& it = gtf_search.first;
                                 it != gtf_search.second; ++it)
                                current_region.push_back(it->second);
                        }
                    }
                    // now clean out the overlaps and generate the final gene
                    // set
                    m_region_list.emplace_back(solve_overlap(current_region));
                    // and give the name
                    m_region_name.push_back(name);
                    m_duplicated_names.insert(name);
                }
                else
                {
                    reporter.report("Duplicated Set: " + token[0]
                                    + ". It will be ignored");
                }
            }
            input.close();
        }
    }
}

void Region::read_background(
    const std::string& background,
    const std::unordered_multimap<std::string, region_bound>& gtf_info,
    const std::unordered_map<std::string, std::set<std::string>>& id_to_name,
    Reporter& reporter)
{
    const std::unordered_map<std::string, int> file_type{
        {"bed", 1}, {"range", 0}, {"gene", 2}};
    // format of the background string should be name:format
    std::vector<std::string> background_info = misc::split(background, ":");
    if (background_info.size() != 2) {
        std::string error =
            "Error: Format of --background should be <File Name>:<File Type>";
        throw std::runtime_error(error);
    }
    // check if we know the format
    auto&& type = file_type.find(background_info[1]);
    if (type == file_type.end()) {
        std::string error = "Error: Undefined file type. Supported formats are "
                            "bed, gene or range";
        throw std::runtime_error(error);
    }
    // now read in the file
    std::ifstream input;
    input.open(background_info[0].c_str());
    if (!input.is_open()) {
        std::string error =
            "Error: Cannot open background file: " + background_info[0];
        throw std::runtime_error(error);
    }
    bool print_warning = false, error = false;
    std::vector<Region::region_bound> current_bound;
    std::string line;
    if (type->second == 0 || type->second == 1) {
        // this is either a range format or a bed file
        // the only difference is range is 1 based and bed is 0 based
        size_t num_line = 0;
        std::vector<std::string> token;
        while (std::getline(input, line)) {
            num_line++;
            misc::trim(line);
            if (line.empty()) continue;
            token = misc::split(line);
            if (token.size() < 3) {
                std::string message =
                    "Error: " + background_info[0]
                    + " contain less than 3 columns, it will be ignored";
                throw std::runtime_error(message);
            }
            if (token.size() <= +BED::STRAND && (m_5prime > 0 || m_3prime > 0)
                && (m_5prime != m_3prime) && !print_warning)
            {
                std::string message = "Warning: You bed file does not contain "
                                      "strand information, we will assume all "
                                      "regions are on the positive strand, "
                                      "e.g. start coordinates always on the 5' "
                                      "end";
                reporter.report(message);
                print_warning = true;
            }
            intptr_t start = 0, end = 0;
            std::string message = "";
            try
            {
                start = misc::convert<intptr_t>(token[+BED::START]);
                if (start >= 0)
                    // That's because bed is 0 based and range format is 1
                    // based. and type for bed is 1 and type for range is 0
                    start = start + type->second;
                else
                {
                    message.append("Error: Negative Start Coordinate at line "
                                   + misc::to_string(num_line) + "!");
                    error = true;
                }
            }
            catch (...)
            {
                message.append("Error: Cannot convert start coordinate! (line: "
                               + misc::to_string(num_line) + ")!");
                error = true;
            }
            try
            {
                // we add 1 to end if it is type range because our end is
                // non-inclusive but for range the end is inclusive
                end = misc::convert<int>(token[+BED::END]) + 1 - type->second;
                if (end < 0) {
                    message.append("Error: Negative End Coordinate at line "
                                   + misc::to_string(num_line) + "!");
                    error = true;
                }
            }
            catch (...)
            {
                message.append("Error: Cannot convert end coordinate! (line: "
                               + std::to_string(num_line) + ")!");
                error = true;
            }
            if (!error || start > end) {
                // don't check if it's already error out
                error = true;
                message.append("Error: Start coordinate should be smaller than "
                               "end coordinate!\n");
                message.append("start: " + std::to_string(start) + "\n");
                message.append("end: " + std::to_string(end) + "\n");
            }
            if (error) break;
            // only include regions that falls into the chromosome of interest
            // the strand location is different depending on the type
            // if it is bed, then we use the STRAND index
            // if not, then we assume the format is CHR START END STRAND
            std::vector<std::string>::size_type strand_index =
                (type->second == 1) ? (+BED::STRAND) : (+BED::END + 1);

            if (token.size() > strand_index) {
                if (token[strand_index] == "-") {
                    start -= m_3prime;
                    if (start < 1) start = m_3prime;
                    end += m_5prime;
                }
                else if (token[strand_index].compare("+") == 0
                         || token[strand_index].compare(".") == 0)
                {
                    start -= m_5prime;
                    if (start < 1) start = 1;
                    end += m_3prime;
                }
                else
                {
                    std::string error = "Error: Undefined strand "
                                        "information. Possibly a malform "
                                        "BED file: "
                                        + token[+BED::STRAND];
                    throw std::runtime_error(error);
                }
            }
            else
            {
                start -= m_5prime;
                if (start < 1) start = 1;
                end += m_3prime;
            }
            region_bound cur_bound;
            cur_bound.chr = get_chrom_code_raw(token[0].c_str());
            cur_bound.start = start;
            cur_bound.end = end;
            // this should help us to avoid problem
            current_bound.push_back(cur_bound);
        }
    }
    else
    {
        // gene list format
        while (std::getline(input, line)) {
            // read in the gene file
            misc::trim(line);
            if (line.empty()) continue;
            auto&& gtf_search = gtf_info.equal_range(line);
            if (gtf_search.first == gtf_search.second) {
                // we cannot find the gene name in the gtf information (which
                // uses gene ID)
                auto&& id_search = id_to_name.find(line);
                if (id_search != id_to_name.end()) {
                    // we found a way to convert the gene name to gene id
                    auto& name = id_search->second;
                    // problem is, one gene name can correspond to multiple gene
                    // id in that case, wee will take all of them
                    for (auto&& translate : name) {
                        auto&& gtf_name_search =
                            gtf_info.equal_range(translate);
                        for (auto&& it = gtf_name_search.first;
                             it != gtf_name_search.second; ++it)
                            // now add in all the genome region
                            current_bound.push_back(it->second);
                    }
                }
            }
            else
            {
                for (auto&& it = gtf_search.first; it != gtf_search.second;
                     ++it)
                    // now add in all the genome region
                    current_bound.push_back(it->second);
            }
        }
    }
    input.close();
    // only indicate there is a background when we read in more than one region
    m_has_background = (current_bound.size() > 0);
    if (!m_has_background) {
        throw std::runtime_error("Error: Background is empty. Please make sure "
                                 "your background file is valid.");
    }
    m_region_list.push_back(solve_overlap(current_bound));
    m_region_name.push_back("Background");
}

// generate the background region using the gtf information
void Region::generate_background(
    const std::unordered_multimap<std::string, Region::region_bound>& gtf_info)
{
    // this is ineffective.But will do for now
    // one way to speed things up is to perform the overlap removal as we loop
    // through
    std::vector<Region::region_bound> temp_storage;
    if (gtf_info.size() != 0) {
        temp_storage.reserve(gtf_info.size());
        for (auto&& gtf : gtf_info) {
            temp_storage.push_back(gtf.second);
        }
        m_has_background = true;
    }
    if (temp_storage.size() == 0) {
        throw std::runtime_error(
            "Error: background generated from the GTF file is empty. Maybe all "
            "gene regions are filtered by feature or you gtf file is empty?");
    }
    m_region_list.push_back(solve_overlap(temp_storage));
    m_region_name.push_back("Background");
}

std::vector<Region::region_bound>
Region::solve_overlap(std::vector<Region::region_bound>& current_region)
{
    // sort the region according to their coordinate
    std::sort(begin(current_region), end(current_region),
              [](region_bound const& t1, region_bound const& t2) {
                  if (t1.chr == t2.chr) {
                      if (t1.start == t2.start) return t1.end < t2.end;
                      return t1.start < t2.start;
                  }
                  else
                      return t1.chr < t2.chr;
              });
    std::vector<Region::region_bound> result;
    result.reserve(current_region.size());
    intptr_t prev_chr = -1;
    intptr_t prev_start = 0;
    intptr_t prev_end = 0;
    // optimize for human for now
    for (auto&& bound : current_region) {
        if (prev_chr == -1) {
            prev_chr = bound.chr;
            prev_start = bound.start;
            prev_end = bound.end;
        }
        else if (prev_chr != bound.chr || bound.start > prev_end)
        {
            // new region
            region_bound cur_bound;
            cur_bound.chr = prev_chr;
            cur_bound.start = prev_start;
            cur_bound.end = prev_end;
            result.push_back(cur_bound);
            prev_chr = bound.chr;
            prev_start = bound.start;
            prev_end = bound.end;
        }
        else if (bound.end > prev_end)
        {
            // same region and the end is further down
            prev_end = bound.end;
        }
    }
    if (prev_chr != -1) {
        // we need to remember to add in the last region
        region_bound cur_bound;
        cur_bound.chr = prev_chr;
        cur_bound.start = prev_start;
        cur_bound.end = prev_end;
        result.push_back(cur_bound);
    }
    result.shrink_to_fit();
    return result;
}

Region::~Region() {}
std::vector<Region::region_bound>::size_type
Region::binary_search_region(const intptr_t chr, const intptr_t loc,
                             std::vector<region_bound>::size_type left,
                             std::vector<region_bound>::size_type right) const
{
    std::vector<region_bound>::size_type midPoint = left + (right - left) / 2;
    while (left < right) {
        midPoint = left + (right - left) / 2;
        auto&& cur = m_region_list.front()[midPoint];
        if (chr > cur.chr) {
            left = midPoint + 1;
        }
        else if (chr == cur.chr)
        {
            if (loc >= cur.end) {
                // if the target location is larger or equal to the end of the
                // range, move up
                left = midPoint + 1;
            }
            else if (loc < cur.start)
            {
                // if the target location is smaller than the start of range,
                // move down
                right = midPoint;
            }
            else
            {
                // we are not bigger than the end and not smaller than the
                // start, therefore we are there
                return midPoint;
            }
        }
        else
        {
            right = midPoint;
        }
    }
    return midPoint;
}

bool Region::check_exclusion(const intptr_t chr, const intptr_t loc)
{
    if (m_region_list.empty() || m_region_list.front().empty()) return false;
    // there is only one region, so we can ignore the for loop
    // problem with check_exclusion is that the input might not be sorted
    // (looking at you, listed input). Therefore, we cannot use our O(1)
    // searching algorithm, but might need something like the binary search
    std::vector<region_bound>::size_type left = 0;
    std::vector<region_bound>::size_type right = m_region_list.front().size();
    std::vector<region_bound>::size_type midPoint = left + (right - left) / 2;
    while (left < right) {
        midPoint = left + (right - left) / 2;
        auto&& cur = m_region_list.front()[midPoint];
        if (chr > cur.chr) {
            left = midPoint + 1;
        }
        else if (chr == cur.chr)
        {
            if (loc >= cur.end) {
                // if the target location is larger or equal to the end of the
                // range, move up
                left = midPoint + 1;
            }
            else if (loc < cur.start)
            {
                // if the target location is smaller than the start of range,
                // move down
                right = midPoint;
            }
            else
            {
                // we are not bigger than the end and not smaller than the
                // start, therefore we are there
                break;
            }
        }
        else
        {
            right = midPoint;
        }
    }

    if (midPoint >= m_region_list.front().size())
        return false;
    else
    {
        auto&& cur = m_region_list.front()[midPoint];
        if (cur.chr == chr && cur.start <= loc && cur.end > loc) return true;
    }
    return false;
}

void Region::update_flag(const intptr_t chr, const std::string& rs,
                         intptr_t loc, std::vector<uintptr_t>& flag)
{
    // while rs is unused here, we will keep it as it is useful for
    // debugging


    // we first make sure the base region is true for all SNPs
    const std::vector<std::string>::size_type region_size =
        m_region_name.size();
    SET_BIT(0, flag.data());

    m_region_snp_count[0]++;
    // if we want to perform competitive p-value calclation and use all SNP
    // as background, we can just add that in
    if (m_genome_wide_background && m_run_perm && region_size > 1) {
        // use everything as background
        SET_BIT(region_size - 1, flag.data());
        m_region_snp_count[region_size - 1]++;
    }

    intptr_t region_chr;
    // indicate if we have a likelihood of changing chr.
    // Situation: SNP on chr1 but out of bound of last chr1 boundary,
    // move to next bound, which is chr2, we don't want it to move at all
    bool moved_chr = false;
    // go through each region (including backgroudn)
    for (std::vector<std::string>::size_type i_region = 1;
         i_region < region_size; ++i_region)
    {
        // find out how many boundaries are there within the current region
        auto&& current_region = m_region_list[i_region];
        const size_t current_region_size = current_region.size();

        // while we still have boundary left
        // allow to go through all boundaries
        // as we assume SNPs are read in the correct order (sorted)
        // we use m_snp_check_index to track the last SNP's boundary index
        // and start our search from there so that we can skip un-necessary
        // comparison
        moved_chr = false;
        while (m_snp_check_index[i_region] < current_region_size) {
            // obtain the current boundary as defined by m_snp_check_index
            auto&& current_bound = current_region[m_snp_check_index[i_region]];
            region_chr = current_bound.chr;

            while (chr > current_bound.chr
                   && m_snp_check_index[i_region] < current_region_size)
            {
                // as long as our chr code for the SNP > the region chr code
                // and we have not exhausted all the region boundary of the
                // current set

                // if we have previously jumped to a specific chromosome, we
                // should not jump again
                if (moved_chr) break;
                // we will move to the next boundary and see if we are now
                // in the correct chr
                m_snp_check_index[i_region]++;
                current_bound = current_region[m_snp_check_index[i_region]];
                region_chr = current_bound.chr;
            }
            // no matter what, we assume we have jumped a chr
            moved_chr = true;
            if (m_snp_check_index[i_region] >= current_region_size
                || chr < region_chr)
                // if the region chr code is larger than the SNP chr or if
                // we have exhaused all boundary, break
                break;
            // when we reach here, the only possible reason for this break
            // is chr < region_chr (So technically, we can remove this if
            // clause)
            if (chr != region_chr) {
                break;
            }
            else if (current_bound.end <= loc)
            {
                // if the region bound is still earlier than or equal to
                // this current SNP's location, we will iterate to next
                // bound (because the end is non-inclusive)
                m_snp_check_index[i_region]++;
            }
            else if (current_bound.start <= loc && current_bound.end > loc)
            {
                // This is the region
                SET_BIT(i_region, flag.data());
                m_region_snp_count[i_region]++;
                break;
            }
            else if (current_bound.start > loc)
                // the current bound has already passed the current region.
                // no need to do anything
                break;
        }
    }
}

void Region::print_region_number(Reporter& reporter) const
{
    std::string message = "";
    if (m_region_name.size() == 1) {
        message = "1 region included";
    }
    else if (m_region_name.size() > 1)
    {
        // -1 to remove the background count, as we are not going to print
        // the background anyway
        message = "A total of " + misc::to_string(m_region_name.size() - 2)
                  + " regions plus the base region are included";
    }
    reporter.report(message);
}

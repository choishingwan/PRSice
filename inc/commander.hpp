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

#ifndef COMMANDER_H
#define COMMANDER_H

#include "enumerators.h"
#include "gzstream.h"
#include "misc.hpp"
#include "storage.hpp"
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <map>
#include <random>
#include <reporter.hpp>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <zlib.h>
#ifdef _WIN32
#include <windows.h>
#endif

const std::string version = "2.2.10";
const std::string date = "2019-09-20";
class Commander
{
public:
    Commander();
    virtual ~Commander();
    bool init(int argc, char* argv[], Reporter& reporter);
    std::string get_base_name() const
    {
        return misc::remove_extension<std::string>(
            misc::base_name<std::string>(m_base_info.file_name));
    }
    const BaseFile& get_base() const { return m_base_info; }
    const CalculatePRS& get_prs_instruction() const { return m_prs_info; }
    const Clumping& get_clump_info() const { return m_clump_info; }
    const GeneSets& get_set() const { return m_prset; }
    const GenoFile& get_target() const { return m_target; }
    const GenoFile& get_reference() const { return m_reference; }
    const Permutations& get_perm() const { return m_perm_info; }
    const Phenotype& get_pheno() const { return m_pheno_info; }
    const PThresholding& get_p_threshold() const { return m_p_thresholds; }
    const QCFiltering& get_base_qc() const { return m_base_filter; }
    const QCFiltering& get_target_qc() const { return m_target_filter; }
    const QCFiltering& get_ref_qc() const { return m_ref_filter; }
    bool use_ref() const
    {
        return !(m_reference.file_list.empty()
                 && m_reference.file_name.empty());
    }
    bool need_ref() const
    {
        return !m_clump_info.no_clump || m_prs_info.use_ref_maf;
    }
    bool use_inter() const { return m_allow_inter; }
    // misc
    std::string delim() const { return m_id_delim; }
    /*!
     * \brief Get output prefix
     * \return the output prefix
     */
    std::string out() const { return m_out_prefix; }
    /*!
     * \brief Get the exclusion_range
     * \return The string containign the exclusion range
     */
    std::string exclusion_range() const { return m_exclusion_range; }
    /*!
     * \brief Return if all score should be printed
     * \return True if all score should be printed
     */
    bool all_scores() const { return m_print_all_scores; }
    /*!
     * \brief Return if non-cumulative PRS should be calculated
     * \return True for non-cumulative PRS calculation
     */

    /*!
     * \brief Return if we should generate the .snp file
     * \return true if we should
     */
    bool print_snp() const { return m_print_snp; }


    /*!
     * \brief Return the maximum memory allowed to use
     * \param detected should be the available memory
     * \return Maximum memory we can use
     */
    unsigned long long max_memory(const unsigned long long detected) const
    {
        return (!m_provided_memory || m_memory > detected) ? detected
                                                           : m_memory;
    }
    unsigned long long memory() const { return m_memory; }

    // prs_snp_filtering
    /*!
     * \brief Return the file name of the exclusion file
     * \return the name of the exclusion file
     */
    std::string exclude_file() const { return m_exclude_file; }
    /*!
     * \brief Return the file name of the extraction file
     * \return the name of the extract file
     */
    std::string extract_file() const { return m_extract_file; }

    /*!
     * \brief Check if we should retain ambiguous SNPs
     * \return true if ambiguous SNPs should be retained
     */
    bool keep_ambig() const { return m_keep_ambig; }


    /*!
     * \brief Check if we want to include non-founders in our analysis
     * \return True if non-founders should be included
     */
    bool nonfounders() const { return m_include_nonfounders; }
    bool enable_mmap() const { return m_enable_mmap; }


protected:
private:
    const std::vector<std::string> supported_types = {"bed", "ped", "bgen"};
    std::string m_id_delim = " ";
    std::string m_out_prefix = "PRSice";
    std::string m_exclusion_range = "";
    std::string m_exclude_file = "";
    std::string m_extract_file = "";
    std::string m_help_message;
    // TODO: might consider using 1e-8 instead
    unsigned long long m_memory = 1e10;
    int m_allow_inter = false;
    int m_enable_mmap = false;
    int m_include_nonfounders = false;
    int m_keep_ambig = false;
    int m_print_all_scores = false;
    int m_print_snp = false;
    int m_user_no_default = false;
    bool m_provided_memory = false;
    bool m_set_delim = false;

    BaseFile m_base_info;
    CalculatePRS m_prs_info;
    Clumping m_clump_info;
    GeneSets m_prset;
    GenoFile m_target;
    GenoFile m_reference;
    Permutations m_perm_info;
    QCFiltering m_target_filter;
    QCFiltering m_ref_filter;
    QCFiltering m_base_filter;
    Phenotype m_pheno_info;
    PThresholding m_p_thresholds;

    std::map<std::string, std::string> m_parameter_log;
    std::string m_error_message = "";
    ////////////////////////////////////////////
    //
    // end of struct definition
    //
    ////////////////////////////////////////////

    /*!
     * \brief Parsing all command line arguments
     * \param argc
     * \param argv
     * \param optString is the string containing the short flags
     * \param longOpts is the struct containing all the options
     * \param reporter is the object to report all messages
     * \return true if we want to continue the program
     */
    bool parse_command(int argc, char* argv[], const char* optString,
                       const struct option longOpts[], Reporter& reporter);
    /*!
     * \brief Print the usage information
     */
    void usage();
    /*!
     * \brief Set the help message
     */
    void set_help_message();


    /*!
     * \brief Function to check if parameters for clumping
     * \param message is the parameter storage
     * \param
     * error_message is the storage for error messages
     * \return  return true if
     * clumping parameters are alright
     */
    bool clump_check();
    /*!
     * \brief Function to check if parameters for reference panel are correct
     * \param message is the parameter storage
     * \param
     * error_message is the storage for error messages
     * \return  return true if
     * reference panel parameters are alright
     */
    bool ref_check();
    /*!
     * \brief Function to check if parameters for covariate are correct
     * \return true if parameters are alright
     */
    bool covariate_check();
    /*!
     * \brief Function to check if parameters for filtering are correct
     * \param error_message is the storage for error messages
     * \return true if parameters are alright
     */
    bool filter_check();
    /*!
     * \brief Function to check if parameters for misc options are correct
     * \param error_message is the storage for error messages
     * \return true if parameters are alright
     */
    bool misc_check();
    /*!
     * \brief Function to check if parameters for prset are correct
     * \param error_message is the storage for error messages
     * \return true if parameters are alright
     */
    bool prset_check();
    /*!
     * \brief Function to check if parameters for prsice thresholding are
     * correct \param error_message is the storage for error messages \return
     * true if parameters are alright
     */
    bool prsice_check();
    /*!
     * \brief Function to check if parameters for target are correct
     * \param error_message is the storage for error messages
     * \return true if parameters are alright
     */
    bool target_check();
    bool pheno_check();


    std::vector<std::string> transform_covariate(const std::string& cov_in);


    /////////////////////////////////////////////////
    /// REFACTORED FUNCTIONS
    /////////////////////////////////////////////////

    int32_t maximum_thread()
    {
        int32_t max_threads = 1;
#if defined(WIN32) || defined(_WIN32) \
    || defined(__WIN32) && !defined(__CYGWIN__)
        // max thread estimation using windows
        SYSTEM_INFO sysinfo;
        GetSystemInfo(&sysinfo);
        max_threads = sysinfo.dwNumberOfProcessors;
        int32_t known_procs = max_threads;
#else
        int32_t known_procs =
            static_cast<int32_t>(sysconf(_SC_NPROCESSORS_ONLN));
        max_threads = (known_procs == -1) ? 1 : known_procs;
#endif
        return max_threads;
    }

    inline void set_string(const std::string& input, const std::string& c,
                           size_t base_index)
    {
        check_duplicate(c);
        m_parameter_log[c] = input;
        m_base_info.column_name[base_index] = input;
        m_base_info.has_column[base_index] = true;
    }
    inline void set_string(const std::string& input, const std::string& c,
                           std::string& target, bool& target_boolean,
                           bool add_quote = false)
    {
        check_duplicate(c);
        if (add_quote) { m_parameter_log[c] = "\"" + input + "\""; }
        else
        {
            m_parameter_log[c] = input;
        }
        target = input;
        target_boolean = true;
    }
    inline void set_string(const std::string& input, const std::string& c,
                           std::string& target)
    {
        bool dumy;
        set_string(input, c, target, dumy);
    }

    template <typename T>
    inline bool load_numeric_vector(const std::string& input,
                                    const std::string& c,
                                    std::vector<T>& target)
    {
        // should always have an input
        if (input.empty()) return false;
        std::string comma = "";
        if (m_parameter_log.find(c) != m_parameter_log.end()) { comma = ","; }
        m_parameter_log[c].append(comma + input);
        if (!input.empty() && input.back() == ',')
        {
            m_error_message.append(
                "Warning: , detected at end of input: " + input
                + ". Have you accidentally included space in "
                  "your input? (Space is not allowed)\n");
        }
        std::vector<std::string> token = misc::split(input, ",");
        try
        {
            for (auto&& bar : token) target.push_back(misc::convert<T>(bar));
        }
        catch (...)
        {
            m_error_message.append("Error: Non numeric argument passed to " + c
                                   + ": " + input + "!\n");
            return false;
        }
        return true;
    }
    size_t number_boolean(const std::string& input, bool& result)
    {
        bool found = false;
        if (misc::hasEnding(input, "T") || misc::hasEnding(input, "TRUE"))
        {
            result = true;
            found = true;
        }
        else if (misc::hasEnding(input, "F") || misc::hasEnding(input, "FALSE"))
        {
            found = true;
            result = false;
        }
        if (!found) { throw std::runtime_error(""); }
        try
        {
            return misc::string_to_size_t(
                input.substr(0, input.length() - 1).c_str());
        }
        catch (const std::runtime_error&)
        {
            throw std::runtime_error("");
        }
    }
    inline bool parse_binary_vector(const std::string& input,
                                    const std::string& c,
                                    std::vector<bool>& target)
    {
        // allow more complex regrex
        // should not come in without something
        if (input.empty()) return false;
        std::string comma = "";
        if (m_parameter_log.find(c) != m_parameter_log.end()) { comma = ","; }
        m_parameter_log[c].append(comma + c);
        if (!input.empty() && input.back() == ',')
        {
            m_error_message.append(
                "Warning: , detected at end of input: " + input
                + ". Have you accidentally included space in "
                  "your input? (Space is not allowed)");
        }
        std::vector<std::string> token = misc::split(input, ",");
        try
        {
            for (auto&& bin : token)
            {
                // check if this is true or false, if, not, try parsing
                std::transform(bin.begin(), bin.end(), bin.begin(), ::toupper);
                if (bin == "T" || bin == "TRUE") { target.push_back(true); }
                else if (bin == "F" || bin == "FALSE")
                {
                    target.push_back(false);
                }
                else
                {
                    try
                    {
                        bool value = false;
                        size_t num_repeat = number_boolean(bin, value);
                        for (size_t i = 0; i < num_repeat; ++i)
                        { target.push_back(value); }
                    }
                    catch (const std::runtime_error&)
                    {
                        m_error_message.append(
                            "Error: Invalid argument passed to " + c + ": "
                            + input
                            + "! Require binary arguments e.g. T/F, "
                              "True/False, or numerical combination e.g. "
                              "10T,3F etc\n");
                        return false;
                    }
                }
            }
        }
        catch (...)
        {
            m_error_message.append(
                "Error: Invalid argument passed to " + c + ": " + input
                + "! Require binary arguments e.g. T/F, True/False\n");
            return false;
        }
        return true;
    }

    inline unsigned long long integer_value(const std::string& str,
                                            const size_t& unit)
    {
        double cur_dist = misc::convert<double>(str) * unit;
        unsigned long long res = static_cast<unsigned long long>(cur_dist);
        if (trunc(cur_dist) != cur_dist)
        {
            throw std::runtime_error("Error: Non-integer value obtained: "
                                     + misc::to_string(str) + " x "
                                     + misc::to_string(unit) + "\n");
        }
        return res;
    }
    inline unsigned long long get_unit(const std::string& in,
                                       const bool bit_thousand = false)
    {
        unsigned long long multiplication = 1;
        const unsigned long long thousand = bit_thousand ? 1024 : 1000;
        if (misc::hasEnding(in, "bp") || misc::hasEnding(in, "b"))
        { multiplication = 1; }
        else if (misc::hasEnding(in, "kb") || misc::hasEnding(in, "k"))
        {
            multiplication = thousand;
        }
        else if (misc::hasEnding(in, "mb") || misc::hasEnding(in, "m"))
        {
            multiplication = thousand * thousand;
        }
        else if (misc::hasEnding(in, "gb") || misc::hasEnding(in, "g"))
        {
            multiplication = thousand * thousand * thousand;
        }
        else if (misc::hasEnding(in, "tb") || misc::hasEnding(in, "t"))
        {
            multiplication = thousand * thousand * thousand * thousand;
        }
        else
        {
            throw std::runtime_error("Error: Undefined input unit: " + in);
        }
        return multiplication;
    }
    inline unsigned long long set_distance(const std::string& input,
                                           const std::string& command,
                                           const size_t default_unit,
                                           bool& error,
                                           const bool for_memory = false)
    {
        check_duplicate(command);
        std::string in = input;
        std::transform(in.begin(), in.end(), in.begin(), ::tolower);
        m_parameter_log[command] = in;
        unsigned long long dist = 0;
        unsigned long long multiplication = 1;
        const size_t thousand = for_memory ? 1024 : 1000;
        try
        {
            try
            {
                dist = integer_value(input, default_unit);
            }
            catch (const std::runtime_error& er)
            {
                m_error_message.append(er.what());
                error = true;
                return ~size_t(0);
            }
            std::string unit = "bp";
            if (default_unit == thousand)
                unit = "kb";
            else if (default_unit == thousand * thousand)
                unit = "mb";
            m_parameter_log[command] = input + unit;
            return dist;
        }
        catch (...)
        {
            try
            {
                multiplication = get_unit(in, for_memory);
            }
            catch (const std::runtime_error& er)
            {
                m_error_message.append(er.what());
            }
        }
        try
        {
            size_t length = 2;
            if (misc::hasEnding(in, "b")) { length = 1; }
            return integer_value(in.substr(in.length() - length),
                                 multiplication);
        }
        catch (const std::runtime_error& er)
        {
            m_error_message.append(er.what());
        }
        m_parameter_log.erase(command);
        error = true;
        return ~size_t(0);
    }

    void check_duplicate(const std::string& c)
    {
        if (m_parameter_log.find(c) != m_parameter_log.end())
        {
            m_error_message.append("Warning: Duplicated argument --" + c
                                   + "\n");
        }
    }

    template <typename Type>
    inline bool set_numeric(const std::string& input, const std::string& c,
                            Type& target, bool& target_boolean)
    {
        check_duplicate(c);
        m_parameter_log[c] = input;
        try
        {
            target = misc::convert<Type>(input);
            target_boolean = true;
            return true;
        }
        catch (...)
        {
            m_error_message.append("Error: Invalid numeric argument passed to "
                                   + c + ": " + input + "!\n");
            return false;
        }
    }

    template <typename Type>
    inline bool set_numeric(const std::string& input, const std::string& c,
                            Type& target)
    {
        bool dumy;
        return set_numeric<Type>(input, c, target, dumy);
    }

    inline void load_string_vector(const std::string& input,
                                   const std::string& c,
                                   std::vector<std::string>& target)
    {
        if (input.empty()) return;
        if (m_parameter_log.find(c) == m_parameter_log.end())
        { m_parameter_log[c] = input; }
        else
        {
            m_parameter_log[c] = "," + input;
        }
        if (!input.empty() && input.back() == ',')
        {
            m_error_message.append(
                "Warning: , detected at end of input: " + input
                + ". Have you accidentally included space in "
                  "your input? (Space is not allowed)\n");
        }
        std::vector<std::string> token = misc::split(input, ",");
        target.insert(target.end(), token.begin(), token.end());
    }

    inline bool set_memory(const std::string& input)
    {
        bool error = false;
        m_memory = set_distance(input, "memory", 1024 * 1024, error, true);
        return !error;
    }

    inline bool set_missing(const std::string& in)
    {
        std::string input = in;
        check_duplicate("missing");
        std::transform(input.begin(), input.end(), input.begin(), ::toupper);
        switch (input.at(0))
        {
        case 'C':
            input = "centre";
            m_prs_info.missing_score = MISSING_SCORE::CENTER;
            break;
        case 'M':
            input = "mean_impute";
            m_prs_info.missing_score = MISSING_SCORE::MEAN_IMPUTE;
            break;
        case 'S':
            input = "set_zero";
            m_prs_info.missing_score = MISSING_SCORE::SET_ZERO;
            break;
        case 'I':
            input = "impute_control";
            m_prs_info.missing_score = MISSING_SCORE::IMPUTE_CONTROL;
            break;
        default:
            m_error_message.append(
                "Error: Unrecognized Missing handling method: " + input
                + "!\n");
            return false;
        }
        m_parameter_log["missing"] = input;
        return true;
    }

    inline bool set_model(const std::string& in)
    {
        std::string input = in;
        std::transform(input.begin(), input.end(), input.begin(), ::toupper);
        check_duplicate("model");
        switch (input.at(0))
        {
        case 'A':
            input = "add";
            m_prs_info.genetic_model = MODEL::ADDITIVE;
            break;
        case 'D':
            input = "dom";
            m_prs_info.genetic_model = MODEL::DOMINANT;
            break;
        case 'R':
            input = "rec";
            m_prs_info.genetic_model = MODEL::RECESSIVE;
            break;
        case 'H':
            input = "het";
            m_prs_info.genetic_model = MODEL::HETEROZYGOUS;
            break;
        default:
            m_error_message.append("Error: Unrecognized model: " + input
                                   + "!\n");
            return false;
        }
        m_parameter_log["model"] = input;
        return true;
    }


    /*!
     * \brief Function parsing string into MISSING enum
     * \param in the input
     * \param message the parameter storage
     * \param error_message the error message storage
     * \return true if successfully parse the input into MISSING enum
     */
    inline bool set_score(const std::string& in)
    {
        std::string input = in;
        std::transform(input.begin(), input.end(), input.begin(), ::toupper);
        check_duplicate("score");
        if (input.length() < 2)
        {
            m_error_message.append("Error: Invalid scoring method: " + input
                                   + "!\n");
        }
        if (input == "AVG")
        {
            input = "avg";
            m_prs_info.scoring_method = SCORING::AVERAGE;
        }
        else if (input == "STD")
        {
            input = "std";
            m_prs_info.scoring_method = SCORING::STANDARDIZE;
        }
        else if (input == "SUM")
        {
            input = "sum";
            m_prs_info.scoring_method = SCORING::SUM;
        }
        else if (input == "CON_STD")
        {
            input = "con_std";
            m_prs_info.scoring_method = SCORING::CONTROL_STD;
        }
        else
        {
            m_error_message.append(
                "Error: Unrecognized scoring method: " + input + "!\n");
            return false;
        }
        m_parameter_log["score"] = input;
        return true;
    }

    bool in_file(const std::vector<std::string>& column_names,
                 const size_t index, const std::string& warning,
                 bool no_default, bool case_sensitive = true,
                 bool print_error = true)
    {
        if ((no_default && !static_cast<bool>(m_base_info.has_column[index]))
            || m_base_info.column_name[index].empty())
        { return false; }
        size_t col_index;
        bool has_col = index_check(m_base_info.column_name[index], column_names,
                                   col_index, case_sensitive);
        if (has_col) { m_base_info.column_index[index] = col_index; }
        else if (m_base_info.has_column[index] && print_error)
        {
            m_error_message.append(warning + ": "
                                   + m_base_info.column_name[index]
                                   + " not found in base file\n");
        }
        m_base_info.has_column[index] = has_col;
        return has_col;
    }


    inline void set_base_info_threshold(const std::vector<std::string>& ref,
                                        bool& error)
    {
        const std::vector<std::string> info =
            misc::split(m_base_info.column_name[+BASE_INDEX::INFO], ",");
        const bool has_input = m_base_info.has_column[+BASE_INDEX::INFO];

        size_t index;
        const bool found = index_check(info[0], ref, index);
        if (found) m_base_info.column_index[+BASE_INDEX::INFO] = index;
        m_base_info.has_column[+BASE_INDEX::INFO] = found;
        if (!found)
        {
            if (has_input)
            {
                m_error_message.append("Warning: INFO field not found in base "
                                       "file, will ignore INFO filtering\n");
            }
            return;
        }
        else if (info.size() != 2) // assume default always valid
        {
            m_error_message.append("Error: Invalid format of "
                                   "--base-info. Should be "
                                   "ColName,Threshold.\n");
            error = true;
            return;
        }
        try
        {
            m_base_filter.info_score = misc::convert<double>(info[1]);
            if ((m_base_filter.info_score < 0 || m_base_filter.info_score > 1))
            {
                if (has_input)
                {
                    m_error_message.append("Error: Base INFO threshold "
                                           "must be within 0 and 1!\n");
                    error = true;
                }
            }
        }
        catch (...)
        {
            if (has_input)
            {
                m_error_message.append(
                    "Error: Invalid argument passed to --base-info: "
                    + m_base_info.column_name[+BASE_INDEX::INFO]
                    + "! Second argument must be numeric\n");
                error = true;
            }
        }
        return;
    }

    // return true if valid
    inline void set_base_maf_filter(const std::vector<std::string>& ref,
                                    bool& error)
    {
        const std::string maf_error =
            "Error: Invalid format of --base-maf. "
            "Should be ColName,Threshold."
            "or ColName,Threshold:ColName,Threshold.\n";
        std::vector<std::string> case_control =
            misc::split(m_base_info.column_name[+BASE_INDEX::MAF], ":");
        const bool print_error = m_base_info.has_column[+BASE_INDEX::MAF];
        if (case_control.size() > 2)
        {
            if (print_error)
            {
                m_error_message.append(maf_error);
                error = true;
            }
            return;
        }
        bool control = true;
        bool found = false;
        BASE_INDEX maf_index = BASE_INDEX::MAF;
        size_t index = 0;
        for (auto&& maf : case_control)
        {
            std::vector<std::string> detail = misc::split(maf, ",");
            found = index_check(detail[0], ref, index);
            m_base_info.has_column[+maf_index] = found;
            if (found) m_base_info.column_index[+maf_index] = index;
            if (!found)
            {
                if (print_error)
                {
                    m_error_message.append(
                        "Warning: MAF field not found in base file. "
                        "Will not perform MAF filtering on the base file\n");
                    error = true;
                }
                return;
            }
            double cur_maf;
            try
            {
                cur_maf = misc::convert<double>(detail[1]);
                if ((cur_maf < 0 || cur_maf > 1))
                {
                    if (print_error)
                    {
                        m_error_message.append("Error: Base MAF threshold must "
                                               "be within 0 and 1!\n");
                        error = true;
                    }
                    return;
                }
            }
            catch (...)
            {
                if (print_error)
                {
                    m_error_message.append(
                        "Error: Invalid argument passed to --base-maf: "
                        + m_base_info.column_name[+BASE_INDEX::MAF]
                        + "! Threshold must be numeric\n");
                    error = true;
                }
                return;
            }
            if (control) { m_base_filter.maf = cur_maf; }
            else
            {
                m_base_filter.maf_case = cur_maf;
            }
            control = false;
        }
        return;
    }
    /*!
     * \brief Get the column index based on file header and the input string
     * \param target is the input string
     * \param ref is the vector containing the column header
     * \return the index to the column containing the column name. -1 if it is
     * not found
     */
    inline bool index_check(const std::string& target,
                            const std::vector<std::string>& ref, size_t& index,
                            bool case_sensitive = true) const
    {
        std::string tmp = "";
        for (size_t i = 0; i < ref.size(); ++i)
        {
            tmp = ref[i];
            if (!case_sensitive)
            {
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
            }
            if (target == tmp)
            {
                index = i;
                return true;
            }
        }
        return false;
    }

    bool get_statistic_column(const std::vector<std::string>& column_names);
    bool base_check();
};

#endif // COMMANDER_H

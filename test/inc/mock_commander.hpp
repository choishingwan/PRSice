#ifndef MOCK_COMMANDER_H
#define MOCK_COMMANDER_H
#include "commander.hpp"

class mockCommander : public Commander
{
public:
    static std::vector<std::string>
    transform_covariate(const std::string& cov_in)
    {
        return Commander::transform_covariate(cov_in);
    }
    bool check_parse_unit_value(const std::string& input, const std::string& c,
                                const size_t default_power, size_t& target,
                                bool memory = false)
    {
        return parse_unit_value(input, c, default_power, target, memory);
    }

    static bool find_first_end_wrapper(const std::string_view& cov,
                                       const size_t idx, size_t& res)
    {
        try
        {
            res = find_first_end(cov, idx);
            return true;
        }
        catch (const std::runtime_error&)
        {
            return false;
        }
    }
    static bool parse_range_wrapper(std::string_view cov,
                                    std::vector<size_t>& res)
    {
        try
        {
            res = parse_range(cov);
            return true;
        }
        catch (std::runtime_error&)
        {
            return false;
        }
    }
    static bool get_range_wrapper(std::string_view cov, size_t start,
                                  size_t end, std::vector<size_t>& res)
    {
        try
        {
            res = get_range(cov, start, end);
            return true;
        }
        catch (std::runtime_error&)
        {
            return false;
        }
    }
    static bool
    update_covariate_ranges_wrapper(std::vector<std::string>& result,
                                    std::vector<size_t> ranges)
    {
        try
        {
            update_covariate_range(ranges, result);
            return true;
        }
        catch (...)
        {
            return false;
        }
    }
    static bool transform_wrapper(const std::string& str,
                                  std::vector<std::string>& result)
    {
        try
        {
            result = transform_covariate(str);
            return true;
        }
        catch (...)
        {
            return false;
        }
    }
    bool parse_command_wrapper(const std::string& command)
    {
        bool early_terminate = false;
        return parse_command_wrapper(command, early_terminate);
    }
    void prepare_header_cov_check_wrapper(
        const std::vector<std::string>& cov_header,
        std::unordered_map<std::string, size_t>& ref_index,
        std::unordered_set<std::string>& included)
    {
        prepare_header_cov_check(cov_header, ref_index, included);
    }
    bool parse_command_wrapper(const std::string& command,
                               bool& early_terminate)
    {
        Reporter reporter("log", 60, true);

        std::vector<std::string> argv_str ;
        if(!command.empty()) argv_str = misc::split("PRSice " + command);
        else argv_str = {"PRSice "};
        std::vector<char*> cstrings;
        cstrings.reserve(argv_str.size());
        for (size_t i = 0; i < argv_str.size(); ++i)
        { cstrings.push_back(const_cast<char*>(argv_str[i].c_str())); }
        int argc = static_cast<int>(argv_str.size());
        try
        {
            early_terminate = false;
            // return false if error
            return !init(argc, &cstrings[0], early_terminate, reporter);
        }
        catch (...)
        {
            // error = false
            return false;
        }
    }

    bool no_default() const { return m_user_no_default; }
    bool target_check_wrapper() { return target_check(); }
    bool prsice_check_wrapper() { return prsice_check(); }
    bool clump_check_wrapper() { return clump_check(); }
    bool ref_check_wrapper() { return ref_check(); }
    bool misc_check_wrapper() { return misc_check(); }
    bool filter_check_wrapper() { return filter_check(); }
    bool prset_check_wrapper() { return prset_check(); }
    bool base_check_wrapper()
    {
        try
        {
            return base_check();
        }
        catch (const std::runtime_error&)
        {
            return false;
        }
    }
    bool base_column_check_wrapper(std::vector<std::string>& column_names)
    {
        return base_column_check(column_names);
    }
    bool pheno_check_wrapper()
    {
        m_ran_base_check = true;
        return pheno_check();
    }
    std::string get_error() const { return m_error_message; }
    int32_t max_thread() { return maximum_thread(); }
    auto get_cov_names_wrap() { return get_cov_names(); }
    size_t
    find_cov_idx_wrap(const std::unordered_set<std::string>& included,
                      const std::unordered_map<std::string, size_t>& ref_index,
                      std::string& missing)
    {
        return find_cov_idx(included, ref_index, missing);
    }
    void reorganize_cov_name_wrap(const std::vector<std::string>& cov_header)
    {
        reorganize_cov_name(cov_header);
    }
    bool process_factor_cov_wrap(
        const std::unordered_set<std::string>& included,
        const std::unordered_map<std::string, size_t>& ref_index,
        const std::unordered_set<std::string>& ori_input)
    {
        return process_factor_cov(included, ref_index, ori_input);
    }
};

#endif // MOCK_COMMANDER_H

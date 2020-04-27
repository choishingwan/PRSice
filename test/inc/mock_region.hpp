#ifndef MOCK_REGION_HPP
#define MOCK_REGION_HPP
#include "genotype.hpp"
#include "region.hpp"

class mock_region : public Region{
public:
    bool test_in_region(const std::string_view & in, const std::vector<std::string> & feature){
        return in_feature(in, feature);
    }
    static bool test_valid_strand(std::string_view input){
        return valid_strand(input);
    }
    std::tuple<std::string, std::string, bool>
    test_get_set_name(const std::string& input)
    {
        return get_set_name(input);
    }
};

#endif // MOCK_REGION_HPP

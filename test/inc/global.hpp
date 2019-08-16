#ifndef GLOBAL_TEST_H
#define GLOBAL_TEST_H
#include "region.hpp"
#include <string>
#include <vector>
extern std::string path;

class FAKE_REGION : public Region
{
public:
    FAKE_REGION(const std::vector<std::string>& bed,
                const std::vector<std::string>& feature,
                const std::vector<std::string>& msigdb,
                const std::vector<std::string>& snp_set,
                const std::string& background, const std::string& gtf,
                const size_t window_5, const size_t window_3,
                const bool genome_wide_background, Reporter* reporter)
    {
        m_bed = bed;
        m_feature = feature;
        m_msigdb = msigdb;
        m_snp_set = snp_set;
        m_background = background;
        m_gtf = gtf;
        m_window_5 = window_5;
        m_window_3 = window_3;
        m_genome_wide_background = genome_wide_background;
        m_reporter = reporter;
    }
};
#endif

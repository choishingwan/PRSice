#ifndef MOCK_GENOTYPE_H
#define MOCK_GENOTYPE_H

#include "genotype.hpp"
#include <memory>
#include "reporter.hpp"

class mockGenotype : public Genotype{
public:

   std::string test_initialize(const GenoFile& geno, const Phenotype& pheno,
                               const std::string& delim, const std::string& type,
                               Reporter* reporter){
       return initialize(geno, pheno, delim, type, reporter);
   }
   std::vector<std::string>
   test_load_genotype_prefix(std::unique_ptr<std::istringstream> in){
       return load_genotype_prefix(std::move(in));
   }
   bool test_ambiguous(const std::string& a, const std::string&b){
       return ambiguous(a, b);
   }
   size_t test_get_rs_column(const std::string& input){
       if(m_reporter == nullptr){
           Reporter reporter("log", 60, true);
           m_reporter = &reporter;
       }
       return get_rs_column(input);
   }
   std::vector<std::string>
   test_load_snp_list(std::unique_ptr<std::istream> input){
       auto res =  load_snp_list(std::move(input));
       std::vector<std::string> result;
       result.insert(result.end(), res.begin(), res.end());
       return result;
   }
   std::string delim()const {return m_delim;}
   std::vector<std::string> genotype_file_names() const { return m_genotype_file_names;}
   std::string keep_file() const { return m_keep_file;}
   std::string remove_file() const { return m_remove_file;}
   std::string sample_file() const { return m_sample_file;}
   bool ignore_fid() const { return m_ignore_fid;}
   // helper
   void set_reporter(Reporter *reporter){
       m_reporter = reporter;
   }
};

#endif // MOCK_GENOTYPE_H

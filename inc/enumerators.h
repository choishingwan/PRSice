#ifndef ENUMERATORS_H
#define ENUMERATORS_H

// From http://stackoverflow.com/a/12927952/1441789
#include <type_traits>
template <typename e>
struct enumeration_traits;

struct enumeration_trait_indexing
{
    static constexpr bool does_index = true;
};

template <typename e>
constexpr typename std::enable_if<enumeration_traits<e>::does_index,
                                  typename std::underlying_type<e>::type>::type
operator+(e val)
{
    return static_cast<typename std::underlying_type<e>::type>(val);
}

template <typename e>
typename std::enable_if<enumeration_traits<e>::does_index, e&>::type
operator++(e& val)
{
    return val = static_cast<e>(+val + 1);
}
// END

enum class BIM
{
    CHR = 0,
    RS,
    CM,
    BP,
    A1,
    A2,
    MAX
};
enum class BASE_INDEX
{
    CHR = 0,
    NONEFFECT,
    BP,
    SE,
    INFO,
    MAF,
    MAF_CASE,
    EFFECT,
    RS,
    P,
    STAT,
    MAX
};

enum class FAM
{
    FID = 0,
    IID,
    FATHER,
    MOTHER,
    SEX,
    PHENOTYPE,
    MAX
};
enum class GTF
{
    CHR = 0,
    SOURCE,
    FEATURE,
    START,
    END,
    SCORE,
    STRAND,
    FRAME,
    ATTRIBUTE,
    MAX
};

enum class BED
{
    CHR = 0,
    START,
    END,
    NAME,
    SCORE,
    STRAND,
    MAX
};
enum class MODEL
{
    ADDITIVE = 0,
    DOMINANT,
    RECESSIVE,
    HETEROZYGOUS
};

enum class MISSING_SCORE
{
    MEAN_IMPUTE = 0,
    IMPUTE_CONTROL,
    SET_ZERO,
    CENTER
};

enum class SCORING
{
    AVERAGE = 0,
    STANDARDIZE,
    CONTROL_STD,
    SUM
};

enum class FILTER_COUNT
{
    DUP_SNP = 0,
    NUM_LINE,
    P_EXCLUDED,
    SELECT,
    REGION,
    AMBIG,
    HAPLOID,
    NOT_CONVERT,
    NEGATIVE,
    INFO,
    CHR,
    MAF,
    MAX
};

template <>
struct enumeration_traits<BASE_INDEX> : enumeration_trait_indexing
{
};
template <>
struct enumeration_traits<FILTER_COUNT> : enumeration_trait_indexing
{
};
template <>
struct enumeration_traits<GTF> : enumeration_trait_indexing
{
};
template <>
struct enumeration_traits<MODEL> : enumeration_trait_indexing
{
};
template <>
struct enumeration_traits<FAM> : enumeration_trait_indexing
{
};
template <>
struct enumeration_traits<BIM> : enumeration_trait_indexing
{
};
template <>
struct enumeration_traits<BED> : enumeration_trait_indexing
{
};
#endif // ENUMERATORS_H

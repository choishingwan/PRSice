#ifndef __PLINK_COMMON_H__
#define __PLINK_COMMON_H__

// This file is part of PLINK 1.90, copyright (C) 2005-2017 Shaun Purcell,
// Christopher Chang.
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


// Resources needed across all plink modules.

#define _FILE_OFFSET_BITS 64
#include <inttypes.h>
#include <iostream>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#define NDEBUG // don't know wha that's for
#include <assert.h>

// Uncomment this to build this without CBLAS/CLAPACK.
// #define NOLAPACK

// Uncomment this to prevent all unstable features from being accessible from
// the command line.
// #define STABLE_BUILD

#define SPECIES_HUMAN 0
#define SPECIES_COW 1
#define SPECIES_DOG 2
#define SPECIES_HORSE 3
#define SPECIES_MOUSE 4
#define SPECIES_RICE 5
#define SPECIES_SHEEP 6
#define SPECIES_UNKNOWN 7
#define SPECIES_DEFAULT SPECIES_HUMAN

#define PROG_NAME_STR "plink"
#define PROG_NAME_CAPS "PLINK"

#ifdef _WIN32
// needed for MEMORYSTATUSEX
#ifndef _WIN64
#define WINVER 0x0500
#endif
#else // Unix
#include <sys/stat.h>
#endif

#ifndef HAVE_NULLPTR
#ifndef __cplusplus
#define nullptr NULL
#else
#if __cplusplus <= 199711L
#ifndef nullptr
#define nullptr NULL
#endif
#endif
#endif
#endif

#ifdef _WIN32
#define PRId64 "I64d"
#define PRIu64 "I64u"
#define fseeko fseeko64
#define ftello ftello64
#include <process.h>
#include <windows.h>
#define pthread_t HANDLE
#define THREAD_RET_TYPE unsigned __stdcall
#define THREAD_RETURN return 0
#define EOLN_STR "\r\n"
#define FOPEN_RB "rb"
#define FOPEN_WB "wb"
#ifdef _WIN64
#define getc_unlocked _fgetc_nolock
#define putc_unlocked _fputc_nolock
#else
#define getc_unlocked getc
#define putc_unlocked putc
#endif
#define uint64_t unsigned long long
#define int64_t long long
#else
#include <pthread.h>
#define THREAD_RET_TYPE void*
#define THREAD_RETURN return nullptr
#ifdef __cplusplus
#ifndef PRId64
#define PRId64 "lld"
#endif
#endif
#define EOLN_STR "\n"
#define FOPEN_RB "r"
#define FOPEN_WB "w"
#ifndef __APPLE__
// argh
// not sure what the right threshold actually is, but this works for now
// (may break on gcc <3.0?  but that shouldn't matter anymore)
// tried defining GCC_VERSION, but that didn't always work
#if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 8)
#define uint64_t unsigned long long
#define int64_t long long
#endif
#endif
#endif

#ifdef _WIN64
#define __LP64__
#define CTZLU __builtin_ctzll
#define CLZLU __builtin_clzll
#else
#define CTZLU __builtin_ctzl
#define CLZLU __builtin_clzl
#ifndef __LP64__
// attempt to patch GCC 6 build failure
#if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 8)
#ifndef uintptr_t
#define uintptr_t unsigned long
#endif
#ifndef intptr_t
#define intptr_t long
#endif
#endif
#endif
#endif

#ifdef __cplusplus
#include <algorithm>
#define HEADER_INLINE inline
#else
#define HEADER_INLINE static inline
#endif

// It would be useful to disable compilation on big-endian platforms, but I
// don't see a decent portable way to do this (see e.g. discussion at
// http://esr.ibiblio.org/?p=5095 ).

#ifdef __LP64__
#ifndef __SSE2__
// It's obviously possible to support this by writing 64-bit non-SSE2 code
// shadowing each SSE2 intrinsic, but this almost certainly isn't worth the
// development/testing effort until regular PLINK 2.0 development is
// complete.  No researcher has ever asked me for this feature.
#error \
    "64-bit builds currently require SSE2.  Try producing a 32-bit build instead."
#endif
#include <emmintrin.h>

#define VECFTYPE __m128
#define VECITYPE __m128i
#define VECDTYPE __m128d

// useful because of its bitwise complement: ~ZEROLU is a word with all 1
// bits, while ~0 is always 32 1 bits.
#define ZEROLU 0LLU

// mainly useful for bitshifts: (ONELU << 32) works in 64-bit builds, while
// (1 << 32) is undefined.  also used to cast some numbers/expressions to
// uintptr_t (e.g. multiplying an int constant by ONELU widens it to 64 bits
// only in 64-bit builds; note that 1LU fails on Win64 while 1LLU doesn't do
// the right thing for 32-bit builds).
#define ONELU 1LLU

#ifdef _WIN32 // i.e. Win64

#ifndef PRIuPTR
#define PRIuPTR PRIu64
#endif
#ifndef PRIdPTR
#define PRIdPTR PRId64
#endif
#define PRIxPTR2 "016I64x"

#else // not _WIN32

#ifndef PRIuPTR
#define PRIuPTR "lu"
#endif
#ifndef PRIdPTR
#define PRIdPTR "ld"
#endif
#define PRIxPTR2 "016lx"

#endif // Win64

#define VEC_BYTES 16

#else // not __LP64__

#define ZEROLU 0LU
#define ONELU 1LU
#ifndef PRIuPTR
#define PRIuPTR "lu"
#endif
#ifndef PRIdPTR
#define PRIdPTR "ld"
#endif
#define PRIxPTR2 "08lx"

// todo: update code so this still works when reduced to 4
#define VEC_BYTES 8

#endif // __LP64__

// use constexpr for these as soon as compiler support is available on all
// platforms
#define FIVEMASK ((~ZEROLU) / 3)
#define AAAAMASK (FIVEMASK * 2)

#define VEC_BYTES_M1 (VEC_BYTES - 1)
#define VEC_BITS (VEC_BYTES * 8)
#define VEC_BITS_M1 (VEC_BITS - 1)


// 64MB of non-workspace memory guaranteed for now.
// Currently also serves as the maximum allele length.
#define NON_BIGSTACK_MIN 67108864

#define PI 3.1415926535897932
#define RECIP_2_32 0.00000000023283064365386962890625
#define RECIP_2_53 0.00000000000000011102230246251565404236316680908203125
// floating point comparison-to-nonzero tolerance, currently 2^{-30}
#define EPSILON 0.000000000931322574615478515625
// less tolerant versions (2^{-35}, 2^{-44}) for some exact calculations
#define SMALLISH_EPSILON 0.00000000002910383045673370361328125
#define SMALL_EPSILON 0.00000000000005684341886080801486968994140625
// at least sqrt(SMALL_EPSILON)
#define BIG_EPSILON 0.000000476837158203125
// 53-bit double precision limit
#define DOUBLE_PREC_LIMIT \
    0.00000000000000011102230246251565404236316680908203125
#define TWO_63 9223372036854775808.0
#define SQRT_HALF 0.70710678118654746

// 2^{-83} bias to give exact tests maximum ability to determine tiny p-values.
// (~2^{-53} is necessary to take advantage of denormalized small numbers, then
// allow tail sum to be up to 2^30.)
#define EXACT_TEST_BIAS \
    0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125

// occasionally used as an infinity substitute that avoids the 32-bit Windows
// performance penalty
// can import from limits.h, we don't bother to include that for now
#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623157e308
#endif

// not quite the same as FLT_MAX since it's a double-precision constant
#define FLT_MAXD 3.4028234663852886e38

#define RET_SUCCESS 0
#define RET_NOMEM 1
#define RET_OPEN_FAIL 2
#define RET_INVALID_FORMAT 3
#define RET_CALC_NOT_YET_SUPPORTED 4
#define RET_INVALID_CMDLINE 5
#define RET_WRITE_FAIL 6
#define RET_READ_FAIL 7
#define RET_HELP 8
#define RET_THREAD_CREATE_FAIL 9
#define RET_ALLELE_MISMATCH 10
#define RET_NULL_CALC 11
#define RET_ALL_SAMPLES_EXCLUDED 12
#define RET_ALL_MARKERS_EXCLUDED 13
#define RET_NETWORK 14
#define LOAD_PHENO_LAST_COL 127

// for 2.0 -> 1.9 backports
#define RET_MALFORMED_INPUT RET_INVALID_FORMAT

#define MISC_AFFECTION_01 1LLU
#define MISC_NONFOUNDERS 2LLU
#define MISC_MAF_SUCC 4LLU
#define MISC_FREQ_COUNTS 8LLU
#define MISC_FREQ_CC 0x10LLU
#define MISC_FREQX 0x20LLU
#define MISC_KEEP_ALLELE_ORDER 0x40LLU
#define MISC_SET_HH_MISSING 0x80LLU
#define MISC_SET_MIXED_MT_MISSING 0x100LLU
#define MISC_KEEP_AUTOCONV 0x200LLU
#define MISC_LOAD_CLUSTER_KEEP_NA 0x400LLU
#define MISC_WRITE_CLUSTER_OMIT_UNASSIGNED 0x800LLU
#define MISC_ALLOW_EXTRA_CHROMS 0x1000LLU
#define MISC_MAKE_FOUNDERS_REQUIRE_2_MISSING 0x2000LLU
#define MISC_MAKE_FOUNDERS_FIRST 0x4000LLU
#define MISC_LASSO_REPORT_ZEROES 0x8000LLU
#define MISC_LASSO_SELECT_COVARS 0x10000LLU
#define MISC_DOUBLE_ID 0x20000LLU
#define MISC_BIALLELIC_ONLY 0x40000LLU
#define MISC_BIALLELIC_ONLY_STRICT 0x80000LLU
#define MISC_BIALLELIC_ONLY_LIST 0x100000LLU
#define MISC_VCF_FILTER 0x200000LLU
#define MISC_GPLINK 0x400000LLU
#define MISC_SNPS_ONLY_JUST_ACGT 0x800000LLU
#define MISC_IMPUTE_SEX 0x1000000LLU
#define MISC_OXFORD_SNPID_CHR 0x2000000LLU
#define MISC_EXTRACT_RANGE 0x4000000LLU
#define MISC_EXCLUDE_RANGE 0x8000000LLU
#define MISC_MERGEX 0x10000000LLU
#define MISC_SET_ME_MISSING 0x20000000LLU
#define MISC_SEXCHECK_YCOUNT 0x40000000LLU
#define MISC_SEXCHECK_YONLY 0x80000000LLU
#define MISC_FAMILY_CLUSTERS 0x100000000LLU
#define MISC_FILL_MISSING_A2 0x200000000LLU
#define MISC_HET_SMALL_SAMPLE 0x400000000LLU
#define MISC_FST_CC 0x800000000LLU
#define MISC_SPLIT_MERGE_NOFAIL 0x1000000000LLU
#define MISC_REAL_REF_ALLELES 0x2000000000LLU
#define MISC_RPLUGIN_DEBUG 0x4000000000LLU
#define MISC_MISSING_GZ 0x8000000000LLU
#define MISC_FREQ_GZ 0x10000000000LLU
#define MISC_HET_GZ 0x20000000000LLU
#define MISC_ALLOW_NO_SAMPLES 0x40000000000LLU
#define MISC_ALLOW_NO_VARS 0x80000000000LLU
#define MISC_VCF_REQUIRE_GT 0x100000000000LLU

// assume for now that .bed must always be accompanied by both .bim and .fam
#define FILTER_ALL_REQ 1LLU
#define FILTER_BIM_REQ 2LLU
#define FILTER_FAM_REQ 4LLU

// ok with --dosage + --map, but not --dosage by itself
#define FILTER_DOSAGEMAP 8LLU

#define FILTER_NODOSAGE 0x10LLU
#define FILTER_NOCNV 0x20LLU
#define FILTER_EXCLUDE_MARKERNAME_SNP 0x40LLU
#define FILTER_BINARY_CASES 0x80LLU
#define FILTER_BINARY_CONTROLS 0x100LLU
#define FILTER_BINARY_FEMALES 0x200LLU
#define FILTER_BINARY_MALES 0x400LLU
#define FILTER_BINARY_FOUNDERS 0x800LLU
#define FILTER_BINARY_NONFOUNDERS 0x1000LLU
#define FILTER_MAKE_FOUNDERS 0x2000LLU
#define FILTER_PRUNE 0x4000LLU
#define FILTER_SNPS_ONLY 0x8000LLU
#define FILTER_TAIL_PHENO 0x10000LLU
#define FILTER_ZERO_CMS 0x20000LLU

#define CALC_RELATIONSHIP 1LLU
#define CALC_IBC 2LLU
#define CALC_DISTANCE 4LLU
#define CALC_PLINK1_DISTANCE_MATRIX 8LLU
#define CALC_PLINK1_IBS_MATRIX 0x10LLU
#define CALC_GDISTANCE_MASK 0x1cLLU
#define CALC_GROUPDIST 0x20LLU
#define CALC_REGRESS_DISTANCE 0x40LLU
#define CALC_UNRELATED_HERITABILITY 0x80LLU
#define CALC_FREQ 0x100LLU
#define CALC_REL_CUTOFF 0x200LLU
#define CALC_WRITE_SNPLIST 0x400LLU
#define CALC_LIST_23_INDELS 0x800LLU
#define CALC_GENOME 0x1000LLU
#define CALC_REGRESS_REL 0x2000LLU
#define CALC_LD_PRUNE 0x4000LLU
#define CALC_REGRESS_PCS 0x8000LLU
#define CALC_REGRESS_PCS_DISTANCE 0x10000LLU
#define CALC_MAKE_BED 0x20000LLU
#define CALC_RECODE 0x40000LLU
#define CALC_MERGE 0x80000LLU
#define CALC_WRITE_COVAR 0x100000LLU
#define CALC_WRITE_CLUSTER 0x200000LLU
#define CALC_MODEL 0x400000LLU
#define CALC_HARDY 0x800000LLU
#define CALC_GXE 0x1000000LLU
#define CALC_IBS_TEST 0x2000000LLU
#define CALC_CLUSTER 0x4000000LLU
#define CALC_HOMOZYG 0x8000000LLU
#define CALC_NEIGHBOR 0x10000000LLU
#define CALC_GLM 0x20000000LLU
#define CALC_MISSING_REPORT 0x40000000LLU
#define CALC_CMH 0x80000000LLU
#define CALC_HOMOG 0x100000000LLU
#define CALC_LASSO 0x200000000LLU
#define CALC_LASSO_LAMBDA 0x400000000LLU
#define CALC_WRITE_SET 0x800000000LLU
#define CALC_LD 0x1000000000LLU
#define CALC_EPI 0x2000000000LLU
#define CALC_TESTMISS 0x4000000000LLU
#define CALC_TESTMISHAP 0x8000000000LLU
#define CALC_SEXCHECK 0x10000000000LLU
#define CALC_CLUMP 0x20000000000LLU
#define CALC_PCA 0x40000000000LLU
#define CALC_BLOCKS 0x80000000000LLU
#define CALC_SCORE 0x100000000000LLU
#define CALC_MENDEL 0x200000000000LLU
#define CALC_HET 0x400000000000LLU
#define CALC_FLIPSCAN 0x800000000000LLU
#define CALC_TDT 0x1000000000000LLU
#define CALC_MAKE_PERM_PHENO 0x2000000000000LLU
#define CALC_QFAM 0x4000000000000LLU
#define CALC_FST 0x8000000000000LLU
#define CALC_SHOW_TAGS 0x10000000000000LLU
#define CALC_MAKE_BIM 0x20000000000000LLU
#define CALC_MAKE_FAM 0x40000000000000LLU
#define CALC_WRITE_VAR_RANGES 0x80000000000000LLU
#define CALC_DUPVAR 0x100000000000000LLU
#define CALC_RPLUGIN 0x200000000000000LLU
#define CALC_DFAM 0x400000000000000LLU
#define CALC_ONLY_BIM                                            \
    (CALC_WRITE_SET | CALC_WRITE_SNPLIST | CALC_WRITE_VAR_RANGES \
     | CALC_LIST_23_INDELS | CALC_MAKE_BIM | CALC_DUPVAR)
#define CALC_ONLY_FAM (CALC_MAKE_PERM_PHENO | CALC_WRITE_COVAR | CALC_MAKE_FAM)
// only room for 6 more basic commands before we need to switch from a single
// uint64_t to uintptr_t*/is_set()/etc.

// necessary to patch heterozygous haploids/female Y chromosome genotypes
// during loading?
#define XMHH_EXISTS 1
#define Y_FIX_NEEDED 2
#define NXMHH_EXISTS 4

#define ALLELE_RECODE 1
#define ALLELE_RECODE_MULTICHAR 2
#define ALLELE_RECODE_ACGT 4

// 0 = non-explicit error
#define VCF_HALF_CALL_ERROR 1
#define VCF_HALF_CALL_MISSING 2
#define VCF_HALF_CALL_HAPLOID 3
#define VCF_HALF_CALL_REFERENCE 4

#define M23_MALE 1
#define M23_FEMALE 2
#define M23_FORCE_MISSING_SEX 4
#define M23_SEX 7

#define MARKER_CMS_OPTIONAL 1
#define MARKER_CMS_FORCED 2

#define UNSORTED_CHROM 1
#define UNSORTED_BP 2
#define UNSORTED_CM 4
#define UNSORTED_SPLIT_CHROM 8

#define ALLOW_NO_SEX 1
#define MUST_HAVE_SEX 2

#define LGEN_REFERENCE 1
#define LGEN_ALLELE_COUNT 2

#define PHENO_ALL 1
#define PHENO_MERGE 2

#define FAM_COL_1 1
#define FAM_COL_34 2
#define FAM_COL_5 4
#define FAM_COL_6 8
#define FAM_COL_13456 15

#define COVAR_KEEP_PHENO_ON_MISSING_COV 1
#define COVAR_NAME 2
#define COVAR_NUMBER 4
#define COVAR_NO_CONST 8
#define COVAR_ALLOW_NONE 0x10

#define DISTANCE_SQ 1
#define DISTANCE_SQ0 2
#define DISTANCE_TRI 3
#define DISTANCE_SHAPEMASK 3
#define DISTANCE_GZ 4
#define DISTANCE_BIN 8
#define DISTANCE_BIN4 0x10
#define DISTANCE_IBS 0x20
#define DISTANCE_1_MINUS_IBS 0x40
#define DISTANCE_ALCT 0x80
#define DISTANCE_TYPEMASK 0xe0
#define DISTANCE_FLAT_MISSING 0x100
#define DISTANCE_CLUSTER 0x200
#define DISTANCE_WTS_NOHEADER 0x400

#define RECODE_01 1
#define RECODE_12 2
#define RECODE_TAB 4
#define RECODE_DELIMX 8
#define RECODE_23 0x10
#define RECODE_A 0x20
#define RECODE_A_TRANSPOSE 0x40
#define RECODE_AD 0x80
#define RECODE_BEAGLE 0x100
#define RECODE_BEAGLE_NOMAP 0x200
#define RECODE_BIMBAM 0x400
#define RECODE_BIMBAM_1CHR 0x800
#define RECODE_COMPOUND 0x1000
#define RECODE_FASTPHASE 0x2000
#define RECODE_FASTPHASE_1CHR 0x4000
#define RECODE_HV 0x8000
#define RECODE_HV_1CHR 0x10000
#define RECODE_LGEN 0x20000
#define RECODE_LGEN_REF 0x40000
#define RECODE_LIST 0x80000
#define RECODE_OXFORD 0x100000
#define RECODE_RLIST 0x200000
#define RECODE_STRUCTURE 0x400000
#define RECODE_TRANSPOSE 0x800000
#define RECODE_PED 0x1000000
#define RECODE_VCF 0x2000000
#define RECODE_TYPEMASK 0x3fffff0
#define RECODE_FID 0x4000000
#define RECODE_IID 0x8000000
#define RECODE_INCLUDE_ALT 0x10000000
#define RECODE_BGZ 0x20000000
#define RECODE_GEN_GZ 0x40000000
#define RECODE_OMIT_NONMALE_Y 0x80000000U

#define GENOME_OUTPUT_GZ 1
#define GENOME_REL_CHECK 2
#define GENOME_OUTPUT_FULL 4
#define GENOME_IBD_UNBOUNDED 8
#define GENOME_NUDGE 0x10
// separate flag to ensure behavior is unchanged under --unbounded
#define GENOME_FILTER_PI_HAT 0x20

#define WRITE_COVAR_PHENO 1
#define WRITE_COVAR_NO_PARENTS 2
#define WRITE_COVAR_NO_SEX 4
#define WRITE_COVAR_FEMALE_2 8
#define WRITE_COVAR_DUMMY 0x10
#define WRITE_COVAR_DUMMY_NO_ROUND 0x20

#define MERGE_MODE_MASK 7
#define MERGE_EQUAL_POS 8
#define MERGE_BINARY 16
#define MERGE_LIST 32

#define SAMPLE_SORT_NONE 1
#define SAMPLE_SORT_NATURAL 2
#define SAMPLE_SORT_ASCII 4
#define SAMPLE_SORT_FILE 8

#define REGRESS_PCS_NORMALIZE_PHENO 1
#define REGRESS_PCS_SEX_SPECIFIC 2
#define REGRESS_PCS_CLIP 4

#define HWE_MIDP 1
#define HWE_THRESH_MIDP 2
#define HWE_THRESH_ALL 4
#define HWE_GZ 8

#define MENDEL_FILTER 1
#define MENDEL_FILTER_VAR_FIRST 2
#define MENDEL_DUOS 4
#define MENDEL_MULTIGEN 8
#define MENDEL_SUMMARIES_ONLY 0x10

#define DUMMY_MISSING_GENO 1
#define DUMMY_MISSING_PHENO 2
#define DUMMY_SCALAR_PHENO 4
#define DUMMY_ACGT 8
#define DUMMY_1234 0x10
#define DUMMY_12 0x20

#define SIMULATE_QT 1
#define SIMULATE_TAGS 2
#define SIMULATE_HAPS 4
#define SIMULATE_ACGT 8
#define SIMULATE_1234 0x10
#define SIMULATE_12 0x20

#define MODEL_ASSOC 1
#define MODEL_FISHER 2
#define MODEL_FISHER_MIDP 4
#define MODEL_PERM 8
#define MODEL_MPERM 0x10
#define MODEL_GENEDROP 0x20
#define MODEL_PERM_COUNT 0x40
#define MODEL_ASSOC_COUNTS 0x80
#define MODEL_ASSOC_FDEPR 0x100
#define MODEL_DMASK 0x1a6
#define MODEL_QT_MEANS 0x200
#define MODEL_PDOM 0x400
#define MODEL_PREC 0x800
#define MODEL_PGEN 0x1000
#define MODEL_PTREND 0x2000
#define MODEL_TRENDONLY 0x4000
#define MODEL_PMASK \
    (MODEL_PDOM | MODEL_PREC | MODEL_PGEN | MODEL_PTREND | MODEL_TRENDONLY)
#define MODEL_LIN 0x8000
#define MODEL_QMASK (MODEL_QT_MEANS | MODEL_LIN)
#define MODEL_SET_TEST 0x10000

#define GLM_LOGISTIC 1
#define GLM_PERM 2
#define GLM_MPERM 4
#define GLM_GENEDROP 8
#define GLM_PERM_COUNT 0x10
#define GLM_GENOTYPIC 0x20
#define GLM_HETHOM 0x40
#define GLM_DOMINANT 0x80
#define GLM_RECESSIVE 0x100
#define GLM_NO_SNP 0x200
#define GLM_HIDE_COVAR 0x400
#define GLM_SEX 0x800
#define GLM_NO_X_SEX 0x1000
#define GLM_INTERACTION 0x2000
#define GLM_STANDARD_BETA 0x4000
#define GLM_BETA 0x8000
#define GLM_TEST_ALL 0x10000
#define GLM_CONDITION_DOMINANT 0x20000
#define GLM_CONDITION_RECESSIVE 0x40000
#define GLM_SET_TEST 0x80000
#define GLM_NO_SNP_EXCL 0x831ea
#define GLM_INTERCEPT 0x100000

#define MPERM_DUMP_BEST 1
#define MPERM_DUMP_ALL 2

// (2^31 - 1000001) / 2
#define APERM_MAX 1073241823

#define ADJUST_GC 2
#define ADJUST_LOG10 4
#define ADJUST_QQ 8
#define ADJUST_LAMBDA 16

#define DUPVAR_REF 1
#define DUPVAR_IDS_ONLY 2
#define DUPVAR_SUPPRESS_FIRST 4

#define CNV_MAKE_MAP 1
#define CNV_MAKE_MAP_LONG 2
#define CNV_CHECK_NO_OVERLAP 4
#define CNV_DEL 8
#define CNV_DUP 0x10
#define CNV_WRITE_FREQ 0x20
#define CNV_UNIQUE 0x40
#define CNV_DROP_NO_SEGMENT 0x80
#define CNV_SAMPLE_PERM 0x100
#define CNV_ENRICHMENT_TEST 0x200
#define CNV_TEST 0x400
#define CNV_TEST_FORCE_1SIDED 0x800
#define CNV_TEST_FORCE_2SIDED 0x1000
#define CNV_TEST_REGION 0x2000
#define CNV_TRACK 0x4000
#define CNV_SEGLIST 0x8000
#define CNV_REPORT_REGIONS 0x10000
#define CNV_VERBOSE_REPORT_REGIONS 0x20000
#define CNV_WRITE 0x40000
#define CNV_EXCLUDE_OFF_BY_1 0x80000

#define CNV_INTERSECT 1
#define CNV_EXCLUDE 2
#define CNV_COUNT 4

#define CNV_OVERLAP 1
#define CNV_OVERLAP_REGION 2
#define CNV_OVERLAP_UNION 3
#define CNV_DISRUPT 4

#define CNV_FREQ_EXCLUDE_ABOVE 1
#define CNV_FREQ_EXCLUDE_BELOW 2
#define CNV_FREQ_EXCLUDE_EXACT 4
#define CNV_FREQ_INCLUDE_EXACT 8
#define CNV_FREQ_FILTER 15
#define CNV_FREQ_OVERLAP 16
#define CNV_FREQ_METHOD2 32

#define SEGMENT_GROUP 1

// default jackknife iterations
#define ITERS_DEFAULT 100000
#define MAX_PCS_DEFAULT 20

#define BIGSTACK_MIN_MB 64
#define BIGSTACK_DEFAULT_MB 2048

#ifdef __LP64__
#define BITCT 64

// unions generally shouldn't be used for reinterpret_cast's job (memcpy is
// the right C-compatible way), but vectors are an exception to this rule.
typedef union {
    VECFTYPE vf;
    VECITYPE vi;
    VECDTYPE vd;
    uintptr_t u8[VEC_BITS / BITCT];
    double d8[VEC_BYTES / sizeof(double)];
    float f4[VEC_BYTES / sizeof(float)];
    uint32_t u4[VEC_BYTES / sizeof(int32_t)];
} __univec;
#else
#define BITCT 32
#endif

#define BITCT2 (BITCT / 2)
#define BYTECT (BITCT / 8)
#define BYTECT4 (BITCT / 32)
#define VEC_WORDS (VEC_BITS / BITCT)
#define VEC_INT32 (VEC_BYTES / 4)

// assumed number of bytes per cache line, for alignment
#define CACHELINE 64

#define CACHELINE_BIT (CACHELINE * 8)
#define CACHELINE_INT32 (CACHELINE / 4)
#define CACHELINE_INT64 (CACHELINE / 8)
#define CACHELINE_WORD (CACHELINE / BYTECT)
#define CACHELINE_DBL (CACHELINE / 8)

// alignment must be a power of 2


#define BITCT_TO_VECCT(val) (((val) + (VEC_BITS - 1)) / VEC_BITS)
#define BITCT_TO_WORDCT(val) (((val) + (BITCT - 1)) / BITCT)
#define BITCT_TO_ALIGNED_WORDCT(val) (VEC_WORDS * BITCT_TO_VECCT(val))

#define QUATERCT_TO_VECCT(val) (((val) + ((VEC_BITS / 2) - 1)) / (VEC_BITS / 2))
#define QUATERCT_TO_WORDCT(val) (((val) + (BITCT2 - 1)) / BITCT2)
#define QUATERCT_TO_ALIGNED_WORDCT(val) (VEC_WORDS * QUATERCT_TO_VECCT(val))

// todo: get rid of (BITCT_TO_WORDCT(x) == QUATERCT_TO_VECCT(x)) and similar
// assumptions, in preparation for AVX2


#define MAXV(aa, bb) (((bb) > (aa)) ? (bb) : (aa))
#define MINV(aa, bb) (((aa) > (bb)) ? (bb) : (aa))

#ifdef _WIN32
// if MAX_THREADS > 65, single WaitForMultipleObjects calls must be
// converted into loops
#define MAX_THREADS 64
#define MAX_THREADS_P1 65
#else
// shouldn't be larger than MODEL_BLOCKSIZE for now
#define MAX_THREADS 512
#define MAX_THREADS_P1 513
#endif

// defined as a macro since type of idx can vary; might want a debug
// compilation mode which performs type-checking, though
#define EXTRACT_2BIT_GENO(ulptr, idx) \
    (((ulptr)[(idx) / BITCT2] >> (2 * ((idx) % BITCT2))) & 3)

// generic maximum line length.  .ped/.vcf/etc. lines can of course be longer
#define MAXLINELEN 131072

// must be at least 2 * MAXLINELEN + 2 to support generic token loader.
#define TEXTBUF_SIZE (2 * MAXLINELEN + 256)

// Maximum length of chromosome, variant, FID, IID, cluster, and set IDs (not
// including terminating null, that's what _P1 is for).  This value supports up
// to 8 IDs per line (maximum so far is 5, for e.g. --hom).
#define MAX_ID_SLEN 16000

#define MAX_ID_BLEN (MAX_ID_SLEN + 1)
#define MAX_ID_SLEN_STR "16000"

// Maximum size of "dynamically" allocated line load buffer.  (This is the
// limit that applies to .vcf and similar files.)  Inconvenient to go higher
// since fgets() takes a int32_t size argument.
#define MAXLINEBUFLEN 0x7fffffc0

// Default --perm-batch-size value in most contexts.  It may actually be better
// to *avoid* a power of two due to the need for transpositions involving this
// stride; see e.g. http://danluu.com/3c-conflict/ ; try 448 instead?  This
// should be tested during PLINK 2.0 development.
#define DEFAULT_PERM_BATCH_SIZE 512

// note that this is NOT foolproof: see e.g.
// http://insanecoding.blogspot.com/2007/11/pathmax-simply-isnt.html .  (This
// is why I haven't bothered with OS-based #ifdefs here.)  But it should be
// good enough in practice.
#define FNAMESIZE 4096

// allow extensions like .model.trend.fisher.set.score.adjusted
#define MAX_POST_EXT 39

// number of types of jackknife values to precompute (x^2, y^2, x, y, xy)
#define JACKKNIFE_VALS_DIST 5
#define JACKKNIFE_VALS_GROUPDIST 3

#ifdef __LP64__
// number of snp-major .bed lines to read at once for distance calc if
// exponent is nonzero.
#define MULTIPLEX_DIST_EXP 64
// number of snp-major .bed lines to read at once for relationship calc
#define MULTIPLEX_REL 60
#else
// N.B. 32-bit version not as carefully tested or optimized, but I'll try to
// make sure it works properly
#define MULTIPLEX_DIST_EXP 28
#define MULTIPLEX_REL 30
#endif

// used to size a few tables
#define EXPECTED_MISSING_FREQ 0.05

// load markers in blocks to enable multithreading and, for quantitative
// phenotypes, PERMORY-style LD exploitation
#define MODEL_BLOCKSIZE 1024
#define MODEL_BLOCKKEEP 64

// string hash table constants, currently only relevant for merge operations
// and annotate()
// (dynamic sizing used for main marker name lookup)

// last prime before 2^19
// size chosen to be likely to fit in L3 cache
#define HASHSIZE 524287
#define HASHSIZE_S 524287

#ifdef __LP64__
#define HASHMEM 4194304
#else
#define HASHMEM 2097152
#endif


// Generic text I/O buffer: any function which reads from/writes to a text file
// or the console may clobber it.  Sized to fit two MAXLINELEN-length lines
// plus a bit extra.
extern char g_textbuf[];

extern const char g_one_char_strs[];
extern const char* g_missing_geno_ptr;
extern const char* g_output_missing_geno_ptr;


extern uintptr_t g_failed_alloc_attempt_size;


// file-scope string constants don't always have the g_ prefix, but multi-file
// constants are always tagged.
extern const char g_errstr_fopen[];
extern const char g_cmdline_format_str[];

extern FILE* g_logfile;

// mostly-safe sprintf buffer.  warning: do NOT put allele codes or
// arbitrary-length lists in here.
extern char g_logbuf[];

extern uint32_t g_debug_on;
extern uint32_t g_log_failed;

// should remove this global: multithreaded functions should use a file-local
// thread_ct which will occasionally be smaller due to job size.
extern uint32_t g_thread_ct;

typedef struct ll_str_struct
{
    struct ll_str_struct* next;
    char ss[];
} Ll_str;


#ifdef STABLE_BUILD
#define UNSTABLE(val)                    \
    sptr = strcpya(&(g_logbuf[9]), val); \
    goto main_unstable_disabled
#else
#define UNSTABLE(val)
#endif


// manually managed, very large double-ended stack
extern unsigned char* g_bigstack_base;
extern unsigned char* g_bigstack_end;


#define END_ALLOC_CHUNK 16
#define END_ALLOC_CHUNK_M1 (END_ALLOC_CHUNK - 1)

// assumes size is divisible by END_ALLOC_CHUNK
// (no value in directly calling this with a constant size parameter: the
// compiler will properly optimize a bigstack_end_alloc() call)

#define bigstack_end_aligned_alloc bigstack_end_alloc


// if we need the digit value, better to use (unsigned char)cc - '0'...


// may as well treat all chars < 32, except tab, as eoln...
// kns = "known non-space" (where tab counts as a space)
/*
HEADER_INLINE int32_t is_eoln_kns(unsigned char ucc) {
  return (ucc < 32);
}
*/


// could assert ucc is not a space/tab
#define is_eoln_kns is_space_or_eoln


// Reads an integer in [1, cap].  Assumes first character is nonspace.  Has the
// overflow detection atoi() lacks.


// Write exactly four digits (padding with zeroes if necessary); useful for
// e.g. floating point encoders.  uii must not be >= 10^4.
char* uitoa_z4(uint32_t uii, char* start);


// These limited-precision converters are usually several times as fast as
// grisu2's descendants; and let's not even speak of sprintf.  (I'm guessing
// that the algorithm cannot be made round-trip-safe without throwing away its
// performance advantage, since we currently multiply by numbers like 1.0e256
// which don't have an exact representation.  But these functions are very,
// very good at what they do.)


// assumes min_width >= 5.


/*
HEADER_INLINE char* ftoa_gx(float dxx, char extra_char, char* start) {
  char* penult = ftoa_g(dxx, start);
  *penult = extra_char;
  return &(penult[1]);
}
*/


void magic_num(uint32_t divisor, uint64_t* multp,
               uint32_t* __restrict pre_shiftp,
               uint32_t* __restrict post_shiftp, uint32_t* __restrict incrp);

// let the compiler worry about the second argument's bit width here
#define SET_BIT(idx, arr) ((arr)[(idx) / BITCT] |= ONELU << ((idx) % BITCT))

#define SET_BIT_DBL(idx, arr) \
    ((arr)[(idx) / BITCT2] |= ONELU << (2 * ((idx) % BITCT2)))

// useful for coercing int32_t loc to unsigned
HEADER_INLINE void set_bit(uint32_t loc, uintptr_t* bitarr)
{
    bitarr[loc / BITCT] |= (ONELU << (loc % BITCT));
}


#define CLEAR_BIT(idx, arr) \
    ((arr)[(idx) / BITCT] &= ~(ONELU << ((idx) % BITCT)))

#define CLEAR_BIT_DBL(idx, arr) \
    ((arr)[(idx) / BITCT2] &= ~(ONELU << (2 * ((idx) % BITCT2))))

HEADER_INLINE void clear_bit(uint32_t loc, uintptr_t* bitarr)
{
    bitarr[loc / BITCT] &= ~(ONELU << (loc % BITCT));
}


#define IS_SET(arr, idx) (((arr)[(idx) / BITCT] >> ((idx) % BITCT)) & 1)

#define IS_SET_DBL(arr, idx) \
    (((arr)[(idx) / BITCT2] >> (2 * ((idx) % BITCT2))) & 1)

// use this instead of IS_SET() for signed 32-bit integers
HEADER_INLINE uint32_t is_set(const uintptr_t* bitarr, uint32_t loc)
{
    return (bitarr[loc / BITCT] >> (loc % BITCT)) & 1;
}


#define IS_NONNULL_AND_SET(arr, idx) ((arr) && IS_SET(arr, idx))


#ifdef __LP64__

// double v indicates that size is a vector count, not a word count.
HEADER_INLINE void fill_vvec_zero(size_t size, VECITYPE* vvec)
{
    size_t ulii;
    for (ulii = 0; ulii < size; ulii++) {
        *vvec++ = _mm_setzero_si128();
    }
}
#else

#endif

HEADER_INLINE void fill_ulong_one(size_t size, uintptr_t* ularr)
{
    size_t ulii;
    for (ulii = 0; ulii < size; ulii++) {
        *ularr++ = ~ZEROLU;
    }
}


HEADER_INLINE void fill_uint_zero(size_t size, uint32_t* uiarr)
{
    size_t ulii;
    for (ulii = 0; ulii < size; ulii++) {
        *uiarr++ = 0;
    }
}


// for hash tables where maximum ID string length is not known in advance.


// okay, time to provide O(c log c)-time instead of O(c^2)-time initialization
// (c = # of chromosomes/contigs).
#define MAX_POSSIBLE_CHROM 65280

// get_id_htable_size(MAX_POSSIBLE_CHROM) (use constexpr once sufficient
// compiler support is available)
#define CHROM_NAME_HTABLE_SIZE 130579

// assumes MAX_POSSIBLE_CHROM is a multiple of 64, otherwise add round-up
#define CHROM_MASK_WORDS (MAX_POSSIBLE_CHROM / BITCT)

// (note that n+1, n+2, n+3, and n+4 are reserved for X/Y/XY/MT)
#define MAX_CHROM_TEXTNUM 95

// get_chrom_code_raw() needs to be modified if this changes
#define MAX_CHROM_TEXTNUM_SLEN 2

#define X_OFFSET 0
#define Y_OFFSET 1
#define XY_OFFSET 2
#define MT_OFFSET 3
#define XYMT_OFFSET_CT 4

#define CHROM_X (MAX_POSSIBLE_CHROM + X_OFFSET)
#define CHROM_Y (MAX_POSSIBLE_CHROM + Y_OFFSET)
#define CHROM_XY (MAX_POSSIBLE_CHROM + XY_OFFSET)
#define CHROM_MT (MAX_POSSIBLE_CHROM + MT_OFFSET)

#ifdef __LP64__
// dog requires 42 bits, and other species require less
#define CHROM_MASK_INITIAL_WORDS 1
#else
#define CHROM_MASK_INITIAL_WORDS 2
#endif


extern const char* g_species_singular;
extern const char* g_species_plural;


// in the usual case where the number of chromosomes/contigs is much less than
// MAX_POSSIBLE_CHROM, this reduces chrom_info's memory consumption and
// improves locality.

#define CHR_OUTPUT_PREFIX 1
#define CHR_OUTPUT_M 2
#define CHR_OUTPUT_MT 4
#define CHR_OUTPUT_0M 8


// does not require null-termination
// only handles 1-99, X, Y, XY, MT, and "chr" prefix
int32_t get_chrom_code_raw(const char* sptr);


// when it's okay to just replace the terminating space/tab with a \0


#ifndef __cplusplus
int32_t llcmp(const void* aa, const void* bb);
#endif


int32_t bsearch_str(const char* id_buf, uintptr_t cur_id_len, const char* lptr,
                    uintptr_t max_id_len, uintptr_t end_idx);

// These ensure the trailing bits are zeroed out.

HEADER_INLINE uint32_t popcount2_long(uintptr_t val)
{
#ifdef __LP64__
    val = (val & 0x3333333333333333LLU) + ((val >> 2) & 0x3333333333333333LLU);
    return (((val + (val >> 4)) & 0x0f0f0f0f0f0f0f0fLLU)
            * 0x0101010101010101LLU)
           >> 56;
#else
    val = (val & 0x33333333) + ((val >> 2) & 0x33333333);
    return (((val + (val >> 4)) & 0x0f0f0f0f) * 0x01010101) >> 24;
#endif
}

HEADER_INLINE uint32_t popcount_long(uintptr_t val)
{
    // the simple version, good enough for all non-time-critical stuff
    return popcount2_long(val - ((val >> 1) & FIVEMASK));
}

// same as is_monomorphic, except it also flags the all-heterozygote case

// uint32_t has_three_genotypes(uintptr_t* lptr, uint32_t sample_ct);

uintptr_t popcount_longs(const uintptr_t* lptr, uintptr_t word_ct);

#define popcount01_longs popcount2_longs


uintptr_t popcount_longs_intersect(const uintptr_t* __restrict lptr1,
                                   const uintptr_t* __restrict lptr2,
                                   uintptr_t word_ct);

#ifdef __LP64__
void count_2freq_dbl_960b(
    const VECITYPE* geno_vvec, const VECITYPE* geno_vvec_end,
    const VECITYPE* __restrict mask1vp, const VECITYPE* __restrict mask2vp,
    uint32_t* __restrict ct1abp, uint32_t* __restrict ct1cp,
    uint32_t* __restrict ct2abp, uint32_t* __restrict ct2cp);

void count_3freq_1920b(const VECITYPE* geno_vvec, const VECITYPE* geno_vvec_end,
                       const VECITYPE* __restrict maskvp,
                       uint32_t* __restrict ctap, uint32_t* __restrict ctbp,
                       uint32_t* __restrict ctcp);
#else
#endif


void fill_all_bits(uintptr_t ct, uintptr_t* bitarr);


void reverse_loadbuf(uintptr_t unfiltered_sample_ct, unsigned char* loadbuf);


HEADER_INLINE uint32_t load_raw(uintptr_t unfiltered_sample_ct4, FILE* bedfile,
                                uintptr_t* rawbuf)
{
    // only use this if all accesses to the data involve
    // 1. some sort of mask, or
    // 2. explicit iteration from 0..(unfiltered_sample_ct-1).
    // otherwise improper trailing bits might cause a segfault, when we should
    // be ignoring them or just issuing a warning.
    // result will be stored within rawbuf
    return (fread(rawbuf, 1, unfiltered_sample_ct4, bedfile)
            < unfiltered_sample_ct4);
}

HEADER_INLINE uintptr_t get_final_mask(uint32_t sample_ct)
{
    uint32_t uii = sample_ct % BITCT2;
    if (uii) {
        return (ONELU << (2 * uii)) - ONELU;
    }
    else
    {
        return ~ZEROLU;
    }
}


// was "collapse_copy_quaterarr_incl", but this should be better way to think
// about it
void copy_quaterarr_nonempty_subset(const uintptr_t* __restrict raw_quaterarr,
                                    const uintptr_t* __restrict subset_mask,
                                    uint32_t raw_quaterarr_size,
                                    uint32_t subset_size,
                                    uintptr_t* __restrict output_quaterarr);

/*
// in-place version of copy_quaterarr_subset (usually destroying original
// data).
// this doesn't seem to provide a meaningful advantage over
// copy_quaterarr_subset in practice, and the latter is more versatile without
// requiring much more memory.
void inplace_quaterarr_proper_subset(const uintptr_t* __restrict subset_mask,
uint32_t orig_quaterarr_size, uint32_t subset_size, uintptr_t* __restrict
main_quaterarr);

HEADER_INLINE void inplace_quaterarr_subset(const uintptr_t* __restrict
subset_mask, uint32_t orig_quaterarr_size, uint32_t subset_size, uintptr_t*
__restrict main_quaterarr) { if (orig_quaterarr_size == subset_size) { return;
  }
  inplace_quaterarr_proper_subset(subset_mask, orig_quaterarr_size, subset_size,
main_quaterarr);
}
*/

uint32_t load_and_collapse_incl(uint32_t unfiltered_sample_ct,
                                uint32_t sample_ct,
                                const uintptr_t* __restrict sample_include,
                                uintptr_t final_mask, uint32_t do_reverse,
                                FILE* bedfile, uintptr_t* __restrict rawbuf,
                                uintptr_t* __restrict mainbuf);

// target_vec := source_vec ANDNOT exclude_vec
// may write an extra word


// excludes (excl_bitarr_1 & excl_bitarr_2).  (union can be excluded by calling
// apply_excl_to_quaterarr_01() twice.)

// initializes output_quatervec bits to 01 iff input_quatervec bits are 01,
// everything else zeroed out


void hh_reset(unsigned char* loadbuf, uintptr_t* sample_include_quaterarr,
              uintptr_t unfiltered_sample_ct);

void hh_reset_y(unsigned char* loadbuf, uintptr_t* sample_include_quaterarr,
                uintptr_t* sample_male_include_quaterarr,
                uintptr_t unfiltered_sample_ct);


uint32_t cubic_real_roots(double coef_a, double coef_b, double coef_c,
                          double* solutions);

void genovec_3freq(const uintptr_t* __restrict geno_vec,
                   const uintptr_t* __restrict include_quatervec,
                   uintptr_t sample_ctl2, uint32_t* __restrict missing_ctp,
                   uint32_t* __restrict het_ctp,
                   uint32_t* __restrict homset_ctp);

#ifdef __LP64__
void count_2freq_dbl_960b(
    const VECITYPE* geno_vvec, const VECITYPE* geno_vvec_end,
    const VECITYPE* __restrict mask1vp, const VECITYPE* __restrict mask2vp,
    uint32_t* __restrict ct1abp, uint32_t* __restrict ct1cp,
    uint32_t* __restrict ct2abp, uint32_t* __restrict ct2cp);

void count_3freq_1920b(const VECITYPE* geno_vvec, const VECITYPE* geno_vvec_end,
                       const VECITYPE* __restrict maskvp,
                       uint32_t* __restrict ctap, uint32_t* __restrict ctbp,
                       uint32_t* __restrict ctcp);
#else
void count_2freq_dbl_24b(const uintptr_t* __restrict geno_vec,
                         const uintptr_t* __restrict mask1p,
                         const uintptr_t* __restrict mask2p,
                         uint32_t* __restrict ct1abp,
                         uint32_t* __restrict ct1cp,
                         uint32_t* __restrict ct2abp,
                         uint32_t* __restrict ct2cp);

void count_3freq_48b(const uintptr_t* __restrict geno_vec,
                     const uintptr_t* __restrict maskp,
                     uint32_t* __restrict ctap, uint32_t* __restrict ctbp,
                     uint32_t* __restrict ctcp);
#endif

void fill_quatervec_55(uint32_t ct, uintptr_t* quatervec);

void vec_datamask(uintptr_t unfiltered_sample_ct, uint32_t matchval,
                  uintptr_t* data_ptr, uintptr_t* mask_ptr,
                  uintptr_t* result_ptr);
uintptr_t popcount2_longs(const uintptr_t* lptr, uintptr_t word_ct);

#define popcount01_longs popcount2_longs

#endif // __PLINK_COMMON_H__

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


#include "plink_common.hpp"


// no leading \n since this is used in LOGPRINTFWW expressions
const char g_errstr_fopen[] = "Error: Failed to open %s.\n";

const char g_cmdline_format_str[] =
    "\n  " PROG_NAME_STR " [input flag(s)...] {command flag(s)...} {other "
    "flag(s)...}\n  " PROG_NAME_STR " --help {flag name(s)...}\n\n";

char g_textbuf[TEXTBUF_SIZE];

// note that \xxx character constants are interpreted in octal.
// technically no need to represent 0-31, but 64 extra bytes of data is
// probably cheaper than the code to subtract 32 everywhere.
const char g_one_char_strs[] =
    "\0\0\1\0\2\0\3\0\4\0\5\0\6\0\7\0\10\0\11\0\12\0\13\0\14\0\15\0\16\0\17\0"
    "\20\0\21\0\22\0\23\0\24\0\25\0\26\0\27\0\30\0\31\0\32\0\33\0\34\0\35\0\36"
    "\0\37\0\40\0\41\0\42\0\43\0\44\0\45\0\46\0\47\0\50\0\51\0\52\0\53\0\54\0"
    "\55\0\56\0\57\0\60\0\61\0\62\0\63\0\64\0\65\0\66\0\67\0\70\0\71\0\72\0\73"
    "\0\74\0\75\0\76\0\77\0\100\0\101\0\102\0\103\0\104\0\105\0\106\0\107\0\110"
    "\0\111\0\112\0\113\0\114\0\115\0\116\0\117\0\120\0\121\0\122\0\123\0\124\0"
    "\125\0\126\0\127\0\130\0\131\0\132\0\133\0\134\0\135\0\136\0\137\0\140\0"
    "\141\0\142\0\143\0\144\0\145\0\146\0\147\0\150\0\151\0\152\0\153\0\154\0"
    "\155\0\156\0\157\0\160\0\161\0\162\0\163\0\164\0\165\0\166\0\167\0\170\0"
    "\171\0\172\0\173\0\174\0\175\0\176\0\177\0\200\0\201\0\202\0\203\0\204\0"
    "\205\0\206\0\207\0\210\0\211\0\212\0\213\0\214\0\215\0\216\0\217\0\220\0"
    "\221\0\222\0\223\0\224\0\225\0\226\0\227\0\230\0\231\0\232\0\233\0\234\0"
    "\235\0\236\0\237\0\240\0\241\0\242\0\243\0\244\0\245\0\246\0\247\0\250\0"
    "\251\0\252\0\253\0\254\0\255\0\256\0\257\0\260\0\261\0\262\0\263\0\264\0"
    "\265\0\266\0\267\0\270\0\271\0\272\0\273\0\274\0\275\0\276\0\277\0\300\0"
    "\301\0\302\0\303\0\304\0\305\0\306\0\307\0\310\0\311\0\312\0\313\0\314\0"
    "\315\0\316\0\317\0\320\0\321\0\322\0\323\0\324\0\325\0\326\0\327\0\330\0"
    "\331\0\332\0\333\0\334\0\335\0\336\0\337\0\340\0\341\0\342\0\343\0\344\0"
    "\345\0\346\0\347\0\350\0\351\0\352\0\353\0\354\0\355\0\356\0\357\0\360\0"
    "\361\0\362\0\363\0\364\0\365\0\366\0\367\0\370\0\371\0\372\0\373\0\374\0"
    "\375\0\376\0\377";
const char* g_missing_geno_ptr = &(g_one_char_strs[96]);
const char* g_output_missing_geno_ptr = &(g_one_char_strs[96]);

uintptr_t g_failed_alloc_attempt_size = 0;


FILE* g_logfile = nullptr;

char g_logbuf[MAXLINELEN * 2];

uint32_t g_debug_on = 0;
uint32_t g_log_failed = 0;
uint32_t g_thread_ct;


// manually managed, very large stack
unsigned char* g_bigstack_base;
unsigned char* g_bigstack_end;


// Okay, time to do banker's rounding when printing doubles.  14 digits of
// precision are used in judging equality to 0.5 (actual precision of doubles
// is 15-17 digits); the intention is to capture all directly loaded or exactly
// computed edge cases (so enough tolerance is needed to survive the internal
// multiplications by powers of 10, etc.), while rounding a negligible number
// of honest-to-god 0.4999999s up and 0.5000001s down.
// To avoid inadvertent printing of an extra digit, there's a deliberate gap
// between the 99.9994999...-type bounds and the largest numbers that would
// actually round down.


void magic_num(uint32_t divisor, uint64_t* multp,
               uint32_t* __restrict pre_shiftp,
               uint32_t* __restrict post_shiftp, uint32_t* __restrict incrp)
{
    // Enables fast integer division by a constant not known until runtime.  See
    // http://ridiculousfish.com/blog/posts/labor-of-division-episode-iii.html .
    // Assumes divisor is not zero, of course.
    // May want to populate a struct instead.
    uint32_t down_multiplier = 0;
    uint32_t down_exponent = 0;
    uint32_t has_magic_down = 0;
    uint32_t quotient;
    uint32_t remainder;
    uint32_t ceil_log_2_d;
    uint32_t exponent;
    uint32_t uii;
    if (divisor & (divisor - 1)) {
        quotient = 0x80000000U / divisor;
        remainder = 0x80000000U - (quotient * divisor);
        ceil_log_2_d = 32 - __builtin_clz(divisor);
        for (exponent = 0;; exponent++) {
            if (remainder >= divisor - remainder) {
                quotient = quotient * 2 + 1;
                remainder = remainder * 2 - divisor;
            }
            else
            {
                quotient = quotient * 2;
                remainder = remainder * 2;
            }
            if ((exponent >= ceil_log_2_d)
                || (divisor - remainder) <= (1U << exponent))
            {
                break;
            }
            if ((!has_magic_down) && (remainder <= (1U << exponent))) {
                has_magic_down = 1;
                down_multiplier = quotient;
                down_exponent = exponent;
            }
        }
        if (exponent < ceil_log_2_d) {
            *multp = quotient + 1;
            *pre_shiftp = 0;
            *post_shiftp = 32 + exponent;
            *incrp = 0;
            return;
        }
        else if (divisor & 1)
        {
            *multp = down_multiplier;
            *pre_shiftp = 0;
            *post_shiftp = 32 + down_exponent;
            *incrp = 1;
            return;
        }
        else
        {
            *pre_shiftp = __builtin_ctz(divisor);
            magic_num(divisor >> (*pre_shiftp), multp, &uii, post_shiftp,
                      incrp);
            return;
        }
    }
    else
    {
        // power of 2
        *multp = 1;
        *pre_shiftp = 0;
        *post_shiftp = __builtin_ctz(divisor);
        *incrp = 0;
    }
}


// hashval computation left to caller since this is frequently used with
// chromosome IDs, where the compiler can optimize the integer modulus
// operation since the hash table size is preset


// Global since species_str() may be called by functions which don't actually
// care about chrom_info.  (chrom_info is really a global variable too, but I
// find it easier to maintain this code when chrom_info dependencies are made
// explicit in the function signatures; in contrast, g_species_singular and
// g_species_plural are just for pretty printing and lend no insight into what
// the functions which reference them are doing.)
const char* g_species_singular = nullptr;
const char* g_species_plural = nullptr;


static inline int32_t single_letter_chrom(uint32_t letter)
{
    letter &= 0xdf;
    if (letter == 'X') {
        return CHROM_X;
    }
    else if (letter == 'Y')
    {
        return CHROM_Y;
    }
    else if (letter == 'M')
    {
        return CHROM_MT;
    }
    else
    {
        return -1;
    }
}

int32_t get_chrom_code_raw(const char* sptr)
{
    // any character <= ' ' is considered a terminator
    // note that char arithmetic tends to be compiled to int32 operations, so we
    // mostly work with ints here
    // assumes MAX_CHROM_TEXTNUM_SLEN == 2
    uint32_t first_char_code = (unsigned char) sptr[0];
    uint32_t second_char_code = (unsigned char) sptr[1];
    if ((first_char_code & 0xdf) == 'C') {
        if (((second_char_code & 0xdf) == 'H')
            && ((((unsigned char) sptr[2]) & 0xdf) == 'R'))
        {
            sptr = &(sptr[3]);
            first_char_code = (unsigned char) sptr[0];
            second_char_code = (unsigned char) sptr[1];
        }
        else
        {
            return -1;
        }
    }
    if (second_char_code > ' ') {
        if (sptr[2] > ' ') {
            return -1;
        }
        const uint32_t first_char_toi = first_char_code - '0';
        if (first_char_toi < 10) {
            const uint32_t second_char_toi = second_char_code - '0';
            if (second_char_toi < 10) {
                return first_char_toi * 10 + second_char_toi;
            }
            else if (!first_char_toi)
            {
                // accept '0X', '0Y', '0M' emitted by Oxford software
                return single_letter_chrom(second_char_code);
            }
        }
        else
        {
            first_char_code &= 0xdf;
            if (first_char_code == 'X') {
                if ((second_char_code == 'Y') || (second_char_code == 'y')) {
                    return CHROM_XY;
                }
            }
            else if (first_char_code == 'M')
            {
                if ((second_char_code == 'T') || (second_char_code == 't')) {
                    return CHROM_MT;
                }
            }
        }
    }
    else
    {
        const uint32_t first_char_toi = first_char_code - '0';
        if (first_char_toi < 10) {
            return first_char_toi;
        }
        else
        {
            return single_letter_chrom(first_char_code);
        }
    }
    return -1;
}


#ifndef __cplusplus
int32_t llcmp(const void* aa, const void* bb)
{
    int64_t diff = *((const int64_t*) aa) - *((const int64_t*) bb);
    if (diff > 0) {
        return 1;
    }
    else if (diff < 0)
    {
        return -1;
    }
    else
    {
        return 0;
    }
}
#endif


#ifdef __LP64__
// Basic SSE2 implementation of Lauradoux/Walisch popcount.
static inline uintptr_t popcount_vecs(const __m128i* vptr, uintptr_t ct)
{
    // popcounts vptr[0..(ct-1)].  Assumes ct is a multiple of 3 (0 ok).
    const __m128i m1 = {FIVEMASK, FIVEMASK};
    const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
    const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
    const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
    uintptr_t tot = 0;
    const __m128i* vend;
    __m128i count1;
    __m128i count2;
    __m128i half1;
    __m128i half2;
    __univec acc;

    while (ct >= 30) {
        ct -= 30;
        vend = &(vptr[30]);
    popcount_vecs_main_loop:
        acc.vi = _mm_setzero_si128();
        do
        {
            count1 = *vptr++;
            count2 = *vptr++;
            half1 = *vptr++;
            half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
            half1 = _mm_and_si128(half1, m1);
            // Two bits can represent values from 0-3, so make each pair in
            // count1 count2 store a partial bitcount covering themselves AND
            // another bit from elsewhere.
            count1 = _mm_sub_epi64(
                count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
            count2 = _mm_sub_epi64(
                count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
            count1 = _mm_add_epi64(count1, half1);
            count2 = _mm_add_epi64(count2, half2);
            // Four bits represent 0-15, so we can safely add four 0-3 partial
            // bitcounts together.
            count1 =
                _mm_add_epi64(_mm_and_si128(count1, m2),
                              _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
            count1 = _mm_add_epi64(
                count1,
                _mm_add_epi64(_mm_and_si128(count2, m2),
                              _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
            // Accumulator stores sixteen 0-255 counts in parallel.
            acc.vi = _mm_add_epi64(
                acc.vi,
                _mm_add_epi64(_mm_and_si128(count1, m4),
                              _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
        } while (vptr < vend);
        acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8),
                               _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
        tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
    }
    if (ct) {
        vend = &(vptr[ct]);
        ct = 0;
        goto popcount_vecs_main_loop;
    }
    return tot;
}


static inline uintptr_t popcount_vecs_intersect(const __m128i* __restrict vptr1,
                                                const __m128i* __restrict vptr2,
                                                uintptr_t ct)
{
    // popcounts vptr1 AND vptr2[0..(ct-1)].  ct is a multiple of 3.
    const __m128i m1 = {FIVEMASK, FIVEMASK};
    const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
    const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
    const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
    uintptr_t tot = 0;
    const __m128i* vend1;
    __m128i count1, count2, half1, half2;
    __univec acc;

    while (ct >= 30) {
        ct -= 30;
        vend1 = &(vptr1[30]);
    popcount_vecs_intersect_main_loop:
        acc.vi = _mm_setzero_si128();
        do
        {
            count1 = _mm_and_si128(*vptr2++, *vptr1++);
            count2 = _mm_and_si128(*vptr2++, *vptr1++);
            half1 = _mm_and_si128(*vptr2++, *vptr1++);
            half2 = _mm_and_si128(_mm_srli_epi64(half1, 1), m1);
            half1 = _mm_and_si128(half1, m1);
            count1 = _mm_sub_epi64(
                count1, _mm_and_si128(_mm_srli_epi64(count1, 1), m1));
            count2 = _mm_sub_epi64(
                count2, _mm_and_si128(_mm_srli_epi64(count2, 1), m1));
            count1 = _mm_add_epi64(count1, half1);
            count2 = _mm_add_epi64(count2, half2);
            count1 =
                _mm_add_epi64(_mm_and_si128(count1, m2),
                              _mm_and_si128(_mm_srli_epi64(count1, 2), m2));
            count1 = _mm_add_epi64(
                count1,
                _mm_add_epi64(_mm_and_si128(count2, m2),
                              _mm_and_si128(_mm_srli_epi64(count2, 2), m2)));
            acc.vi = _mm_add_epi64(
                acc.vi,
                _mm_add_epi64(_mm_and_si128(count1, m4),
                              _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
        } while (vptr1 < vend1);
        acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8),
                               _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
        tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
    }
    if (ct) {
        vend1 = &(vptr1[ct]);
        ct = 0;
        goto popcount_vecs_intersect_main_loop;
    }
    return tot;
}
#endif

uintptr_t popcount_longs(const uintptr_t* lptr, uintptr_t word_ct)
{
    // Efficiently popcounts lptr[0..(word_ct - 1)].  In the 64-bit case, lptr[]
    // must be 16-byte aligned.
    // The popcount_longs_nzbase() wrapper takes care of starting from a later
    // index.
    uintptr_t tot = 0;
    const uintptr_t* lptr_end = &(lptr[word_ct]);
#ifdef __LP64__
    uintptr_t six_ct;
    const __m128i* vptr;
    vptr = (const __m128i*) lptr;
    six_ct = word_ct / 6;
    tot += popcount_vecs(vptr, six_ct * 3);
    lptr = &(lptr[six_ct * 6]);
#else
    // The humble 16-bit lookup table actually beats
    // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
    // on my development machine by a hair.
    // However, if we take the hint from Lauradoux/Walisch and postpone the
    // multiply and right shift, this is no longer true.  Ah well.
    const uintptr_t* lptr_six_end;
    uintptr_t tmp_stor;
    uintptr_t loader;
    uintptr_t ulii;
    uintptr_t uljj;
    lptr_six_end = &(lptr[word_ct - (word_ct % 6)]);
    while (lptr < lptr_six_end) {
        loader = *lptr++;
        ulii = loader - ((loader >> 1) & FIVEMASK);
        loader = *lptr++;
        uljj = loader - ((loader >> 1) & FIVEMASK);
        loader = *lptr++;
        ulii += (loader >> 1) & FIVEMASK;
        uljj += loader & FIVEMASK;
        ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
        ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
        tmp_stor = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

        loader = *lptr++;
        ulii = loader - ((loader >> 1) & FIVEMASK);
        loader = *lptr++;
        uljj = loader - ((loader >> 1) & FIVEMASK);
        loader = *lptr++;
        ulii += (loader >> 1) & FIVEMASK;
        uljj += loader & FIVEMASK;
        ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
        ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
        tmp_stor += (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

        // Each 8-bit slot stores a number in 0..48.  Multiplying by 0x01010101
        // is equivalent to the left-shifts and adds we need to sum those four
        // 8-bit numbers in the high-order slot.
        tot += (tmp_stor * 0x01010101) >> 24;
    }
#endif
    while (lptr < lptr_end) {
        tot += popcount_long(*lptr++);
    }
    return tot;
}


uintptr_t popcount_longs_intersect(const uintptr_t* __restrict lptr1,
                                   const uintptr_t* __restrict lptr2,
                                   uintptr_t word_ct)
{
    uintptr_t tot = 0;
    const uintptr_t* lptr1_end = &(lptr1[word_ct]);
#ifdef __LP64__
    uintptr_t six_ct = word_ct / 6;
    tot += popcount_vecs_intersect((const __m128i*) lptr1,
                                   (const __m128i*) lptr2, six_ct * 3);
    lptr1 = &(lptr1[six_ct * 6]);
    lptr2 = &(lptr2[six_ct * 6]);
#else
    const uintptr_t* lptr1_six_end;
    uintptr_t tmp_stor;
    uintptr_t loader;
    uintptr_t ulii;
    uintptr_t uljj;
    lptr1_six_end = &(lptr1[word_ct - (word_ct % 6)]);
    while (lptr1 < lptr1_six_end) {
        loader = (*lptr1++) & (*lptr2++);
        ulii = loader - ((loader >> 1) & FIVEMASK);
        loader = (*lptr1++) & (*lptr2++);
        uljj = loader - ((loader >> 1) & FIVEMASK);
        loader = (*lptr1++) & (*lptr2++);
        ulii += (loader >> 1) & FIVEMASK;
        uljj += loader & FIVEMASK;
        ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
        ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
        tmp_stor = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

        loader = (*lptr1++) & (*lptr2++);
        ulii = loader - ((loader >> 1) & FIVEMASK);
        loader = (*lptr1++) & (*lptr2++);
        uljj = loader - ((loader >> 1) & FIVEMASK);
        loader = (*lptr1++) & (*lptr2++);
        ulii += (loader >> 1) & FIVEMASK;
        uljj += loader & FIVEMASK;
        ulii = (ulii & 0x33333333) + ((ulii >> 2) & 0x33333333);
        ulii += (uljj & 0x33333333) + ((uljj >> 2) & 0x33333333);
        tmp_stor += (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);

        // Each 8-bit slot stores a number in 0..48.  Multiplying by 0x01010101
        // is equivalent to the left-shifts and adds we need to sum those four
        // 8-bit numbers in the high-order slot.
        tot += (tmp_stor * 0x01010101) >> 24;
    }
#endif
    while (lptr1 < lptr1_end) {
        tot += popcount_long((*lptr1++) & (*lptr2++));
    }
    return tot;
}


void fill_all_bits(uintptr_t ct, uintptr_t* bitarr)
{
    // leaves bits beyond the end unset
    // ok for ct == 0
    uintptr_t quotient = ct / BITCT;
    uintptr_t remainder = ct % BITCT;
    fill_ulong_one(quotient, bitarr);
    if (remainder) {
        bitarr[quotient] = (ONELU << remainder) - ONELU;
    }
}

void reverse_loadbuf(uintptr_t unfiltered_sample_ct, unsigned char* loadbuf)
{
    // unfiltered_sample_ct can be zero
    uintptr_t sample_bidx = 0;
    unsigned char* loadbuf_end = &(loadbuf[(unfiltered_sample_ct + 3) / 4]);
    unsigned char ucc;
    unsigned char ucc2;
    uintptr_t unfiltered_sample_ctd;
    uint32_t* loadbuf_alias32;
    uint32_t uii;
    uint32_t ujj;
#ifdef __LP64__
    const __m128i m1 = {FIVEMASK, FIVEMASK};
    __m128i* loadbuf_alias;
    __m128i vii;
    __m128i vjj;
    // todo: use this vector loop even when loadbuf is unaligned, so stuff like
    // recode_load_to() is faster
    if (!(((uintptr_t) loadbuf) & 15)) {
        loadbuf_alias = (__m128i*) loadbuf;
        unfiltered_sample_ctd = unfiltered_sample_ct / 64;
        for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
            vii = *loadbuf_alias;
            // we want to exchange 00 and 11, and leave 01/10 untouched.  So
            // make vjj := 11 iff vii is 00/11, and vjj := 00 otherwise; then
            // xor.
            vjj = _mm_andnot_si128(_mm_xor_si128(vii, _mm_srli_epi64(vii, 1)),
                                   m1);
            vjj = _mm_or_si128(vjj, _mm_slli_epi64(vjj, 1));
            *loadbuf_alias++ = _mm_xor_si128(vii, vjj);
        }
        loadbuf = (unsigned char*) loadbuf_alias;
    }
    else if (!(((uintptr_t) loadbuf) & 3))
    {
        loadbuf_alias32 = (uint32_t*) loadbuf;
        unfiltered_sample_ctd = unfiltered_sample_ct / BITCT2;
        for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
            uii = *loadbuf_alias32;
            ujj = 0x55555555 & (~(uii ^ (uii >> 1)));
            ujj *= 3;
            *loadbuf_alias32++ = uii ^ ujj;
        }
        loadbuf = (unsigned char*) loadbuf_alias32;
    }
#else
    if (!(((uintptr_t) loadbuf) & 3)) {
        loadbuf_alias32 = (uint32_t*) loadbuf;
        unfiltered_sample_ctd = unfiltered_sample_ct / BITCT2;
        for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
            uii = *loadbuf_alias32;
            ujj = 0x55555555 & (~(uii ^ (uii >> 1)));
            ujj *= 3;
            *loadbuf_alias32++ = uii ^ ujj;
        }
        loadbuf = (unsigned char*) loadbuf_alias32;
    }
#endif
    for (; loadbuf < loadbuf_end;) {
        ucc = *loadbuf;
        ucc2 = 0x55 & (~(ucc ^ (ucc >> 1)));
        ucc2 *= 3;
        *loadbuf++ = ucc ^ ucc2;
    }
    uii = unfiltered_sample_ct & 3;
    if (uii) {
        loadbuf[-1] &= (0xff >> (8 - 2 * uii));
    }
}

void copy_quaterarr_nonempty_subset(const uintptr_t* __restrict raw_quaterarr,
                                    const uintptr_t* __restrict subset_mask,
                                    uint32_t raw_quaterarr_size,
                                    uint32_t subset_size,
                                    uintptr_t* __restrict output_quaterarr)
{
    // in plink 2.0, we probably want (0-based) bit raw_quaterarr_size of
    // subset_mask to be always allocated and unset.  This removes a few special
    // cases re: iterating past the end of arrays.
    assert(subset_size);
    assert(raw_quaterarr_size >= subset_size);
    uintptr_t cur_output_word = 0;
    uintptr_t* output_quaterarr_last =
        &(output_quaterarr[subset_size / BITCT2]);
    const uint32_t word_write_halfshift_end = subset_size % BITCT2;
    uint32_t word_write_halfshift = 0;
    // if < 2/3-filled, use sparse copy algorithm
    if (subset_size * (3 * ONELU) < raw_quaterarr_size * (2 * ONELU)) {
        uint32_t subset_mask_widx = 0;
        while (1) {
            const uintptr_t cur_include_word = subset_mask[subset_mask_widx];
            if (cur_include_word) {
                uint32_t wordhalf_idx = 0;
#ifdef __LP64__
                uint32_t cur_include_halfword = (uint32_t) cur_include_word;
#else
                uint32_t cur_include_halfword = (uint16_t) cur_include_word;
#endif
                while (1) {
                    if (cur_include_halfword) {
                        uintptr_t raw_quaterarr_word =
                            raw_quaterarr[subset_mask_widx * 2 + wordhalf_idx];
                        do
                        {
                            uint32_t rqa_idx_lowbits =
                                __builtin_ctz(cur_include_halfword);
                            cur_output_word |=
                                ((raw_quaterarr_word >> (rqa_idx_lowbits * 2))
                                 & 3)
                                << (word_write_halfshift * 2);
                            if (++word_write_halfshift == BITCT2) {
                                *output_quaterarr++ = cur_output_word;
                                word_write_halfshift = 0;
                                cur_output_word = 0;
                            }
                            cur_include_halfword &= cur_include_halfword - 1;
                        } while (cur_include_halfword);
                    }
                    if (wordhalf_idx) {
                        break;
                    }
                    wordhalf_idx++;
#ifdef __LP64__
                    cur_include_halfword = cur_include_word >> 32;
#else
                    cur_include_halfword = cur_include_word >> 16;
#endif
                }
                if (output_quaterarr == output_quaterarr_last) {
                    if (word_write_halfshift == word_write_halfshift_end) {
                        if (word_write_halfshift_end) {
                            *output_quaterarr_last = cur_output_word;
                        }
                        return;
                    }
                }
            }
            subset_mask_widx++;
        }
    }
    // blocked copy
    while (1) {
        const uintptr_t cur_include_word = *subset_mask++;
        uint32_t wordhalf_idx = 0;
#ifdef __LP64__
        uintptr_t cur_include_halfword = (uint32_t) cur_include_word;
#else
        uint32_t cur_include_halfword = (uint16_t) cur_include_word;
#endif
        while (1) {
            uintptr_t raw_quaterarr_word = *raw_quaterarr++;
            while (cur_include_halfword) {
                uint32_t rqa_idx_lowbits = CTZLU(cur_include_halfword);
                uintptr_t halfword_invshifted =
                    (~cur_include_halfword) >> rqa_idx_lowbits;
                uintptr_t raw_quaterarr_curblock_unmasked =
                    raw_quaterarr_word >> (rqa_idx_lowbits * 2);
                uint32_t rqa_block_len = CTZLU(halfword_invshifted);
                uint32_t block_len_limit = BITCT2 - word_write_halfshift;
                cur_output_word |= raw_quaterarr_curblock_unmasked
                                   << (2 * word_write_halfshift);
                if (rqa_block_len < block_len_limit) {
                    word_write_halfshift += rqa_block_len;
                    cur_output_word &=
                        (ONELU << (word_write_halfshift * 2)) - ONELU;
                }
                else
                {
                    // no need to mask, extra bits vanish off the high end
                    *output_quaterarr++ = cur_output_word;
                    word_write_halfshift = rqa_block_len - block_len_limit;
                    if (word_write_halfshift) {
                        cur_output_word =
                            (raw_quaterarr_curblock_unmasked
                             >> (2 * block_len_limit))
                            & ((ONELU << (2 * word_write_halfshift)) - ONELU);
                    }
                    else
                    {
                        // avoid potential right-shift-64
                        cur_output_word = 0;
                    }
                }
                cur_include_halfword &=
                    (~(ONELU << (rqa_block_len + rqa_idx_lowbits))) + ONELU;
            }
            if (wordhalf_idx) {
                break;
            }
            wordhalf_idx++;
#ifdef __LP64__
            cur_include_halfword = cur_include_word >> 32;
#else
            cur_include_halfword = cur_include_word >> 16;
#endif
        }
        if (output_quaterarr == output_quaterarr_last) {
            if (word_write_halfshift == word_write_halfshift_end) {
                if (word_write_halfshift_end) {
                    *output_quaterarr_last = cur_output_word;
                }
                return;
            }
        }
    }
}

/*
void inplace_quaterarr_proper_subset(const uintptr_t* __restrict subset_mask,
uint32_t orig_quaterarr_size, uint32_t subset_size, uintptr_t* __restrict
main_quaterarr) { assert(orig_quaterarr_size > subset_size);
  // worthwhile to special-case this since we get to entirely skip
  // reading/writing these words
  if (!(~subset_mask[0])) {
    const uintptr_t* subset_mask_initial = subset_mask;
    // guaranteed to terminate since orig_quaterarr_size > subset_size.
    do {
      subset_mask++;
    } while (!(~subset_mask[0]));
    const uint32_t quaterarr_word_skip_ct = 2 * ((uintptr_t)(subset_mask -
subset_mask_initial)); main_quaterarr =
&(main_quaterarr[quaterarr_word_skip_ct]); const uint32_t item_skip_ct =
quaterarr_word_skip_ct * BITCT2; orig_quaterarr_size -= item_skip_ct;
    subset_size -= item_skip_ct;
  }
  uintptr_t cur_output_word = 0;
  uintptr_t* main_quaterarr_writer = main_quaterarr;
  uintptr_t* main_quaterarr_write_last = &(main_quaterarr[subset_size /
BITCT2]); const uint32_t word_write_halfshift_end = subset_size % BITCT2;
  uint32_t word_write_halfshift = 0;
  // if <= 2/3-filled, use sparse copy algorithm
  if (subset_size * (3 * ONELU) <= orig_quaterarr_size * (2 * ONELU)) {
    uint32_t subset_mask_widx = 0;
    while (1) {
      const uintptr_t cur_include_word = subset_mask[subset_mask_widx];
      if (cur_include_word) {
    uint32_t wordhalf_idx = 0;
#ifdef __LP64__
    uint32_t cur_include_halfword = (uint32_t)cur_include_word;
#else
    uint32_t cur_include_halfword = (uint16_t)cur_include_word;
#endif
    while (1) {
      if (cur_include_halfword) {
        uintptr_t orig_quaterarr_word = main_quaterarr[subset_mask_widx * 2 +
wordhalf_idx]; do { uint32_t rqa_idx_lowbits =
__builtin_ctz(cur_include_halfword); cur_output_word |= ((orig_quaterarr_word >>
(rqa_idx_lowbits * 2)) & 3) << (word_write_halfshift * 2); if
(++word_write_halfshift == BITCT2) { *main_quaterarr_writer++ = cur_output_word;
        word_write_halfshift = 0;
        cur_output_word = 0;
          }
          cur_include_halfword &= cur_include_halfword - 1;
        } while (cur_include_halfword);
      }
      if (wordhalf_idx) {
        break;
      }
      wordhalf_idx++;
#ifdef __LP64__
      cur_include_halfword = cur_include_word >> 32;
#else
      cur_include_halfword = cur_include_word >> 16;
#endif
    }
    if (main_quaterarr_writer == main_quaterarr_write_last) {
      if (word_write_halfshift == word_write_halfshift_end) {
            if (word_write_halfshift_end) {
          *main_quaterarr_writer = cur_output_word;
        }
        return;
      }
    }
      }
      subset_mask_widx++;
    }
  }
  // blocked copy
  while (1) {
    const uintptr_t cur_include_word = *subset_mask++;
    uint32_t wordhalf_idx = 0;
#ifdef __LP64__
    uintptr_t cur_include_halfword = (uint32_t)cur_include_word;
#else
    uint32_t cur_include_halfword = (uint16_t)cur_include_word;
#endif
    while (1) {
      uintptr_t orig_quaterarr_word = *main_quaterarr++;
      while (cur_include_halfword) {
    uint32_t rqa_idx_lowbits = CTZLU(cur_include_halfword);
    uintptr_t halfword_invshifted = (~cur_include_halfword) >> rqa_idx_lowbits;
    uintptr_t orig_quaterarr_curblock_unmasked = orig_quaterarr_word >>
(rqa_idx_lowbits * 2); uint32_t rqa_block_len = CTZLU(halfword_invshifted);
    uint32_t block_len_limit = BITCT2 - word_write_halfshift;
    cur_output_word |= orig_quaterarr_curblock_unmasked << (2 *
word_write_halfshift); if (rqa_block_len < block_len_limit) {
      word_write_halfshift += rqa_block_len;
      cur_output_word &= (ONELU << (word_write_halfshift * 2)) - ONELU;
    } else {
      // no need to mask, extra bits vanish off the high end

      *main_quaterarr_writer++ = cur_output_word;
      word_write_halfshift = rqa_block_len - block_len_limit;
      if (word_write_halfshift) {
        cur_output_word = (orig_quaterarr_curblock_unmasked >> (2 *
block_len_limit)) & ((ONELU << (2 * word_write_halfshift)) - ONELU); } else {
        cur_output_word = 0;
      }
    }
    cur_include_halfword &= (~(ONELU << (rqa_block_len + rqa_idx_lowbits))) +
ONELU;
      }
      if (wordhalf_idx) {
    break;
      }
      wordhalf_idx++;
#ifdef __LP64__
      cur_include_halfword = cur_include_word >> 32;
#else
      cur_include_halfword = cur_include_word >> 16;
#endif
    }
    if (main_quaterarr_writer == main_quaterarr_write_last) {
      if (word_write_halfshift == word_write_halfshift_end) {
    if (word_write_halfshift_end) {
      *main_quaterarr_writer = cur_output_word;
    }
    return;
      }
    }
  }
}
*/

uint32_t load_and_collapse_incl(uint32_t unfiltered_sample_ct,
                                uint32_t sample_ct,
                                const uintptr_t* __restrict sample_include,
                                uintptr_t final_mask, uint32_t do_reverse,
                                FILE* bedfile, uintptr_t* __restrict rawbuf,
                                uintptr_t* __restrict mainbuf)
{
    assert(unfiltered_sample_ct);
    uint32_t unfiltered_sample_ct4 = (unfiltered_sample_ct + 3) / 4;
    if (unfiltered_sample_ct == sample_ct) {
        rawbuf = mainbuf;
    }
    if (load_raw(unfiltered_sample_ct4, bedfile, rawbuf)) {
        return RET_READ_FAIL;
    }
    if (unfiltered_sample_ct != sample_ct) {
        copy_quaterarr_nonempty_subset(
            rawbuf, sample_include, unfiltered_sample_ct, sample_ct, mainbuf);
    }
    else
    {
        mainbuf[(unfiltered_sample_ct - 1) / BITCT2] &= final_mask;
    }
    if (do_reverse) {
        reverse_loadbuf(sample_ct, (unsigned char*) mainbuf);
    }
    // mainbuf should contains the information
    return 0;
}


void hh_reset(unsigned char* loadbuf, uintptr_t* sample_include_quaterarr,
              uintptr_t unfiltered_sample_ct)
{
    uintptr_t sample_bidx = 0;
    unsigned char* loadbuf_end = &(loadbuf[(unfiltered_sample_ct + 3) / 4]);
    unsigned char* iicp;
    unsigned char ucc;
    unsigned char ucc2;
    uintptr_t unfiltered_sample_ctd;
    uint32_t* loadbuf_alias32;
    uint32_t uii;
    uint32_t ujj;
#ifdef __LP64__
    uint32_t* sample_include_quaterarr_alias32;
    __m128i* loadbuf_alias;
    __m128i* iivp;
    __m128i vii;
    __m128i vjj;
    if (!(((uintptr_t) loadbuf) & 15)) {
        loadbuf_alias = (__m128i*) loadbuf;
        iivp = (__m128i*) sample_include_quaterarr;
        unfiltered_sample_ctd = unfiltered_sample_ct / 64;
        for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
            vii = *loadbuf_alias;
            vjj = _mm_and_si128(_mm_andnot_si128(vii, _mm_srli_epi64(vii, 1)),
                                *iivp++);
            *loadbuf_alias++ = _mm_sub_epi64(vii, vjj);
        }
        loadbuf = (unsigned char*) loadbuf_alias;
        iicp = (unsigned char*) iivp;
    }
    else if (!(((uintptr_t) loadbuf) & 3))
    {
        loadbuf_alias32 = (uint32_t*) loadbuf;
        sample_include_quaterarr_alias32 = (uint32_t*) sample_include_quaterarr;
        unfiltered_sample_ctd = unfiltered_sample_ct / BITCT2;
        for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
            uii = *loadbuf_alias32;
            ujj = ((uii >> 1) & (~uii)) & (*sample_include_quaterarr_alias32++);
            *loadbuf_alias32++ = uii - ujj;
        }
        loadbuf = (unsigned char*) loadbuf_alias32;
        iicp = (unsigned char*) sample_include_quaterarr_alias32;
    }
    else
    {
        iicp = (unsigned char*) sample_include_quaterarr;
    }
#else
    if (!(((uintptr_t) loadbuf) & 3)) {
        loadbuf_alias32 = (uint32_t*) loadbuf;
        unfiltered_sample_ctd = unfiltered_sample_ct / BITCT2;
        for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
            uii = *loadbuf_alias32;
            ujj = ((uii >> 1) & (~uii)) & (*sample_include_quaterarr++);
            *loadbuf_alias32++ = uii - ujj;
        }
        loadbuf = (unsigned char*) loadbuf_alias32;
    }
    iicp = (unsigned char*) sample_include_quaterarr;
#endif
    for (; loadbuf < loadbuf_end;) {
        ucc = *loadbuf;
        ucc2 = ((ucc >> 1) & (~ucc)) & (*iicp++);
        *loadbuf++ = ucc - ucc2;
    }
}

void hh_reset_y(unsigned char* loadbuf, uintptr_t* sample_include_quaterarr,
                uintptr_t* sample_male_include_quaterarr,
                uintptr_t unfiltered_sample_ct)
{
    uintptr_t sample_bidx = 0;
    unsigned char* loadbuf_end = &(loadbuf[(unfiltered_sample_ct + 3) / 4]);
    unsigned char* iicp;
    unsigned char* imicp;
    unsigned char ucc;
    unsigned char ucc2;
    unsigned char ucc3;
    uintptr_t unfiltered_sample_ctd;
    uint32_t* loadbuf_alias32;
    uint32_t uii;
    uint32_t ujj;
    uint32_t ukk;
#ifdef __LP64__
    const __m128i m1 = {FIVEMASK, FIVEMASK};
    uint32_t* sample_include_quaterarr_alias32;
    uint32_t* sample_male_include_quaterarr_alias32;
    __m128i* loadbuf_alias;
    __m128i* iivp;
    __m128i* imivp;
    __m128i vii;
    __m128i vjj;
    __m128i vkk;
    if (!(((uintptr_t) loadbuf) & 15)) {
        loadbuf_alias = (__m128i*) loadbuf;
        iivp = (__m128i*) sample_include_quaterarr;
        imivp = (__m128i*) sample_male_include_quaterarr;
        unfiltered_sample_ctd = unfiltered_sample_ct / 64;
        for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
            // sample_include_quaterarr & ~sample_male_include_quaterarr: force
            // to 01 sample_male_include_quaterarr: convert 10 to 01, keep
            // everything else
            vii = *imivp++;
            vjj = *iivp++;
            vkk = _mm_and_si128(*loadbuf_alias,
                                _mm_or_si128(vii, _mm_slli_epi64(vii, 1)));
            *loadbuf_alias++ = _mm_or_si128(
                _mm_andnot_si128(vii, vjj),
                _mm_sub_epi64(
                    vkk,
                    _mm_and_si128(_mm_andnot_si128(vkk, _mm_srli_epi64(vkk, 1)),
                                  m1)));
        }
        loadbuf = (unsigned char*) loadbuf_alias;
        iicp = (unsigned char*) iivp;
        imicp = (unsigned char*) imivp;
    }
    else if (!(((uintptr_t) loadbuf) & 3))
    {
        loadbuf_alias32 = (uint32_t*) loadbuf;
        sample_include_quaterarr_alias32 = (uint32_t*) sample_include_quaterarr;
        sample_male_include_quaterarr_alias32 =
            (uint32_t*) sample_male_include_quaterarr;
        unfiltered_sample_ctd = unfiltered_sample_ct / 16;
        for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
            uii = *sample_male_include_quaterarr_alias32++;
            ujj = *sample_include_quaterarr_alias32++;
            ukk = (*loadbuf_alias32) & (uii * 3);
            *loadbuf_alias32++ =
                ((~uii) & ujj) | (ukk - ((~ukk) & (ukk >> 1) & 0x55555555));
        }
        loadbuf = (unsigned char*) loadbuf_alias32;
        iicp = (unsigned char*) sample_include_quaterarr_alias32;
        imicp = (unsigned char*) sample_male_include_quaterarr_alias32;
    }
    else
    {
        iicp = (unsigned char*) sample_include_quaterarr;
        imicp = (unsigned char*) sample_male_include_quaterarr;
    }
#else
    if (!(((uintptr_t) loadbuf) & 3)) {
        loadbuf_alias32 = (uint32_t*) loadbuf;
        unfiltered_sample_ctd = unfiltered_sample_ct / 16;
        for (; sample_bidx < unfiltered_sample_ctd; sample_bidx++) {
            uii = *sample_male_include_quaterarr++;
            ujj = *sample_include_quaterarr++;
            ukk = (*loadbuf_alias32) & (uii * 3);
            *loadbuf_alias32++ =
                ((~uii) & ujj) | (ukk - ((~ukk) & (ukk >> 1) & 0x55555555));
        }
        loadbuf = (unsigned char*) loadbuf_alias32;
    }
    iicp = (unsigned char*) sample_include_quaterarr;
    imicp = (unsigned char*) sample_male_include_quaterarr;
#endif
    for (; loadbuf < loadbuf_end;) {
        ucc = *imicp++;
        ucc2 = *iicp++;
        ucc3 = (*loadbuf) & (ucc * 3);
        *loadbuf++ = ((~ucc) & ucc2) | (ucc3 - ((~ucc3) & (ucc3 >> 1) & 0x55));
    }
}

uint32_t cubic_real_roots(double coef_a, double coef_b, double coef_c,
                          double* solutions)
{
    // Analytically finds all real roots of x^3 + ax^2 + bx + c, saving them in
    // solutions[] (sorted from smallest to largest), and returning the count.
    // Multiple roots are only returned/counted once.
    // Additional research into numerical stability may be in order here.
    double a2 = coef_a * coef_a;
    double qq = (a2 - 3 * coef_b) * (1.0 / 9.0);
    double rr =
        (2 * a2 * coef_a - 9 * coef_a * coef_b + 27 * coef_c) * (1.0 / 54.0);
    double r2 = rr * rr;
    double q3 = qq * qq * qq;
    double adiv3 = coef_a * (1.0 / 3.0);
    double sq;
    double dxx;
    if (r2 < q3) {
        // three real roots
        sq = sqrt(qq);
        dxx = acos(rr / (qq * sq)) * (1.0 / 3.0);
        sq *= -2;
        solutions[0] = sq * cos(dxx) - adiv3;
        solutions[1] = sq * cos(dxx + (2.0 * PI / 3.0)) - adiv3;
        solutions[2] = sq * cos(dxx - (2.0 * PI / 3.0)) - adiv3;
        // now sort and check for within-epsilon equality
        if (solutions[0] > solutions[1]) {
            dxx = solutions[0];
            solutions[0] = solutions[1];
            if (dxx > solutions[2]) {
                solutions[1] = solutions[2];
                solutions[2] = dxx;
            }
            else
            {
                solutions[1] = dxx;
            }
            if (solutions[0] > solutions[1]) {
                dxx = solutions[0];
                solutions[0] = solutions[1];
                solutions[1] = dxx;
            }
        }
        else if (solutions[1] > solutions[2])
        {
            dxx = solutions[1];
            solutions[1] = solutions[2];
            solutions[2] = dxx;
        }
        if (solutions[1] - solutions[0] < EPSILON) {
            solutions[1] = solutions[2];
            return (solutions[1] - solutions[0] < EPSILON) ? 1 : 2;
        }
        return (solutions[2] - solutions[1] < EPSILON) ? 2 : 3;
    }
    dxx = -pow(fabs(rr) + sqrt(r2 - q3), 1.0 / 3.0);
    if (dxx == 0.0) {
        solutions[0] = -adiv3;
        return 1;
    }
    if (rr < 0.0) {
        dxx = -dxx;
    }
    sq = qq / dxx;
    solutions[0] = dxx + sq - adiv3;
    // use of regular epsilon here has actually burned us
    if (fabs(dxx - sq) >= (EPSILON * 8)) {
        return 1;
    }
    if (dxx >= 0.0) {
        solutions[1] = solutions[0];
        solutions[0] = -dxx - adiv3;
    }
    else
    {
        solutions[1] = -dxx - adiv3;
    }
    return 2;
}


void genovec_3freq(const uintptr_t* __restrict geno_vec,
                   const uintptr_t* __restrict include_quatervec,
                   uintptr_t sample_ctl2, uint32_t* __restrict missing_ctp,
                   uint32_t* __restrict het_ctp,
                   uint32_t* __restrict homset_ctp)
{
    // generic routine for getting all counts.
    const uintptr_t* geno_vec_end = &(geno_vec[sample_ctl2]);
    uintptr_t loader;
    uintptr_t loader2;
    uintptr_t loader3;
    uint32_t acc_even = 0;
    uint32_t acc_odd = 0;
    uint32_t acc_and = 0;
#ifdef __LP64__
    uintptr_t cur_decr = 120;
    const uintptr_t* geno_vec_12x_end;
    sample_ctl2 -= sample_ctl2 % 12;
    while (sample_ctl2 >= 120) {
    genovec_3freq_loop:
        geno_vec_12x_end = &(geno_vec[cur_decr]);
        count_3freq_1920b(
            (const __m128i*) geno_vec, (const __m128i*) geno_vec_12x_end,
            (const __m128i*) include_quatervec, &acc_even, &acc_odd, &acc_and);
        geno_vec = geno_vec_12x_end;
        include_quatervec = &(include_quatervec[cur_decr]);
        sample_ctl2 -= cur_decr;
    }
    if (sample_ctl2) {
        cur_decr = sample_ctl2;
        goto genovec_3freq_loop;
    }
#else
    const uintptr_t* geno_vec_twelve_end =
        &(geno_vec[sample_ctl2 - (sample_ctl2 % 12)]);
    while (geno_vec < geno_vec_twelve_end) {
        count_3freq_48b(geno_vec, include_quatervec, &acc_even, &acc_odd,
                        &acc_and);
        geno_vec = &(geno_vec[12]);
        include_quatervec = &(include_quatervec[12]);
    }
#endif
    while (geno_vec < geno_vec_end) {
        loader = *geno_vec++;
        loader2 = *include_quatervec++;
        loader3 = loader2 & (loader >> 1);
        acc_even += popcount2_long(loader & loader2);
        acc_odd += popcount2_long(loader3);
        acc_and += popcount2_long(loader & loader3);
    }
    *missing_ctp = acc_even - acc_and;
    *het_ctp = acc_odd - acc_and;
    *homset_ctp = acc_and;
}


#ifdef __LP64__
void count_2freq_dbl_960b(
    const VECITYPE* geno_vvec, const VECITYPE* geno_vvec_end,
    const VECITYPE* __restrict mask1vp, const VECITYPE* __restrict mask2vp,
    uint32_t* __restrict ct1abp, uint32_t* __restrict ct1cp,
    uint32_t* __restrict ct2abp, uint32_t* __restrict ct2cp)
{
    const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
    const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
    __m128i loader;
    __m128i loader2;
    __m128i loader3;
    __m128i to_ct1_ab;
    __m128i to_ct_abtmp;
    __m128i to_ct1_c;
    __m128i to_ct2_ab;
    __m128i to_ct2_c;
    __univec acc1_ab;
    __univec acc1_c;
    __univec acc2_ab;
    __univec acc2_c;

    acc1_ab.vi = _mm_setzero_si128();
    acc1_c.vi = _mm_setzero_si128();
    acc2_ab.vi = _mm_setzero_si128();
    acc2_c.vi = _mm_setzero_si128();
    do
    {
        loader = *geno_vvec++;
        loader2 = *mask1vp++;
        loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
        loader2 = _mm_and_si128(loader2, loader);
        to_ct1_ab = _mm_add_epi64(loader3, loader2);
        to_ct1_c = _mm_andnot_si128(loader3, loader2);
        loader2 = *mask2vp++;
        loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
        loader2 = _mm_and_si128(loader2, loader);
        to_ct2_ab = _mm_add_epi64(loader3, loader2);
        to_ct2_c = _mm_andnot_si128(loader3, loader2);
        to_ct1_ab =
            _mm_add_epi64(_mm_and_si128(to_ct1_ab, m2),
                          _mm_and_si128(_mm_srli_epi64(to_ct1_ab, 2), m2));
        to_ct2_ab =
            _mm_add_epi64(_mm_and_si128(to_ct2_ab, m2),
                          _mm_and_si128(_mm_srli_epi64(to_ct2_ab, 2), m2));

        loader = *geno_vvec++;
        loader2 = *mask1vp++;
        loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
        loader2 = _mm_and_si128(loader2, loader);
        to_ct_abtmp = _mm_add_epi64(loader3, loader2);
        to_ct1_c = _mm_add_epi64(to_ct1_c, _mm_andnot_si128(loader3, loader2));
        to_ct1_ab = _mm_add_epi64(
            to_ct1_ab,
            _mm_add_epi64(_mm_and_si128(to_ct_abtmp, m2),
                          _mm_and_si128(_mm_srli_epi64(to_ct_abtmp, 2), m2)));
        loader2 = *mask2vp++;
        loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
        loader2 = _mm_and_si128(loader2, loader);
        to_ct_abtmp = _mm_add_epi64(loader3, loader2);
        to_ct2_c = _mm_add_epi64(to_ct2_c, _mm_andnot_si128(loader3, loader2));
        to_ct2_ab = _mm_add_epi64(
            to_ct2_ab,
            _mm_add_epi64(_mm_and_si128(to_ct_abtmp, m2),
                          _mm_and_si128(_mm_srli_epi64(to_ct_abtmp, 2), m2)));

        loader = *geno_vvec++;
        loader2 = *mask1vp++;
        loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
        loader2 = _mm_and_si128(loader2, loader);
        to_ct_abtmp = _mm_add_epi64(loader3, loader2);
        to_ct1_c = _mm_add_epi64(to_ct1_c, _mm_andnot_si128(loader3, loader2));
        to_ct1_ab = _mm_add_epi64(
            to_ct1_ab,
            _mm_add_epi64(_mm_and_si128(to_ct_abtmp, m2),
                          _mm_and_si128(_mm_srli_epi64(to_ct_abtmp, 2), m2)));
        loader2 = *mask2vp++;
        loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
        loader2 = _mm_and_si128(loader2, loader);
        to_ct_abtmp = _mm_add_epi64(loader3, loader2);
        to_ct2_c = _mm_add_epi64(to_ct2_c, _mm_andnot_si128(loader3, loader2));
        to_ct2_ab = _mm_add_epi64(
            to_ct2_ab,
            _mm_add_epi64(_mm_and_si128(to_ct_abtmp, m2),
                          _mm_and_si128(_mm_srli_epi64(to_ct_abtmp, 2), m2)));

        to_ct1_c =
            _mm_add_epi64(_mm_and_si128(to_ct1_c, m2),
                          _mm_and_si128(_mm_srli_epi64(to_ct1_c, 2), m2));
        to_ct2_c =
            _mm_add_epi64(_mm_and_si128(to_ct2_c, m2),
                          _mm_and_si128(_mm_srli_epi64(to_ct2_c, 2), m2));

        acc1_ab.vi = _mm_add_epi64(
            acc1_ab.vi,
            _mm_add_epi64(_mm_and_si128(to_ct1_ab, m4),
                          _mm_and_si128(_mm_srli_epi64(to_ct1_ab, 4), m4)));
        acc1_c.vi = _mm_add_epi64(
            acc1_c.vi,
            _mm_add_epi64(_mm_and_si128(to_ct1_c, m4),
                          _mm_and_si128(_mm_srli_epi64(to_ct1_c, 4), m4)));
        acc2_ab.vi = _mm_add_epi64(
            acc2_ab.vi,
            _mm_add_epi64(_mm_and_si128(to_ct2_ab, m4),
                          _mm_and_si128(_mm_srli_epi64(to_ct2_ab, 4), m4)));
        acc2_c.vi = _mm_add_epi64(
            acc2_c.vi,
            _mm_add_epi64(_mm_and_si128(to_ct2_c, m4),
                          _mm_and_si128(_mm_srli_epi64(to_ct2_c, 4), m4)));
    } while (geno_vvec < geno_vvec_end);
    const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
    acc1_ab.vi =
        _mm_add_epi64(_mm_and_si128(acc1_ab.vi, m8),
                      _mm_and_si128(_mm_srli_epi64(acc1_ab.vi, 8), m8));
    acc1_c.vi = _mm_and_si128(
        _mm_add_epi64(acc1_c.vi, _mm_srli_epi64(acc1_c.vi, 8)), m8);
    acc2_ab.vi =
        _mm_add_epi64(_mm_and_si128(acc2_ab.vi, m8),
                      _mm_and_si128(_mm_srli_epi64(acc2_ab.vi, 8), m8));
    acc2_c.vi = _mm_and_si128(
        _mm_add_epi64(acc2_c.vi, _mm_srli_epi64(acc2_c.vi, 8)), m8);
    *ct1abp += ((acc1_ab.u8[0] + acc1_ab.u8[1]) * 0x1000100010001LLU) >> 48;
    *ct1cp += ((acc1_c.u8[0] + acc1_c.u8[1]) * 0x1000100010001LLU) >> 48;
    *ct2abp += ((acc2_ab.u8[0] + acc2_ab.u8[1]) * 0x1000100010001LLU) >> 48;
    *ct2cp += ((acc2_c.u8[0] + acc2_c.u8[1]) * 0x1000100010001LLU) >> 48;
}

void count_3freq_1920b(const VECITYPE* geno_vvec, const VECITYPE* geno_vvec_end,
                       const VECITYPE* __restrict maskvp,
                       uint32_t* __restrict even_ctp,
                       uint32_t* __restrict odd_ctp,
                       uint32_t* __restrict homset_ctp)
{
    const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
    const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
    __m128i loader;
    __m128i loader2;
    __m128i loader3;
    __m128i even1;
    __m128i odd1;
    __m128i homset1;
    __m128i even2;
    __m128i odd2;
    __m128i homset2;
    __univec acc_even;
    __univec acc_odd;
    __univec acc_homset;

    acc_even.vi = _mm_setzero_si128();
    acc_odd.vi = _mm_setzero_si128();
    acc_homset.vi = _mm_setzero_si128();
    do
    {
        loader = *geno_vvec++;
        loader2 = *maskvp++;
        odd1 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
        even1 = _mm_and_si128(loader2, loader);
        homset1 = _mm_and_si128(odd1, loader);
        loader = *geno_vvec++;
        loader2 = *maskvp++;
        loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
        even1 = _mm_add_epi64(even1, _mm_and_si128(loader2, loader));
        odd1 = _mm_add_epi64(odd1, loader3);
        homset1 = _mm_add_epi64(homset1, _mm_and_si128(loader3, loader));
        loader = *geno_vvec++;
        loader2 = *maskvp++;
        loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
        even1 = _mm_add_epi64(even1, _mm_and_si128(loader2, loader));
        odd1 = _mm_add_epi64(odd1, loader3);
        homset1 = _mm_add_epi64(homset1, _mm_and_si128(loader3, loader));

        even1 = _mm_add_epi64(_mm_and_si128(even1, m2),
                              _mm_and_si128(_mm_srli_epi64(even1, 2), m2));
        odd1 = _mm_add_epi64(_mm_and_si128(odd1, m2),
                             _mm_and_si128(_mm_srli_epi64(odd1, 2), m2));
        homset1 = _mm_add_epi64(_mm_and_si128(homset1, m2),
                                _mm_and_si128(_mm_srli_epi64(homset1, 2), m2));

        loader = *geno_vvec++;
        loader2 = *maskvp++;
        odd2 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
        even2 = _mm_and_si128(loader2, loader);
        homset2 = _mm_and_si128(odd2, loader);
        loader = *geno_vvec++;
        loader2 = *maskvp++;
        loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
        even2 = _mm_add_epi64(even2, _mm_and_si128(loader2, loader));
        odd2 = _mm_add_epi64(odd2, loader3);
        homset2 = _mm_add_epi64(homset2, _mm_and_si128(loader3, loader));
        loader = *geno_vvec++;
        loader2 = *maskvp++;
        loader3 = _mm_and_si128(loader2, _mm_srli_epi64(loader, 1));
        even2 = _mm_add_epi64(even2, _mm_and_si128(loader2, loader));
        odd2 = _mm_add_epi64(odd2, loader3);
        homset2 = _mm_add_epi64(homset2, _mm_and_si128(loader3, loader));

        even1 = _mm_add_epi64(
            even1, _mm_add_epi64(_mm_and_si128(even2, m2),
                                 _mm_and_si128(_mm_srli_epi64(even2, 2), m2)));
        odd1 = _mm_add_epi64(
            odd1, _mm_add_epi64(_mm_and_si128(odd2, m2),
                                _mm_and_si128(_mm_srli_epi64(odd2, 2), m2)));
        homset1 = _mm_add_epi64(
            homset1,
            _mm_add_epi64(_mm_and_si128(homset2, m2),
                          _mm_and_si128(_mm_srli_epi64(homset2, 2), m2)));

        acc_even.vi = _mm_add_epi64(
            acc_even.vi,
            _mm_add_epi64(_mm_and_si128(even1, m4),
                          _mm_and_si128(_mm_srli_epi64(even1, 4), m4)));
        acc_odd.vi = _mm_add_epi64(
            acc_odd.vi,
            _mm_add_epi64(_mm_and_si128(odd1, m4),
                          _mm_and_si128(_mm_srli_epi64(odd1, 4), m4)));
        acc_homset.vi = _mm_add_epi64(
            acc_homset.vi,
            _mm_add_epi64(_mm_and_si128(homset1, m4),
                          _mm_and_si128(_mm_srli_epi64(homset1, 4), m4)));
    } while (geno_vvec < geno_vvec_end);
    const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
    acc_even.vi =
        _mm_add_epi64(_mm_and_si128(acc_even.vi, m8),
                      _mm_and_si128(_mm_srli_epi64(acc_even.vi, 8), m8));
    acc_odd.vi =
        _mm_add_epi64(_mm_and_si128(acc_odd.vi, m8),
                      _mm_and_si128(_mm_srli_epi64(acc_odd.vi, 8), m8));
    acc_homset.vi =
        _mm_add_epi64(_mm_and_si128(acc_homset.vi, m8),
                      _mm_and_si128(_mm_srli_epi64(acc_homset.vi, 8), m8));
    *even_ctp += ((acc_even.u8[0] + acc_even.u8[1]) * 0x1000100010001LLU) >> 48;
    *odd_ctp += ((acc_odd.u8[0] + acc_odd.u8[1]) * 0x1000100010001LLU) >> 48;
    *homset_ctp +=
        ((acc_homset.u8[0] + acc_homset.u8[1]) * 0x1000100010001LLU) >> 48;
}
#else
void count_2freq_dbl_24b(const uintptr_t* __restrict geno_vec,
                         const uintptr_t* __restrict mask1p,
                         const uintptr_t* __restrict mask2p,
                         uint32_t* __restrict ct1abp,
                         uint32_t* __restrict ct1cp,
                         uint32_t* __restrict ct2abp,
                         uint32_t* __restrict ct2cp)
{
    uintptr_t loader = *geno_vec++;
    uintptr_t loader2 = *mask1p++;
    uintptr_t loader3 = (loader >> 1) & loader2;
    uintptr_t to_ct1_ab;
    uintptr_t to_ct1_c;
    uintptr_t to_ct2_ab;
    uintptr_t to_ct2_c;
    uintptr_t to_ct_abtmp;
    uintptr_t partial1_ab;
    uintptr_t partial1_c;
    uintptr_t partial2_ab;
    uintptr_t partial2_c;
    loader2 &= loader;
    to_ct1_ab = loader2 + loader3;
    to_ct1_c = loader2 & (~loader3);
    loader2 = *mask2p++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    to_ct2_ab = loader2 + loader3;
    to_ct2_c = loader2 & (~loader3);

    to_ct1_ab = (to_ct1_ab & 0x33333333) + ((to_ct1_ab >> 2) & 0x33333333);
    to_ct2_ab = (to_ct2_ab & 0x33333333) + ((to_ct2_ab >> 2) & 0x33333333);

    loader = *geno_vec++;
    loader2 = *mask1p++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    to_ct_abtmp = loader2 + loader3;
    to_ct1_c += loader2 & (~loader3);
    to_ct1_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);
    loader2 = *mask2p++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    to_ct_abtmp = loader2 + loader3;
    to_ct2_c += loader2 & (~loader3);
    to_ct2_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);

    loader = *geno_vec++;
    loader2 = *mask1p++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    to_ct_abtmp = loader2 + loader3;
    to_ct1_c += loader2 & (~loader3);
    to_ct1_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);
    loader2 = *mask2p++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    to_ct_abtmp = loader2 + loader3;
    to_ct2_c += loader2 & (~loader3);
    to_ct2_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);

    partial1_ab = (to_ct1_ab & 0x0f0f0f0f) + ((to_ct1_ab >> 4) & 0x0f0f0f0f);
    partial1_c = (to_ct1_c & 0x33333333) + ((to_ct1_c >> 2) & 0x33333333);
    partial2_ab = (to_ct2_ab & 0x0f0f0f0f) + ((to_ct2_ab >> 4) & 0x0f0f0f0f);
    partial2_c = (to_ct2_c & 0x33333333) + ((to_ct2_c >> 2) & 0x33333333);

    loader = *geno_vec++;
    loader2 = *mask1p++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    to_ct1_ab = loader2 + loader3;
    to_ct1_c = loader2 & (~loader3);
    loader2 = *mask2p++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    to_ct2_ab = loader2 + loader3;
    to_ct2_c = loader2 & (~loader3);

    to_ct1_ab = (to_ct1_ab & 0x33333333) + ((to_ct1_ab >> 2) & 0x33333333);
    to_ct2_ab = (to_ct2_ab & 0x33333333) + ((to_ct2_ab >> 2) & 0x33333333);

    loader = *geno_vec++;
    loader2 = *mask1p++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    to_ct_abtmp = loader2 + loader3;
    to_ct1_c += loader2 & (~loader3);
    to_ct1_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);
    loader2 = *mask2p++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    to_ct_abtmp = loader2 + loader3;
    to_ct2_c += loader2 & (~loader3);
    to_ct2_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);

    loader = *geno_vec++;
    loader2 = *mask1p++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    to_ct_abtmp = loader2 + loader3;
    to_ct1_c += loader2 & (~loader3);
    to_ct1_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);
    loader2 = *mask2p++;
    loader3 = (loader >> 1) & loader2;
    loader2 &= loader;
    to_ct_abtmp = loader2 + loader3;
    to_ct2_c += loader2 & (~loader3);
    to_ct2_ab += (to_ct_abtmp & 0x33333333) + ((to_ct_abtmp >> 2) & 0x33333333);

    partial1_ab += (to_ct1_ab & 0x0f0f0f0f) + ((to_ct1_ab >> 4) & 0x0f0f0f0f);
    partial1_c += (to_ct1_c & 0x33333333) + ((to_ct1_c >> 2) & 0x33333333);
    partial2_ab += (to_ct2_ab & 0x0f0f0f0f) + ((to_ct2_ab >> 4) & 0x0f0f0f0f);
    partial2_c += (to_ct2_c & 0x33333333) + ((to_ct2_c >> 2) & 0x33333333);

    partial1_c = (partial1_c & 0x0f0f0f0f) + ((partial1_c >> 4) & 0x0f0f0f0f);
    partial2_c = (partial2_c & 0x0f0f0f0f) + ((partial2_c >> 4) & 0x0f0f0f0f);

    *ct1abp += (partial1_ab * 0x01010101) >> 24;
    *ct1cp += (partial1_c * 0x01010101) >> 24;
    *ct2abp += (partial2_ab * 0x01010101) >> 24;
    *ct2cp += (partial2_c * 0x01010101) >> 24;
}

void count_3freq_48b(const uintptr_t* __restrict geno_vec,
                     const uintptr_t* __restrict maskp,
                     uint32_t* __restrict ctap, uint32_t* __restrict ctbp,
                     uint32_t* __restrict ctcp)
{
    uintptr_t loader = *geno_vec++;
    uintptr_t loader2 = *maskp++;
    uint32_t to_ct_a1 = loader & loader2;
    uint32_t to_ct_b1 = (loader >> 1) & loader2;
    uint32_t to_ct_c1 = loader & to_ct_b1;
    uintptr_t loader3;
    uint32_t to_ct_a2;
    uint32_t to_ct_b2;
    uint32_t to_ct_c2;
    uintptr_t partial_a;
    uintptr_t partial_b;
    uintptr_t partial_c;
    loader = *geno_vec++;
    loader2 = *maskp++;
    loader3 = (loader >> 1) & loader2;
    to_ct_a1 += loader & loader2;
    to_ct_b1 += loader3;
    to_ct_c1 += loader & loader3;
    loader = *geno_vec++;
    loader2 = *maskp++;
    loader3 = (loader >> 1) & loader2;
    to_ct_a1 += loader & loader2;
    to_ct_b1 += loader3;
    to_ct_c1 += loader & loader3;

    loader = *geno_vec++;
    loader2 = *maskp++;
    to_ct_a2 = loader & loader2;
    to_ct_b2 = (loader >> 1) & loader2;
    to_ct_c2 = loader & to_ct_b2;
    loader = *geno_vec++;
    loader2 = *maskp++;
    loader3 = (loader >> 1) & loader2;
    to_ct_a2 += loader & loader2;
    to_ct_b2 += loader3;
    to_ct_c2 += loader & loader3;
    loader = *geno_vec++;
    loader2 = *maskp++;
    loader3 = (loader >> 1) & loader2;
    to_ct_a2 += loader & loader2;
    to_ct_b2 += loader3;
    to_ct_c2 += loader & loader3;

    to_ct_a1 = (to_ct_a1 & 0x33333333) + ((to_ct_a1 >> 2) & 0x33333333);
    to_ct_a1 += (to_ct_a2 & 0x33333333) + ((to_ct_a2 >> 2) & 0x33333333);
    partial_a = (to_ct_a1 & 0x0f0f0f0f) + ((to_ct_a1 >> 4) & 0x0f0f0f0f);
    to_ct_b1 = (to_ct_b1 & 0x33333333) + ((to_ct_b1 >> 2) & 0x33333333);
    to_ct_b1 += (to_ct_b2 & 0x33333333) + ((to_ct_b2 >> 2) & 0x33333333);
    partial_b = (to_ct_b1 & 0x0f0f0f0f) + ((to_ct_b1 >> 4) & 0x0f0f0f0f);
    to_ct_c1 = (to_ct_c1 & 0x33333333) + ((to_ct_c1 >> 2) & 0x33333333);
    to_ct_c1 += (to_ct_c2 & 0x33333333) + ((to_ct_c2 >> 2) & 0x33333333);
    partial_c = (to_ct_c1 & 0x0f0f0f0f) + ((to_ct_c1 >> 4) & 0x0f0f0f0f);

    loader = *geno_vec++;
    loader2 = *maskp++;
    to_ct_a1 = loader & loader2;
    to_ct_b1 = (loader >> 1) & loader2;
    to_ct_c1 = loader & to_ct_b1;
    loader = *geno_vec++;
    loader2 = *maskp++;
    loader3 = (loader >> 1) & loader2;
    to_ct_a1 += loader & loader2;
    to_ct_b1 += loader3;
    to_ct_c1 += loader & loader3;
    loader = *geno_vec++;
    loader2 = *maskp++;
    loader3 = (loader >> 1) & loader2;
    to_ct_a1 += loader & loader2;
    to_ct_b1 += loader3;
    to_ct_c1 += loader & loader3;

    loader = *geno_vec++;
    loader2 = *maskp++;
    to_ct_a2 = loader & loader2;
    to_ct_b2 = (loader >> 1) & loader2;
    to_ct_c2 = loader & to_ct_b2;
    loader = *geno_vec++;
    loader2 = *maskp++;
    loader3 = (loader >> 1) & loader2;
    to_ct_a2 += loader & loader2;
    to_ct_b2 += loader3;
    to_ct_c2 += loader & loader3;
    loader = *geno_vec;
    loader2 = *maskp;
    loader3 = (loader >> 1) & loader2;
    to_ct_a2 += loader & loader2;
    to_ct_b2 += loader3;
    to_ct_c2 += loader & loader3;

    to_ct_a1 = (to_ct_a1 & 0x33333333) + ((to_ct_a1 >> 2) & 0x33333333);
    to_ct_a1 += (to_ct_a2 & 0x33333333) + ((to_ct_a2 >> 2) & 0x33333333);
    partial_a += (to_ct_a1 & 0x0f0f0f0f) + ((to_ct_a1 >> 4) & 0x0f0f0f0f);
    to_ct_b1 = (to_ct_b1 & 0x33333333) + ((to_ct_b1 >> 2) & 0x33333333);
    to_ct_b1 += (to_ct_b2 & 0x33333333) + ((to_ct_b2 >> 2) & 0x33333333);
    partial_b += (to_ct_b1 & 0x0f0f0f0f) + ((to_ct_b1 >> 4) & 0x0f0f0f0f);
    to_ct_c1 = (to_ct_c1 & 0x33333333) + ((to_ct_c1 >> 2) & 0x33333333);
    to_ct_c1 += (to_ct_c2 & 0x33333333) + ((to_ct_c2 >> 2) & 0x33333333);
    partial_c += (to_ct_c1 & 0x0f0f0f0f) + ((to_ct_c1 >> 4) & 0x0f0f0f0f);

    *ctap += (partial_a * 0x01010101) >> 24;
    *ctbp += (partial_b * 0x01010101) >> 24;
    *ctcp += (partial_c * 0x01010101) >> 24;
}
#endif

void fill_quatervec_55(uint32_t ct, uintptr_t* quatervec)
{
    uint32_t rem = ct & (BITCT - 1);
#ifdef __LP64__
    const __m128i m1 = {FIVEMASK, FIVEMASK};
    __m128i* vecp = (__m128i*) quatervec;
    __m128i* vec_end = (__m128i*) (&(quatervec[2 * (ct / BITCT)]));
    uintptr_t* second_to_last;
    while (vecp < vec_end) {
        *vecp++ = m1;
    }
    if (rem) {
        second_to_last = (uintptr_t*) vecp;
        if (rem > BITCT2) {
            second_to_last[0] = FIVEMASK;
            second_to_last[1] = FIVEMASK >> ((BITCT - rem) * 2);
        }
        else
        {
            second_to_last[0] = FIVEMASK >> ((BITCT2 - rem) * 2);
            second_to_last[1] = 0;
        }
    }
#else
    uintptr_t* vec_end = &(quatervec[2 * (ct / BITCT)]);
    while (quatervec < vec_end) {
        *quatervec++ = FIVEMASK;
    }
    if (rem) {
        if (rem > BITCT2) {
            quatervec[0] = FIVEMASK;
            quatervec[1] = FIVEMASK >> ((BITCT - rem) * 2);
        }
        else
        {
            quatervec[0] = FIVEMASK >> ((BITCT2 - rem) * 2);
            quatervec[1] = 0;
        }
    }
#endif
}


void vec_datamask(uintptr_t unfiltered_sample_ct, uint32_t matchval,
                  uintptr_t* data_ptr, uintptr_t* mask_ptr,
                  uintptr_t* result_ptr)
{
    // vec_ptr assumed to be standard 00/01 bit vector
    // sets result_vec bits to 01 iff data_ptr bits are equal to matchval and
    // vec_ptr bit is set, 00 otherwise.
    // currently assumes matchval is not 1.
    assert(unfiltered_sample_ct);
#ifdef __LP64__
    __m128i* data_read = (__m128i*) data_ptr;
    __m128i* mask_read = (__m128i*) mask_ptr;
    __m128i* data_read_end =
        &(data_read[QUATERCT_TO_VECCT(unfiltered_sample_ct)]);
    __m128i* writer = (__m128i*) result_ptr;
    __m128i loader;
#else
    uintptr_t* data_read_end =
        &(data_ptr[QUATERCT_TO_ALIGNED_WORDCT(unfiltered_sample_ct)]);
    uintptr_t loader;
#endif
    if (matchval) {
        if (matchval == 2) {
#ifdef __LP64__
            do
            {
                loader = *data_read++;
                *writer++ = _mm_and_si128(
                    _mm_andnot_si128(loader, _mm_srli_epi64(loader, 1)),
                    *mask_read++);
            } while (data_read < data_read_end);
#else
            do
            {
                loader = *data_ptr++;
                *result_ptr++ = (~loader) & (loader >> 1) & (*mask_ptr++);
            } while (data_ptr < data_read_end);
#endif
        }
        else
        {
#ifdef __LP64__
            do
            {
                loader = *data_read++;
                *writer++ = _mm_and_si128(
                    _mm_and_si128(loader, _mm_srli_epi64(loader, 1)),
                    *mask_read++);
            } while (data_read < data_read_end);
#else
            do
            {
                loader = *data_ptr++;
                *result_ptr++ = loader & (loader >> 1) & (*mask_ptr++);
            } while (data_ptr < data_read_end);
#endif
        }
    }
    else
    {
#ifdef __LP64__
        do
        {
            loader = *data_read++;
            *writer++ = _mm_andnot_si128(
                _mm_or_si128(loader, _mm_srli_epi64(loader, 1)), *mask_read++);
        } while (data_read < data_read_end);
#else
        do
        {
            loader = *data_ptr++;
            *result_ptr++ = (~(loader | (loader >> 1))) & (*mask_ptr++);
        } while (data_ptr < data_read_end);
#endif
    }
}


static inline uintptr_t popcount2_vecs(const __m128i* vptr, uintptr_t ct)
{
    // assumes ct is a multiple of 6.
    const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
    const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
    const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
    uintptr_t tot = 0;
    const __m128i* vend;
    __m128i loader1;
    __m128i loader2;
    __m128i count1;
    __m128i count2;
    __univec acc;

    while (ct >= 30) {
        ct -= 30;
        vend = &(vptr[30]);
    popcount2_vecs_main_loop:
        acc.vi = _mm_setzero_si128();
        do
        {
            loader1 = *vptr++;
            loader2 = *vptr++;
            count1 =
                _mm_add_epi64(_mm_and_si128(loader1, m2),
                              _mm_and_si128(_mm_srli_epi64(loader1, 2), m2));
            count2 =
                _mm_add_epi64(_mm_and_si128(loader2, m2),
                              _mm_and_si128(_mm_srli_epi64(loader2, 2), m2));

            loader1 = *vptr++;
            loader2 = *vptr++;
            count1 = _mm_add_epi64(
                count1,
                _mm_add_epi64(_mm_and_si128(loader1, m2),
                              _mm_and_si128(_mm_srli_epi64(loader1, 2), m2)));
            count2 = _mm_add_epi64(
                count2,
                _mm_add_epi64(_mm_and_si128(loader2, m2),
                              _mm_and_si128(_mm_srli_epi64(loader2, 2), m2)));

            loader1 = *vptr++;
            loader2 = *vptr++;
            count1 = _mm_add_epi64(
                count1,
                _mm_add_epi64(_mm_and_si128(loader1, m2),
                              _mm_and_si128(_mm_srli_epi64(loader1, 2), m2)));
            count2 = _mm_add_epi64(
                count2,
                _mm_add_epi64(_mm_and_si128(loader2, m2),
                              _mm_and_si128(_mm_srli_epi64(loader2, 2), m2)));

            acc.vi = _mm_add_epi64(
                acc.vi,
                _mm_add_epi64(_mm_and_si128(count1, m4),
                              _mm_and_si128(_mm_srli_epi64(count1, 4), m4)));
            acc.vi = _mm_add_epi64(
                acc.vi,
                _mm_add_epi64(_mm_and_si128(count2, m4),
                              _mm_and_si128(_mm_srli_epi64(count2, 4), m4)));
        } while (vptr < vend);
        acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8),
                               _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
        tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
    }
    if (ct) {
        vend = &(vptr[ct]);
        ct = 0;
        goto popcount2_vecs_main_loop;
    }
    return tot;
}

uintptr_t popcount2_longs(const uintptr_t* lptr, uintptr_t word_ct)
{
    // treats lptr[] as an array of two-bit instead of one-bit numbers
    uintptr_t tot = 0;
    const uintptr_t* lptr_end = &(lptr[word_ct]);
#ifdef __LP64__
    uintptr_t twelve_ct;
    const __m128i* vptr;
    vptr = (const __m128i*) lptr;
    twelve_ct = word_ct / 12;
    tot += popcount2_vecs(vptr, twelve_ct * 6);
    lptr = &(lptr[twelve_ct * 12]);
#else
    const uintptr_t* lptr_six_end;
    uintptr_t loader1;
    uintptr_t loader2;
    uintptr_t ulii;
    uintptr_t uljj;
    lptr_six_end = &(lptr[word_ct - (word_ct % 6)]);
    while (lptr < lptr_six_end) {
        loader1 = *lptr++;
        loader2 = *lptr++;
        ulii = (loader1 & 0x33333333) + ((loader1 >> 2) & 0x33333333);
        uljj = (loader2 & 0x33333333) + ((loader2 >> 2) & 0x33333333);
        loader1 = *lptr++;
        loader2 = *lptr++;
        ulii += (loader1 & 0x33333333) + ((loader1 >> 2) & 0x33333333);
        uljj += (loader2 & 0x33333333) + ((loader2 >> 2) & 0x33333333);
        loader1 = *lptr++;
        loader2 = *lptr++;
        ulii += (loader1 & 0x33333333) + ((loader1 >> 2) & 0x33333333);
        uljj += (loader2 & 0x33333333) + ((loader2 >> 2) & 0x33333333);
        ulii = (ulii & 0x0f0f0f0f) + ((ulii >> 4) & 0x0f0f0f0f);
        ulii += (uljj & 0x0f0f0f0f) + ((uljj >> 4) & 0x0f0f0f0f);

        // Each 8-bit slot stores a number in 0..48.  Multiplying by 0x01010101
        // is equivalent to the left-shifts and adds we need to sum those four
        // 8-bit numbers in the high-order slot.
        tot += (ulii * 0x01010101) >> 24;
    }
#endif
    while (lptr < lptr_end) {
        tot += popcount2_long(*lptr++);
    }
    return tot;
}


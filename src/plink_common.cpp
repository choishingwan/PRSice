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

uint32_t aligned_malloc(uintptr_t size, uintptr_t** aligned_pp)
{
#if defined __LP64__ && !defined __APPLE__
    // Avoid random segfaults on 64-bit machines which have 8-byte- instead of
    // 16-byte-aligned malloc().  (Slightly different code is needed if malloc()
    // does not even guarantee 8-byte alignment.)
    uintptr_t* malloc_ptr = (uintptr_t*) malloc(size + VEC_BYTES);
    if (!malloc_ptr) {
        g_failed_alloc_attempt_size = size + VEC_BYTES;
        return 1;
    }
    *aligned_pp = (uintptr_t*) ((((uintptr_t) malloc_ptr) + VEC_BYTES)
                                & (~(VEC_BYTES_M1 * ONELU)));
    (*aligned_pp)[-1] = (uintptr_t) malloc_ptr;
#else
    // no SSE2 concerns here
    *aligned_pp = (uintptr_t*) malloc(size);
    if (!(*aligned_pp)) {
        g_failed_alloc_attempt_size = size;
        return 1;
    }
#endif
    return 0;
}

void aligned_free(uintptr_t* aligned_pp)
{
#if defined __LP64__ && !defined __APPLE__
    free((uintptr_t*) (aligned_pp[-1]));
#else
    free(aligned_pp);
#endif
}

void logstr(const char* ss)
{
    if (!g_debug_on) {
        fputs(ss, g_logfile);
        if (ferror(g_logfile)) {
            putc_unlocked('\n', stdout);
            fflush(stdout);
            fprintf(stderr,
                    "Warning: Logging failure on:\n%s\nFurther logging will "
                    "not be attempted in this run.\n",
                    ss);
            g_log_failed = 1;
        }
    }
    else
    {
        if (g_log_failed) {
            fflush(stdout);
            fputs(ss, stderr);
        }
        else
        {
            fputs(ss, g_logfile);
            if (ferror(g_logfile)) {
                putc_unlocked('\n', stdout);
                fflush(stdout);
                fprintf(stderr,
                        "Error: Debug logging failure.  Dumping to stderr:\n%s",
                        ss);
                g_log_failed = 1;
            }
            else
            {
                fflush(g_logfile);
            }
        }
    }
}

void logprint(const char* ss)
{
    logstr(ss);
    fputs(ss, stdout);
}

void logerrprint(const char* ss)
{
    logstr(ss);
    fflush(stdout);
    fputs(ss, stderr);
}


void logerrprintb()
{
    logstr(g_logbuf);
    fflush(stdout);
    fputs(g_logbuf, stderr);
}

void wordwrap(uint32_t suffix_len, char* ss)
{
    // Input: A null-terminated string with no intermediate newlines.  If
    //        suffix_len is zero, there should be a terminating \n; otherwise,
    //        the last character should be a space.
    // Effect: Spaces are replaced with newlines in a manner that plays well
    // with
    //         80 column terminal windows.  (Multi-space blocks are never
    //         collapsed.)
    char* token_start = ss;
    char* line_end = &(ss[79]);
    char* token_end;
    while (1) {
        while (*token_start == ' ') {
            token_start++;
        }
        if (token_start > line_end) {
            do
            {
                *line_end = '\n';
                line_end = &(line_end[80]);
            } while (token_start > line_end);
        }
        token_end = strchr(token_start, ' ');
        if (!token_end) {
            if (&(token_start[79]) == line_end) {
                return;
            }
            token_end = strchr(token_start, '\0');
            if (!suffix_len) {
                if (token_end <= &(line_end[1])) {
                    // okay if end-of-string is one past the end, because
                    // function assumes last character is \n in suffix_len == 0
                    // case
                    assert(token_end[-1] == '\n');
                    return;
                }
            }
            else
            {
                if (&(token_end[suffix_len]) <= line_end) {
                    return;
                }
                // because of terminal space assumption, token_start actually
                // points to the end of the string
                assert(token_start[-1] == ' ');
            }
            token_start[-1] = '\n';
            return;
        }
        if (token_end > line_end) {
            if (&(token_start[79]) != line_end) {
                token_start[-1] = '\n';
                line_end = &(token_start[79]);
                if (token_end > line_end) {
                    // single really long token, can't do anything beyond
                    // putting it on its own line
                    *token_end = '\n';
                    line_end = &(token_end[80]);
                }
            }
            else
            {
                // single really long token, *and* previous token was either
                // nonexistent or long
                *token_end = '\n';
                line_end = &(token_end[80]);
            }
        }
        token_start = &(token_end[1]);
    }
}

void wordwrapb(uint32_t suffix_len) { wordwrap(suffix_len, g_logbuf); }


// manually managed, very large stack
unsigned char* g_bigstack_base;
unsigned char* g_bigstack_end;

unsigned char* bigstack_alloc(uintptr_t size)
{
    unsigned char* alloc_ptr;
    size = round_up_pow2(size, CACHELINE);
    if (bigstack_left() < size) {
        g_failed_alloc_attempt_size = size;
        return nullptr;
    }
    alloc_ptr = g_bigstack_base;
    g_bigstack_base += size;
    return alloc_ptr;
}


#ifdef __LP64__
static inline uint32_t scan_uint_capped_finish(const char* ss, uint64_t cap,
                                               uint32_t* valp)
{
    uint64_t val = *valp;
    while (1) {
        // a little bit of unrolling seems to help
        const uint64_t cur_digit = (uint64_t)((unsigned char) (*ss++)) - 48;
        if (cur_digit >= 10) {
            break;
        }
        // val = val * 10 + cur_digit;
        const uint64_t cur_digit2 = (uint64_t)((unsigned char) (*ss++)) - 48;
        if (cur_digit2 >= 10) {
            val = val * 10 + cur_digit;
            if (val > cap) {
                return 1;
            }
            break;
        }
        val = val * 100 + cur_digit * 10 + cur_digit2;
        if (val > cap) {
            return 1;
        }
    }
    *valp = val;
    return 0;
}

uint32_t scan_posint_capped(const char* ss, uint64_t cap, uint32_t* valp)
{
    // '0' has ascii code 48
    *valp = (uint32_t)((unsigned char) (*ss++)) - 48;
    if (*valp >= 10) {
        // permit leading '+' (ascii 43), but not '++' or '+-'
        if (*valp != 0xfffffffbU) {
            return 1;
        }
        *valp = (uint32_t)((unsigned char) (*ss++)) - 48;
        if (*valp >= 10) {
            return 1;
        }
    }
    while (!(*valp)) {
        *valp = (uint32_t)((unsigned char) (*ss++)) - 48;
        if ((*valp) >= 10) {
            return 1;
        }
    }
    return scan_uint_capped_finish(ss, cap, valp);
}

uint32_t scan_uint_capped(const char* ss, uint64_t cap, uint32_t* valp)
{
    // Reads an integer in [0, cap].  Assumes first character is nonspace.
    uint32_t val = (uint32_t)((unsigned char) (*ss++)) - 48;
    if (val >= 10) {
        if (val != 0xfffffffbU) {
            // '-' has ascii code 45, so unsigned 45 - 48 = 0xfffffffdU
            if ((val != 0xfffffffdU) || (*ss != '0')) {
                return 1;
            }
            // accept "-0", "-00", etc.
            while (*(++ss) == '0')
                ;
            *valp = 0;
            return ((uint32_t)((unsigned char) (*ss)) - 48) < 10;
        }
        // accept leading '+'
        val = (uint32_t)((unsigned char) (*ss++)) - 48;
        if (val >= 10) {
            return 1;
        }
    }
    *valp = val;
    return scan_uint_capped_finish(ss, cap, valp);
}

uint32_t scan_int_abs_bounded(const char* ss, uint64_t bound, int32_t* valp)
{
    // Reads an integer in [-bound, bound].  Assumes first character is
    // nonspace.
    *valp = (uint32_t)((unsigned char) (*ss++)) - 48;
    int32_t sign = 1;
    if (((uint32_t) *valp) >= 10) {
        if (*valp == -3) {
            sign = -1;
        }
        else if (*valp != -5)
        {
            return 1;
        }
        *valp = (uint32_t)((unsigned char) (*ss++)) - 48;
        if (((uint32_t) *valp) >= 10) {
            return 1;
        }
    }
    if (scan_uint_capped_finish(ss, bound, (uint32_t*) valp)) {
        return 1;
    }
    *valp *= sign;
    return 0;
}
#else // not __LP64__
uint32_t scan_posint_capped32(const char* ss, uint32_t cap_div_10,
                              uint32_t cap_mod_10, uint32_t* valp)
{
    // '0' has ascii code 48
    uint32_t val = (uint32_t)((unsigned char) (*ss++)) - 48;
    if (val >= 10) {
        if (val != 0xfffffffbU) {
            return 1;
        }
        val = (uint32_t)((unsigned char) (*ss++)) - 48;
        if (val >= 10) {
            return 1;
        }
    }
    while (!val) {
        val = (uint32_t)((unsigned char) (*ss++)) - 48;
        if (val >= 10) {
            return 1;
        }
    }
    while (1) {
        const uint32_t cur_digit = (uint32_t)((unsigned char) (*ss++)) - 48;
        if (cur_digit >= 10) {
            *valp = val;
            return 0;
        }
        // avoid integer overflow in middle of computation
        if ((val >= cap_div_10)
            && ((val > cap_div_10) || (cur_digit > cap_mod_10)))
        {
            return 1;
        }
        val = val * 10 + cur_digit;
    }
}

uint32_t scan_uint_capped32(const char* ss, uint32_t cap_div_10,
                            uint32_t cap_mod_10, uint32_t* valp)
{
    // Reads an integer in [0, cap].  Assumes first character is nonspace.
    uint32_t val = (uint32_t)((unsigned char) (*ss++)) - 48;
    if (val >= 10) {
        if (val != 0xfffffffbU) {
            if ((val != 0xfffffffdU) || (*ss != '0')) {
                return 1;
            }
            while (*(++ss) == '0')
                ;
            *valp = 0;
            return ((uint32_t)((unsigned char) (*ss)) - 48) < 10;
        }
        val = (uint32_t)((unsigned char) (*ss++)) - 48;
        if (val >= 10) {
            return 1;
        }
    }
    while (1) {
        const uint32_t cur_digit = (uint32_t)((unsigned char) (*ss++)) - 48;
        if (cur_digit >= 10) {
            *valp = val;
            return 0;
        }
        if ((val >= cap_div_10)
            && ((val > cap_div_10) || (cur_digit > cap_mod_10)))
        {
            return 1;
        }
        val = val * 10 + cur_digit;
    }
}

uint32_t scan_int_abs_bounded32(const char* ss, uint32_t bound_div_10,
                                uint32_t bound_mod_10, int32_t* valp)
{
    // Reads an integer in [-bound, bound].  Assumes first character is
    // nonspace.
    uint32_t val = (uint32_t)((unsigned char) (*ss++)) - 48;
    int32_t sign = 1;
    if (val >= 10) {
        if (val == 0xfffffffdU) {
            sign = -1;
        }
        else if (val != 0xfffffffbU)
        {
            return 1;
        }
        val = (uint32_t)((unsigned char) (*ss++)) - 48;
        if (val >= 10) {
            return 1;
        }
    }
    while (1) {
        const uint32_t cur_digit = (uint32_t)((unsigned char) (*ss++)) - 48;
        if (cur_digit >= 10) {
            *valp = sign * ((int32_t) val);
            return 0;
        }
        if ((val >= bound_div_10)
            && ((val > bound_div_10) || (cur_digit > bound_mod_10)))
        {
            return 1;
        }
        val = val * 10 + cur_digit;
    }
}
#endif


char* next_token_mult(char* sptr, uint32_t ct)
{
    assert(ct);
    if (!sptr) {
        return nullptr;
    }
    unsigned char ucc = *sptr;
    do
    {
        while (ucc > 32) {
            ucc = *(++sptr);
        }
        while ((ucc == ' ') || (ucc == '\t')) {
            ucc = *(++sptr);
        }
        if (ucc <= 32) {
            return nullptr;
        }
    } while (--ct);
    return sptr;
}

// number-to-string encoders

static const char digit2_table[200] = {
    '0', '0', '0', '1', '0', '2', '0', '3', '0', '4', '0', '5', '0', '6', '0',
    '7', '0', '8', '0', '9', '1', '0', '1', '1', '1', '2', '1', '3', '1', '4',
    '1', '5', '1', '6', '1', '7', '1', '8', '1', '9', '2', '0', '2', '1', '2',
    '2', '2', '3', '2', '4', '2', '5', '2', '6', '2', '7', '2', '8', '2', '9',
    '3', '0', '3', '1', '3', '2', '3', '3', '3', '4', '3', '5', '3', '6', '3',
    '7', '3', '8', '3', '9', '4', '0', '4', '1', '4', '2', '4', '3', '4', '4',
    '4', '5', '4', '6', '4', '7', '4', '8', '4', '9', '5', '0', '5', '1', '5',
    '2', '5', '3', '5', '4', '5', '5', '5', '6', '5', '7', '5', '8', '5', '9',
    '6', '0', '6', '1', '6', '2', '6', '3', '6', '4', '6', '5', '6', '6', '6',
    '7', '6', '8', '6', '9', '7', '0', '7', '1', '7', '2', '7', '3', '7', '4',
    '7', '5', '7', '6', '7', '7', '7', '8', '7', '9', '8', '0', '8', '1', '8',
    '2', '8', '3', '8', '4', '8', '5', '8', '6', '8', '7', '8', '8', '8', '9',
    '9', '0', '9', '1', '9', '2', '9', '3', '9', '4', '9', '5', '9', '6', '9',
    '7', '9', '8', '9', '9'};

char* uint32toa(uint32_t uii, char* start)
{
    // Memory-efficient fast integer writer.  (You can do a bit better sometimes
    // by using a larger lookup table, but on average I doubt that pays off.)
    // Returns a pointer to the end of the integer (not null-terminated).
    uint32_t quotient;
    if (uii < 1000) {
        if (uii < 10) {
            *start++ = '0' + uii;
            return start;
        }
        if (uii < 100) {
            goto uint32toa_2;
        }
        quotient = uii / 100;
        *start++ = '0' + quotient;
    }
    else
    {
        if (uii < 10000000) {
            if (uii >= 100000) {
                if (uii < 1000000) {
                    goto uint32toa_6;
                }
                quotient = uii / 1000000;
                *start++ = '0' + quotient;
                goto uint32toa_6b;
            }
            if (uii < 10000) {
                goto uint32toa_4;
            }
            quotient = uii / 10000;
            *start++ = '0' + quotient;
        }
        else
        {
            if (uii >= 100000000) {
                quotient = uii / 100000000;
                if (uii >= 1000000000) {
                    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
                }
                else
                {
                    *start++ = '0' + quotient;
                }
                uii -= 100000000 * quotient;
            }
            quotient = uii / 1000000;
            start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        uint32toa_6b:
            uii -= 1000000 * quotient;
        uint32toa_6:
            quotient = uii / 10000;
            start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        }
        uii -= 10000 * quotient;
    uint32toa_4:
        // could make a uitoa_z4() call here, but that's slightly slower
        quotient = uii / 100;
        start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    }
    uii -= 100 * quotient;
uint32toa_2:
    return memcpya(start, &(digit2_table[uii * 2]), 2);
}

char* int32toa(int32_t ii, char* start)
{
    uint32_t uii = ii;
    if (ii < 0) {
        // -INT_MIN is undefined, but negating the unsigned int equivalent works
        *start++ = '-';
        uii = -uii;
    }
    return uint32toa(uii, start);
}

char* uitoa_z4(uint32_t uii, char* start)
{
    uint32_t quotient = uii / 100;
    assert(quotient < 100);
    uii -= 100 * quotient;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    return memcpya(start, &(digit2_table[uii * 2]), 2);
}

char* uitoa_z6(uint32_t uii, char* start)
{
    uint32_t quotient = uii / 10000;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    return uitoa_z4(uii - 10000 * quotient, start);
}

char* uitoa_z8(uint32_t uii, char* start)
{
    uint32_t quotient = uii / 1000000;
    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    return uitoa_z6(uii - 1000000 * quotient, start);
}

char* uint32toa_w4(uint32_t uii, char* start)
{
    uint32_t quotient;
    if (uii < 1000) {
        if (uii < 10) {
            // assumes little-endian
            *((uint32_t*) start) = 0x30202020 + (uii << 24);
            return &(start[4]);
        }
        if (uii < 100) {
            memset(start, 32, 2);
        }
        else
        {
            quotient = uii / 100;
            *start++ = ' ';
            *start++ = '0' + quotient;
            uii -= quotient * 100;
        }
        return memcpya(start, &(digit2_table[uii * 2]), 2);
    }
    else
    {
        // presumably the field width is 4 for a reason; don't bother optimizing
        // this
        return uint32toa(uii, start);
    }
}

char* uint32toa_w6(uint32_t uii, char* start)
{
    uint32_t quotient;
    if (uii < 1000) {
        if (uii < 10) {
            start = memseta(start, 32, 5);
            *start++ = '0' + uii;
            return start;
        }
        if (uii < 100) {
            start = memseta(start, 32, 4);
            goto uint32toa_w6_2;
        }
        quotient = uii / 100;
        // the little-endian trick doesn't seem to help here.  possibly relevant
        // differences from uint32toa_w4() and _w8(): sequential dependence on
        // quotient, need to interpret pointer as a char* again
        start = memseta(start, 32, 3);
        *start++ = '0' + quotient;
    }
    else
    {
        if (uii < 10000000) {
            if (uii >= 100000) {
                if (uii < 1000000) {
                    goto uint32toa_w6_6;
                }
                quotient = uii / 1000000;
                *start++ = '0' + quotient;
                goto uint32toa_w6_6b;
            }
            else if (uii >= 10000)
            {
                *start++ = ' ';
                quotient = uii / 10000;
                *start++ = '0' + quotient;
            }
            else
            {
                start = memseta(start, 32, 2);
                goto uint32toa_w6_4;
            }
        }
        else
        {
            if (uii >= 100000000) {
                quotient = uii / 100000000;
                if (uii >= 1000000000) {
                    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
                }
                else
                {
                    *start++ = '0' + quotient;
                }
                uii -= 100000000 * quotient;
            }
            quotient = uii / 1000000;
            start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        uint32toa_w6_6b:
            uii -= 1000000 * quotient;
        uint32toa_w6_6:
            quotient = uii / 10000;
            start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        }
        uii -= 10000 * quotient;
    uint32toa_w6_4:
        quotient = uii / 100;
        start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    }
    uii -= 100 * quotient;
uint32toa_w6_2:
    return memcpya(start, &(digit2_table[uii * 2]), 2);
}

char* uint32toa_w7(uint32_t uii, char* start)
{
    uint32_t quotient;
    if (uii < 1000) {
        if (uii < 10) {
            start = memseta(start, 32, 6);
            *start++ = '0' + uii;
            return start;
        }
        if (uii < 100) {
            start = memseta(start, 32, 5);
            goto uint32toa_w7_2;
        }
        quotient = uii / 100;
        start = memseta(start, 32, 4);
        *start++ = '0' + quotient;
    }
    else
    {
        if (uii < 10000000) {
            if (uii >= 100000) {
                if (uii >= 1000000) {
                    quotient = uii / 1000000;
                    *start++ = '0' + quotient;
                    goto uint32toa_w7_6b;
                }
                *start++ = ' ';
                goto uint32toa_w7_6;
            }
            else if (uii >= 10000)
            {
                start = memseta(start, 32, 2);
                quotient = uii / 10000;
                *start++ = '0' + quotient;
            }
            else
            {
                start = memseta(start, 32, 3);
                goto uint32toa_w7_4;
            }
        }
        else
        {
            if (uii >= 100000000) {
                quotient = uii / 100000000;
                if (uii >= 1000000000) {
                    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
                }
                else
                {
                    *start++ = '0' + quotient;
                }
                uii -= 100000000 * quotient;
            }
            quotient = uii / 1000000;
            start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        uint32toa_w7_6b:
            uii -= 1000000 * quotient;
        uint32toa_w7_6:
            quotient = uii / 10000;
            start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        }
        uii -= 10000 * quotient;
    uint32toa_w7_4:
        quotient = uii / 100;
        start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    }
    uii -= 100 * quotient;
uint32toa_w7_2:
    return memcpya(start, &(digit2_table[uii * 2]), 2);
}

char* uint32toa_w8(uint32_t uii, char* start)
{
    uint32_t quotient;
    if (uii < 1000) {
        if (uii < 10) {
#ifdef __LP64__
            *((uintptr_t*) start) =
                0x3020202020202020LLU + (((uintptr_t) uii) << 56);
            return &(start[8]);
#else
            start = memseta(start, 32, 7);
            *start++ = '0' + uii;
            return start;
#endif
        }
        if (uii < 100) {
            start = memseta(start, 32, 6);
            goto uint32toa_w8_2;
        }
        quotient = uii / 100;
        start = memseta(start, 32, 5);
        *start++ = '0' + quotient;
    }
    else
    {
        if (uii < 10000000) {
            if (uii >= 100000) {
                if (uii < 1000000) {
                    start = memseta(start, 32, 2);
                    goto uint32toa_w8_6;
                }
                quotient = uii / 1000000;
                *start = ' ';
                start[1] = '0' + quotient;
                start += 2;
                goto uint32toa_w8_6b;
            }
            else if (uii < 10000)
            {
                start = memseta(start, 32, 4);
                goto uint32toa_w8_4;
            }
            memset(start, 32, 3);
            quotient = uii / 10000;
            start[3] = '0' + quotient;
            start += 4;
        }
        else
        {
            if (uii >= 100000000) {
                quotient = uii / 100000000;
                if (uii >= 1000000000) {
                    start = memcpya(start, &(digit2_table[quotient * 2]), 2);
                }
                else
                {
                    *start++ = '0' + quotient;
                }
                uii -= 100000000 * quotient;
            }
            quotient = uii / 1000000;
            start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        uint32toa_w8_6b:
            uii -= 1000000 * quotient;
        uint32toa_w8_6:
            quotient = uii / 10000;
            start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        }
        uii -= 10000 * quotient;
    uint32toa_w8_4:
        quotient = uii / 100;
        start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    }
    uii -= 100 * quotient;
uint32toa_w8_2:
    return memcpya(start, &(digit2_table[uii * 2]), 2);
}

char* uint32toa_w10(uint32_t uii, char* start)
{
    // if we decide to reduce code size and optimize only one field width, this
    // should be it
    uint32_t quotient;
    if (uii < 1000) {
        if (uii < 10) {
            start = memseta(start, 32, 9);
            *start++ = '0' + uii;
            return start;
        }
        if (uii < 100) {
            start = memseta(start, 32, 8);
            goto uint32toa_w10_2;
        }
        quotient = uii / 100;
        start = memseta(start, 32, 7);
        *start++ = '0' + quotient;
    }
    else
    {
        if (uii < 10000000) {
            if (uii >= 100000) {
                if (uii < 1000000) {
                    start = memseta(start, 32, 4);
                    goto uint32toa_w10_6;
                }
                quotient = uii / 1000000;
                memset(start, 32, 3);
                start[3] = '0' + quotient;
                start += 4;
                goto uint32toa_w10_6b;
            }
            else if (uii < 10000)
            {
                start = memseta(start, 32, 6);
                goto uint32toa_w10_4;
            }
            memset(start, 32, 5);
            quotient = uii / 10000;
            start[5] = '0' + quotient;
            start += 6;
        }
        else
        {
            if (uii >= 100000000) {
                quotient = uii / 100000000;
                if (uii >= 1000000000) {
                    memcpy(start, &(digit2_table[quotient * 2]), 2);
                }
                else
                {
                    *start = ' ';
                    start[1] = '0' + quotient;
                }
                uii -= 100000000 * quotient;
            }
            else
            {
                memset(start, 32, 2);
            }
            quotient = uii / 1000000;
            memcpy(&(start[2]), &(digit2_table[quotient * 2]), 2);
            start += 4;
        uint32toa_w10_6b:
            uii -= 1000000 * quotient;
        uint32toa_w10_6:
            quotient = uii / 10000;
            start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        }
        uii -= 10000 * quotient;
    uint32toa_w10_4:
        quotient = uii / 100;
        start = memcpya(start, &(digit2_table[quotient * 2]), 2);
    }
    uii -= 100 * quotient;
uint32toa_w10_2:
    return memcpya(start, &(digit2_table[uii * 2]), 2);
}


static inline char* uitoa_trunc3(uint32_t uii, char* start)
{
    *start++ = '0' + (uii / 100);
    uii %= 100;
    if (!uii) {
        return start;
    }
    memcpy(start, &(digit2_table[uii * 2]), 2);
    if (start[1] != '0') {
        return &(start[2]);
    }
    return &(start[1]);
}

static inline char* uitoa_trunc4(uint32_t uii, char* start)
{
    uint32_t quotient = uii / 100;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    uii -= 100 * quotient;
    if (uii) {
        start += 2;
        memcpy(start, &(digit2_table[uii * 2]), 2);
    }
    if (start[1] != '0') {
        return &(start[2]);
    }
    return &(start[1]);
}

static inline char* uitoa_trunc6(uint32_t uii, char* start)
{
    uint32_t quotient = uii / 10000;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    uii -= 10000 * quotient;
    if (uii) {
        quotient = uii / 100;
        start += 2;
        memcpy(start, &(digit2_table[quotient * 2]), 2);
        uii -= 100 * quotient;
        if (uii) {
            start += 2;
            memcpy(start, &(digit2_table[uii * 2]), 2);
        }
    }
    if (start[1] != '0') {
        return &(start[2]);
    }
    return &(start[1]);
}

static inline char* uitoa_trunc8(uint32_t uii, char* start)
{
    uint32_t quotient = uii / 1000000;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    uii -= 1000000 * quotient;
    if (uii) {
        quotient = uii / 10000;
        start += 2;
        memcpy(start, &(digit2_table[quotient * 2]), 2);
        uii -= 10000 * quotient;
        if (uii) {
            quotient = uii / 100;
            start += 2;
            memcpy(start, &(digit2_table[quotient * 2]), 2);
            uii -= 100 * quotient;
            if (uii) {
                start += 2;
                memcpy(start, &(digit2_table[uii * 2]), 2);
            }
        }
    }
    if (start[1] != '0') {
        return &(start[2]);
    }
    return &(start[1]);
}


static inline char* qrtoa_1p2(uint32_t quotient, uint32_t remainder,
                              char* start)
{
    *start++ = '0' + quotient;
    if (!remainder) {
        return start;
    }
    *start++ = '.';
    memcpy(start, &(digit2_table[remainder * 2]), 2);
    if (start[1] != '0') {
        return &(start[2]);
    }
    return &(start[1]);
}

static inline char* qrtoa_1p3(uint32_t quotient, uint32_t remainder,
                              char* start)
{
    // quotient = (int32_t)dxx;
    // remainder = ((int32_t)(dxx * 1000)) - (quotient * 1000);
    *start++ = '0' + quotient;
    if (!remainder) {
        return start;
    }
    *start++ = '.';
    quotient = remainder / 10;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    remainder -= 10 * quotient;
    if (remainder) {
        start[2] = '0' + remainder;
        return &(start[3]);
    }
    if (start[1] != '0') {
        return &(start[2]);
    }
    return &(start[1]);
}

static inline char* qrtoa_1p5(uint32_t quotient, uint32_t remainder,
                              char* start)
{
    *start++ = '0' + quotient;
    if (!remainder) {
        return start;
    }
    *start++ = '.';
    quotient = remainder / 1000;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    remainder -= 1000 * quotient;
    if (remainder) {
        quotient = remainder / 10;
        start += 2;
        memcpy(start, &(digit2_table[quotient * 2]), 2);
        remainder -= 10 * quotient;
        if (remainder) {
            start[2] = '0' + remainder;
            return &(start[3]);
        }
    }
    if (start[1] != '0') {
        return &(start[2]);
    }
    return &(start[1]);
}

static inline char* qrtoa_1p7(uint32_t quotient, uint32_t remainder,
                              char* start)
{
    *start++ = '0' + quotient;
    if (!remainder) {
        return start;
    }
    *start++ = '.';
    quotient = remainder / 100000;
    memcpy(start, &(digit2_table[quotient * 2]), 2);
    remainder -= 100000 * quotient;
    if (remainder) {
        quotient = remainder / 1000;
        start += 2;
        memcpy(start, &(digit2_table[quotient * 2]), 2);
        remainder -= 1000 * quotient;
        if (remainder) {
            quotient = remainder / 10;
            start += 2;
            memcpy(start, &(digit2_table[quotient * 2]), 2);
            remainder -= 10 * quotient;
            if (remainder) {
                start[2] = '0' + remainder;
                return &(start[3]);
            }
        }
    }
    if (start[1] != '0') {
        return &(start[2]);
    }
    return &(start[1]);
}

// Okay, time to do banker's rounding when printing doubles.  14 digits of
// precision are used in judging equality to 0.5 (actual precision of doubles
// is 15-17 digits); the intention is to capture all directly loaded or exactly
// computed edge cases (so enough tolerance is needed to survive the internal
// multiplications by powers of 10, etc.), while rounding a negligible number
// of honest-to-god 0.4999999s up and 0.5000001s down.
// To avoid inadvertent printing of an extra digit, there's a deliberate gap
// between the 99.9994999...-type bounds and the largest numbers that would
// actually round down.
static const double banker_round5[] = {0.499995, 0.500005};
static const double banker_round6[] = {0.4999995, 0.5000005};
static const double banker_round7[] = {0.49999995, 0.50000005};
static const double banker_round8[] = {0.499999995, 0.500000005};
static const double banker_round9[] = {0.4999999995, 0.5000000005};
static const double banker_round10[] = {0.49999999995, 0.50000000005};
static const double banker_round11[] = {0.499999999995, 0.500000000005};
static const double banker_round12[] = {0.4999999999995, 0.5000000000005};

static inline uint32_t double_bround(double dxx, const double* banker_round)
{
    uint32_t result = (int32_t) dxx;
    return result
           + (int32_t)((dxx - ((int32_t) result)) + banker_round[result & 1]);
}

// These are separate functions so the compiler can optimize the integer
// divisions.
static inline void double_bround1(double dxx, const double* banker_round,
                                  uint32_t* quotientp, uint32_t* remainderp)
{
    dxx *= 10;
    uint32_t remainder = (int32_t) dxx;
    remainder +=
        (int32_t)((dxx - ((int32_t) remainder)) + banker_round[remainder & 1]);
    *quotientp = remainder / 10;
    *remainderp = remainder - (*quotientp) * 10;
}

static inline void double_bround2(double dxx, const double* banker_round,
                                  uint32_t* quotientp, uint32_t* remainderp)
{
    dxx *= 100;
    uint32_t remainder = (int32_t) dxx;
    remainder +=
        (int32_t)((dxx - ((int32_t) remainder)) + banker_round[remainder & 1]);
    *quotientp = remainder / 100;
    *remainderp = remainder - (*quotientp) * 100;
}

static inline void double_bround3(double dxx, const double* banker_round,
                                  uint32_t* quotientp, uint32_t* remainderp)
{
    dxx *= 1000;
    uint32_t remainder = (int32_t) dxx;
    remainder +=
        (int32_t)((dxx - ((int32_t) remainder)) + banker_round[remainder & 1]);
    *quotientp = remainder / 1000;
    *remainderp = remainder - (*quotientp) * 1000;
}

static inline void double_bround4(double dxx, const double* banker_round,
                                  uint32_t* quotientp, uint32_t* remainderp)
{
    dxx *= 10000;
    uint32_t remainder = (int32_t) dxx;
    remainder +=
        (int32_t)((dxx - ((int32_t) remainder)) + banker_round[remainder & 1]);
    *quotientp = remainder / 10000;
    *remainderp = remainder - (*quotientp) * 10000;
}

static inline void double_bround5(double dxx, const double* banker_round,
                                  uint32_t* quotientp, uint32_t* remainderp)
{
    dxx *= 100000;
    uint32_t remainder = (int32_t) dxx;
    remainder +=
        (int32_t)((dxx - ((int32_t) remainder)) + banker_round[remainder & 1]);
    *quotientp = remainder / 100000;
    *remainderp = remainder - (*quotientp) * 100000;
}

static inline void double_bround6(double dxx, const double* banker_round,
                                  uint32_t* quotientp, uint32_t* remainderp)
{
    dxx *= 1000000;
    uint32_t remainder = (int32_t) dxx;
    remainder +=
        (int32_t)((dxx - ((int32_t) remainder)) + banker_round[remainder & 1]);
    *quotientp = remainder / 1000000;
    *remainderp = remainder - (*quotientp) * 1000000;
}

static inline void double_bround7(double dxx, const double* banker_round,
                                  uint32_t* quotientp, uint32_t* remainderp)
{
    dxx *= 10000000;
    uint32_t remainder = (int32_t) dxx;
    remainder +=
        (int32_t)((dxx - ((int32_t) remainder)) + banker_round[remainder & 1]);
    *quotientp = remainder / 10000000;
    *remainderp = remainder - (*quotientp) * 10000000;
}

char* dtoa_so6(double dxx, char* start)
{
    // 6 sig fig number, 0.999995 <= dxx < 999999.5
    // 'so' = "significand only"
    // Just hardcoding all six cases, in the absence of a better approach...
    uint32_t uii;
    uint32_t quotient;
    uint32_t remainder;
    if (dxx < 99.999949999999) {
        if (dxx < 9.9999949999999) {
            double_bround5(dxx, banker_round8, &quotient, &remainder);
            return qrtoa_1p5(quotient, remainder, start);
        }
        double_bround4(dxx, banker_round8, &quotient, &remainder);
        start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        if (!remainder) {
            return start;
        }
        *start++ = '.';
        quotient = remainder / 100;
        memcpy(start, &(digit2_table[quotient * 2]), 2);
        remainder -= 100 * quotient;
        if (remainder) {
            start += 2;
        dtoa_so6_pretail:
            memcpy(start, &(digit2_table[remainder * 2]), 2);
        }
    dtoa_so6_tail:
        if (start[1] != '0') {
            return &(start[2]);
        }
        return &(start[1]);
    }
    else if (dxx < 9999.9949999999)
    {
        if (dxx < 999.99949999999) {
            double_bround3(dxx, banker_round8, &uii, &remainder);
            quotient = uii / 100;
            *start++ = '0' + quotient;
            quotient = uii - 100 * quotient;
            start = memcpya(start, &(digit2_table[quotient * 2]), 2);
            if (!remainder) {
                return start;
            }
            *start++ = '.';
            quotient = remainder / 10;
            memcpy(start, &(digit2_table[quotient * 2]), 2);
            remainder -= quotient * 10;
            if (!remainder) {
                goto dtoa_so6_tail;
            }
            start[2] = '0' + remainder;
            return &(start[3]);
        }
        double_bround2(dxx, banker_round8, &uii, &remainder);
        quotient = uii / 100;
        start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        quotient = uii - (100 * quotient);
        start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        if (!remainder) {
            return start;
        }
        *start++ = '.';
        goto dtoa_so6_pretail;
    }
    else if (dxx < 99999.949999999)
    {
        double_bround1(dxx, banker_round8, &uii, &remainder);
        quotient = uii / 10000;
        *start = '0' + quotient;
        uii -= 10000 * quotient;
        quotient = uii / 100;
        start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
        uii = uii - 100 * quotient;
        start = memcpya(start, &(digit2_table[uii * 2]), 2);
        if (!remainder) {
            return start;
        }
        *start++ = '.';
        *start = '0' + remainder;
        return &(start[1]);
    }
    else
    {
        return uitoa_z6(double_bround(dxx, banker_round8), start);
    }
}

// Briefly had banker's rounding for floats, but then I realized that the only
// float-printing function calls are --make-grm related, they all request 6-7
// digits of precision, and at that point it's impossible to distinguish exact
// 0.5-matches in the remainder.  So we just have generic rounding functions
// here, with similar interfaces to the double-rounding functions to minimize
// the need for separate reasoning about this code.
static inline uint32_t float_round(float fxx)
{
    return (uint32_t)((int32_t)(fxx + 0.5));
}

static inline void float_round6(float fxx, uint32_t* quotientp,
                                uint32_t* remainderp)
{
    uint32_t remainder = float_round(fxx * 1000000);
    *quotientp = remainder / 1000000;
    *remainderp = remainder - (*quotientp) * 1000000;
}

char* dtoa_so3(double dxx, char* start)
{
    // 3 sig fig number, 0.995 <= dxx < 999.5
    uint32_t quotient;
    uint32_t remainder;
    if (dxx < 99.949999999999) {
        if (dxx < 9.9949999999999) {
            double_bround2(dxx, banker_round11, &quotient, &remainder);
            return qrtoa_1p2(quotient, remainder, start);
        }
        double_bround1(dxx, banker_round11, &quotient, &remainder);
        start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        if (!remainder) {
            return start;
        }
        *start++ = '.';
    }
    else
    {
        quotient = double_bround(dxx, banker_round11);
        start = memcpya(start, &(digit2_table[(quotient / 10) * 2]), 2);
        remainder = quotient % 10;
    }
    *start = '0' + remainder;
    return &(start[1]);
}

char* dtoa_so4(double dxx, char* start)
{
    // 4 sig fig number, 0.9995 <= dxx < 9999.5
    uint32_t uii;
    uint32_t quotient;
    uint32_t remainder;
    if (dxx < 99.994999999999) {
        if (dxx < 9.9994999999999) {
            double_bround3(dxx, banker_round10, &quotient, &remainder);
            return qrtoa_1p3(quotient, remainder, start);
        }
        double_bround2(dxx, banker_round10, &quotient, &remainder);
        start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        if (!remainder) {
            return start;
        }
        *start++ = '.';
        memcpy(start, &(digit2_table[remainder * 2]), 2);
        if (start[1] != '0') {
            return &(start[2]);
        }
        return &(start[1]);
    }
    else if (dxx < 999.94999999999)
    {
        double_bround1(dxx, banker_round10, &uii, &remainder);
        quotient = uii / 100;
        *start = '0' + quotient;
        quotient = uii - 100 * quotient;
        start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
        if (!remainder) {
            return start;
        }
        *start = '.';
        start[1] = '0' + remainder;
        return &(start[2]);
    }
    else
    {
        uitoa_z4(double_bround(dxx, banker_round10), start);
        return &(start[4]);
    }
}

char* dtoa_so8(double dxx, char* start)
{
    // 8 sig fig number, 0.99999995 <= dxx < 99999999.5
    uint32_t uii;
    uint32_t quotient;
    uint32_t remainder;
    if (dxx < 99.999999499999) {
        if (dxx < 9.9999999499999) {
            double_bround7(dxx, banker_round6, &quotient, &remainder);
            return qrtoa_1p7(quotient, remainder, start);
        }
        double_bround6(dxx, banker_round6, &quotient, &remainder);
        start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        if (!remainder) {
            return start;
        }
        *start++ = '.';
        quotient = remainder / 10000;
        memcpy(start, &(digit2_table[quotient * 2]), 2);
        remainder -= 10000 * quotient;
        if (remainder) {
            start += 2;
        dtoa_so8_pretail4:
            quotient = remainder / 100;
            memcpy(start, &(digit2_table[quotient * 2]), 2);
            remainder -= 100 * quotient;
            if (remainder) {
                start += 2;
            dtoa_so8_pretail2:
                memcpy(start, &(digit2_table[remainder * 2]), 2);
            }
        }
    dtoa_so8_tail:
        if (start[1] != '0') {
            return &(start[2]);
        }
        return &(start[1]);
    }
    else if (dxx < 9999.9999499999)
    {
        if (dxx < 999.99999499999) {
            double_bround5(dxx, banker_round6, &uii, &remainder);
            quotient = uii / 100;
            *start++ = '0' + quotient;
            quotient = uii - 100 * quotient;
            start = memcpya(start, &(digit2_table[quotient * 2]), 2);
            if (!remainder) {
                return start;
            }
            *start++ = '.';
            quotient = remainder / 1000;
            memcpy(start, &(digit2_table[quotient * 2]), 2);
            remainder -= quotient * 1000;
            if (!remainder) {
                goto dtoa_so8_tail;
            }
            start += 2;
        dtoa_so8_pretail3:
            quotient = remainder / 10;
            memcpy(start, &(digit2_table[quotient * 2]), 2);
            remainder -= quotient * 10;
            if (!remainder) {
                goto dtoa_so8_tail;
            }
            start[2] = '0' + remainder;
            return &(start[3]);
        }
        double_bround4(dxx, banker_round6, &uii, &remainder);
        quotient = uii / 100;
        start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        quotient = uii - (100 * quotient);
        start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        if (!remainder) {
            return start;
        }
        *start++ = '.';
        goto dtoa_so8_pretail4;
    }
    else if (dxx < 999999.99499999)
    {
        if (dxx < 99999.999499999) {
            double_bround3(dxx, banker_round6, &uii, &remainder);
            quotient = uii / 10000;
            *start = '0' + quotient;
            uii -= 10000 * quotient;
            quotient = uii / 100;
            start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
            uii -= 100 * quotient;
            start = memcpya(start, &(digit2_table[uii * 2]), 2);
            if (!remainder) {
                return start;
            }
            *start++ = '.';
            goto dtoa_so8_pretail3;
        }
        double_bround2(dxx, banker_round6, &uii, &remainder);
        quotient = uii / 10000;
        start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        uii -= 10000 * quotient;
        quotient = uii / 100;
        start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        uii -= 100 * quotient;
        start = memcpya(start, &(digit2_table[uii * 2]), 2);
        if (!remainder) {
            return start;
        }
        *start++ = '.';
        goto dtoa_so8_pretail2;
    }
    else if (dxx < 9999999.9499999)
    {
        double_bround1(dxx, banker_round6, &uii, &remainder);
        quotient = uii / 1000000;
        *start = '0' + quotient;
        uii -= 1000000 * quotient;
        quotient = uii / 10000;
        start = memcpya(&(start[1]), &(digit2_table[quotient * 2]), 2);
        uii -= 10000 * quotient;
        quotient = uii / 100;
        start = memcpya(start, &(digit2_table[quotient * 2]), 2);
        uii -= 100 * quotient;
        start = memcpya(start, &(digit2_table[uii * 2]), 2);
        if (!remainder) {
            return start;
        }
        *start = '.';
        start[1] = '0' + remainder;
        return &(start[2]);
    }
    else
    {
        return uitoa_z8(double_bround(dxx, banker_round6), start);
    }
}

char* dtoa_e(double dxx, char* start)
{
    uint32_t xp10 = 0;
    uint32_t quotient;
    uint32_t remainder;
    char sign;
    if (dxx != dxx) {
        // do this first to avoid generating exception
        return memcpyl3a(start, "nan");
    }
    else if (dxx < 0)
    {
        *start++ = '-';
        dxx = -dxx;
    }
    if (dxx >= 9.9999994999999e-1) {
        if (dxx >= 9.9999994999999e7) {
            if (dxx >= 9.9999994999999e127) {
                if (dxx == INFINITY) {
                    return memcpyl3a(start, "inf");
                }
                else if (dxx >= 9.9999994999999e255)
                {
                    dxx *= 1.0e-256;
                    xp10 |= 256;
                }
                else
                {
                    dxx *= 1.0e-128;
                    xp10 |= 128;
                }
            }
            if (dxx >= 9.9999994999999e63) {
                dxx *= 1.0e-64;
                xp10 |= 64;
            }
            if (dxx >= 9.9999994999999e31) {
                dxx *= 1.0e-32;
                xp10 |= 32;
            }
            if (dxx >= 9.9999994999999e15) {
                dxx *= 1.0e-16;
                xp10 |= 16;
            }
            if (dxx >= 9.9999994999999e7) {
                dxx *= 1.0e-8;
                xp10 |= 8;
            }
        }
        if (dxx >= 9.9999994999999e3) {
            dxx *= 1.0e-4;
            xp10 |= 4;
        }
        if (dxx >= 9.9999994999999e1) {
            dxx *= 1.0e-2;
            xp10 |= 2;
        }
        if (dxx >= 9.9999994999999) {
            dxx *= 1.0e-1;
            xp10++;
        }
        sign = '+';
    }
    else
    {
        if (dxx < 9.9999994999999e-8) {
            // general case
            if (dxx < 9.9999994999999e-128) {
                if (dxx == 0.0) {
                    return memcpya(start, "0.000000e+00", 12);
                }
                if (dxx < 9.9999994999999e-256) {
                    dxx *= 1.0e256;
                    xp10 |= 256;
                }
                else
                {
                    dxx *= 1.0e128;
                    xp10 |= 128;
                }
            }
            if (dxx < 9.9999994999999e-64) {
                dxx *= 1.0e64;
                xp10 |= 64;
            }
            if (dxx < 9.9999994999999e-32) {
                dxx *= 1.0e32;
                xp10 |= 32;
            }
            if (dxx < 9.9999994999999e-16) {
                dxx *= 1.0e16;
                xp10 |= 16;
            }
            if (dxx < 9.9999994999999e-8) {
                dxx *= 100000000;
                xp10 |= 8;
            }
        }
        if (dxx < 9.999994999999e-4) {
            dxx *= 10000;
            xp10 |= 4;
        }
        if (dxx < 9.9999994999999e-2) {
            dxx *= 100;
            xp10 |= 2;
        }
        if (dxx < 9.9999994999999e-1) {
            dxx *= 10;
            xp10++;
        }
        sign = '-';
    }
    double_bround6(dxx, banker_round7, &quotient, &remainder);
    *start++ = '0' + quotient;
    *start++ = '.';
    start = uitoa_z6(remainder, start);
    *start++ = 'e';
    *start++ = sign;
    if (xp10 >= 100) {
        quotient = xp10 / 100;
        *start++ = '0' + quotient;
        xp10 -= quotient * 100;
    }
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
}

char* ftoa_e(float fxx, char* start)
{
    uint32_t xp10 = 0;
    uint32_t quotient;
    uint32_t remainder;
    char sign;
    if (fxx != fxx) {
        // do this first to avoid generating exception
        return memcpyl3a(start, "nan");
    }
    else if (fxx < 0)
    {
        *start++ = '-';
        fxx = -fxx;
    }
    if (fxx >= 9.9999995e-1) {
        if (fxx >= 9.9999995e15) {
            if (fxx == INFINITY) {
                return memcpyl3a(start, "inf");
            }
            else if (fxx >= 9.9999995e31)
            {
                fxx *= 1.0e-32;
                xp10 |= 32;
            }
            else
            {
                fxx *= 1.0e-16;
                xp10 |= 16;
            }
        }
        if (fxx >= 9.9999995e7) {
            fxx *= 1.0e-8;
            xp10 |= 8;
        }
        if (fxx >= 9.9999995e3) {
            fxx *= 1.0e-4;
            xp10 |= 4;
        }
        if (fxx >= 9.9999995e1) {
            fxx *= 1.0e-2;
            xp10 |= 2;
        }
        if (fxx >= 9.9999995) {
            fxx *= 1.0e-1;
            xp10++;
        }
        sign = '+';
    }
    else
    {
        if (fxx < 9.9999995e-16) {
            if (fxx == 0.0) {
                return memcpya(start, "0.000000e+00", 12);
            }
            else if (fxx < 9.9999995e-32)
            {
                fxx *= 1.0e32;
                xp10 |= 32;
            }
            else
            {
                fxx *= 1.0e16;
                xp10 |= 16;
            }
        }
        if (fxx < 9.9999995e-8) {
            fxx *= 100000000;
            xp10 |= 8;
        }
        if (fxx < 9.9999995e-4) {
            fxx *= 10000;
            xp10 |= 4;
        }
        if (fxx < 9.9999995e-2) {
            fxx *= 100;
            xp10 |= 2;
        }
        if (fxx < 9.9999995e-1) {
            fxx *= 10;
            xp10++;
        }
        sign = '-';
    }
    float_round6(fxx, &quotient, &remainder);
    *start++ = '0' + quotient;
    *start++ = '.';
    start = uitoa_z6(remainder, start);
    *start++ = 'e';
    *start++ = sign;
    return memcpya(start, &(digit2_table[xp10 * 2]), 2);
}

char* dtoa_f_w9p6(double dxx, char* start)
{
    uint32_t quotient;
    uint32_t remainder;
    if (dxx != dxx) {
        return memcpya(start, "      nan", 9);
    }
    else if (dxx < 9.9999994999999)
    {
        if (dxx < 0) {
            *start++ = '-';
            dxx = -dxx;
            if (dxx >= 9.9999994999999) {
                goto dtoa_f_w9p6_10;
            }
        }
        else
        {
            *start++ = ' ';
        }
        double_bround6(dxx, banker_round7, &quotient, &remainder);
        *start++ = '0' + quotient;
    dtoa_f_w9p6_dec:
        *start++ = '.';
        return uitoa_z6(remainder, start);
    }
dtoa_f_w9p6_10:
    if (dxx < 999.99999949999) {
        double_bround6(dxx,
                       (dxx < 99.999999499999) ? banker_round6 : banker_round5,
                       &quotient, &remainder);
        start = uint32toa(quotient, start);
        goto dtoa_f_w9p6_dec;
    }
    if (dxx == INFINITY) {
        return memcpya(start, "      inf", 9);
    }
    start += sprintf(start, "%.6f", dxx);
    return start;
}

char* dtoa_f_w7p4(double dxx, char* start)
{
    const double* br_ptr;
    uint32_t quotient;
    uint32_t remainder;
    if (dxx != dxx) {
        return memcpya(start, "    nan", 7);
    }
    else if (dxx < 9.9999499999999)
    {
        if (dxx < 0) {
            *start++ = '-';
            dxx = -dxx;
            if (dxx >= 9.9999499999999) {
                goto dtoa_f_w7p4_10;
            }
        }
        else
        {
            *start++ = ' ';
        }
        double_bround4(dxx, banker_round9, &quotient, &remainder);
        *start++ = '0' + quotient;
    dtoa_f_w7p4_dec:
        *start++ = '.';
        quotient = remainder / 100;
        remainder -= 100 * quotient;
        return memcpya(memcpya(start, &(digit2_table[quotient * 2]), 2),
                       &(digit2_table[remainder * 2]), 2);
    }
dtoa_f_w7p4_10:
    if (dxx < 99999.999949999) {
        if (dxx < 999.99994999999) {
            if (dxx < 99.999949999999) {
                br_ptr = banker_round8;
            }
            else
            {
                br_ptr = banker_round7;
            }
        }
        else if (dxx < 9999.9999499999)
        {
            br_ptr = banker_round6;
        }
        else
        {
            br_ptr = banker_round5;
        }
        double_bround4(dxx, br_ptr, &quotient, &remainder);
        start = uint32toa(quotient, start);
        goto dtoa_f_w7p4_dec;
    }
    if (dxx == INFINITY) {
        return memcpya(start, "    inf", 7);
    }
    start += sprintf(start, "%.4f", dxx);
    return start;
}

char* dtoa_g(double dxx, char* start)
{
    uint32_t xp10 = 0;
    uint32_t quotient;
    uint32_t remainder;
    if (dxx != dxx) {
        return memcpyl3a(start, "nan");
    }
    else if (dxx < 0)
    {
        *start++ = '-';
        dxx = -dxx;
    }
    if (dxx < 9.9999949999999e-5) {
        // 6 sig fig exponential notation, small
        if (dxx < 9.9999949999999e-16) {
            if (dxx < 9.9999949999999e-128) {
                if (dxx == 0.0) {
                    *start = '0';
                    return &(start[1]);
                }
                else if (dxx < 9.9999949999999e-256)
                {
                    dxx *= 1.0e256;
                    xp10 |= 256;
                }
                else
                {
                    dxx *= 1.0e128;
                    xp10 |= 128;
                }
            }
            if (dxx < 9.9999949999999e-64) {
                dxx *= 1.0e64;
                xp10 |= 64;
            }
            if (dxx < 9.9999949999999e-32) {
                dxx *= 1.0e32;
                xp10 |= 32;
            }
            if (dxx < 9.9999949999999e-16) {
                dxx *= 1.0e16;
                xp10 |= 16;
            }
        }
        if (dxx < 9.9999949999999e-8) {
            dxx *= 100000000;
            xp10 |= 8;
        }
        if (dxx < 9.9999949999999e-4) {
            dxx *= 10000;
            xp10 |= 4;
        }
        if (dxx < 9.9999949999999e-2) {
            dxx *= 100;
            xp10 |= 2;
        }
        if (dxx < 9.9999949999999e-1) {
            dxx *= 10;
            xp10++;
        }
        double_bround5(dxx, banker_round8, &quotient, &remainder);
        start = memcpya(qrtoa_1p5(quotient, remainder, start), "e-", 2);
        if (xp10 >= 100) {
            quotient = xp10 / 100;
            *start++ = '0' + quotient;
            xp10 -= 100 * quotient;
        }
        return memcpya(start, &(digit2_table[xp10 * 2]), 2);
    }
    else if (dxx >= 999999.49999999)
    {
        // 6 sig fig exponential notation, large
        if (dxx >= 9.9999949999999e15) {
            if (dxx >= 9.9999949999999e127) {
                if (dxx == INFINITY) {
                    return memcpyl3a(start, "inf");
                }
                else if (dxx >= 9.9999949999999e255)
                {
                    dxx *= 1.0e-256;
                    xp10 |= 256;
                }
                else
                {
                    dxx *= 1.0e-128;
                    xp10 |= 128;
                }
            }
            if (dxx >= 9.9999949999999e63) {
                dxx *= 1.0e-64;
                xp10 |= 64;
            }
            if (dxx >= 9.9999949999999e31) {
                dxx *= 1.0e-32;
                xp10 |= 32;
            }
            if (dxx >= 9.9999949999999e15) {
                dxx *= 1.0e-16;
                xp10 |= 16;
            }
        }
        if (dxx >= 9.9999949999999e7) {
            dxx *= 1.0e-8;
            xp10 |= 8;
        }
        if (dxx >= 9.9999949999999e3) {
            dxx *= 1.0e-4;
            xp10 |= 4;
        }
        if (dxx >= 9.9999949999999e1) {
            dxx *= 1.0e-2;
            xp10 |= 2;
        }
        if (dxx >= 9.9999949999999e0) {
            dxx *= 1.0e-1;
            xp10++;
        }
        double_bround5(dxx, banker_round8, &quotient, &remainder);
        start = memcpya(qrtoa_1p5(quotient, remainder, start), "e+", 2);
        if (xp10 >= 100) {
            quotient = xp10 / 100;
            *start++ = '0' + quotient;
            xp10 -= 100 * quotient;
        }
        return memcpya(start, &(digit2_table[xp10 * 2]), 2);
    }
    else if (dxx >= 0.99999949999999)
    {
        return dtoa_so6(dxx, start);
    }
    else
    {
        // 6 sig fig decimal, no less than ~0.0001
        start = memcpya(start, "0.", 2);
        if (dxx < 9.9999949999999e-3) {
            dxx *= 100;
            start = memcpya(start, "00", 2);
        }
        if (dxx < 9.9999949999999e-2) {
            dxx *= 10;
            *start++ = '0';
        }
        return uitoa_trunc6(double_bround(dxx * 1000000, banker_round8), start);
    }
}

char* dtoa_g_wxp3(double dxx, uint32_t min_width, char* start)
{
    assert(min_width >= 5);
    uint32_t xp10 = 0;
    char wbuf[16];
    char* wpos = wbuf;
    uint32_t quotient;
    uint32_t remainder;
    if (dxx != dxx) {
        memcpy(memseta(start, 32, min_width - 4), " nan", 4);
        return &(start[min_width]);
    }
    else if (dxx < 0)
    {
        *wpos++ = '-';
        dxx = -dxx;
    }
    if (dxx < 9.9949999999999e-5) {
        // 3 sig fig exponential notation, small
        if (dxx < 9.9949999999999e-16) {
            if (dxx < 9.9949999999999e-128) {
                if (dxx == 0.0) {
                    memset(start, 32, min_width - 1);
                    start[min_width - 1] = '0';
                    return &(start[min_width]);
                }
                else if (dxx < 9.9949999999999e-256)
                {
                    dxx *= 1.0e256;
                    xp10 |= 256;
                }
                else
                {
                    dxx *= 1.0e128;
                    xp10 |= 128;
                }
            }
            if (dxx < 9.9949999999999e-64) {
                dxx *= 1.0e64;
                xp10 |= 64;
            }
            if (dxx < 9.9949999999999e-32) {
                dxx *= 1.0e32;
                xp10 |= 32;
            }
            if (dxx < 9.9949999999999e-16) {
                dxx *= 1.0e16;
                xp10 |= 16;
            }
        }
        if (dxx < 9.9949999999999e-8) {
            dxx *= 100000000;
            xp10 |= 8;
        }
        if (dxx < 9.9949999999999e-4) {
            dxx *= 10000;
            xp10 |= 4;
        }
        if (dxx < 9.9949999999999e-2) {
            dxx *= 100;
            xp10 |= 2;
        }
        if (dxx < 9.9949999999999e-1) {
            dxx *= 10;
            xp10++;
        }
        double_bround2(dxx, banker_round11, &quotient, &remainder);
        wpos = qrtoa_1p2(quotient, remainder, wpos);
        remainder = wpos - wbuf;
        if (xp10 >= 100) {
            if (remainder < min_width - 5) {
                memcpy(memseta(start, 32, min_width - 5 - remainder), wbuf,
                       remainder);
                start = &(start[min_width - 5]);
            }
            else
            {
                start = memcpya(start, wbuf, remainder);
            }
            quotient = xp10 / 100;
            start = memcpyax(start, "e-", 2, '0' + quotient);
            xp10 -= 100 * quotient;
        }
        else
        {
            if (remainder < min_width - 4) {
                memcpy(memseta(start, 32, min_width - 4 - remainder), wbuf,
                       remainder);
                start = &(start[min_width - 4]);
            }
            else
            {
                start = memcpya(start, wbuf, remainder);
            }
            start = memcpya(start, "e-", 2);
        }
        return memcpya(start, &(digit2_table[xp10 * 2]), 2);
    }
    else if (dxx >= 999.49999999999)
    {
        // 3 sig fig exponential notation, large
        if (dxx >= 9.9949999999999e15) {
            if (dxx >= 9.9949999999999e127) {
                if (dxx == INFINITY) {
                    start = memseta(start, 32, min_width - 4);
                    if (wpos == wbuf) {
                        return memcpya(start, " inf", 4);
                    }
                    else
                    {
                        return memcpya(start, "-inf", 4);
                    }
                }
                else if (dxx >= 9.9949999999999e255)
                {
                    dxx *= 1.0e-256;
                    xp10 |= 256;
                }
                else
                {
                    dxx *= 1.0e-128;
                    xp10 |= 128;
                }
            }
            if (dxx >= 9.9949999999999e63) {
                dxx *= 1.0e-64;
                xp10 |= 64;
            }
            if (dxx >= 9.9949999999999e31) {
                dxx *= 1.0e-32;
                xp10 |= 32;
            }
            if (dxx >= 9.9949999999999e15) {
                dxx *= 1.0e-16;
                xp10 |= 16;
            }
        }
        if (dxx >= 9.9949999999999e7) {
            dxx *= 1.0e-8;
            xp10 |= 8;
        }
        if (dxx >= 9.9949999999999e3) {
            dxx *= 1.0e-4;
            xp10 |= 4;
        }
        if (dxx >= 9.9949999999999e1) {
            dxx *= 1.0e-2;
            xp10 |= 2;
        }
        if (dxx >= 9.9949999999999e0) {
            dxx *= 1.0e-1;
            xp10++;
        }
        double_bround2(dxx, banker_round11, &quotient, &remainder);
        wpos = qrtoa_1p2(quotient, remainder, wpos);
        remainder = wpos - wbuf;
        if (xp10 >= 100) {
            if (remainder < min_width - 5) {
                memcpy(memseta(start, 32, min_width - 5 - remainder), wbuf,
                       remainder);
                start = &(start[min_width - 5]);
            }
            else
            {
                start = memcpya(start, wbuf, remainder);
            }
            quotient = xp10 / 100;
            start = memcpyax(start, "e+", 2, '0' + quotient);
            xp10 -= 100 * quotient;
        }
        else
        {
            if (remainder < min_width - 4) {
                memcpy(memseta(start, 32, min_width - 4 - remainder), wbuf,
                       remainder);
                start = &(start[min_width - 4]);
            }
            else
            {
                start = memcpya(start, wbuf, remainder);
            }
            start = memcpya(start, "e+", 2);
        }
        return memcpya(start, &(digit2_table[xp10 * 2]), 2);
    }
    else
    {
        if (dxx >= 0.99949999999999) {
            wpos = dtoa_so3(dxx, wpos);
        }
        else
        {
            // 3 sig fig decimal, no less than ~0.001
            wpos = memcpya(wpos, "0.", 2);
            if (dxx < 9.9949999999999e-3) {
                dxx *= 100;
                wpos = memcpya(wpos, "00", 2);
            }
            if (dxx < 9.9949999999999e-2) {
                dxx *= 10;
                *wpos++ = '0';
            }
            wpos =
                uitoa_trunc3(double_bround(dxx * 1000, banker_round11), wpos);
        }
        remainder = wpos - wbuf;
        if (remainder < min_width) {
            memcpy(memseta(start, 32, min_width - remainder), wbuf, remainder);
            return &(start[min_width]);
        }
        else
        {
            return memcpya(start, wbuf, remainder);
        }
    }
}

char* dtoa_g_wxp4(double dxx, uint32_t min_width, char* start)
{
    uint32_t xp10 = 0;
    char wbuf[16];
    char* wpos = wbuf;
    uint32_t quotient;
    uint32_t remainder;
    if (dxx != dxx) {
        if (min_width > 3) {
            start = memseta(start, 32, min_width - 3);
        }
        return memcpyl3a(start, "nan");
    }
    else if (dxx < 0)
    {
        *wpos++ = '-';
        dxx = -dxx;
    }
    if (dxx < 9.9994999999999e-5) {
        // 4 sig fig exponential notation, small
        if (dxx < 9.9994999999999e-16) {
            if (dxx < 9.9994999999999e-128) {
                if (dxx == 0.0) {
                    memset(start, 32, min_width - 1);
                    start[min_width - 1] = '0';
                    return &(start[min_width]);
                }
                else if (dxx < 9.9994999999999e-256)
                {
                    dxx *= 1.0e256;
                    xp10 |= 256;
                }
                else
                {
                    dxx *= 1.0e128;
                    xp10 |= 128;
                }
            }
            if (dxx < 9.9994999999999e-64) {
                dxx *= 1.0e64;
                xp10 |= 64;
            }
            if (dxx < 9.9994999999999e-32) {
                dxx *= 1.0e32;
                xp10 |= 32;
            }
            if (dxx < 9.9994999999999e-16) {
                dxx *= 1.0e16;
                xp10 |= 16;
            }
        }
        if (dxx < 9.9994999999999e-8) {
            dxx *= 100000000;
            xp10 |= 8;
        }
        if (dxx < 9.9994999999999e-4) {
            dxx *= 10000;
            xp10 |= 4;
        }
        if (dxx < 9.9994999999999e-2) {
            dxx *= 100;
            xp10 |= 2;
        }
        if (dxx < 9.9994999999999e-1) {
            dxx *= 10;
            xp10++;
        }
        double_bround3(dxx, banker_round10, &quotient, &remainder);
        wpos = qrtoa_1p3(quotient, remainder, wpos);
        remainder = wpos - wbuf;
        if (xp10 >= 100) {
            if (remainder + 5 < min_width) {
                memcpy(memseta(start, 32, min_width - (remainder + 5)), wbuf,
                       remainder);
                start = &(start[min_width - 5]);
            }
            else
            {
                start = memcpya(start, wbuf, remainder);
            }
            quotient = xp10 / 100;
            start = memcpyax(start, "e-", 2, '0' + quotient);
            xp10 -= 100 * quotient;
        }
        else
        {
            if (remainder + 4 < min_width) {
                memcpy(memseta(start, 32, min_width - (remainder + 4)), wbuf,
                       remainder);
                start = &(start[min_width - 4]);
            }
            else
            {
                start = memcpya(start, wbuf, remainder);
            }
            start = memcpya(start, "e-", 2);
        }
        return memcpya(start, &(digit2_table[xp10 * 2]), 2);
    }
    else if (dxx >= 9999.4999999999)
    {
        // 4 sig fig exponential notation, large
        if (dxx >= 9.9994999999999e15) {
            if (dxx >= 9.9994999999999e127) {
                if (dxx == INFINITY) {
                    if (min_width > 4) {
                        start = memseta(start, 32, min_width - 4);
                    }
                    if (wpos == wbuf) {
                        return memcpya(start, " inf", 4);
                    }
                    else
                    {
                        return memcpya(start, "-inf", 4);
                    }
                }
                else if (dxx >= 9.9994999999999e255)
                {
                    dxx *= 1.0e-256;
                    xp10 |= 256;
                }
                else
                {
                    dxx *= 1.0e-128;
                    xp10 |= 128;
                }
            }
            if (dxx >= 9.9994999999999e63) {
                dxx *= 1.0e-64;
                xp10 |= 64;
            }
            if (dxx >= 9.9994999999999e31) {
                dxx *= 1.0e-32;
                xp10 |= 32;
            }
            if (dxx >= 9.9994999999999e15) {
                dxx *= 1.0e-16;
                xp10 |= 16;
            }
        }
        if (dxx >= 9.9994999999999e7) {
            dxx *= 1.0e-8;
            xp10 |= 8;
        }
        if (dxx >= 9.9994999999999e3) {
            dxx *= 1.0e-4;
            xp10 |= 4;
        }
        if (dxx >= 9.9994999999999e1) {
            dxx *= 1.0e-2;
            xp10 |= 2;
        }
        if (dxx >= 9.9994999999999e0) {
            dxx *= 1.0e-1;
            xp10++;
        }
        double_bround3(dxx, banker_round10, &quotient, &remainder);
        wpos = qrtoa_1p3(quotient, remainder, wpos);
        remainder = wpos - wbuf;
        if (xp10 >= 100) {
            if (remainder + 5 < min_width) {
                memcpy(memseta(start, 32, min_width - (remainder + 5)), wbuf,
                       remainder);
                start = &(start[min_width - 5]);
            }
            else
            {
                start = memcpya(start, wbuf, remainder);
            }
            quotient = xp10 / 100;
            start = memcpyax(start, "e+", 2, '0' + quotient);
            xp10 -= 100 * quotient;
        }
        else
        {
            if (remainder + 4 < min_width) {
                memcpy(memseta(start, 32, min_width - (remainder + 4)), wbuf,
                       remainder);
                start = &(start[min_width - 4]);
            }
            else
            {
                start = memcpya(start, wbuf, remainder);
            }
            start = memcpya(start, "e+", 2);
        }
        return memcpya(start, &(digit2_table[xp10 * 2]), 2);
    }
    else
    {
        if (dxx >= 0.99994999999999) {
            wpos = dtoa_so4(dxx, wpos);
        }
        else
        {
            // 4 sig fig decimal, no less than ~0.0001
            wpos = memcpya(wpos, "0.", 2);
            if (dxx < 9.9994999999999e-3) {
                dxx *= 100;
                wpos = memcpya(wpos, "00", 2);
            }
            if (dxx < 9.9994999999999e-2) {
                dxx *= 10;
                *wpos++ = '0';
            }
            wpos =
                uitoa_trunc4(double_bround(dxx * 10000, banker_round10), wpos);
        }
        remainder = wpos - wbuf;
        if (remainder < min_width) {
            memcpy(memseta(start, 32, min_width - remainder), wbuf, remainder);
            return &(start[min_width]);
        }
        else
        {
            return memcpya(start, wbuf, remainder);
        }
    }
}

char* dtoa_g_wxp8(double dxx, uint32_t min_width, char* start)
{
    uint32_t xp10 = 0;
    char wbuf[16];
    char* wpos = wbuf;
    uint32_t quotient;
    uint32_t remainder;
    if (dxx != dxx) {
        if (min_width > 3) {
            start = memseta(start, 32, min_width - 3);
        }
        return memcpyl3a(start, "nan");
    }
    else if (dxx < 0)
    {
        *wpos++ = '-';
        dxx = -dxx;
    }
    if (dxx < 9.9999999499999e-5) {
        // 8 sig fig exponential notation, small
        if (dxx < 9.9999999499999e-16) {
            if (dxx < 9.9999999499999e-128) {
                if (dxx == 0.0) {
                    memset(start, 32, min_width - 1);
                    start[min_width - 1] = '0';
                    return &(start[min_width]);
                }
                else if (dxx < 9.9999999499999e-256)
                {
                    dxx *= 1.0e256;
                    xp10 |= 256;
                }
                else
                {
                    dxx *= 1.0e128;
                    xp10 |= 128;
                }
            }
            if (dxx < 9.9999999499999e-64) {
                dxx *= 1.0e64;
                xp10 |= 64;
            }
            if (dxx < 9.9999999499999e-32) {
                dxx *= 1.0e32;
                xp10 |= 32;
            }
            if (dxx < 9.9999999499999e-16) {
                dxx *= 1.0e16;
                xp10 |= 16;
            }
        }
        if (dxx < 9.9999999499999e-8) {
            dxx *= 100000000;
            xp10 |= 8;
        }
        if (dxx < 9.9999999499999e-4) {
            dxx *= 10000;
            xp10 |= 4;
        }
        if (dxx < 9.9999999499999e-2) {
            dxx *= 100;
            xp10 |= 2;
        }
        if (dxx < 9.9999999499999e-1) {
            dxx *= 10;
            xp10++;
        }
        double_bround7(dxx, banker_round6, &quotient, &remainder);
        wpos = qrtoa_1p7(quotient, remainder, wpos);
        remainder = wpos - wbuf;
        if (xp10 >= 100) {
            if (remainder + 5 < min_width) {
                memcpy(memseta(start, 32, min_width - (remainder + 5)), wbuf,
                       remainder);
                start = &(start[min_width - 5]);
            }
            else
            {
                start = memcpya(start, wbuf, remainder);
            }
            quotient = xp10 / 100;
            start = memcpyax(start, "e-", 2, '0' + quotient);
            xp10 -= 100 * quotient;
        }
        else
        {
            if (remainder + 4 < min_width) {
                memcpy(memseta(start, 32, min_width - (remainder + 4)), wbuf,
                       remainder);
                start = &(start[min_width - 4]);
            }
            else
            {
                start = memcpya(start, wbuf, remainder);
            }
            start = memcpya(start, "e-", 2);
        }
        return memcpya(start, &(digit2_table[xp10 * 2]), 2);
    }
    else if (dxx >= 99999999.499999)
    {
        // 8 sig fig exponential notation, large
        if (dxx >= 9.9999999499999e15) {
            if (dxx >= 9.9999999499999e127) {
                if (dxx == INFINITY) {
                    if (min_width > 4) {
                        start = memseta(start, 32, min_width - 4);
                    }
                    if (wpos == wbuf) {
                        return memcpya(start, " inf", 4);
                    }
                    else
                    {
                        return memcpya(start, "-inf", 4);
                    }
                }
                else if (dxx >= 9.9999999499999e255)
                {
                    dxx *= 1.0e-256;
                    xp10 |= 256;
                }
                else
                {
                    dxx *= 1.0e-128;
                    xp10 |= 128;
                }
            }
            if (dxx >= 9.9999999499999e63) {
                dxx *= 1.0e-64;
                xp10 |= 64;
            }
            if (dxx >= 9.9999999499999e31) {
                dxx *= 1.0e-32;
                xp10 |= 32;
            }
            if (dxx >= 9.9999999499999e15) {
                dxx *= 1.0e-16;
                xp10 |= 16;
            }
        }
        if (dxx >= 9.9999999499999e7) {
            dxx *= 1.0e-8;
            xp10 |= 8;
        }
        if (dxx >= 9.9999999499999e3) {
            dxx *= 1.0e-4;
            xp10 |= 4;
        }
        if (dxx >= 9.9999999499999e1) {
            dxx *= 1.0e-2;
            xp10 |= 2;
        }
        if (dxx >= 9.9999999499999e0) {
            dxx *= 1.0e-1;
            xp10++;
        }
        double_bround7(dxx, banker_round6, &quotient, &remainder);
        wpos = qrtoa_1p7(quotient, remainder, wpos);
        remainder = wpos - wbuf;
        if (xp10 >= 100) {
            if (remainder + 5 < min_width) {
                memcpy(memseta(start, 32, min_width - (remainder + 5)), wbuf,
                       remainder);
                start = &(start[min_width - 5]);
            }
            else
            {
                start = memcpya(start, wbuf, remainder);
            }
            quotient = xp10 / 100;
            start = memcpyax(start, "e+", 2, '0' + quotient);
            xp10 -= 100 * quotient;
        }
        else
        {
            if (remainder + 4 < min_width) {
                memcpy(memseta(start, 32, min_width - (remainder + 4)), wbuf,
                       remainder);
                start = &(start[min_width - 4]);
            }
            else
            {
                start = memcpya(start, wbuf, remainder);
            }
            start = memcpya(start, "e+", 2);
        }
        return memcpya(start, &(digit2_table[xp10 * 2]), 2);
    }
    else
    {
        if (dxx >= 0.99999999499999) {
            wpos = dtoa_so8(dxx, wpos);
        }
        else
        {
            // 8 sig fig decimal, no less than ~0.0001
            wpos = memcpya(wpos, "0.", 2);
            if (dxx < 9.9999999499999e-3) {
                dxx *= 100;
                wpos = memcpya(wpos, "00", 2);
            }
            if (dxx < 9.9999999499999e-2) {
                dxx *= 10;
                *wpos++ = '0';
            }
            wpos = uitoa_trunc8(double_bround(dxx * 100000000, banker_round6),
                                wpos);
        }
        remainder = wpos - wbuf;
        if (remainder < min_width) {
            memcpy(memseta(start, 32, min_width - remainder), wbuf, remainder);
            return &(start[min_width]);
        }
        else
        {
            return memcpya(start, wbuf, remainder);
        }
    }
}

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


uint32_t next_unset_unsafe(const uintptr_t* bitarr, uint32_t loc)
{
    const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
    uintptr_t ulii = (~(*bitarr_ptr)) >> (loc % BITCT);
    if (ulii) {
        return loc + CTZLU(ulii);
    }
    do
    {
        ulii = *(++bitarr_ptr);
    } while (ulii == ~ZEROLU);
    return ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + CTZLU(~ulii);
}

#ifdef __LP64__
uintptr_t next_unset_ul_unsafe(const uintptr_t* bitarr, uintptr_t loc)
{
    const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
    uintptr_t ulii = (~(*bitarr_ptr)) >> (loc % BITCT);
    if (ulii) {
        return loc + CTZLU(ulii);
    }
    do
    {
        ulii = *(++bitarr_ptr);
    } while (ulii == ~ZEROLU);
    return (((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + CTZLU(~ulii));
}
#endif

uint32_t next_unset(const uintptr_t* bitarr, uint32_t loc, uint32_t ceil)
{
    // safe version.
    assert(ceil >= 1);
    const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
    uintptr_t ulii = (~(*bitarr_ptr)) >> (loc % BITCT);
    const uintptr_t* bitarr_last;
    if (ulii) {
        loc += CTZLU(ulii);
        return MINV(loc, ceil);
    }
    bitarr_last = &(bitarr[(ceil - 1) / BITCT]);
    do
    {
        if (bitarr_ptr >= bitarr_last) {
            return ceil;
        }
        ulii = *(++bitarr_ptr);
    } while (ulii == ~ZEROLU);
    loc = ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + CTZLU(~ulii);
    return MINV(loc, ceil);
}

#ifdef __LP64__
uintptr_t next_unset_ul(const uintptr_t* bitarr, uintptr_t loc, uintptr_t ceil)
{
    const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
    uintptr_t ulii = (~(*bitarr_ptr)) >> (loc % BITCT);
    const uintptr_t* bitarr_last;
    if (ulii) {
        ulii = loc + CTZLU(ulii);
        return MINV(ulii, ceil);
    }
    bitarr_last = &(bitarr[(ceil - 1) / BITCT]);
    do
    {
        if (bitarr_ptr >= bitarr_last) {
            return ceil;
        }
        ulii = *(++bitarr_ptr);
    } while (ulii == ~ZEROLU);
    ulii = ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + CTZLU(~ulii);
    return MINV(ulii, ceil);
}
#endif

uint32_t next_set_unsafe(const uintptr_t* bitarr, uint32_t loc)
{
    const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
    uintptr_t ulii = (*bitarr_ptr) >> (loc % BITCT);
    if (ulii) {
        return loc + CTZLU(ulii);
    }
    do
    {
        ulii = *(++bitarr_ptr);
    } while (!ulii);
    return ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + CTZLU(ulii);
}

#ifdef __LP64__
uintptr_t next_set_ul_unsafe(const uintptr_t* bitarr, uintptr_t loc)
{
    const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
    uintptr_t ulii = (*bitarr_ptr) >> (loc % BITCT);
    if (ulii) {
        return loc + CTZLU(ulii);
    }
    do
    {
        ulii = *(++bitarr_ptr);
    } while (!ulii);
    return ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + CTZLU(ulii);
}
#endif

uint32_t next_set(const uintptr_t* bitarr, uint32_t loc, uint32_t ceil)
{
    const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
    uintptr_t ulii = (*bitarr_ptr) >> (loc % BITCT);
    const uintptr_t* bitarr_last;
    uint32_t rval;
    if (ulii) {
        rval = loc + CTZLU(ulii);
        return MINV(rval, ceil);
    }
    bitarr_last = &(bitarr[(ceil - 1) / BITCT]);
    do
    {
        if (bitarr_ptr >= bitarr_last) {
            return ceil;
        }
        ulii = *(++bitarr_ptr);
    } while (!ulii);
    rval = ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + CTZLU(ulii);
    return MINV(rval, ceil);
}

#ifdef __LP64__
uintptr_t next_set_ul(const uintptr_t* bitarr, uintptr_t loc, uintptr_t ceil)
{
    const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
    uintptr_t ulii = (*bitarr_ptr) >> (loc % BITCT);
    const uintptr_t* bitarr_last;
    if (ulii) {
        ulii = loc + CTZLU(ulii);
        return MINV(ulii, ceil);
    }
    bitarr_last = &(bitarr[(ceil - 1) / BITCT]);
    do
    {
        if (bitarr_ptr >= bitarr_last) {
            return ceil;
        }
        ulii = *(++bitarr_ptr);
    } while (!ulii);
    ulii = ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + CTZLU(ulii);
    return MINV(ulii, ceil);
}
#endif

uint32_t prev_unset_unsafe(const uintptr_t* bitarr, uint32_t loc)
{
    // unlike the next_{un}set family, this always returns a STRICTLY earlier
    // position
    const uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
    uint32_t remainder = loc % BITCT;
    uintptr_t ulii;
    if (remainder) {
        ulii = (~(*bitarr_ptr)) & ((ONELU << remainder) - ONELU);
        if (ulii) {
            return (loc | (BITCT - 1)) - CLZLU(ulii);
        }
    }
    do
    {
        ulii = ~(*(--bitarr_ptr));
    } while (!ulii);
    return ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + BITCT - 1 - CLZLU(ulii);
}

/*
uint32_t prev_unset(uintptr_t* bitarr, uint32_t loc, uint32_t floor) {
  uintptr_t* bitarr_ptr = &(bitarr[loc / BITCT]);
  uint32_t remainder = loc % BITCT;
  uintptr_t* bitarr_first;
  uintptr_t ulii;
  if (remainder) {
    ulii = (~(*bitarr_ptr)) & ((ONELU << remainder) - ONELU);
    if (ulii) {
      loc = (loc | (BITCT - 1)) - CLZLU(ulii);
      return MAXV(loc, floor);
    }
  }
  bitarr_first = &(bitarr[floor / BITCT]);
  do {
    if (bitarr_ptr == bitarr_first) {
      return floor;
    }
    ulii = ~(*(--bitarr_ptr));
  } while (!ulii);
  loc = ((uintptr_t)(bitarr_ptr - bitarr)) * BITCT + BITCT - 1 - CLZLU(ulii);
  return MAXV(loc, floor);
}
*/


int32_t bigstack_calloc_ui(uintptr_t ct, uint32_t** uip_ptr)
{
    *uip_ptr = (uint32_t*) bigstack_alloc(ct * sizeof(int32_t));
    if (!(*uip_ptr)) {
        return 1;
    }
    fill_uint_zero(ct, *uip_ptr);
    return 0;
}

// MurmurHash3, from
// https://code.google.com/p/smhasher/source/browse/trunk/MurmurHash3.cpp

static inline uint32_t rotl32(uint32_t x, int8_t r)
{
    return (x << r) | (x >> (32 - r));
}

static inline uint32_t getblock32(const uint32_t* p, int i) { return p[i]; }

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

static inline uint32_t fmix32(uint32_t h)
{
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;

    return h;
}

uint32_t murmurhash3_32(const void* key, uint32_t len)
{
    const uint8_t* data = (const uint8_t*) key;
    const int32_t nblocks = len / 4;

    uint32_t h1 = 0;
    // uint32_t h1 = seed;

    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;

    //----------
    // body

    const uint32_t* blocks = (const uint32_t*) (data + nblocks * 4);

    int32_t i;
    uint32_t k1;
    for (i = -nblocks; i; i++) {
        k1 = getblock32(blocks, i);

        k1 *= c1;
        k1 = rotl32(k1, 15);
        k1 *= c2;

        h1 ^= k1;
        h1 = rotl32(h1, 13);
        h1 = h1 * 5 + 0xe6546b64;
    }

    //----------
    // tail

    const uint8_t* tail = (const uint8_t*) (data + nblocks * 4);

    k1 = 0;

    switch (len & 3)
    {
    case 3:
        k1 ^= tail[2] << 16;
        // fall through
    case 2:
        k1 ^= tail[1] << 8;
        // fall through
    case 1:
        k1 ^= tail[0];
        k1 *= c1;
        k1 = rotl32(k1, 15);
        k1 *= c2;
        h1 ^= k1;
    }

    //----------
    // finalization

    h1 ^= len;

    return fmix32(h1);
}

uint32_t is_composite6(uintptr_t num)
{
    // assumes num is congruent to 1 or 5 mod 6.
    // can speed this up by ~50% by hardcoding avoidance of multiples of 5/7,
    // but this isn't currently a bottleneck so I'll keep this simple
    uintptr_t divisor = 5;
    while (divisor * divisor <= num) {
        if (!(num % divisor)) {
            return 1;
        }
        divisor += 2;
        if (!(num % divisor)) {
            return 1;
        }
        divisor += 4;
    }
    return 0;
}

uintptr_t geqprime(uintptr_t floor)
{
    // assumes floor is odd and greater than 1.  Returns 5 if floor = 3,
    // otherwise returns the first prime >= floor.
    uintptr_t ulii = floor % 3;
    if (!ulii) {
        floor += 2;
    }
    else if (ulii == 1)
    {
        goto geqprime_1mod6;
    }
    while (is_composite6(floor)) {
        floor += 2;
    geqprime_1mod6:
        if (!is_composite6(floor)) {
            return floor;
        }
        floor += 4;
    }
    return floor;
}

int32_t populate_id_htable(uintptr_t unfiltered_ct,
                           const uintptr_t* exclude_arr, uintptr_t item_ct,
                           const char* item_ids, uintptr_t max_id_len,
                           uint32_t store_dups, uint32_t id_htable_size,
                           uint32_t* id_htable)
{
    // While unique IDs are normally assumed (and enforced) here, --extract and
    // --exclude are an exception, since we want to be able to e.g. exclude all
    // variants named '.'.  Since there could be millions of them, ordinary
    // O(n^2) hash table duplicate resolution is unacceptably slow; instead, we
    // allocate additional linked lists past the end of id_htable to track all
    // unfiltered indexes of duplicate names.  (This requires the
    // alloc_and_populate_id_htable interface; bigstack_end_alloc doesn't work
    // there.)
    uintptr_t item_uidx = 0;
    uint32_t extra_alloc = 0;
    uint32_t prev_llidx = 0;
    // needs to be synced with extract_exclude_flag_norange()
    uint32_t* extra_alloc_base = (uint32_t*) g_bigstack_base;
    uint32_t item_idx = 0;
    const char* sptr;
    uintptr_t prev_uidx;
    uintptr_t cur_bigstack_left;
    uint32_t max_extra_alloc;
    uint32_t slen;
    uint32_t hashval;
    uint32_t next_incr;
    uint32_t top_diff;
    uint32_t hash_result;
    uint32_t cur_dup;
    fill_uint_one(id_htable_size, id_htable);
    if (!store_dups) {
        for (; item_idx < item_ct; item_uidx++, item_idx++) {
            next_unset_ul_unsafe_ck(exclude_arr, &item_uidx);
            sptr = &(item_ids[item_uidx * max_id_len]);
            slen = strlen(sptr);
            hashval = murmurhash3_32(sptr, slen) % id_htable_size;
            next_incr = 1;
            while (1) {
                hash_result = id_htable[hashval];
                if (hash_result == 0xffffffffU) {
                    id_htable[hashval] = item_uidx;
                    break;
                }
                else if (!memcmp(sptr, &(item_ids[hash_result * max_id_len]),
                                 slen + 1))
                {
                    // could add an allow_dups parameter which controls whether
                    // this is an error
                    LOGERRPRINTFWW("Error: Duplicate ID '%s'.\n", sptr);
                    return RET_INVALID_FORMAT;
                }
                // defend against overflow
                top_diff = id_htable_size - hashval;
                if (top_diff > next_incr) {
                    hashval += next_incr;
                }
                else
                {
                    hashval = next_incr - top_diff;
                }
                next_incr += 2; // quadratic probing
            }
        }
    }
    else
    {
        cur_bigstack_left = bigstack_left();
#ifdef __LP64__
        if (cur_bigstack_left >= 0x400000000LLU) {
            max_extra_alloc = 0xfffffffeU;
        }
        else
        {
            max_extra_alloc = cur_bigstack_left / sizeof(int32_t);
        }
#else
        max_extra_alloc = cur_bigstack_left / sizeof(int32_t);
#endif
        for (; item_idx < item_ct; item_uidx++, item_idx++) {
            next_unset_ul_unsafe_ck(exclude_arr, &item_uidx);
            sptr = &(item_ids[item_uidx * max_id_len]);
            slen = strlen(sptr);
            hashval = murmurhash3_32(sptr, slen) % id_htable_size;
            next_incr = 1;
            while (1) {
                hash_result = id_htable[hashval];
                if (hash_result == 0xffffffffU) {
                    id_htable[hashval] = item_uidx;
                    break;
                }
                else
                {
                    cur_dup = hash_result >> 31;
                    if (cur_dup) {
                        prev_llidx = hash_result << 1;
                        prev_uidx = extra_alloc_base[prev_llidx];
                    }
                    else
                    {
                        prev_uidx = hash_result;
                    }
                    if (!memcmp(sptr, &(item_ids[prev_uidx * max_id_len]),
                                slen + 1))
                    {
                        if (extra_alloc + 4 > max_extra_alloc) {
                            return RET_NOMEM;
                        }
                        // point to linked list entry instead
                        if (!cur_dup) {
                            extra_alloc_base[extra_alloc] = hash_result;
                            extra_alloc_base[extra_alloc + 1] =
                                0xffffffffU; // list end
                            prev_llidx = extra_alloc;
                            extra_alloc += 2;
                        }
                        extra_alloc_base[extra_alloc] = item_uidx;
                        extra_alloc_base[extra_alloc + 1] = prev_llidx;
                        id_htable[hashval] = 0x80000000U | (extra_alloc >> 1);
                        extra_alloc += 2;
                        break; // bugfix
                    }
                }
                top_diff = id_htable_size - hashval;
                if (top_diff > next_incr) {
                    hashval += next_incr;
                }
                else
                {
                    hashval = next_incr - top_diff;
                }
                next_incr += 2;
            }
        }
        if (extra_alloc) {
            bigstack_alloc(extra_alloc * sizeof(int32_t));
        }
    }
    return 0;
}


// hashval computation left to caller since this is frequently used with
// chromosome IDs, where the compiler can optimize the integer modulus
// operation since the hash table size is preset
uint32_t unklen_id_htable_find(const char* cur_id, const char* const* item_ids,
                               const uint32_t* id_htable, uint32_t hashval,
                               uint32_t id_htable_size)
{
    // returns 0xffffffffU on failure
    uint32_t next_incr = 1;
    while (1) {
        const uint32_t hash_result = id_htable[hashval];
        if (hash_result == 0xffffffffU) {
            return 0xffffffffU;
        }
        const char* htable_entry = item_ids[hash_result];
        if (!strcmp(cur_id, htable_entry)) {
            return hash_result;
        }
        const uint32_t top_diff = id_htable_size - hashval;
        if (top_diff > next_incr) {
            hashval += next_incr;
        }
        else
        {
            hashval = next_incr - top_diff;
        }
        next_incr += 2;
    }
}

static inline uint32_t nonstd_chrom_name_htable_find(
    const char* chrom_name, const char* const* nonstd_names,
    const uint32_t* nonstd_id_htable, uint32_t name_slen)
{
    const uint32_t hashval =
        murmurhash3_32(chrom_name, name_slen) % CHROM_NAME_HTABLE_SIZE;
    return unklen_id_htable_find(chrom_name, nonstd_names, nonstd_id_htable,
                                 hashval, CHROM_NAME_HTABLE_SIZE);
}


// Global since species_str() may be called by functions which don't actually
// care about chrom_info.  (chrom_info is really a global variable too, but I
// find it easier to maintain this code when chrom_info dependencies are made
// explicit in the function signatures; in contrast, g_species_singular and
// g_species_plural are just for pretty printing and lend no insight into what
// the functions which reference them are doing.)
const char* g_species_singular = nullptr;
const char* g_species_plural = nullptr;

int32_t init_chrom_info(Chrom_info* chrom_info_ptr)
{
    // "constructor".  initializes with maximum capacity.  doesn't use bigstack.
    // chrom_mask, haploid_mask: bits
    // chrom_file_order, chrom_idx_to_foidx: int32s
    // chrom_fo_vidx_start: int32s, with an extra trailing element
    // nonstd_names: intptr_ts
    // nonstd_id_htable: CHROM_NAME_HTABLE_SIZE int32s

    assert(!(MAX_POSSIBLE_CHROM % VEC_BYTES));
    const uintptr_t vecs_required =
        2 * BITCT_TO_VECCT(MAX_POSSIBLE_CHROM)
        + 3 * (MAX_POSSIBLE_CHROM / VEC_INT32) + 1
        + (MAX_POSSIBLE_CHROM / VEC_WORDS)
        + (CHROM_NAME_HTABLE_SIZE + (VEC_INT32 - 1)) / VEC_INT32;

    // needed for proper cleanup
    chrom_info_ptr->name_ct = 0;
    chrom_info_ptr->incl_excl_name_stack = nullptr;
    if (aligned_malloc(vecs_required * VEC_BYTES,
                       &(chrom_info_ptr->chrom_mask)))
    {
        return RET_NOMEM;
    }
    uintptr_t* alloc_iter =
        &(chrom_info_ptr
              ->chrom_mask[BITCT_TO_VECCT(MAX_POSSIBLE_CHROM) * VEC_WORDS]);
    chrom_info_ptr->haploid_mask = alloc_iter;
    alloc_iter = &(alloc_iter[BITCT_TO_VECCT(MAX_POSSIBLE_CHROM) * VEC_WORDS]);
    chrom_info_ptr->chrom_file_order = (uint32_t*) alloc_iter;
    alloc_iter = &(alloc_iter[(MAX_POSSIBLE_CHROM / VEC_INT32) * VEC_WORDS]);
    chrom_info_ptr->chrom_fo_vidx_start = (uint32_t*) alloc_iter;
    alloc_iter =
        &(alloc_iter[((MAX_POSSIBLE_CHROM / VEC_INT32) + 1) * VEC_WORDS]);
    chrom_info_ptr->chrom_idx_to_foidx = (uint32_t*) alloc_iter;
    alloc_iter = &(alloc_iter[(MAX_POSSIBLE_CHROM / VEC_INT32) * VEC_WORDS]);
    chrom_info_ptr->nonstd_names = (char**) alloc_iter;
    alloc_iter = &(alloc_iter[MAX_POSSIBLE_CHROM]);
    chrom_info_ptr->nonstd_id_htable = (uint32_t*) alloc_iter;
    // alloc_iter = &(alloc_iter[((CHROM_NAME_HTABLE_SIZE + (VEC_INT32 - 1)) /
    // VEC_INT32) * VEC_WORDS]); postpone nonstd_id_htable initialization until
    // first nonstandard ID is loaded fill_uint_one(CHROM_NAME_HTABLE_SIZE,
    // chrom_info_ptr->nonstd_id_htable);
    return 0;
}

// if these are defined within init_species(), they may not persist after
// function exit
static const char species_singular_constants[][7] = {
    "person", "cow", "dog", "horse", "mouse", "plant", "sheep", "sample"};
static const char species_plural_constants[][8] = {
    "people", "cattle", "dogs", "horses", "mice", "plants", "sheep", "samples"};

void init_species(uint32_t species_code, Chrom_info* chrom_info_ptr)
{
    // human: 22, X, Y, XY, MT
    // cow: 29, X, Y, MT
    // dog: 38, X, Y, XY, MT
    // horse: 31, X, Y
    // mouse: 19, X, Y
    // rice: 12
    // sheep: 26, X, Y
    const int32_t species_xymt_codes[] = {
        23, 24, 25, 26, 30, 31, -1, 33, 39, 40, 41, 42, 32, 33,
        -1, -1, 20, 21, -1, -1, -1, -1, -1, -1, 27, 28, -1, -1};
    const uint32_t species_autosome_ct[] = {22, 29, 38, 31, 19, 12, 26};
    const uint32_t species_max_code[] = {26, 33, 42, 33, 21, 12, 28};
    fill_ulong_zero(CHROM_MASK_WORDS, chrom_info_ptr->chrom_mask);
    chrom_info_ptr->output_encoding = 0;
    chrom_info_ptr->zero_extra_chroms = 0;
    chrom_info_ptr->species = species_code;
    chrom_info_ptr->is_include_stack = 0;
    g_species_singular = species_singular_constants[species_code];
    g_species_plural = species_plural_constants[species_code];
    if (species_code != SPECIES_UNKNOWN) {
        // these are assumed to be already initialized in the SPECIES_UNKNOWN
        // case

        // bugfix: haploid_mask was being cleared in --chr-set case
        fill_ulong_zero(CHROM_MASK_WORDS, chrom_info_ptr->haploid_mask);
        memcpy(chrom_info_ptr->xymt_codes,
               &(species_xymt_codes[species_code * XYMT_OFFSET_CT]),
               XYMT_OFFSET_CT * sizeof(int32_t));
        chrom_info_ptr->autosome_ct = species_autosome_ct[species_code];
        chrom_info_ptr->max_code = species_max_code[species_code];
        switch (species_code)
        {
        case SPECIES_HUMAN: chrom_info_ptr->haploid_mask[0] = 0x1800000; break;
        case SPECIES_COW: chrom_info_ptr->haploid_mask[0] = 0xc0000000LU; break;
        case SPECIES_DOG:
#ifdef __LP64__
            chrom_info_ptr->haploid_mask[0] = 0x18000000000LLU;
#else
            chrom_info_ptr->haploid_mask[1] = 0x180;
#endif
            break;
        case SPECIES_HORSE:
#ifdef __LP64__
            chrom_info_ptr->haploid_mask[0] = 0x300000000LLU;
#else
            chrom_info_ptr->haploid_mask[1] = 3;
#endif
            break;
        case SPECIES_MOUSE: chrom_info_ptr->haploid_mask[0] = 0x300000; break;
        case SPECIES_RICE: chrom_info_ptr->haploid_mask[0] = 0x1fff; break;
        case SPECIES_SHEEP: chrom_info_ptr->haploid_mask[0] = 0x18000000; break;
        }
    }
    fill_uint_one(chrom_info_ptr->max_code + 1,
                  chrom_info_ptr->chrom_idx_to_foidx);
}

void init_default_chrom_mask(Chrom_info* chrom_info_ptr)
{
    if (chrom_info_ptr->species != SPECIES_UNKNOWN) {
        fill_all_bits(chrom_info_ptr->max_code + 1, chrom_info_ptr->chrom_mask);
    }
    else
    {
        fill_all_bits(chrom_info_ptr->autosome_ct + 1,
                      chrom_info_ptr->chrom_mask);
        // --chr-set support
        for (uint32_t xymt_idx = 0; xymt_idx < XYMT_OFFSET_CT; ++xymt_idx) {
            int32_t cur_code = chrom_info_ptr->xymt_codes[xymt_idx];
            if (cur_code != -1) {
                set_bit(chrom_info_ptr->xymt_codes[xymt_idx],
                        chrom_info_ptr->chrom_mask);
            }
        }
    }
}


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

int32_t get_chrom_code(const char* chrom_name, const Chrom_info* chrom_info_ptr,
                       uint32_t name_slen)
{
    // requires chrom_name to be null-terminated
    // in practice, name_slen will usually already be known, may as well avoid
    // redundant strlen() calls even though this uglifies the interface
    // does not perform exhaustive error-checking
    // -1 = --allow-extra-chr ok, -2 = total fail
    const int32_t chrom_code_raw = get_chrom_code_raw(chrom_name);
    if (((const uint32_t) chrom_code_raw) <= chrom_info_ptr->max_code) {
        return chrom_code_raw;
    }
    if (chrom_code_raw != -1) {
        if (chrom_code_raw >= MAX_POSSIBLE_CHROM) {
            return chrom_info_ptr
                ->xymt_codes[chrom_code_raw - MAX_POSSIBLE_CHROM];
        }
        return -2;
    }
    if (!chrom_info_ptr->name_ct) {
        return -1;
    }
    // 0xffffffffU gets casted to -1
    return (int32_t) nonstd_chrom_name_htable_find(
        chrom_name, (const char* const*) chrom_info_ptr->nonstd_names,
        chrom_info_ptr->nonstd_id_htable, name_slen);
}

uint32_t get_variant_chrom_fo_idx(const Chrom_info* chrom_info_ptr,
                                  uintptr_t variant_uidx)
{
    const uint32_t* variant_binsearch = chrom_info_ptr->chrom_fo_vidx_start;
    uint32_t chrom_fo_min = 0;
    uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
    while (chrom_ct - chrom_fo_min > 1) {
        const uint32_t chrom_fo_cur = (chrom_ct + chrom_fo_min) / 2;
        if (variant_binsearch[chrom_fo_cur] > variant_uidx) {
            chrom_ct = chrom_fo_cur;
        }
        else
        {
            chrom_fo_min = chrom_fo_cur;
        }
    }
    return chrom_fo_min;
}

void chrom_error(const char* chrom_name, const char* file_descrip,
                 const Chrom_info* chrom_info_ptr, uintptr_t line_idx,
                 int32_t error_code)
{
    // assumes chrom_name is null-terminated
    const int32_t raw_code = get_chrom_code_raw(chrom_name);
    logprint("\n");
    if (line_idx) {
        LOGERRPRINTFWW("Error: Invalid chromosome code '%s' on line %" PRIuPTR
                       " of %s.\n",
                       chrom_name, line_idx, file_descrip);
    }
    else
    {
        LOGERRPRINTFWW("Error: Invalid chromosome code '%s' in %s.\n",
                       chrom_name, file_descrip);
    }
    if ((raw_code > ((int32_t) chrom_info_ptr->max_code))
        && ((raw_code <= MAX_CHROM_TEXTNUM + XYMT_OFFSET_CT)
            || (raw_code >= MAX_POSSIBLE_CHROM)))
    {
        if (chrom_info_ptr->species != SPECIES_UNKNOWN) {
            if (chrom_info_ptr->species == SPECIES_HUMAN) {
                logerrprint("(This is disallowed for humans.  Check if the "
                            "problem is with your data, or if\nyou forgot to "
                            "define a different chromosome set with e.g. "
                            "--chr-set.).\n");
            }
            else
            {
                logerrprint("(This is disallowed by the PLINK 1.07 species "
                            "flag you used.  You can\ntemporarily work around "
                            "this restriction with --chr-set; contact the "
                            "developers\nif you want the flag to be "
                            "permanently redefined.)\n");
            }
        }
        else
        {
            logerrprint("(This is disallowed by your --chr-set/--autosome-num "
                        "parameters.  Check if the\nproblem is with your data, "
                        "or your command line.)\n");
        }
    }
    else if (error_code == -1)
    {
        logerrprint("(Use --allow-extra-chr to force it to be accepted.)\n");
    }
}

int32_t try_to_add_chrom_name(const char* chrom_name, const char* file_descrip,
                              uintptr_t line_idx, uint32_t name_slen,
                              uint32_t allow_extra_chroms,
                              int32_t* chrom_idx_ptr,
                              Chrom_info* chrom_info_ptr)
{
    // assumes chrom_name is nonstandard (i.e. not "2", "chr2", "chrX", etc.)
    // requires chrom_name to be null-terminated
    // assumes chrom_idx currently has the return value of get_chrom_code()
    if ((!allow_extra_chroms) || ((*chrom_idx_ptr) == -2)) {
        chrom_error(chrom_name, file_descrip, chrom_info_ptr, line_idx,
                    *chrom_idx_ptr);
        return RET_MALFORMED_INPUT;
    }

    // quasi-bugfix: remove redundant hash table check

    if (chrom_name[0] == '#') {
        // redundant with some of the comment-skipping loaders, but this isn't
        // performance-critical
        logprint("\n");
        logerrprint("Error: Chromosome/contig names may not begin with '#'.\n");
        return RET_MALFORMED_INPUT;
    }
    if (name_slen > MAX_ID_SLEN) {
        logprint("\n");
        if (line_idx) {
            LOGERRPRINTFWW("Error: Line %" PRIuPTR " of %s has an excessively "
                           "long chromosome/contig "
                           "name. (The " PROG_NAME_CAPS
                           " limit is " MAX_ID_SLEN_STR " characters.)\n",
                           line_idx, file_descrip);
        }
        else
        {
            LOGERRPRINTFWW("Error: Excessively long chromosome/contig name in "
                           "%s. (The " PROG_NAME_CAPS
                           " limit is " MAX_ID_SLEN_STR " characters.)\n",
                           file_descrip);
        }
        return RET_MALFORMED_INPUT;
    }
    const uint32_t max_code_p1 = chrom_info_ptr->max_code + 1;
    const uint32_t name_ct = chrom_info_ptr->name_ct;
    const uint32_t chrom_code_end = max_code_p1 + name_ct;
    if (chrom_code_end == MAX_POSSIBLE_CHROM) {
        logprint("\n");
        logerrprint(
            "Error: Too many distinct nonstandard chromosome/contig names.\n");
        return RET_MALFORMED_INPUT;
    }
    if (!name_ct) {
        // lazy initialization
        fill_uint_one(CHROM_NAME_HTABLE_SIZE, chrom_info_ptr->nonstd_id_htable);
    }
    char** nonstd_names = chrom_info_ptr->nonstd_names;
    nonstd_names[chrom_code_end] = (char*) malloc(name_slen + 1);
    if (!nonstd_names[chrom_code_end]) {
        return RET_NOMEM;
    }
    Ll_str* name_stack_ptr = chrom_info_ptr->incl_excl_name_stack;
    uint32_t in_name_stack = 0;
    while (name_stack_ptr) {
        // there shouldn't be many of these, so sorting is unimportant
        if (!strcmp(chrom_name, name_stack_ptr->ss)) {
            in_name_stack = 1;
            break;
        }
        name_stack_ptr = name_stack_ptr->next;
    }
    if ((in_name_stack && chrom_info_ptr->is_include_stack)
        || ((!in_name_stack) && (!chrom_info_ptr->is_include_stack)))
    {
        SET_BIT(chrom_code_end, chrom_info_ptr->chrom_mask);
        if (chrom_info_ptr->haploid_mask[0] & 1) {
            SET_BIT(chrom_code_end, chrom_info_ptr->haploid_mask);
        }
    }
    memcpy(nonstd_names[chrom_code_end], chrom_name, name_slen + 1);
    *chrom_idx_ptr = (int32_t) chrom_code_end;
    chrom_info_ptr->name_ct = name_ct + 1;
    uint32_t* id_htable = chrom_info_ptr->nonstd_id_htable;
    uint32_t hashval =
        murmurhash3_32(chrom_name, name_slen) % CHROM_NAME_HTABLE_SIZE;
    uint32_t next_incr = 1;
    while (1) {
        if (id_htable[hashval] == 0xffffffffU) {
            id_htable[hashval] = chrom_code_end;
            return 0;
        }
        // no overflow danger here
        hashval += next_incr;
        if (hashval >= CHROM_NAME_HTABLE_SIZE) {
            hashval -= CHROM_NAME_HTABLE_SIZE;
        }
        next_incr += 2; // quadratic probing
    }
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


int32_t bsearch_str(const char* id_buf, uintptr_t cur_id_len, const char* lptr,
                    uintptr_t max_id_len, uintptr_t end_idx)
{
    // does not assume null-terminated id_buf, or nonempty array.
    // N.B. max_id_len includes null terminator as usual, while cur_id_len does
    // NOT.
    uintptr_t start_idx = 0;
    uintptr_t mid_idx;
    int32_t ii;
    if (cur_id_len >= max_id_len) {
        return -1;
    }
    while (start_idx < end_idx) {
        mid_idx = (start_idx + end_idx) / 2;
        ii = memcmp(id_buf, &(lptr[mid_idx * max_id_len]), cur_id_len);
        if (ii > 0) {
            start_idx = mid_idx + 1;
        }
        else if ((ii < 0) || lptr[mid_idx * max_id_len + cur_id_len])
        {
            end_idx = mid_idx;
        }
        else
        {
            return ((uint32_t) mid_idx);
        }
    }
    return -1;
}


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

static inline uintptr_t
popcount_vecs_exclude(const __m128i* __restrict vptr,
                      const __m128i* __restrict exclude_ptr, uintptr_t ct)
{
    // popcounts vptr ANDNOT exclude_ptr[0..(ct-1)].  ct is a multiple of 3.
    const __m128i m1 = {FIVEMASK, FIVEMASK};
    const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
    const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
    const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
    uintptr_t tot = 0;
    const __m128i* vend;
    __m128i count1, count2, half1, half2;
    __univec acc;

    while (ct >= 30) {
        ct -= 30;
        vend = &(vptr[30]);
    popcount_vecs_exclude_main_loop:
        acc.vi = _mm_setzero_si128();
        do
        {
            // nots the FIRST value
            count1 = _mm_andnot_si128(*exclude_ptr++, *vptr++);
            count2 = _mm_andnot_si128(*exclude_ptr++, *vptr++);
            half1 = _mm_andnot_si128(*exclude_ptr++, *vptr++);
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
        } while (vptr < vend);
        acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8),
                               _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
        tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
    }
    if (ct) {
        vend = &(vptr[ct]);
        ct = 0;
        goto popcount_vecs_exclude_main_loop;
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

uintptr_t popcount_bit_idx(const uintptr_t* lptr, uintptr_t start_idx,
                           uintptr_t end_idx)
{
    uintptr_t start_idxl = start_idx / BITCT;
    uintptr_t start_idxlr = start_idx & (BITCT - 1);
    uintptr_t end_idxl = end_idx / BITCT;
    uintptr_t end_idxlr = end_idx & (BITCT - 1);
    uintptr_t ct = 0;
    if (start_idxl == end_idxl) {
        return popcount_long(lptr[start_idxl]
                             & ((ONELU << end_idxlr) - (ONELU << start_idxlr)));
    }
    if (start_idxlr) {
        ct = popcount_long(lptr[start_idxl++] >> start_idxlr);
    }
    if (end_idxl > start_idxl) {
        ct += popcount_longs_nzbase(lptr, start_idxl, end_idxl);
    }
    if (end_idxlr) {
        ct += popcount_long(lptr[end_idxl] & ((ONELU << end_idxlr) - ONELU));
    }
    return ct;
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
#endif

#ifdef __LP64__
void count_set_freq_60v(const __m128i* vptr, const __m128i* vend,
                        const __m128i* __restrict include_vec,
                        uint32_t* __restrict set_ctp,
                        uint32_t* __restrict missing_ctp)
{
    const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
    const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
    const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
    __m128i loader;
    __m128i loader2;
    __m128i loader3;
    __m128i odds;
    __m128i evens;
    __m128i missings;
    __univec acc;
    __univec accm;
    acc.vi = _mm_setzero_si128();
    accm.vi = _mm_setzero_si128();
    do
    {
        loader = *vptr++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader3 = *include_vec++;
        odds = _mm_and_si128(loader2, loader3);
        evens = _mm_and_si128(odds, loader);
        missings = _mm_and_si128(loader, _mm_andnot_si128(loader2, loader3));

        loader = *vptr++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader3 = *include_vec++;
        odds = _mm_add_epi64(odds, _mm_and_si128(loader2, loader3));
        loader3 = _mm_and_si128(loader, loader3);
        evens = _mm_add_epi64(evens, _mm_and_si128(loader2, loader3));
        missings = _mm_add_epi64(missings, _mm_andnot_si128(loader2, loader3));

        loader = *vptr++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader3 = *include_vec++;
        odds = _mm_add_epi64(odds, _mm_and_si128(loader2, loader3));
        loader3 = _mm_and_si128(loader, loader3);
        evens = _mm_add_epi64(evens, _mm_and_si128(loader2, loader3));
        missings = _mm_add_epi64(missings, _mm_andnot_si128(loader2, loader3));

        odds = _mm_add_epi64(_mm_and_si128(odds, m2),
                             _mm_and_si128(_mm_srli_epi64(odds, 2), m2));
        missings =
            _mm_add_epi64(_mm_and_si128(missings, m2),
                          _mm_and_si128(_mm_srli_epi64(missings, 2), m2));
        odds = _mm_add_epi64(
            odds, _mm_add_epi64(_mm_and_si128(evens, m2),
                                _mm_and_si128(_mm_srli_epi64(evens, 2), m2)));

        // each 4-bit value here <= 6, so safe to add before m4 mask
        accm.vi = _mm_add_epi64(
            accm.vi,
            _mm_and_si128(_mm_add_epi64(missings, _mm_srli_epi64(missings, 4)),
                          m4));

        acc.vi = _mm_add_epi64(
            acc.vi, _mm_add_epi64(_mm_and_si128(odds, m4),
                                  _mm_and_si128(_mm_srli_epi64(odds, 4), m4)));
    } while (vptr < vend);
    // and each 8-bit value here <= 120
    accm.vi =
        _mm_and_si128(_mm_add_epi64(accm.vi, _mm_srli_epi64(accm.vi, 8)), m8);

    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8),
                           _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
    *set_ctp += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
    *missing_ctp += ((accm.u8[0] + accm.u8[1]) * 0x1000100010001LLU) >> 48;
}

void count_set_freq_hap_120v(const __m128i* vptr, const __m128i* vend,
                             const __m128i* __restrict include_vec,
                             uint32_t* __restrict set_ctp,
                             uint32_t* __restrict missing_ctp)
{
    const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
    const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
    const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
    __univec acc;
    __univec accm;
    __m128i loader;
    __m128i loader2;
    __m128i loader3;
    __m128i partial;
    __m128i partialm;
    __m128i partial2;
    __m128i partial2m;
    acc.vi = _mm_setzero_si128();
    accm.vi = _mm_setzero_si128();
    do
    {
        loader = *vptr++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader3 = *include_vec++;
        partial = _mm_and_si128(loader3, _mm_and_si128(loader, loader2));
        partialm = _mm_and_si128(loader3, _mm_xor_si128(loader, loader2));
        loader = *vptr++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader3 = *include_vec++;
        partial = _mm_add_epi64(
            partial, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
        partialm = _mm_add_epi64(
            partialm, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
        loader = *vptr++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader3 = *include_vec++;
        partial = _mm_add_epi64(
            partial, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
        partialm = _mm_add_epi64(
            partialm, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
        partial2 = _mm_add_epi64(_mm_and_si128(partial, m2),
                                 _mm_and_si128(_mm_srli_epi64(partial, 2), m2));
        partial2m =
            _mm_add_epi64(_mm_and_si128(partialm, m2),
                          _mm_and_si128(_mm_srli_epi64(partialm, 2), m2));

        loader = *vptr++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader3 = *include_vec++;
        partial = _mm_and_si128(loader3, _mm_and_si128(loader, loader2));
        partialm = _mm_and_si128(loader3, _mm_xor_si128(loader, loader2));
        loader = *vptr++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader3 = *include_vec++;
        partial = _mm_add_epi64(
            partial, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
        partialm = _mm_add_epi64(
            partialm, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
        loader = *vptr++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader3 = *include_vec++;
        partial = _mm_add_epi64(
            partial, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
        partialm = _mm_add_epi64(
            partialm, _mm_and_si128(loader3, _mm_and_si128(loader, loader2)));
        partial2 = _mm_add_epi64(
            partial2,
            _mm_add_epi64(_mm_and_si128(partial, m2),
                          _mm_and_si128(_mm_srli_epi64(partial, 2), m2)));
        partial2m = _mm_add_epi64(
            partial2m,
            _mm_add_epi64(_mm_and_si128(partialm, m2),
                          _mm_and_si128(_mm_srli_epi64(partialm, 2), m2)));
        acc.vi = _mm_add_epi64(
            acc.vi,
            _mm_add_epi64(_mm_and_si128(partial2, m4),
                          _mm_and_si128(_mm_srli_epi64(partial2, 4), m4)));
        accm.vi = _mm_add_epi64(
            accm.vi,
            _mm_add_epi64(_mm_and_si128(partial2m, m4),
                          _mm_and_si128(_mm_srli_epi64(partial2m, 4), m4)));
    } while (vptr < vend);
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8),
                           _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
    accm.vi = _mm_add_epi64(_mm_and_si128(accm.vi, m8),
                            _mm_and_si128(_mm_srli_epi64(accm.vi, 8), m8));
    *set_ctp += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
    *missing_ctp += ((accm.u8[0] + accm.u8[1]) * 0x1000100010001LLU) >> 48;
}

void count_set_freq_x_60v(const __m128i* vptr, const __m128i* vend,
                          const __m128i* __restrict include_vec,
                          const __m128i* __restrict male_vec,
                          uint32_t* __restrict set_ctp,
                          uint32_t* __restrict missing_ctp)
{
    const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
    const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
    const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
    __m128i loader;
    __m128i loader2;
    __m128i loader3;
    __m128i loader4;
    __m128i set_odds;
    __m128i set_evens;
    __m128i missings_nm;
    __m128i missings_m;
    __m128i males;
    __univec acc;
    __univec accm;
    acc.vi = _mm_setzero_si128();
    accm.vi = _mm_setzero_si128();
    do
    {
        loader = *vptr++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader3 = *include_vec++;
        loader4 = _mm_andnot_si128(*male_vec, loader3);
        set_evens =
            _mm_and_si128(loader, loader4); // subtract missings_nm later
        set_odds = _mm_and_si128(loader2, loader4);
        missings_nm = _mm_andnot_si128(loader2, set_evens);
        males = _mm_and_si128(loader3, *male_vec++);
        set_evens = _mm_or_si128(
            set_evens, _mm_and_si128(_mm_and_si128(loader, loader2), males));
        missings_m = _mm_and_si128(_mm_xor_si128(loader, loader2), males);

        loader = *vptr++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader3 = *include_vec++;
        loader4 = _mm_andnot_si128(*male_vec, loader3);
        set_odds = _mm_add_epi64(set_odds, _mm_and_si128(loader2, loader4));
        loader4 = _mm_and_si128(loader, loader4);
        set_evens = _mm_add_epi64(set_evens, loader4);
        missings_nm =
            _mm_add_epi64(missings_nm, _mm_andnot_si128(loader2, loader4));
        loader4 = _mm_and_si128(loader3, *male_vec++);
        set_evens = _mm_add_epi64(
            set_evens, _mm_and_si128(_mm_and_si128(loader, loader2), loader4));
        missings_m = _mm_add_epi64(
            missings_m, _mm_and_si128(_mm_xor_si128(loader, loader2), loader4));
        males = _mm_add_epi64(males, loader4);

        loader = *vptr++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader3 = *include_vec++;
        loader4 = _mm_andnot_si128(*male_vec, loader3);
        set_odds = _mm_add_epi64(set_odds, _mm_and_si128(loader2, loader4));
        loader4 = _mm_and_si128(loader, loader4);
        set_evens = _mm_add_epi64(set_evens, loader4);
        missings_nm =
            _mm_add_epi64(missings_nm, _mm_andnot_si128(loader2, loader4));
        loader4 = _mm_and_si128(loader3, *male_vec++);
        set_evens = _mm_add_epi64(
            set_evens, _mm_and_si128(_mm_and_si128(loader, loader2), loader4));
        missings_m = _mm_add_epi64(
            missings_m, _mm_and_si128(_mm_xor_si128(loader, loader2), loader4));
        males = _mm_add_epi64(males, loader4);

        set_evens = _mm_sub_epi64(set_evens, missings_nm);
        missings_nm = _mm_slli_epi64(
            _mm_add_epi64(_mm_and_si128(missings_nm, m2),
                          _mm_and_si128(_mm_srli_epi64(missings_nm, 2), m2)),
            1);
        set_odds =
            _mm_add_epi64(_mm_and_si128(set_odds, m2),
                          _mm_and_si128(_mm_srli_epi64(set_odds, 2), m2));
        missings_nm = _mm_add_epi64(
            missings_nm,
            _mm_add_epi64(_mm_and_si128(missings_m, m2),
                          _mm_and_si128(_mm_srli_epi64(missings_m, 2), m2)));
        set_odds = _mm_add_epi64(
            set_odds,
            _mm_add_epi64(_mm_and_si128(set_evens, m2),
                          _mm_and_si128(_mm_srli_epi64(set_evens, 2), m2)));
        missings_nm = _mm_add_epi64(
            missings_nm,
            _mm_add_epi64(_mm_and_si128(males, m2),
                          _mm_and_si128(_mm_srli_epi64(males, 2), m2)));
        acc.vi = _mm_add_epi64(
            acc.vi,
            _mm_add_epi64(_mm_and_si128(set_odds, m4),
                          _mm_and_si128(_mm_srli_epi64(set_odds, 4), m4)));
        accm.vi = _mm_add_epi64(
            accm.vi,
            _mm_add_epi64(_mm_and_si128(missings_nm, m4),
                          _mm_and_si128(_mm_srli_epi64(missings_nm, 4), m4)));
    } while (vptr < vend);
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8),
                           _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
    accm.vi = _mm_add_epi64(_mm_and_si128(accm.vi, m8),
                            _mm_and_si128(_mm_srli_epi64(accm.vi, 8), m8));
    *set_ctp += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
    *missing_ctp += ((accm.u8[0] + accm.u8[1]) * 0x1000100010001LLU) >> 48;
}

void count_set_freq_y_120v(const __m128i* vptr, const __m128i* vend,
                           const __m128i* __restrict include_vec,
                           const __m128i* __restrict nonmale_vec,
                           uint32_t* __restrict set_ctp,
                           uint32_t* __restrict missing_ctp)
{
    const __m128i m2 = {0x3333333333333333LLU, 0x3333333333333333LLU};
    const __m128i m4 = {0x0f0f0f0f0f0f0f0fLLU, 0x0f0f0f0f0f0f0f0fLLU};
    const __m128i m8 = {0x00ff00ff00ff00ffLLU, 0x00ff00ff00ff00ffLLU};
    __m128i loader;
    __m128i loader2;
    __m128i loader3;
    __m128i loader4;
    __m128i sets1;
    __m128i missings1;
    __m128i sets2;
    __m128i missings2;
    __univec acc;
    __univec accm;
    acc.vi = _mm_setzero_si128();
    accm.vi = _mm_setzero_si128();
    do
    {
        loader = *vptr++;
        loader3 = *include_vec++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader4 = *nonmale_vec++;
        sets1 = _mm_and_si128(_mm_andnot_si128(loader4, loader3),
                              _mm_and_si128(loader, loader2));
        missings1 = _mm_and_si128(
            loader3, _mm_or_si128(loader4, _mm_xor_si128(loader, loader2)));

        loader = *vptr++;
        loader3 = *include_vec++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader4 = *nonmale_vec++;
        sets1 = _mm_add_epi64(sets1,
                              _mm_and_si128(_mm_andnot_si128(loader4, loader3),
                                            _mm_and_si128(loader, loader2)));
        missings1 = _mm_add_epi64(
            missings1,
            _mm_and_si128(
                loader3,
                _mm_or_si128(loader4, _mm_xor_si128(loader, loader2))));

        loader = *vptr++;
        loader3 = *include_vec++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader4 = *nonmale_vec++;
        sets1 = _mm_add_epi64(sets1,
                              _mm_and_si128(_mm_andnot_si128(loader4, loader3),
                                            _mm_and_si128(loader, loader2)));
        missings1 = _mm_add_epi64(
            missings1,
            _mm_and_si128(
                loader3,
                _mm_or_si128(loader4, _mm_xor_si128(loader, loader2))));
        sets1 = _mm_add_epi64(_mm_and_si128(sets1, m2),
                              _mm_and_si128(_mm_srli_epi64(sets1, 2), m2));
        missings1 =
            _mm_add_epi64(_mm_and_si128(missings1, m2),
                          _mm_and_si128(_mm_srli_epi64(missings1, 2), m2));

        loader = *vptr++;
        loader3 = *include_vec++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader4 = *nonmale_vec++;
        sets2 = _mm_and_si128(_mm_andnot_si128(loader4, loader3),
                              _mm_and_si128(loader, loader2));
        missings2 = _mm_and_si128(
            loader3, _mm_or_si128(loader4, _mm_xor_si128(loader, loader2)));

        loader = *vptr++;
        loader3 = *include_vec++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader4 = *nonmale_vec++;
        sets2 = _mm_add_epi64(sets2,
                              _mm_and_si128(_mm_andnot_si128(loader4, loader3),
                                            _mm_and_si128(loader, loader2)));
        missings2 = _mm_add_epi64(
            missings2,
            _mm_and_si128(
                loader3,
                _mm_or_si128(loader4, _mm_xor_si128(loader, loader2))));

        loader = *vptr++;
        loader3 = *include_vec++;
        loader2 = _mm_srli_epi64(loader, 1);
        loader4 = *nonmale_vec++;
        sets2 = _mm_add_epi64(sets2,
                              _mm_and_si128(_mm_andnot_si128(loader4, loader3),
                                            _mm_and_si128(loader, loader2)));
        missings2 = _mm_add_epi64(
            missings2,
            _mm_and_si128(
                loader3,
                _mm_or_si128(loader4, _mm_xor_si128(loader, loader2))));
        sets1 = _mm_add_epi64(
            sets1, _mm_add_epi64(_mm_and_si128(sets2, m2),
                                 _mm_and_si128(_mm_srli_epi64(sets2, 2), m2)));
        missings1 = _mm_add_epi64(
            missings1,
            _mm_add_epi64(_mm_and_si128(missings2, m2),
                          _mm_and_si128(_mm_srli_epi64(missings2, 2), m2)));
        acc.vi = _mm_add_epi64(
            acc.vi, _mm_add_epi64(_mm_and_si128(sets1, m4),
                                  _mm_and_si128(_mm_srli_epi64(sets1, 4), m4)));
        accm.vi = _mm_add_epi64(
            accm.vi,
            _mm_add_epi64(_mm_and_si128(missings1, m4),
                          _mm_and_si128(_mm_srli_epi64(missings1, 4), m4)));
    } while (vptr < vend);
    acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8),
                           _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
    accm.vi = _mm_add_epi64(_mm_and_si128(accm.vi, m8),
                            _mm_and_si128(_mm_srli_epi64(accm.vi, 8), m8));
    *set_ctp += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
    *missing_ctp += ((accm.u8[0] + accm.u8[1]) * 0x1000100010001LLU) >> 48;
}

uintptr_t count_01_vecs(const __m128i* vptr, uintptr_t vct)
{
    // counts number of aligned 01s (i.e. PLINK missing genotypes) in
    // [vptr, vend).  Assumes number of words in interval is a multiple of 12.
    const __m128i m1 = {FIVEMASK, FIVEMASK};
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

    while (vct >= 60) {
        vct -= 60;
        vend = &(vptr[60]);
    count_01_vecs_main_loop:
        acc.vi = _mm_setzero_si128();
        do
        {
            loader1 = *vptr++;
            loader2 = *vptr++;
            count1 = _mm_and_si128(
                _mm_andnot_si128(_mm_srli_epi64(loader1, 1), loader1), m1);
            count2 = _mm_and_si128(
                _mm_andnot_si128(_mm_srli_epi64(loader2, 1), loader2), m1);
            loader1 = *vptr++;
            loader2 = *vptr++;
            count1 = _mm_add_epi64(
                count1,
                _mm_and_si128(
                    _mm_andnot_si128(_mm_srli_epi64(loader1, 1), loader1), m1));
            count2 = _mm_add_epi64(
                count2,
                _mm_and_si128(
                    _mm_andnot_si128(_mm_srli_epi64(loader2, 1), loader2), m1));
            loader1 = *vptr++;
            loader2 = *vptr++;
            count1 = _mm_add_epi64(
                count1,
                _mm_and_si128(
                    _mm_andnot_si128(_mm_srli_epi64(loader1, 1), loader1), m1));
            count2 = _mm_add_epi64(
                count2,
                _mm_and_si128(
                    _mm_andnot_si128(_mm_srli_epi64(loader2, 1), loader2), m1));
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
        } while (vptr < vend);
        acc.vi = _mm_add_epi64(_mm_and_si128(acc.vi, m8),
                               _mm_and_si128(_mm_srli_epi64(acc.vi, 8), m8));
        tot += ((acc.u8[0] + acc.u8[1]) * 0x1000100010001LLU) >> 48;
    }
    if (vct) {
        vend = &(vptr[vct]);
        vct = 0;
        goto count_01_vecs_main_loop;
    }
    return tot;
}

#else

#endif

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

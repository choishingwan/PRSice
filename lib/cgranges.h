/* The MIT License

   Copyright (c) 2019 Dana-Farber Cancer Institute

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/
#ifndef CRANGES_H
#define CRANGES_H

#include <stdint.h>

typedef struct {    // a contig
	char *name;     // name of the contig
	int32_t len;    // max length seen in data
	int32_t root_k;
	int64_t n, off; // sum of lengths of previous contigs
} cr_ctg_t;

typedef struct {    // an interval
	uint64_t x;     // prior to cr_index(), x = ctg_id<<32|start_pos; after: x = start_pos<<32|end_pos
	uint32_t y:31, rev:1;
	int32_t label;  // NOT used
} cr_intv_t;

typedef struct {
	int64_t n_r, m_r;     // number and max number of intervals
	cr_intv_t *r;         // list of intervals (of size _n_r_)
	int32_t n_ctg, m_ctg; // number and max number of contigs
	cr_ctg_t *ctg;        // list of contigs (of size _n_ctg_)
	void *hc;             // dictionary for converting contig names to integers
} cgranges_t;

#ifdef __cplusplus
extern "C" {
#endif

// retrieve start and end positions from a cr_intv_t object
static inline int32_t cr_st(const cr_intv_t *r) { return (int32_t)(r->x>>32); }
static inline int32_t cr_en(const cr_intv_t *r) { return (int32_t)r->x; }
static inline int32_t cr_start(const cgranges_t *cr, int64_t i) { return cr_st(&cr->r[i]); }
static inline int32_t cr_end(const cgranges_t *cr, int64_t i) { return cr_en(&cr->r[i]); }
static inline int32_t cr_label(const cgranges_t *cr, int64_t i) { return cr->r[i].label; }

// Initialize
cgranges_t *cr_init(void);

// Deallocate
void cr_destroy(cgranges_t *cr);

// Add an interval
cr_intv_t *cr_add(cgranges_t *cr, const char *ctg, int32_t st, int32_t en, int32_t label_int);

// Sort and index intervals
void cr_index(cgranges_t *cr);

int64_t cr_overlap(const cgranges_t *cr, const char *ctg, int32_t st, int32_t en, int64_t **b_, int64_t *m_b_);

// Add a contig and length. Call this for desired contig ordering. _len_ can be 0.
int32_t cr_add_ctg(cgranges_t *cr, const char *ctg, int32_t len);

// Get the contig ID given its name
int32_t cr_get_ctg(const cgranges_t *cr, const char *ctg);

#ifdef __cplusplus
}
#endif

#endif

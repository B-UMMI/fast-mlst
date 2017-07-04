/*-
 * Copyright (c) 2017, Alexandre P. Francisco <aplf@ist.utl.pt>
 * Copyright (c) 2017, Cátia vaz <cvaz@cc.isel.ipl.pt>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sautils.h"

static int
array_cmp(int32_t *u, int k, int32_t *v, int l)
{
    int i = 0;
    while (i < k && i < l && u[i] == v[i])
        i ++;
    if (i >= k || i >= l)
        return 0;
    return u[i] - v[i];
}

static int
sa_search_low(int *SA, int32_t *s, int n, int32_t *p, int m)
{
    int lo = 0, hi = n - 1, mid;
    while (lo <= hi) {
        mid = lo + (hi - lo) / 2;
        int cmp = array_cmp(p, m, s + SA[mid], m);
        if (cmp <= 0) hi = mid - 1;
        else lo = mid + 1;
    }
    return lo;
}

static int
sa_search_high(int *SA, int32_t *s, int n, int32_t *p, int m)
{
    int lo = 0, hi = n - 1, mid;
    while (lo <= hi) {
        mid = lo + (hi - lo) / 2;
        int cmp = array_cmp(p, m, s + SA[mid], m);
        if (cmp < 0) hi = mid - 1;
        else lo = mid + 1;
    }
    return lo;
}

static int
hamming_distance(int32_t *s, int32_t *r, int m, int k, int p)
{
    int x = 0, l, b = m/k;
    /* Extend right */
    for (l = p + b; l < m && x < k; l++)
        if (s[l] != r[l]) x++;
    /* Extend left */
    for (l = p - 1; l >= 0 && x < k; l--)
        if (s[l] != r[l]) x++;
    return x;
}

int
solve_query(int32_t *s, int32_t *sa, int32_t *q, int d, int m, int k,
    int32_pair_t *rv)
{
    int j, ltk, nv;

    int *filter = malloc(sizeof(int)*d);
    memset(filter, 0xff, sizeof(int)*d);

    ltk = nv = 0;
    for (int ik = 0; ik < m+1 - m/k; ik += m/k) {
        int r = sa_search_low(sa, s, d*(m+1), q + ik, m/k);
        int t = sa_search_high(sa, s, d*(m+1), q + ik, m/k);
        for(; r < t && nv < d; r++) {
            j = sa[r]/(m+1);
            if (sa[r]%(m+1) == ik && filter[j] != 1) {
                filter[j] = 1;
                nv ++;
                int x = hamming_distance(q, s + j*(m+1), m, k, ik);
                if (x < k) {
                    rv[ltk].id = j;
                    rv[ltk].n = x;
                    ltk++;
                }
            }
        }
    }
    free(filter);
    fprintf(stderr, "#hits: %d (%d)\n", ltk, nv);
    return ltk;
}

int
int32_pair_cmp(const void *p, const void *q)
{
    int32_pair_t * ip = (int32_pair_t *) p;
    int32_pair_t * iq = (int32_pair_t *) q;

    if (ip->n != iq->n)
        return ip->n - iq->n;

    return ip->id - iq->id;
}

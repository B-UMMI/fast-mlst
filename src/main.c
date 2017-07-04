/*-
 * Copyright (c) 2017, Alexandre P. Francisco <aplf@ist.utl.pt>
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

/* Fast querying MLST data. Based on ideas discussed in:
 * M Crochemore, AP Francisco, SP Pissis, and C Vaz: Towards distance-based
 * phylogenetic inference in average-case linear-time. In WABI'2017.
 */

#define _POSIX_C_SOURCE 2

#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <assert.h>

#include "main.h"
#include "sautils.h"

#define cpuTime() (clock()*1e-6)

static int32_t *lidx;
static int32_t *profiles;
static int32_t *sa;
static int32_t n, n_al, n_ST;

int
main(int argc, char * argv[])
{
    char iname[132] = { 0 }, lname[132] = { 0 }, mode = 'q', *lblock;
    int32_t opt = -1, k = -1, sigma = -1, *isa = NULL, wn = 0, fd, lfd, *mblock,
        *q, nr, i;
    int32_pair_t *r;
    struct stat sb;
    FILE *fptr = NULL, *lptr = NULL;

    /* Command line options :
     *
     *  i - index name
     *  b - build index
     *  q - query index
     *
     * Option 'i' requires an argument, a string. Option 'q' requires also an
     * argument, the maximum error allowed. Both the query and profiles to index
     * should be provided through stdin.
     */
    while ((opt = getopt(argc, argv, "i:q:b")) != -1) {
        switch (opt) {
        case 'i':
            strncpy(iname, optarg, 127);
            strncpy(lname, optarg, 127);
            break;
        case 'b':
            mode = 'b'; 
            break;
        case 'q':
            mode = 'q';
            k = atoi(optarg);
            break;
        default: /* 'h' and invalid options.  */
            usage(argv[0]);
            exit(EXIT_FAILURE);
        }
    }
            
    if (optind > argc || iname[0] == 0) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    strcat(iname, ".idx");
    strcat(lname, ".ids");

    if (mode == 'q') {
        fd = open(iname, O_RDONLY);
        fstat(fd, &sb);
        mblock = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
        if (mblock == MAP_FAILED) {
            perror("Error loading index");
            return EXIT_FAILURE;
        }

        n_ST = mblock[0];
        n_al = mblock[1];
        n = n_ST * (n_al + 1);
        profiles = mblock + 2;
        sa = profiles + (n + 1);
        lidx = sa + (n + 1);

        lfd = open(lname, O_RDONLY);
        fstat(lfd, &sb);
        lblock = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, lfd, 0);
        if (lblock == MAP_FAILED) {
            perror("Error loading index");
            return EXIT_FAILURE;
        }

        /* Read query */
        q = malloc(sizeof(int32_t)*(n_al + 1));
        wn = read_query(stdin, q, n_al);
        if (wn != n_al) {
            fprintf(stderr, "ERROR while loading query, giving up...\n");
            return EXIT_FAILURE;
        }

        r = malloc(sizeof(int32_pair_t)*n_ST);
        nr = solve_query(profiles, sa, q, n_ST, n_al, k+1, r);
        qsort(r, nr, sizeof(int32_pair_t), int32_pair_cmp);

        for (i = 0; i < nr; i++)
            printf("%s\t%d\n", lblock + lidx[r[i].id], r[i].n);

        free(r);
        free(q);
        close(fd);
        close(lfd);
        return EXIT_SUCCESS;
    }

    /* Load data... */
    fprintf(stderr, "[%f] Loading data...\n", cpuTime());

    lptr = fopen(lname, "w");
    sigma = load_STs(stdin, lptr);
    fclose(lptr);

    if (sigma < 0) {
        fprintf(stderr, "ERROR while loading data, giving up...\n");
        return EXIT_FAILURE;
    }

    /* Build/update index... */
    fprintf(stderr, "[%f] Constructing SA and ISA...\n", cpuTime());
    n = n_ST * (n_al + 1);
    sa = malloc(sizeof(int32_t)*(n+1));
    memset(sa, 0xff, sizeof(int32_t)*(n+1));
    isa = malloc(sizeof(int32_t)*(n+1));
    memset(isa, 0xff, sizeof(int32_t)*(n+1));

    memcpy(isa, profiles, sizeof(int32_t)*n);
    suffixsort(isa, sa, n, sigma+1, -1);
    /* We do not need the ISA! */
    free(isa);

    fprintf(stderr, "[%f] Writing index...\n", cpuTime());
    fptr = fopen(iname,"wb");
    if (fptr == NULL) {
        perror("Error opening index");
        return EXIT_FAILURE;
    }
    wn = fwrite(&n_ST, sizeof(n_ST), 1, fptr);
    wn += fwrite(&n_al, sizeof(n_al), 1, fptr);
    wn += fwrite(profiles, sizeof(int32_t), (n + 1), fptr);
    wn += fwrite(sa, sizeof(int32_t), (n + 1), fptr);
    wn += fwrite(lidx, sizeof(int32_t), n_ST, fptr);
    fclose(fptr);

    free(lidx);
    free(profiles);
    free(sa);

    if (wn != 2 + 2*(n+1) + n_ST) {
        fprintf(stderr,
            "An error occured while writing the index, exiting...\n");
        return EXIT_FAILURE;
    }

    fprintf(stderr, "[%f] done!\n", cpuTime());
    return EXIT_SUCCESS;
}

int
read_query(FILE * fd, int32_t *q, int32_t l) {
    char *buffer = NULL, *tok;
    int32_t bsize = 0, i;

    readline(fd, &buffer, &bsize);

    tok = strtok(buffer, "\t\n, ");

    if (tok == NULL)
            return -1;

    /* Let us get the ST_id for this line. */
    /* id = atoi(tok); unused */
 
    for (i = 0; i < l && (tok = strtok(NULL, "\t\n, ")) != NULL; i++)
        q[i] = atoi(tok) + 1;

    q[i] = 0;

    free(buffer);

    return i;
}

int
load_STs(FILE *fd, FILE *lfd)
{
    char *buffer = NULL;
    int32_t bsize = 0, max_ST = 0;
    int32_t i = 0, sigma = 0, pl = 0;
    char *tok;

    max_ST = n_ST = n_al = 0;

    while (readline(fd, &buffer, &bsize) != EOF) {
    
        /* If the number of alleles is unknown, then let us find it. */
        if (n_al == 0) {

            size_t tmpsz = strlen(buffer) + 1;
            char * tmpbf = malloc(sizeof(char)*tmpsz);
            strncpy(tmpbf, buffer, tmpsz);
            tok = strtok(tmpbf, "\t\n, ");

            /* Let us count tokens. */
            for (i = 0; (tok = strtok(NULL, "\t\n, ")) != NULL; i++);
            n_al = i;
            fprintf(stderr, "%d alleles\n", n_al);

            free(tmpbf);
        }

        tok = strtok(buffer, "\t\n, ");

        if (tok == NULL)
            return -1;
    
        /* Checking free space... */
        if (n_ST >= max_ST) {
            if (max_ST != 0) max_ST <<= 1;
            else max_ST = 1024;

            profiles = realloc(profiles, sizeof(int32_t)*(max_ST*(n_al + 1) + 1));
            lidx = realloc(lidx, sizeof(int32_t)*max_ST);
        }

        /* Let us get the ST_id for this line. */
        /* id = atoi(tok); unused */
        lidx[n_ST] = 0;
        if (n_ST > 0)
            lidx[n_ST] = lidx[n_ST-1] + pl;
        fprintf(lfd, "%s%c", tok, 0);
        pl = strlen(tok) + 1;

        /* Parse and store the profile. */
        for (i = 0; (tok = strtok(NULL, "\t\n, ")) != NULL; i++) {
            if (i >= n_al)
                return -1;
            int val = atoi(tok) + 1;
            if (val > sigma)
                sigma = val;
            profiles[n_ST*(n_al+1) + i] = val;
        }

        if (i != n_al)
            return -1;
                
        profiles[n_ST*(n_al+1) + i] = 0;

        n_ST ++;
    }

    profiles[n_ST*(n_al+1)] = -1;
    free(buffer);
    fprintf(stderr, "%d profiles\n", n_ST);

    return sigma;
}

int
readline(FILE *fd, char **bf, int *bz)
{
    int rt, c, i = 0;

    while ((c = getc(fd)) != '\n' && c != EOF) {
        if (i >= *bz) {
            if (*bz != 0) *bz <<= 1;
            else *bz = 1024; 
            *bf = realloc(*bf, sizeof(char)*(*bz + 1));
        }

        (*bf)[i++] = c;
    }

    rt = c;
    (*bf)[i] = '\0';

    return rt;
}

int
st_diff(int i, int j)
{
    int diff, k;

    for (diff = k = 0; k < n_al; k++)
        if (profiles[i*(n_al+1) + k] != profiles[j*(n_al+1) + k])
           diff++;

    return diff;
}

void 
usage(char * cmd)
{

    fprintf(stderr, "\nUsage: %s [OPTION]...\n", cmd);
    fprintf(stderr, "\n\
Indexes and queries a set of STs. The input must be a tab separated vector of\n\
integers, with the same length per line.\n\
\n\
Options:\n\
");
    fprintf(stderr, "  -i INAME   The name for the index (mandatory).\n");
    fprintf(stderr, "  -q K       List matches with at most K errors.\n");
    fprintf(stderr, "  -b         The index should be (re)built.\n");
    fprintf(stderr, "\n");

}


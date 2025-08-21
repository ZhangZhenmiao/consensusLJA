#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE 4096
#define MAX_ID   256
#define INIT_INTERVAL_CAP 16
#define HASH_SIZE 200003  // prime for hash table

typedef struct {
    int start, end;
    double ident_bases;
    int aln_len;
} Interval;

typedef struct RefGroup {
    char sseqid[MAX_ID];
    int qlen;
    Interval* intervals;
    int count, capacity;
    struct RefGroup* next;
} RefGroup;

typedef struct Query {
    char qseqid[MAX_ID];
    RefGroup* refs;       // linked list of references
    struct Query* next;   // hash collision
} Query;

static unsigned hash_str(const char* s) {
    unsigned h = 5381; int c;
    while ((c = *s++)) h = ((h << 5) + h) + c;
    return h;
}

// hash table of queries
static Query* queries[HASH_SIZE];

static Query* get_query(const char* q) {
    unsigned idx = hash_str(q) % HASH_SIZE;
    Query* qptr = queries[idx];
    while (qptr) {
        if (!strcmp(qptr->qseqid, q)) return qptr;
        qptr = qptr->next;
    }
    qptr = calloc(1, sizeof(Query));
    strncpy(qptr->qseqid, q, MAX_ID - 1);
    qptr->refs = NULL;
    qptr->next = queries[idx];
    queries[idx] = qptr;
    return qptr;
}

static RefGroup* get_ref(Query* q, const char* s, int qlen) {
    for (RefGroup* r = q->refs; r; r = r->next)
        if (!strcmp(r->sseqid, s)) return r;

    RefGroup* r = calloc(1, sizeof(RefGroup));
    strncpy(r->sseqid, s, MAX_ID - 1);
    r->qlen = qlen;
    r->intervals = NULL;
    r->count = r->capacity = 0;
    r->next = q->refs;
    q->refs = r;
    return r;
}

static void add_interval(RefGroup* r, int start, int end, double pident) {
    if (start > end) { int tmp = start; start = end; end = tmp; }
    int aln_len = end - start + 1;
    double ident_bases = aln_len * (pident / 100.0);
    if (r->count == r->capacity) {
        r->capacity = r->capacity ? r->capacity * 2 : INIT_INTERVAL_CAP;
        r->intervals = realloc(r->intervals, r->capacity * sizeof(Interval));
    }
    r->intervals[r->count].start = start;
    r->intervals[r->count].end = end;
    r->intervals[r->count].aln_len = aln_len;
    r->intervals[r->count].ident_bases = ident_bases;
    r->count++;
}

static int cmp_interval(const void* a, const void* b) {
    return ((Interval*)a)->start - ((Interval*)b)->start;
}

// safer merging: only add non-overlapping parts
static void compute_merged(RefGroup* r, int* merged_aln, double* merged_ident) {
    if (r->count == 0) {
        *merged_aln = 0;
        *merged_ident = 0.0;
        return;
    }

    qsort(r->intervals, r->count, sizeof(Interval), cmp_interval);

    int cur_s = r->intervals[0].start;
    int cur_e = r->intervals[0].end;
    double cur_ident = r->intervals[0].ident_bases;

    *merged_aln = 0;
    *merged_ident = 0.0;

    for (int i = 1; i < r->count; i++) {
        Interval* iv = &r->intervals[i];

        if (iv->start <= cur_e) {
            // overlap: only add the non-overlapping tail
            int nonoverlap_start = cur_e + 1;
            if (iv->end >= nonoverlap_start) {
                int add_len = iv->end - nonoverlap_start + 1;
                double ident_rate = iv->ident_bases / iv->aln_len;
                cur_ident += ident_rate * add_len;
                cur_e = iv->end;
            }
        }
        else {
            // finish previous merged block
            *merged_aln += cur_e - cur_s + 1;
            *merged_ident += cur_ident;

            // start new block
            cur_s = iv->start;
            cur_e = iv->end;
            cur_ident = iv->ident_bases;
        }
    }

    // add the last block
    *merged_aln += cur_e - cur_s + 1;
    *merged_ident += cur_ident;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s blast6.txt summary.tsv\n", argv[0]);
        return 1;
    }

    FILE* in = fopen(argv[1], "r");
    if (!in) { perror("input"); return 1; }
    FILE* out = fopen(argv[2], "w");
    if (!out) { perror("output"); return 1; }

    fprintf(out, "qseqid\tsseqid\tqlen\taligned_bases\tcoverage(%%)\tidentity(%%)\n");
    fprintf(out, "# Note: Best reference per query = largest coverage; ties broken by higher identity\n");

    char line[MAX_LINE];
    while (fgets(line, sizeof(line), in)) {
        char qseqid[MAX_ID], sseqid[MAX_ID];
        double pident; int length, qlen, qstart, qend, sstart, send;
        double evalue, bitscore;

        int n = sscanf(line, "%255s\t%255s\t%lf\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf",
            qseqid, sseqid, &pident, &length, &qlen, &qstart, &qend, &sstart, &send, &evalue, &bitscore);
        if (n < 11) continue;
        if (pident < 70.0 || evalue > 1e-3) continue;  // skip low-quality hits

        Query* q = get_query(qseqid);
        RefGroup* r = get_ref(q, sseqid, qlen);
        add_interval(r, qstart, qend, pident);
    }
    fclose(in);

    // compute best reference per query
    for (int i = 0; i < HASH_SIZE; i++) {
        for (Query* q = queries[i]; q; q = q->next) {
            RefGroup* best_r = NULL;
            int best_aln = -1;
            double best_ident = -1.0;

            for (RefGroup* r = q->refs; r; r = r->next) {
                int merged_aln; double merged_ident;
                compute_merged(r, &merged_aln, &merged_ident);

                if (merged_aln > best_aln ||
                    (merged_aln == best_aln && merged_ident > best_ident)) {
                    best_aln = merged_aln;
                    best_ident = merged_ident;
                    best_r = r;
                }
            }

            if (best_r) {
                int merged_aln; double merged_ident;
                compute_merged(best_r, &merged_aln, &merged_ident);

                double coverage = best_r->qlen > 0 ? (double)merged_aln / best_r->qlen : 0.0;
                double identity = merged_aln > 0 ? merged_ident / merged_aln : 0.0;

                // clamp
                if (identity < 0.0) identity = 0.0;
                if (identity > 1.0) identity = 1.0;

                fprintf(out, "%s\t%s\t%d\t%d\t%.2f\t%.2f\n",
                    q->qseqid, best_r->sseqid, best_r->qlen, merged_aln,
                    coverage * 100.0, identity * 100.0);
            }
        }
    }

    // free memory
    for (int i = 0; i < HASH_SIZE; i++) {
        for (Query* q = queries[i]; q;) {
            for (RefGroup* r = q->refs; r;) {
                RefGroup* rnext = r->next;
                free(r->intervals);
                free(r);
                r = rnext;
            }
            Query* qnext = q->next;
            free(q);
            q = qnext;
        }
    }
    fclose(out);
    return 0;
}

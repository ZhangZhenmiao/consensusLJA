#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

// Read a full line from gzFile (unlimited length)
static ssize_t gz_getline_dyn(gzFile fp, char** line, size_t* cap) {
    if (*line == NULL || *cap == 0) {
        *cap = 1024;
        *line = malloc(*cap);
    }

    size_t len = 0;
    int c;

    while ((c = gzgetc(fp)) != -1) {
        if (len + 1 >= *cap) {
            *cap *= 2;
            *line = realloc(*line, *cap);
        }
        if (c == '\n') break;
        (*line)[len++] = (char)c;
    }

    if (c == -1 && len == 0)
        return -1;  // EOF

    (*line)[len] = '\0';
    return (ssize_t)len;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input.fastq.gz> <output.fastq>\n", argv[0]);
        return 1;
    }

    gzFile fin = gzopen(argv[1], "rb");
    if (!fin) {
        fprintf(stderr, "Error opening input.\n");
        return 1;
    }

    FILE* fout = fopen(argv[2], "w");
    if (!fout) {
        fprintf(stderr, "Error opening output.\n");
        gzclose(fin);
        return 1;
    }

    char* h = NULL, * s = NULL, * p = NULL, * q = NULL;
    size_t cap_h = 0, cap_s = 0, cap_p = 0, cap_q = 0;

    while (1) {
        if (gz_getline_dyn(fin, &h, &cap_h) < 0) break;
        if (gz_getline_dyn(fin, &s, &cap_s) < 0) break;
        if (gz_getline_dyn(fin, &p, &cap_p) < 0) break;
        if (gz_getline_dyn(fin, &q, &cap_q) < 0) break;

        // Trim after first whitespace
        for (size_t i = 0; h[i]; i++) {
            if (h[i] == ' ' || h[i] == '\t') {
                h[i] = '\0';
                break;
            }
        }

        fprintf(fout, "%s\n%s\n%s\n%s\n", h, s, p, q);
    }

    free(h); free(s); free(p); free(q);
    gzclose(fin);
    fclose(fout);

    return 0;
}

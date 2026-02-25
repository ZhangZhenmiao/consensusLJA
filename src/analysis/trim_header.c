#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdbool.h>

typedef struct {
    bool is_gz;
    FILE *fp;
    gzFile gzp;
} auto_file;

/* ---------- helpers ---------- */

static bool has_gz_suffix(const char *s) {
    size_t n = strlen(s);
    return (n > 3 && strcmp(s + n - 3, ".gz") == 0);
}

static auto_file auto_open(const char *path, const char *mode) {
    auto_file f = {0};
    f.is_gz = has_gz_suffix(path);
    if (f.is_gz)
        f.gzp = gzopen(path, mode);
    else
        f.fp = fopen(path, mode);
    return f;
}

static void auto_close(auto_file *f) {
    if (f->is_gz)
        gzclose(f->gzp);
    else
        fclose(f->fp);
}

/* ---------- unified getline ---------- */

static ssize_t auto_getline(auto_file *f, char **line, size_t *cap) {
    if (*line == NULL || *cap == 0) {
        *cap = 1024;
        *line = malloc(*cap);
    }

    size_t len = 0;
    int c;

    while (1) {
        if (f->is_gz)
            c = gzgetc(f->gzp);
        else
            c = fgetc(f->fp);

        if (c == EOF || c == -1) break;

        if (len + 1 >= *cap) {
            *cap *= 2;
            *line = realloc(*line, *cap);
        }

        if (c == '\n') break;
        (*line)[len++] = (char)c;
    }

    if (c == EOF || c == -1) {
        if (len == 0) return -1;
    }

    (*line)[len] = '\0';
    return (ssize_t)len;
}

/* ---------- unified write ---------- */

static void auto_puts(auto_file *f, const char *s) {
    if (f->is_gz) {
        gzputs(f->gzp, s);
        gzputc(f->gzp, '\n');
    } else {
        fputs(s, f->fp);
        fputc('\n', f->fp);
    }
}

/* ---------- main ---------- */

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr,
                "Usage: %s <input.fastq[.gz]> <output.fastq[.gz]>\n",
                argv[0]);
        return 1;
    }

    auto_file fin  = auto_open(argv[1], "rb");
    auto_file fout = auto_open(argv[2], "wb");

    if ((!fin.is_gz && !fin.fp) || (fin.is_gz && !fin.gzp)) {
        fprintf(stderr, "Error opening input\n");
        return 1;
    }
    if ((!fout.is_gz && !fout.fp) || (fout.is_gz && !fout.gzp)) {
        fprintf(stderr, "Error opening output\n");
        auto_close(&fin);
        return 1;
    }

    char *h = NULL, *s = NULL, *p = NULL, *q = NULL;
    size_t ch = 0, cs = 0, cp = 0, cq = 0;

    while (1) {
        if (auto_getline(&fin, &h, &ch) < 0) break;
        if (auto_getline(&fin, &s, &cs) < 0) break;
        if (auto_getline(&fin, &p, &cp) < 0) break;
        if (auto_getline(&fin, &q, &cq) < 0) break;

        /* trim header at first whitespace */
        for (size_t i = 0; h[i]; i++) {
            if (h[i] == ' ' || h[i] == '\t') {
                h[i] = '\0';
                break;
            }
        }

        auto_puts(&fout, h);
        auto_puts(&fout, s);
        auto_puts(&fout, p);
        auto_puts(&fout, q);
    }

    free(h); free(s); free(p); free(q);
    auto_close(&fin);
    auto_close(&fout);

    return 0;
}
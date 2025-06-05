/* Simple tool that prints the i-th line of every input line read. */

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/* 64 KiB ought to be enough for anyone.
If lines are longer than this, the output will not be correct. */
static char line[65536];

bool parse_long(const char *s, long long *i) {
    char *end = NULL;
    if (*s == '\0') return false;
    *i = strtoll(s, &end, 10);
    return *end == '\0';
}

int main(int argc, char *argv[]) {
    long long interval = 0;
    bool opt_n = false;
    int argi = 1;
    if (argi < argc && strcmp(argv[argi], "-n") == 0) {
        opt_n = true;
        ++argi;
    }
    if (argc - argi != 1 || !parse_long(argv[argi], &interval) || interval < 1) {
        printf("Usage: print-every [-n] <i>\n"
            "Prints every <i>-th line.\n"
            "With option -n, prefixes each line with the original line number.\n");
        return 1;
    }
    long long number = 0;
    long long next = interval;
    while (fgets(line, sizeof(line), stdin)) {
        ++number;
        if (--next == 0) {
            if (opt_n) {
                printf("%20lld  %s", number, line);
            } else {
                fputs(line, stdout);
            }                
            next = interval;
        }
    }
}

/*

    ttcdt-huff - Huffman encoding library.

    ttcdt <dev@triptico.com>

    This software is released into the public domain.

*/

#include <stdio.h>
#include <string.h>

#include "duffman.h"

/* total number of tests and oks */
int tests = 0;
int oks = 0;
int verbose = 1;

/* failed tests messages */
char *failed_msgs[5000];
int i_failed_msgs = 0;

/** code **/

void _do_test(char *str, int ok, int src_line)
{
    char tmp[1024];

    sprintf(tmp, "stress.c:%d: (#%d) %s: %s\n",
        src_line, tests, str, ok ? "OK!" : "*** Failed ***"
    );

    if (!ok || verbose)
        printf("%s", tmp);

    tests++;

    if (ok)
        oks++;
    else
        failed_msgs[i_failed_msgs++] = strdup(tmp);
}

#define do_test(str, ok) _do_test(str, ok, __LINE__)


unsigned char udata[30000];
int uz;
unsigned char cdata[30000];
int cz;
unsigned char buf[30000];
int bz;

void test_1(char *id)
{
    const unsigned char *ptr;

    memset(cdata, '\0', sizeof(cdata));
    memset(buf, '\0', sizeof(buf));

    ptr = duffman_compress(udata, uz, cdata);
    cz = ptr - cdata;

    duffman_size(cdata, &bz);

    if (verbose)
        printf("test: %s -- uz: %d, cz: %d, bz: %d\n", id, uz, cz, bz);

    do_test("Data size stored in compressed block", bz == uz);

    ptr = duffman_decompress(cdata, buf);

    do_test("Decompression without errors", ptr != NULL);

    if (ptr != NULL)
        do_test("Pre-compressed un-compressed comparison", memcmp(udata, buf, uz) == 0);
}


int main(int argc, char *argv[])
{
    FILE *f;
    int n;

    f = fopen("carcosa.txt", "r");
    uz = fread(udata, 1, sizeof(udata), f);
    fclose(f);

    test_1("carcosa.txt");

    for (n = 0; n < 20000; n++)
        udata[n] = 'A' + (n % ('Z' - 'A'));
    uz = 20000;

    test_1("Letter sequence");

    printf("\n*** Total tests passed: %d/%d\n", oks, tests);

    if (oks == tests)
        printf("*** ALL TESTS PASSED\n");
    else {
        int n;

        printf("*** %d %s\n", tests - oks, "TESTS ---FAILED---");

        printf("\nFailed tests:\n\n");
        for (n = 0; n < i_failed_msgs; n++)
            printf("%s", failed_msgs[n]);
    }

    return oks == tests ? 0 : 1;
}

/* duffman - Huffman encoding software by grunfink - public domain */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "duffman.h"

void usage(char *argv0)
{
    printf("%s - Huffman file [de]compressor\n", argv0);

    printf("Usage:\n");
    printf("  %s -C     Compress STDIN to STDOUT\n",   argv0);
    printf("  %s -D     Decompress STDIN to STDOUT\n", argv0);
}


#define CHUNK_SIZE 16384

int compress(FILE *i, FILE *o)
{
    int ret = 0;
    int z;
    uint8_t bi[CHUNK_SIZE];
    uint8_t bo[CHUNK_SIZE * 2];

    while ((z = fread(bi, 1, CHUNK_SIZE, i))) {
        uint8_t *ptr;
        int cz;

        ptr = duffman_compress(bi, z, bo);
        cz = ptr - bo;

        if (cz >= z) {
            /* non-compressed block */
            z = -z;
            fwrite(&z, sizeof(z), 1, o);
            fwrite(bi, 1, -z, o);
        }
        else {
            fwrite(&cz, sizeof(cz), 1, o);
            fwrite(bo, 1, cz, o);
        }
    }

    return ret;
}


int decompress(FILE *i, FILE *o)
{
    int ret = 0;
    int z;
    uint8_t bi[CHUNK_SIZE];
    uint8_t bo[CHUNK_SIZE];

    while (fread(&z, sizeof(z), 1, i)) {
        if (z < 0) {
            /* non-compressed block */
            fread(bi, 1, -z, i);
            fwrite(bi, 1, -z, o);
        }
        else {
            int dz;

            fread(bi, 1, z, i);

            duffman_size(bi, &dz);

            if (dz > CHUNK_SIZE) {
                fprintf(stderr, "error: corrupted stream\n");
                ret = 4;
                break;
            }

            if (duffman_decompress(bi, bo) == NULL) {
                fprintf(stderr, "error: corrupted stream\n");
                ret = 4;
                break;
            }

            fwrite(bo, 1, dz, o);
        }
    }

    return ret;
}


int main(int argc, char *argv[])
{
    int ret = 0;

    if (argc == 1) {
        usage(argv[0]);
        ret = 1;
    }
    else
    if (strcmp(argv[1], "-C") == 0) {
        ret = compress(stdin, stdout);
    }
    else
    if (strcmp(argv[1], "-D") == 0) {
        ret = decompress(stdin, stdout);
    }
    else {
        usage(argv[0]);
        ret = 2;
    }

    return ret;
}

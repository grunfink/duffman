/* duffman - Huffman encoding software by grunfink - OH YEAH - public domain */

/* Yes, I know actual code in .h files is noob, so STFU */

#ifndef DUFFMAN_H
#define DUFFMAN_H

#include <string.h>
#include <stdint.h>

#define DUFFMAN_VERSION "2.00-wrk"

struct node {
    int freq;       /* frequency */
    int next;       /* next (in sort order) */
    int branch[2];  /* branches (0: left, 1: right) */
    int chr;        /* char */
};

#define NUM_NODES 512


/** utility functions **/

static uint8_t *write_bits(uint8_t *obuf, int *mask, int count, int v)
/* write @count bits of @v into @obuf, in reverse order */
{
    while (count--) {
        if (v & 0x1)
            *obuf |= *mask;
        else
            *obuf &= ~*mask;

        v >>= 1;
        (*mask) >>= 1;

        if (*mask == 0) {
            *mask = 0x80;
            obuf++;
        }
    }

    return obuf;
}


static const uint8_t *read_bits(const uint8_t *ibuf, int *mask, int count, int *v)
/* read @count bits from @ibuf into @v, in reverse order */
{
    int om = 0x1;

    while (count--) {
        if (*ibuf & *mask)
            *v |= om;
        else
            *v &= ~om;

        (*mask) >>= 1;
        om <<= 1;

        if (*mask == 0) {
            *mask = 0x80;
            ibuf++;
        }
    }

    return ibuf;
}


/** compression **/

static int insert_node(int i, int *z, int f, int l, int r, int c, struct node *tree)
/* inserts a node, sorted by frequency */
{
    int n, p;

    /* loop until a value with a higher frequency is found */
    for (n = i, p = -1; n != -1 && tree[n].freq < f; p = n, n = tree[n].next);

    /* insert in that position */
    if (p != -1)
        tree[p].next = *z;
    else
        i = *z;

    tree[*z].next      = n;
    tree[*z].freq      = f;
    tree[*z].branch[0] = l;
    tree[*z].branch[1] = r;
    tree[*z].chr       = c;

    *z = *z + 1;

    return i;
}


static int build_tree_from_data(const uint8_t *ibuf, int uz, struct node *tree)
/* builds a tree from an input buffer */
{
    int n, m, r;
    int freqs[256];
    int tz = 0;

    /* reset */
    memset(tree,  '\0', sizeof(struct node) * NUM_NODES);
    memset(freqs, '\0', sizeof(freqs));

    /* count frequency of symbols in data */
    for (n = 0; n < uz; n++)
        freqs[ibuf[n]]++;

    /* add all symbols as leaf nodes */
    m = -1;
    for (n = 0; n < 256; n++) {
        if (freqs[n])
            m = insert_node(m, &tz, freqs[n], -1, -1, n, tree);
    }

    /* build the internal nodes:
       they are created by iterating the tree from the lowest
       frequency nodes. An internal node has two leaf nodes
       and has a frequency that is the sum of both */
    for (r = m; tree[r].next != -1; r = tree[n].next) {
        n = tree[r].next;
        m = insert_node(m, &tz, tree[r].freq + tree[n].freq, r, n, '\0', tree);
    }

    return r;
}


static void build_symbols(const struct node *tree, int r, int b, int v, int *n_bits, int *values)
/* builds the bits and values from a tree (recursive) */
{
    if (tree[r].branch[0] == -1 && tree[r].branch[1] == -1) {
        /* found leaf node */
        n_bits[tree[r].chr] = b;
        values[tree[r].chr] = v;
    }
    else {
        build_symbols(tree, tree[r].branch[0], b + 1, v,            n_bits, values);
        build_symbols(tree, tree[r].branch[1], b + 1, v | (1 << b), n_bits, values);
    }
}


static uint8_t *compress_tree(const struct node *tree, uint8_t *obuf)
/* compresses a tree into @obuf. Returns a pointer to the next byte in @obuf */
{
    int n, c, xc;
    int mask = 0x80;
    uint8_t *cptr;

    /* save the position where the count of leaf nodes will be stored */
    cptr = obuf;
    obuf++;

    n = xc = 0;
    while (tree[n].branch[0] == -1 && tree[n].branch[1] == -1) {
        c = 0;

        /* count how many entries in the tree are consecutive */
        while (tree[n].branch[0] == -1 && tree[n].branch[1] == -1 && xc == tree[n].chr) {
            c++;
            xc = tree[n].chr + 1;
            n++;
        }

        if (c) {
            /* store a run-length count (bit: 0) */
            obuf = write_bits(obuf, &mask, 1, 0);
            obuf = write_bits(obuf, &mask, 8, c);
        }

        /* store all unexpected values as verbatim symbols (bit: 1) */
        while (tree[n].branch[0] == -1 && tree[n].branch[1] == -1 && xc != tree[n].chr) {
            obuf = write_bits(obuf, &mask, 1, 1);
            obuf = write_bits(obuf, &mask, 8, tree[n].chr);
            xc = tree[n].chr + 1;

            n++;
        }
    }

    /* store the count */
    /* note: 0 means 256 nodes (truncation helps us) */
    *cptr = n;

    /* align to byte */
    if (mask != 0x80) {
        obuf++;
        mask = 0x80;
    }

    /* save the position where count of internal nodes will be stored */
    cptr = obuf;
    obuf++;

    /* store only branches */
    for (c = 0; tree[n].branch[0] || tree[n].branch[1]; n++, c++) {
        obuf = write_bits(obuf, &mask, 9, tree[n].branch[0]);
        obuf = write_bits(obuf, &mask, 9, tree[n].branch[1]);
    }

    /* store the count */
    *cptr = c;

    /* align to byte */
    if (mask != 0x80)
        obuf++;

    return obuf;
}


static uint8_t *compress_stream(const uint8_t *ibuf, int uz,
                                uint8_t *obuf, int *n_bits, int *values)
/* compresses a stream of bytes in @ibuf to Huffman symbols written onto @obuf.
   Returns the pointer to the next byte of @obuf */
{
    int n;
    int mask = 0x80;

    for (n = 0; n < uz; n++)
        obuf = write_bits(obuf, &mask, n_bits[ibuf[n]], values[ibuf[n]]);

    if (mask != 0x80)
        obuf++;

    return obuf;
}


/** decompression **/

const uint8_t *decompress_tree(const uint8_t *ibuf, int *r, struct node *tree)
/* decompresses a compressed tree from inside @ibuf into a usable tree.
   The root node will be stored into @r.
   The tree will not have frequency information,
   but that is not necessary for decompression */
{
    int n, c, xc;
    int mask = 0x80;

    /* reset the tree */
    memset(tree, '\0', sizeof(struct node) * NUM_NODES);

    /* get the count of leaf nodes */
    c = *ibuf;
    ibuf++;

    /* '0' elements means there is an entry for every char */
    if (c == 0)
        c = 256;

    n = xc = 0;
    while (n < c) {
        int p, v;

        /* read prefix and value */
        p = v = 0;
        ibuf = read_bits(ibuf, &mask, 1, &p);
        ibuf = read_bits(ibuf, &mask, 8, &v);

        if (p == 0) {
            /* prefix is 0: run-length sequence of consecutive values */
            if (v == 0)
                v = 256;

            while (v) {
                tree[n].branch[0] = tree[n].branch[1] = -1;
                tree[n].chr = xc;
                xc++;
                n++;
                v--;
            }
        }
        else {
            /* prefix is 1: as-is char */
            tree[n].branch[0] = tree[n].branch[1] = -1;
            tree[n].chr = v;
            xc = v + 1;
            n++;
        }
    }

    if (mask != 0x80) {
        ibuf++;
        mask = 0x80;
    }

    /* get the count of internal nodes */
    c += *ibuf;
    ibuf++;

    for (; n < c; n++) {
        ibuf = read_bits(ibuf, &mask, 9, &tree[n].branch[0]);
        ibuf = read_bits(ibuf, &mask, 9, &tree[n].branch[1]);
    }

    if (mask != 0x80)
        ibuf++;

    /* the root node is always the last one */
    *r = n - 1;

    return ibuf;
}


const uint8_t *decompress_stream(const struct node *tree, int r,
                                 const uint8_t *ibuf, int uz, uint8_t *obuf)
/* decompresses a compressed stream of @uz Huffman symbols.
   Returns the pointer to the next byte of @ibuf */
{
    int n = 0;
    int mask = 0x80;

    for (n = 0; n < uz; n++) {
        int c = -1;
        int nr = r;

        do {
            if (nr < 0 || nr >= NUM_NODES)
                goto end;           /* corrupted stream */
            else
            if (tree[nr].branch[0] == -1 && tree[nr].branch[1] == -1)
                c = tree[nr].chr;   /* symbol found */
            else {
                int b;

                /* pick the branch */
                ibuf = read_bits(ibuf, &mask, 1, &b);
                nr = tree[nr].branch[b & 0x01];
            }
        } while (c == -1);

        obuf[n] = c;
    }

    if (mask != 0x80)
        ibuf++;

end:
    return ibuf;
}


/** interface **/

/**
 * duffman_compress - Compresses a block of data.
 * @ibuf: input buffer
 * @uz: data size in bytes
 * @obuf: output buffer
 *
 * Compresses the @ibuf block of @uz bytes into the buffer
 * pointed by @obuf. @uz must be non-zero. @obuf must
 * be at least @uz size.
 *
 * Returns the pointer to the next byte in @ibuf.
 */
uint8_t *duffman_compress(const uint8_t *ibuf, int uz, uint8_t *obuf)
/* compresses @uz bytes from @ibuf into @obuf.
   Returns the pointer to the next byte of @obuf */
{
    struct node tree[NUM_NODES];
    int r;
    int n_bits[256];
    int values[256];
    int mask = 0x80;

    /* build a tree using this data buffer */
    r = build_tree_from_data(ibuf, uz, tree);

    /* build bits and values */
    build_symbols(tree, r, 0, 0, n_bits, values);

    /* store the number of bytes the decompressed data contains */
    obuf = write_bits(obuf, &mask, 24, uz);

    /* store the tree in compressed form */
    obuf = compress_tree(tree, obuf);

    /* compress the data stream */
    return compress_stream(ibuf, uz, obuf, n_bits, values);
}


/**
 * duffman_size - Gets the size of stored data.
 * @ibuf: input buffer
 * @uz: pointer to store the uncompressed data size.
 *
 * Gets the number of bytes @ibuf will expand to
 * after decompression.
 *
 * Returns the pointer to the next byte in @ibuf.
 */
const uint8_t *duffman_size(const uint8_t *ibuf, int *uz)
/* returns the number of bytes @ibuf will expand to */
{
    int mask = 0x80;

    return read_bits(ibuf, &mask, 24, uz);
}


/**
 * duffman_decompress - Decompresses a block of compressed data.
 * @ibuf: input buffer
 * @obuf: output buffer
 *
 * Decompresses the compressed data block in @ibuf into the
 * buffer pointed by @obuf. The buffer must have enough
 * size for the uncompressed block (see duffman_size()).
 *
 * Returns the pointer to the next byte in @ibuf.
 */
const uint8_t *duffman_decompress(const uint8_t *ibuf, uint8_t *obuf)
/* decompresses @ib into @ob. Returns the pointer to the next byte of @ib */
{
    struct node tree[NUM_NODES];
    int r, uz = 0;

    /* take the expected data size */
    ibuf = duffman_size(ibuf, &uz);

    /* decompress the tree, getting also the root node */
    ibuf = decompress_tree(ibuf, &r, tree);

    /* decompress the stream */
    return decompress_stream(tree, r, ibuf, uz, obuf);
}

#endif /* DUFFMAN_H */

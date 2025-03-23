#include <stdlib.h>
#include <string.h>
 
#include <htree.h>
#include <huff.h>
#include <pqueue.h>

struct meta {
  uint8_t lbs;
  uint8_t pqs;

  struct hnode* ss;
  struct pqueue* pq;
  struct hnode** cpy;
};

static struct meta* meta_new();
static void meta_del(struct meta*);
static size_t meta_size(struct meta*);
static int meta_get(FILE*, struct meta*);
static int meta_put(FILE*, struct meta*);

int encode(FILE* src, FILE* dst) {
  struct meta* meta = meta_new();

  for (uint8_t c = getc(src); !feof(src); c = getc(src))
    if (meta->ss[c].sym == 255) {
      meta->ss[c].sym = c;
      meta->ss[c].frq = 1;
      pq_ins(meta->pq, &meta->ss[c]);
    } else
      pq_incr(meta->pq, meta->ss[c].idx);

  meta->pqs = meta->pq->size;
  meta->cpy = malloc(sizeof(struct hnode) * meta->pqs);

  for (int i = 0; i < meta->pqs; ++i)
    meta->cpy[i] = meta->pq->heap[i];

  struct hnode* root = ht_new(meta->pq);

  fseek(src, 0, SEEK_SET);
  fseek(dst, meta_size(meta), SEEK_SET);

  uint8_t cache = 0;

  for (uint8_t c = getc(src); !feof(src); c = getc(src))
    for (int i = 0; i < meta->ss[c].len; ++i) {
      uint64_t blk = i < 64 ? meta->ss[c].code[0] : meta->ss[c].code[1];
      size_t pos = i < 64 ? i : i - 64;

      if (blk & (1UL << (63 - pos)))
        cache = cache | (0b00000001 << (7 - meta->lbs));

      if (++(meta->lbs) == 8) {
        putc(cache, dst);
        meta->lbs = 0;
        cache = 0;
      }
    }

  if (meta->lbs != 0)
    putc(cache, dst);
  else
    meta->lbs = 8;

  fseek(dst, 0, SEEK_SET);
  meta_put(dst, meta);

  ht_del(root);
  meta_del(meta);

  return 0;
}

int decode(FILE* src, FILE* dst) {
  struct meta* meta = meta_new();

  if (meta_get(src, meta) != 0)
    return -1;

  struct hnode* root = ht_new(meta->pq);
  struct hnode* next = root;

  uint8_t b = getc(src);

  while (!feof(src)) {
    uint8_t n = getc(src);

    if (ferror(src))
      return -1;

    for (int j = 0; j < (feof(src) ? meta->lbs : 8); ++j) {
      if (b & (1 << (7 - j)))
        next = next->rgt;
      else
        next = next->lft;

      if (next == NULL)
        return -1;
      if (next->sym != 255) {
        putc(next->sym, dst);
        next = root;
      }
    }

    b = n;
  }

  ht_del(root);
  meta_del(meta);

  return 0;
}

static struct meta* meta_new() {
  struct meta* meta = malloc(sizeof(struct meta));

  meta->lbs = 0;
  meta->pqs = 0;
  meta->ss = malloc(sizeof(struct hnode) * (UINT8_MAX));
  meta->pq = pq_new(UINT8_MAX);
  meta->cpy = NULL;
  meta->pq->cmp = ht_cmp;
  meta->pq->incr = ht_incr;
  meta->pq->repl = ht_repl;

  for (int i = 0; i < UINT8_MAX; ++i) {
    meta->ss[i].sym = 255;
    meta->ss[i].frq = 0;
    meta->ss[i].idx = 0;
    meta->ss[i].code[0] = 0;
    meta->ss[i].code[1] = 0;
    meta->ss[i].len = 0;
    meta->ss[i].lft = NULL;
    meta->ss[i].rgt = NULL;
  }

  return meta;
}

static void meta_del(struct meta* meta) {
  if (meta->cpy)
    free(meta->cpy);
  pq_del(meta->pq);
  free(meta->ss);
  free(meta);
}

static inline size_t meta_size(struct meta* meta) {
  return sizeof(meta->lbs) + sizeof(meta->pqs) + meta->pqs * 9;
}

static int meta_get(FILE* src, struct meta* meta) {
  meta->lbs = getc(src);
  if (ferror(src))
    return -1;

  meta->pq->size = getc(src);
  if (ferror(src))
    return -1;

  for (int i = 0; i < meta->pq->size; ++i) {
    uint8_t s = getc(src);
    if (ferror(src))
      return -1;

    meta->ss[s].sym = s;
    meta->ss[s].idx = i;

    fread(&meta->ss[s].frq, sizeof(meta->ss[s].frq), 1, src);
    if (ferror(src))
      return -1;

    meta->pq->heap[i] = &meta->ss[s];
  }

  return 0;
}

static int meta_put(FILE* dst, struct meta* meta) {
  putc(meta->lbs, dst);
  putc(meta->pqs, dst);

  for (int i = 0; i < meta->pqs; ++i) {
    putc(meta->cpy[i]->sym, dst);
    fwrite(&meta->cpy[i]->frq, sizeof(meta->cpy[i]->frq), 1, dst);
  }

  return 0;
}

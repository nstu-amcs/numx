#ifndef HTREE_H
#define HTREE_H

#include <stdint.h>

#include <pqueue.h>

struct hnode {
  uint8_t sym;

  size_t frq;
  size_t idx;

  uint64_t code[2];
  size_t len;

  struct hnode *lft, *rgt;
};

struct hnode* ht_new(struct pqueue*);
void ht_del(struct hnode*);

int ht_cmp(void*, void*);
void ht_incr(void*);
void ht_repl(void*, size_t);

#endif  // HTREE

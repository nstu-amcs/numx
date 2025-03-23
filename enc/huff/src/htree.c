#include <stdlib.h>

#include <htree.h>

static void extend(struct hnode* a, int d);

struct hnode* ht_new(struct pqueue* pq) {
  if (pq->size == 0)
    return NULL;
  if (pq->size == 1) {
    struct hnode* root = malloc(sizeof(struct hnode));
    root->sym = 255;
    root->lft = pq_extr(pq);
    extend(root->lft, 0);
    return root;
  }

  while (pq->size != 1) {
    struct hnode* a = pq_extr(pq);
    struct hnode* b = pq_extr(pq);
    struct hnode* r = malloc(sizeof(struct hnode));

    r->sym = 255;
    r->frq = a->frq + b->frq;
    r->lft = a;
    r->rgt = b;

    extend(a, 0);
    extend(b, 1);
    pq_ins(pq, r);
  }

  return pq_extr(pq);
}

void ht_del(struct hnode* root) {
  if (root == NULL)
    return;
  if (root->lft)
    ht_del(root->lft);
  if (root->rgt)
    ht_del(root->rgt);
  if (root->sym == -1)
    free(root);
}

int ht_cmp(void* a, void* b) {
  if (a == NULL)
    return -1;
  if (b == NULL)
    return 1;

  struct hnode* n1 = a;
  struct hnode* n2 = b;

  if (n1->frq > n2->frq)
    return 1;
  else if (n1->frq < n2->frq)
    return -1;
  else
    return 0;
}

void ht_incr(void* a) {
  if (a == NULL)
    return;
  struct hnode* n = a;
  n->frq++;
}

void ht_repl(void* a, size_t i) {
  if (a == NULL)
    return;
  struct hnode* n = a;
  n->idx = i;
}

static void extend(struct hnode* a, int d) {
  (a->len)++;
  a->code[1] >>= 1;

  if (a->code[0] & 1UL)
    a->code[1] |= (1UL << 63);

  a->code[0] >>= 1;

  if (d == 1)
    a->code[0] |= (1UL << 63);
  if (a->lft)
    extend(a->lft, d);
  if (a->rgt)
    extend(a->rgt, d);
}

#include <stdlib.h>

#include <pqueue.h>

static void swap(struct pqueue* pq, size_t a, size_t b);
static void fixd(struct pqueue* pq, size_t i);
static void fixu(struct pqueue* pq, size_t i);

struct pqueue* pq_new(size_t cap) {
  struct pqueue* pq = malloc(sizeof(struct pqueue));

  pq->size = 0;
  pq->heap = malloc(cap * sizeof(void*));

  return pq;
}

void pq_del(struct pqueue* pq) {
  if (pq == NULL)
    return;
  free(pq->heap);
  free(pq);
}

void pq_incr(struct pqueue* pq, size_t i) {
  pq->incr(pq->heap[i - 1]);
  fixd(pq, i);
}

void pq_ins(struct pqueue* pq, void* n) {
  pq->repl(n, ++pq->size);
  pq->heap[pq->size - 1] = n;
  fixu(pq, pq->size);
}

void* pq_extr(struct pqueue* pq) {
  void** heap = pq->heap - 1;
  void* min = heap[1];

  pq->repl(heap[pq->size], 1);
  heap[1] = heap[pq->size--];
  fixd(pq, 1);

  return min;
}

static void swap(struct pqueue* pq, size_t a, size_t b) {
  void* tmp = pq->heap[a - 1];

  pq->repl(pq->heap[b - 1], a);
  pq->heap[a - 1] = pq->heap[b - 1];
  pq->repl(tmp, b);
  pq->heap[b - 1] = tmp;
}

static void fixd(struct pqueue* pq, size_t i) {
  void** heap = pq->heap - 1;
  size_t l = i << 1, r = (i << 1) + 1, s = i;

  if (l <= pq->size && pq->cmp(heap[s], heap[l]) == 1)
    s = l;
  if (r <= pq->size && pq->cmp(heap[s], heap[r]) == 1)
    s = r;
  if (s != i) {
    swap(pq, i, s);
    fixd(pq, s);
  }
}

static void fixu(struct pqueue* pq, size_t i) {
  void** heap = pq->heap - 1;
  size_t p = i >> 1;

  if (p > 0 && pq->cmp(heap[i], heap[p]) == -1) {
    swap(pq, i, p);
    fixu(pq, p);
  }
}

#ifndef PQUEUE_H
#define PQUEUE_H

#include <stddef.h>

struct pqueue {
  void** heap;
  size_t size;

  int (*cmp)(void*, void*);
  void (*incr)(void*);
  void (*repl)(void*, size_t);
};

struct pqueue* pq_new(size_t cap);
void pq_del(struct pqueue*);
void pq_ins(struct pqueue*, void*);
void pq_incr(struct pqueue*, size_t);
void* pq_extr(struct pqueue*);

#endif  // PQUEUE_H

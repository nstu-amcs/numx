#ifndef HUFF_H
#define HUFF_H

#include <stdio.h>

int encode(FILE* src, FILE* dst);
int decode(FILE* src, FILE* dst);

#endif  // HUFF_H

#ifndef NUMX_COM_CAP_H
#define NUMX_COM_CAP_H

#define PUB
#define PRV static

typedef struct vcap
{
    void (*run)(void *, int, ...);
    void *ctx;
} vcap;

typedef struct icap
{
    int (*run)(void *, int, ...);
    void *ctx;
} icap;

typedef struct ucap
{
    int (*run)(void *, int, ...);
    void *ctx;
} ucap;

typedef struct dcap
{
    double (*run)(void *, int, ...);
    void *ctx;
} dcap;

#endif // NUMX_COM_CAP_H

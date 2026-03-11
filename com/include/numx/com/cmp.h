#ifndef NUMX_COM_CMP_H
#define NUMX_COM_CMP_H

#define min(a, b) ((a) < (b) ? (a) : (b))
#define fabs(f) (((f) < 0) ? -(f) : (f))
#define isclose(a, b, tol) (fabs((a) - (b)) < (tol))
#define less(a, b, tol) ((b) - (a) > (tol))

int iasc(void *, int, ...);
int idsc(void *, int, ...);

int dasc(void *, int, ...);
int ddsc(void *, int, ...);

#endif // NUMX_COM_CMP_H

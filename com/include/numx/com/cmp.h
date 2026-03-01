#ifndef NUMX_COM_CMP_H
#define NUMX_COM_CMP_H

#define fabs(f) (((f) < 0) ? -(f) : (f))
#define isclose(a, b, tol) (fabs((a) - (b)) < (tol))

int iasc(void *, int, ...);
int idsc(void *, int, ...);

int dasc(void *, int, ...);
int ddsc(void *, int, ...);

#endif // NUMX_COM_CMP_H

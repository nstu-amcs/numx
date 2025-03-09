#include <numx/vec/iss.h>
#include <stdarg.h>

void iss_itr_cap(void* ctx, int n, ...) {
  va_list arg;
  va_start(arg, n);

  ((struct iss_itr*)ctx)->num = va_arg(arg, int);
  ((struct iss_itr*)ctx)->err = va_arg(arg, double);

  va_end(arg);
}
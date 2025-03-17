#include <math.h>
#include <munit.h>
#include <numx/pde/sse.h>
#include <stdarg.h>
#include <stdio.h>

// clang-format off

double fabs(double v) {
  return v > 0 ? v : -v;
}

// double sx(void* ctx, int n, ...) {
//   (void)ctx;
//   (void)n;
// 
//   va_list arg;
//   va_start(arg, n);
// 
//   double xb = va_arg(arg, double);
//   double xe = va_arg(arg, double);
//   double xs = va_arg(arg, double);
//   double x0 = va_arg(arg, double);
// 
//   (void)xb;
//   (void)xe;
// 
//   va_end(arg);
// 
//   if (xs == 0)
//     return 0.2;
// 
//   if (x0 < 0.5)
//     return xs < 0.02 ? xs : xs / 2;
// 
//   return xs * 2;
// }

double sx(void* ctx, int n, ...) {
  (void)ctx;
  (void)n;

  return 0.0635;
}

double tgt(struct vtx* v) {
  return sin(v->x + v->y);
}

double lam(double n, struct vtx* v) {
  (void)n;
  (void)v;

  return 5;
}

double gam(double n, struct vtx* v) {
  (void)n;
  (void)v;

  return 0.4;
}

double ext(double n, struct vtx* v) {
  (void)n;

  return 10.4 * sin(v->x + v->y);
}

double dir(double n, struct vtx* v) {
  (void)n;

  return tgt(v);
}

double neup(double n, struct vtx* v) {
  (void)n;
  (void)v;

  return cos(v->x + v->y);
}

double neun(double n, struct vtx* v) {
  (void)n;
  (void)v;

  return -cos(v->x + v->y);
}


static fun pay[] = {
  &lam,
  &gam,
  &ext,
  &dir,
  &neun,
  &neup,
};

MunitResult test_fdm(const MunitParameter pps[], void* ctx) {
  (void)pps;
  (void)ctx;

  struct obj obj;
  struct vec res;
  struct iss_itr itr;

  struct fun_cut dat;

  fun_cut_new(&dat);
  fun_cut_cov(&dat, pay, sizeof(pay) / sizeof(fun));

  obj_new(&obj);
  obj_get(&obj, fopen("obj/qud-let.obj", "r"), (struct obj_get_ops){
      .get_vtx_ctx = 0,
      .get_hxd_ctx = 0,
      .get_seg_ctx = &((struct icap){.ctx = &dat, .call = &fctx_get}),
      .get_qud_ctx = &((struct icap){.ctx = &dat, .call = &ectx_get}),
  });
  obj_gen(&obj, (struct obj_gen_ops){
      .sx = &((struct dcap){.ctx = 0, .call = &sx}),
      .sy = &((struct dcap){.ctx = 0, .call = &sx}),
      .sz = 0,
      .eps = 0.01,
  });

  vec_new(&res, obj.ax.len * obj.ay.len * obj.az.len);

  pde_sse_fdm_slv(&obj, &res, (struct sse_fdm_ops){
    .ops = {
      .iss.mod = ISS_RLX,
      .iss.ops.rlx = {
        .omg = 1.2,
        .ops.eps = 1e-10,
        .ops.max = 10000,
        .ops.itr.ctx = &itr,
        .ops.itr.call = &iss_itr_cap
      }
    }
  });

  int sp = obj_srh(&obj, 0.5, 0.5, 0);
  double err = fabs(res.dat[sp] - tgt(&((struct vtx){.x = 0.5, .y = 0.5, .z = 0})));

  printf("\n\t-> val = %.7e\n\t-> err = %.7e\n\t-> iss/num = %d\n\t-> iss/err = %.7e\n", res.dat[sp], err, itr.num, itr.err);

  return MUNIT_OK;
}

static MunitTest test[] = {
  { "/fdm", test_fdm, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
  { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }
};

static const MunitSuite suite = {
  "/pde/sse", test, NULL, 1, MUNIT_SUITE_OPTION_NONE,
};

int main(int argc, char** argv) {
  return munit_suite_main(&suite, NULL, argc, argv);
}

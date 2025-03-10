#include <math.h>
#include <munit.h>
#include <numx/pde/sse.h>
#include <stdio.h>

// clang-format off

double fabs(double v) {
  return v > 0 ? v : -v;
}

double sq_sx(void* ctx, int n, ...) {
  (void)ctx;
  (void)n;

  return 0.1;
}

double sq_dir(struct vtx* v) {
  return v->x + v->y;
}

double sq_neu(struct vtx* v) {
  (void)v;

  return 1;
}

double sq_neu_dn(struct vtx* v) {
  (void)v;

  return -1;
}

double sq_ext(struct vtx* v) {

  return 0.4 * (v->x + v->y);
}

static vfun dat[] = {
  &sq_dir, // 0
  &sq_ext, // 1
  &sq_neu, // 2
  &sq_neu_dn, // 3
};

struct test {
  const char* name;

  struct {
    double tx;
    double ty;
    double tz;
    double tv;

    double eps;

    struct obj_gen_ops ops;
  } pps;
};

struct test tests[] = {
  { "sq-1", { 
    .tx = 0.5, .ty = 0.5, .tz = 0, .tv = 0.84147, .eps = 0.5,
    .ops.sx = &((struct dcap){.call = &sq_sx, .ctx = 0}),
    .ops.sy = &((struct dcap){.call = &sq_sx, .ctx = 0}),
    .ops.sz = 0,
    .ops.eps = 0.00005,
    .ops.with_vtx = 0,
    .ops.with_seg = 0,
    .ops.with_qud = 0,
    .ops.with_hxd = 0,
  }},
  { "sq-2", { 
    .tx = 0.5, .ty = 0.5, .tz = 0, .tv = 0.84147, .eps = 0.5,
    .ops.sx = &((struct dcap){.call = &sq_sx, .ctx = 0}),
    .ops.sy = &((struct dcap){.call = &sq_sx, .ctx = 0}),
    .ops.sz = 0,
    .ops.eps = 0.00005,
    .ops.with_vtx = 0,
    .ops.with_seg = 0,
    .ops.with_qud = 0,
    .ops.with_hxd = 0,
  }},
  { "sq-3", { 
    .tx = 0.5, .ty = 0.5, .tz = 0, .tv = 0.84147, .eps = 0.5,
    .ops.sx = &((struct dcap){.call = &sq_sx, .ctx = 0}),
    .ops.sy = &((struct dcap){.call = &sq_sx, .ctx = 0}),
    .ops.sz = 0,
    .ops.eps = 0.00005,
    .ops.with_vtx = 0,
    .ops.with_seg = 0,
    .ops.with_qud = 0,
    .ops.with_hxd = 0,
  }},
  { "t-1", { 
    .tx = 0, .ty = 8, .tz = 0, .tv = 8, .eps = 0.5,
    .ops.sx = &((struct dcap){.call = &sq_sx, .ctx = 0}),
    .ops.sy = &((struct dcap){.call = &sq_sx, .ctx = 0}),
    .ops.sz = 0,
    .ops.eps = 0.00005,
    .ops.with_vtx = 0,
    .ops.with_seg = 0,
    .ops.with_qud = 0,
    .ops.with_hxd = 0,
  }},
};

struct ctx {
  struct test* test;

  struct obj* obj;
  struct vec* res;

  struct itr {
    int k;
    double res;
  } itr;
};

void cback(void* ctx, int n, ...);

void* test_sup(const MunitParameter pps[], void* non) {
  (void)non;

  int n = atoi(pps[0].value);

  struct test* t = &tests[n];
  struct obj* o = malloc(sizeof(struct obj));

  struct ocut dcut;

  ocut_new(&dcut);
  ocut_cov(&dcut, dat, sizeof(dat) / sizeof(double (*)(struct vtx*)));
  
  char fname[64] = "obj/";

  strcat(fname, t->name);
  strcat(fname, ".obj");

  FILE* f = fopen(fname, "r");

  obj_get(o, f, (struct obj_get_ops){
    .get_vtx_ctx = 0,
    .get_hxd_ctx = 0,
    .get_seg_ctx = &((struct icap){
      .ctx = &dcut,
      .call = &fctx_get,
    }),
    .get_qud_ctx = &((struct icap){
      .ctx = &dcut,
      .call = &ectx_get,
    })
  });

  fclose(f);
  obj_gen(o, t->pps.ops);

  struct ctx* ctx = malloc(sizeof(struct ctx));
  
  ctx->test = t;
  ctx->obj = o;
  ctx->res = malloc(sizeof(struct vec));

  return ctx;
}

void test_tdn(void* ctx) {
  struct ctx* c = (struct ctx*)ctx;

  obj_cls(c->obj);
  vec_cls(c->res);

  free(c->obj);
  free(c->res);
}

MunitResult test_fdm(const MunitParameter pps[], void* ctx) {
  (void)pps;

  struct ctx* c = (struct ctx*)ctx;
  struct iss_itr itr;

  munit_assert_int(0, ==, pde_sse_fdm_slv(c->obj, c->res, (struct sse_fdm_ops){
    .iss_type = ISS_RLX,
    .iss_ops.rlx = {
      .omg = 1.2,
      .ops.eps = 1e-6,
      .ops.max = 10000,
      .ops.itr.ctx = &itr,
      .ops.itr.call = &iss_itr_cap
    }
  }));

  int sp = obj_srh(c->obj, 
    c->test->pps.tx,
    c->test->pps.ty,
    c->test->pps.tz
  );

  double err = fabs(c->res->dat[sp] - c->test->pps.tv);

  printf("\n[%s]\n\t-> n: %d\n\t-> target: %.7e\n\t-> result: %.7e\n\t-> error: %.7e\n\t-> iss/num: %d\n\t-> iss/err: %.7e\n", 
    c->test->name, 
    c->obj->ax.len * c->obj->ay.len,
    c->test->pps.tv,
    c->res->dat[sp], 
    err,
    itr.num,
    itr.err
  );

  return MUNIT_OK;
}

static char* cases[] = { "3" };

static MunitParameterEnum pps[] = {
  { "case", cases },
  { NULL, NULL },
};

static MunitTest test[] = {
  { "/fdm", test_fdm, test_sup, test_tdn, MUNIT_TEST_OPTION_NONE, pps},
  { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }
};

static const MunitSuite suite = {
  "/pde/sse", test, NULL, 1, MUNIT_SUITE_OPTION_NONE,
};

int main(int argc, char** argv) {
  return munit_suite_main(&suite, NULL, argc, argv);
}

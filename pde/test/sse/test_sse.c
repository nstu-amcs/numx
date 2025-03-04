#include <munit.h>
#include <numx/pde/sse.h>
#include <stdio.h>

// clang-format off

double dir(struct vtx* v) {
  (void)v;

  return 0;
}

double neu(struct vtx* v) {
  (void)v;

  return 1;
}

double ext(struct vtx* v) {
  return 0.4 * (v->x + v->y);
}

double step(void* ctx, int n, ...) {
  (void)ctx;
  (void)n;

  return 0.01;
}

__attribute_maybe_unused__
static MunitResult test_ctx(const MunitParameter pps[], void* ctx) {
  (void)pps;
  (void)ctx;

  struct ocut dat;

  ocut_new(&dat);
  ocut_add(&dat, &dir);
  ocut_add(&dat, &neu);
  ocut_add(&dat, &ext);

  FILE *f = fopen("data/case1.obj", "r");
  struct obj obj;

  munit_assert_int(0, ==, obj_new(&obj));
  munit_assert_int(0, ==, obj_get(&obj, f, (struct obj_get_ops){
    .get_vtx_ctx = 0,
    .get_hxd_ctx = 0,
    .get_seg_ctx = &((struct icap){
      .ctx = &dat,
      .call = &fctx_get,
    }),
    .get_qud_ctx = &((struct icap){
      .ctx = &dat,
      .call = &ectx_get,
    })
  }));

  munit_assert_int(DIR, ==, ((struct fctx*)obj.seg.dat[0].ctx)->cnd.type);
  munit_assert_int(DIR, ==, ((struct fctx*)obj.seg.dat[1].ctx)->cnd.type);
  munit_assert_int(NEU, ==, ((struct fctx*)obj.seg.dat[2].ctx)->cnd.type);
  munit_assert_int(DIR, ==, ((struct fctx*)obj.seg.dat[3].ctx)->cnd.type);

  munit_assert(&dir == ((struct fctx*)obj.seg.dat[0].ctx)->cnd.pps.dir.tmp);
  munit_assert(&dir == ((struct fctx*)obj.seg.dat[1].ctx)->cnd.pps.dir.tmp);
  munit_assert(&neu == ((struct fctx*)obj.seg.dat[2].ctx)->cnd.pps.neu.tta);
  munit_assert(&dir == ((struct fctx*)obj.seg.dat[3].ctx)->cnd.pps.dir.tmp);

  munit_assert_double_equal(5.0, ((struct ectx*)obj.qud.dat[0].ctx)->lam, 15);
  munit_assert_double_equal(0.4, ((struct ectx*)obj.qud.dat[0].ctx)->gam, 15);

  munit_assert(&ext == ((struct ectx*)obj.qud.dat[0].ctx)->ext);

  munit_assert_int(0, ==, obj_cls(&obj));

  ocut_cls(&dat);
  fclose(f);

  return MUNIT_OK;
}

MunitResult test_fdm_case1(const MunitParameter pps[], void* ctx) {
  (void)pps;
  (void)ctx;

  struct ocut dat;
  struct obj obj;

  ocut_new(&dat);
  ocut_add(&dat, &dir);
  ocut_add(&dat, &neu);
  ocut_add(&dat, &ext);

  FILE *f = fopen("data/case1.obj", "r");

  munit_assert_int(0, ==, obj_new(&obj));
  munit_assert_int(0, ==, obj_get(&obj, f, (struct obj_get_ops){
    .get_vtx_ctx = 0,
    .get_hxd_ctx = 0,
    .get_seg_ctx = &((struct icap){
      .ctx = &dat,
      .call = &fctx_get,
    }),
    .get_qud_ctx = &((struct icap){
      .ctx = &dat,
      .call = &ectx_get,
    })
  }));

  fclose(f);

  struct vec x;
  struct dcap cap = {
    .ctx = 0,
    .call = &step
  };
  
  munit_assert_int(0, ==, obj_gen(&obj, (struct obj_gen_ops){
    .with_vtx = false,
    .with_seg = false,
    .with_qud = false,
    .with_hxd = false,

    .sx = &cap,
    .sy = &cap,
    .sz = 0
  }));

  munit_assert_int(0, ==, pde_sse_fdm_slv(&obj, &x));

  printf("Cortex: %lf\n", x.dat[x.n / 2 + 1]);

  munit_assert_int(0, ==, obj_cls(&obj));
  munit_assert_int(0, ==, ocut_cls(&dat));


  return MUNIT_OK;
}

MunitResult test_fdm_case2(const MunitParameter pps[], void* ctx) {
  (void)pps;
  (void)ctx;

  struct ocut dat;
  struct obj obj;

  ocut_new(&dat);
  ocut_add(&dat, &dir);
  ocut_add(&dat, &ext);

  FILE *f = fopen("data/case2.obj", "r");

  munit_assert_int(0, ==, obj_new(&obj));
  munit_assert_int(0, ==, obj_get(&obj, f, (struct obj_get_ops){
    .get_vtx_ctx = 0,
    .get_hxd_ctx = 0,
    .get_seg_ctx = &((struct icap){
      .ctx = &dat,
      .call = &fctx_get,
    }),
    .get_qud_ctx = &((struct icap){
      .ctx = &dat,
      .call = &ectx_get,
    })
  }));

  fclose(f);

  struct vec x;
  struct dcap cap = {
    .ctx = 0,
    .call = &step
  };
  
  munit_assert_int(0, ==, obj_gen(&obj, (struct obj_gen_ops){
    .with_vtx = false,
    .with_seg = false,
    .with_qud = false,
    .with_hxd = false,

    .sx = &cap,
    .sy = &cap,
    .sz = 0
  }));

  munit_assert_int(0, ==, pde_sse_fdm_slv(&obj, &x));

  printf("Cortex: %lf\n", x.dat[x.n / 2 + 1]);

  munit_assert_int(0, ==, obj_cls(&obj));
  munit_assert_int(0, ==, ocut_cls(&dat));

  return MUNIT_OK;
}

static MunitTest tests[] = {
  { "/ctx", test_ctx, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
  { "/fdm/case1", test_fdm_case1, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
  { "/fdm/case2", test_fdm_case2, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
  { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }
};

static const MunitSuite suite = {
  "/pde/sse", tests, NULL, 1, MUNIT_SUITE_OPTION_NONE,
};

int main(int argc, char** argv) {
  return munit_suite_main(&suite, NULL, argc, argv);
}

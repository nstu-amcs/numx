#include <munit.h>
#include <numx/pde/geo.h>
#include <stdio.h>

// clang-format off

static MunitResult test_get(const MunitParameter pps[], void* ctx) {
  (void)pps;
  (void)ctx;

  FILE *f = fopen("data/case1.obj", "r");

  struct obj obj;

  munit_assert_int(0, ==, obj_new(&obj));
  munit_assert_int(0, ==, obj_get(&obj, f, (struct obj_get_ops){ 0, 0, 0, 0 }));

  fclose(f);

  munit_assert_int(4, ==, obj.vtx.len);
  munit_assert_int(4, ==, obj.seg.len);
  munit_assert_int(1, ==, obj.qud.len);
  munit_assert_int(0, ==, obj.hxd.len);

  munit_assert_int(0, ==, obj_cls(&obj));

  return MUNIT_OK;
}

static MunitResult test_gen_idle(const MunitParameter pps[], void* ctx) {
  (void)pps;
  (void)ctx;

  FILE *f = fopen("data/case1.obj", "r");

  struct obj obj;

  munit_assert_int(0, ==, obj_new(&obj));
  munit_assert_int(0, ==, obj_get(&obj, f, (struct obj_get_ops){ 0, 0, 0, 0 }));

  fclose(f);

  munit_assert_int(0, ==, obj_gen(&obj, (struct obj_gen_ops){
    .sx = 0,
    .sy = 0,
    .sz = 0
  }));

  munit_assert_int(2, ==, obj.ax.len);
  munit_assert_int(2, ==, obj.ay.len);
  munit_assert_int(1, ==, obj.az.len);

  munit_assert_double_equal(0.0, obj.ax.dat[0], 15);
  munit_assert_double_equal(1.0, obj.ax.dat[1], 15);
  munit_assert_double_equal(0.0, obj.ay.dat[0], 15);
  munit_assert_double_equal(1.0, obj.ay.dat[1], 15);
  munit_assert_double_equal(0.0, obj.az.dat[0], 15);

  munit_assert_int(0, ==, obj.seg.dat[0].vtx[0]);
  munit_assert_int(1, ==, obj.seg.dat[0].vtx[1]);
  munit_assert_int(1, ==, obj.seg.dat[1].vtx[0]);
  munit_assert_int(3, ==, obj.seg.dat[1].vtx[1]);
  munit_assert_int(2, ==, obj.seg.dat[2].vtx[0]);
  munit_assert_int(3, ==, obj.seg.dat[2].vtx[1]);
  munit_assert_int(0, ==, obj.seg.dat[3].vtx[0]);
  munit_assert_int(2, ==, obj.seg.dat[3].vtx[1]);

  munit_assert_int(0, ==, obj.qud.dat[0].vtx[0]);
  munit_assert_int(1, ==, obj.qud.dat[0].vtx[1]);
  munit_assert_int(2, ==, obj.qud.dat[0].vtx[2]);
  munit_assert_int(3, ==, obj.qud.dat[0].vtx[3]);

  munit_assert_int(0, ==, obj_cls(&obj));

  return MUNIT_OK;
}

double step(void* ctx, int n, ...) {
  (void)ctx;
  (void)n;

  return 0.5;
}

static MunitResult test_gen_step(const MunitParameter pps[], void* ctx) {
  (void)pps;
  (void)ctx;

  FILE *f = fopen("data/case1.obj", "r");

  struct obj obj;

  munit_assert_int(0, ==, obj_new(&obj));
  munit_assert_int(0, ==, obj_get(&obj, f, (struct obj_get_ops){ 0, 0, 0, 0 }));

  fclose(f);

  struct dcap s = cap(&step, 0);

  munit_assert_int(0, ==, obj_gen(&obj, (struct obj_gen_ops){
    .sx = &s,
    .sy = &s,
    .sz = &s
  }));

  munit_assert_int(3, ==, obj.ax.len);
  munit_assert_int(3, ==, obj.ay.len);
  munit_assert_int(1, ==, obj.az.len);

  munit_assert_double_equal(0.0, obj.ax.dat[0], 15);
  munit_assert_double_equal(0.5, obj.ax.dat[1], 15);
  munit_assert_double_equal(1.0, obj.ax.dat[2], 15);
  munit_assert_double_equal(0.0, obj.ay.dat[0], 15);
  munit_assert_double_equal(0.5, obj.ay.dat[1], 15);
  munit_assert_double_equal(1.0, obj.ay.dat[2], 15);
  munit_assert_double_equal(0.0, obj.az.dat[0], 15);

  munit_assert_int(0, ==, obj.seg.dat[0].vtx[0]);
  munit_assert_int(2, ==, obj.seg.dat[0].vtx[1]);

  munit_assert_int(2, ==, obj.seg.dat[1].vtx[0]);
  munit_assert_int(8, ==, obj.seg.dat[1].vtx[1]);

  munit_assert_int(6, ==, obj.seg.dat[2].vtx[0]);
  munit_assert_int(8, ==, obj.seg.dat[2].vtx[1]);

  munit_assert_int(0, ==, obj.seg.dat[3].vtx[0]);
  munit_assert_int(6, ==, obj.seg.dat[3].vtx[1]);

  munit_assert_int(0, ==, obj.qud.dat[0].vtx[0]);
  munit_assert_int(2, ==, obj.qud.dat[0].vtx[1]);
  munit_assert_int(6, ==, obj.qud.dat[0].vtx[2]);
  munit_assert_int(8, ==, obj.qud.dat[0].vtx[3]);

  munit_assert_int(0, ==, obj_cls(&obj));

  return MUNIT_OK;
}

static MunitTest tests[] = {
  { "/get", test_get, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
  { "/gen/idle", test_gen_idle, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
  { "/gen/step", test_gen_step, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
  { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }
};

static const MunitSuite suite = {
  "/pde/geo", tests, NULL, 1, MUNIT_SUITE_OPTION_NONE,
};

int main(int argc, char** argv) {
  return munit_suite_main(&suite, NULL, argc, argv);
}

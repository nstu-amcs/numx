#include <munit.h>
#include <numx/pde/geo.h>
#include <stdio.h>

// clang-format off

static MunitResult test_square(const MunitParameter pps[], void* ctx) {
  (void)pps;
  (void)ctx;

  FILE *f = fopen("data/case1.obj", "r");

  struct obj obj;

  munit_assert_int(0, ==, obj_new(&obj));
  munit_assert_int(0, ==, obj_get(&obj, f, (struct obj_get_ops){ 0, 0, 0 }));

  fclose(f);

  munit_assert_int(4, ==, obj.vtx.len);
  munit_assert_int(4, ==, obj.seg.len);
  munit_assert_int(1, ==, obj.qud.len);
  munit_assert_int(0, ==, obj.hxd.len);

  munit_assert_int(0, ==, obj_cls(&obj));

  return MUNIT_OK;
}

static MunitTest tests[] = {
  { "/square", test_square, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
  { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }
};

static const MunitSuite suite = {
  "/pde/geo", tests, NULL, 1, MUNIT_SUITE_OPTION_NONE,
};

int main(int argc, char** argv) {
  return munit_suite_main(&suite, NULL, argc, argv);
}
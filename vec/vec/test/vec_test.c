#include <munit.h>
#include <numx/vec/vec.h>
#include <stdio.h>
#include <string.h>

static MunitResult test_new(const MunitParameter[], void*) {
  struct vec v;

  vec_new_ini(&v, 3, 1.0, 2.0, 3.0);

  munit_assert_double_equal(1.0, v.dat[0], 15);
  munit_assert_double_equal(2.0, v.dat[1], 15);
  munit_assert_double_equal(3.0, v.dat[2], 15);

  vec_cls(&v);

  return MUNIT_OK;
}

static MunitTest tests[] = {
    {"/new", test_new, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
    {NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
};

static const MunitSuite suite = {"/vec/vec", tests, NULL, 1, MUNIT_SUITE_OPTION_NONE};

int main(int argc, char** argv) { return munit_suite_main(&suite, NULL, argc, argv); }

#include <munit.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

#include <numx/non/apx.h>

static MunitResult test_cub(const MunitParameter p[], void* c) {
  (void)p;
  (void)c;

  struct imtx km;
  struct vec xv;

  if (vec_new(&xv, 3)) {
    munit_logf(MUNIT_LOG_ERROR, "vec: %s", strerror(errno));
    return MUNIT_FAIL;
  }

  xv.dat[0] = -6;
  xv.dat[1] = -4;
  xv.dat[2] = -2;

  if (imtx_new(&km, (struct imtx_pps){4, 3})) {
    munit_logf(MUNIT_LOG_ERROR, "mtx: %s", strerror(errno));
    return MUNIT_FAIL;
  }

  km.dat[0][0] = 2;
  km.dat[0][1] = -2;
  km.dat[0][2] = 1;

  if (apx_cub(&xv, &km)) {
    munit_logf(MUNIT_LOG_ERROR, "cub: %s", strerror(errno));
    return MUNIT_FAIL;
  }

  FILE* f = fopen("result.dat", "w+");

  fprintf(f, "a: ");

  for (uint i = 0; i < 3; ++i)
    fprintf(f, "%lf ", km.dat[0][i]);

  fprintf(f, "\nb: ");

  for (uint i = 0; i < 3; ++i)
    fprintf(f, "%lf ", km.dat[1][i]);

  fprintf(f, "\nc: ");

  for (uint i = 0; i < 3; ++i)
    fprintf(f, "%lf ", km.dat[2][i]);

  fprintf(f, "\nd: ");

  for (uint i = 0; i < 3; ++i)
    fprintf(f, "%lf ", km.dat[3][i]);

  fclose(f);

  return MUNIT_OK;
}

static MunitTest cub_tests[] = {
  {"/simple", test_cub, 0, 0, MUNIT_TEST_OPTION_NONE, 0},
  {NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL}};

static MunitSuite suites[] = {
  {"/cub", cub_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
  {NULL, NULL, NULL, 1, MUNIT_SUITE_OPTION_NONE},
};

static const MunitSuite suite = {"non/apx", NULL, suites, 1, MUNIT_SUITE_OPTION_NONE};

int main(int argc, char** argv) {
  return munit_suite_main(&suite, NULL, argc, argv);
}

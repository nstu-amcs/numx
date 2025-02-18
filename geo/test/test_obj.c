#include <munit.h>
#include <numx/geo/obj.h>
#include <stdio.h>
#include <string.h>

static MunitResult test_obj_get(const MunitParameter[], void*) {
  struct obj o;

  munit_assert_int(0, ==, obj_new(&o));
  munit_assert_int(0, ==, obj_get(&o, fopen("data/test.obj", "r")));

  munit_assert_int(8, ==, o.v.len);
  munit_assert_int(6, ==, o.q.len);
  munit_assert_int(1, ==, o.h.len);

  munit_assert_double_equal(0.0, ((struct vtx*)o.v.dat[0])->x, 15);
  munit_assert_double_equal(0.0, ((struct vtx*)o.v.dat[0])->y, 15);
  munit_assert_double_equal(0.0, ((struct vtx*)o.v.dat[0])->z, 15);

  munit_assert_double_equal(1.0, ((struct vtx*)o.v.dat[1])->x, 15);
  munit_assert_double_equal(0.0, ((struct vtx*)o.v.dat[1])->y, 15);
  munit_assert_double_equal(0.0, ((struct vtx*)o.v.dat[1])->z, 15);

  munit_assert_double_equal(0.0, ((struct vtx*)o.v.dat[2])->x, 15);
  munit_assert_double_equal(1.0, ((struct vtx*)o.v.dat[2])->y, 15);
  munit_assert_double_equal(0.0, ((struct vtx*)o.v.dat[2])->z, 15);

  munit_assert_double_equal(1.0, ((struct vtx*)o.v.dat[3])->x, 15);
  munit_assert_double_equal(1.0, ((struct vtx*)o.v.dat[3])->y, 15);
  munit_assert_double_equal(0.0, ((struct vtx*)o.v.dat[3])->z, 15);

  munit_assert_double_equal(0.0, ((struct vtx*)o.v.dat[4])->x, 15);
  munit_assert_double_equal(0.0, ((struct vtx*)o.v.dat[4])->y, 15);
  munit_assert_double_equal(1.0, ((struct vtx*)o.v.dat[4])->z, 15);

  munit_assert_double_equal(1.0, ((struct vtx*)o.v.dat[5])->x, 15);
  munit_assert_double_equal(0.0, ((struct vtx*)o.v.dat[5])->y, 15);
  munit_assert_double_equal(1.0, ((struct vtx*)o.v.dat[5])->z, 15);

  munit_assert_double_equal(0.0, ((struct vtx*)o.v.dat[6])->x, 15);
  munit_assert_double_equal(1.0, ((struct vtx*)o.v.dat[6])->y, 15);
  munit_assert_double_equal(1.0, ((struct vtx*)o.v.dat[6])->z, 15);

  munit_assert_double_equal(1.0, ((struct vtx*)o.v.dat[7])->x, 15);
  munit_assert_double_equal(1.0, ((struct vtx*)o.v.dat[7])->y, 15);
  munit_assert_double_equal(1.0, ((struct vtx*)o.v.dat[7])->z, 15);

  munit_assert_int(1, ==, ((struct qud*)o.q.dat[0])->vtx[0]);
  munit_assert_int(3, ==, ((struct qud*)o.q.dat[0])->vtx[1]);
  munit_assert_int(5, ==, ((struct qud*)o.q.dat[0])->vtx[2]);
  munit_assert_int(7, ==, ((struct qud*)o.q.dat[0])->vtx[3]);

  munit_assert_int(3, ==, ((struct qud*)o.q.dat[1])->vtx[0]);
  munit_assert_int(4, ==, ((struct qud*)o.q.dat[1])->vtx[1]);
  munit_assert_int(7, ==, ((struct qud*)o.q.dat[1])->vtx[2]);
  munit_assert_int(8, ==, ((struct qud*)o.q.dat[1])->vtx[3]);

  munit_assert_int(5, ==, ((struct qud*)o.q.dat[2])->vtx[0]);
  munit_assert_int(6, ==, ((struct qud*)o.q.dat[2])->vtx[1]);
  munit_assert_int(7, ==, ((struct qud*)o.q.dat[2])->vtx[2]);
  munit_assert_int(8, ==, ((struct qud*)o.q.dat[2])->vtx[3]);

  munit_assert_int(2, ==, ((struct qud*)o.q.dat[3])->vtx[0]);
  munit_assert_int(4, ==, ((struct qud*)o.q.dat[3])->vtx[1]);
  munit_assert_int(6, ==, ((struct qud*)o.q.dat[3])->vtx[2]);
  munit_assert_int(8, ==, ((struct qud*)o.q.dat[3])->vtx[3]);

  munit_assert_int(1, ==, ((struct qud*)o.q.dat[4])->vtx[0]);
  munit_assert_int(2, ==, ((struct qud*)o.q.dat[4])->vtx[1]);
  munit_assert_int(5, ==, ((struct qud*)o.q.dat[4])->vtx[2]);
  munit_assert_int(6, ==, ((struct qud*)o.q.dat[4])->vtx[3]);

  munit_assert_int(1, ==, ((struct qud*)o.q.dat[5])->vtx[0]);
  munit_assert_int(2, ==, ((struct qud*)o.q.dat[5])->vtx[1]);
  munit_assert_int(3, ==, ((struct qud*)o.q.dat[5])->vtx[2]);
  munit_assert_int(4, ==, ((struct qud*)o.q.dat[5])->vtx[3]);

  munit_assert_int(1, ==, ((struct hxd*)o.h.dat[0])->vtx[0]);
  munit_assert_int(2, ==, ((struct hxd*)o.h.dat[0])->vtx[1]);
  munit_assert_int(3, ==, ((struct hxd*)o.h.dat[0])->vtx[2]);
  munit_assert_int(4, ==, ((struct hxd*)o.h.dat[0])->vtx[3]);
  munit_assert_int(5, ==, ((struct hxd*)o.h.dat[0])->vtx[4]);
  munit_assert_int(6, ==, ((struct hxd*)o.h.dat[0])->vtx[5]);
  munit_assert_int(7, ==, ((struct hxd*)o.h.dat[0])->vtx[6]);
  munit_assert_int(8, ==, ((struct hxd*)o.h.dat[0])->vtx[7]);

  munit_assert_int(0, ==, obj_cls(&o));

  return MUNIT_OK;
}

static MunitTest tests[] = {
    {"/get", test_obj_get, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
    {NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
};

static const MunitSuite suite = {"/geo/obj", tests, NULL, 1, MUNIT_SUITE_OPTION_NONE};

int main(int argc, char** argv) { return munit_suite_main(&suite, NULL, argc, argv); }

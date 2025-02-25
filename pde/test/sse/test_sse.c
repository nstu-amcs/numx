#include <munit.h>
#include <numx/pde/sse.h>
#include <stdio.h>

// clang-format off

static void* test_fde_sup(const MunitParameter pps[], void* ctx) {
  (void)ctx;

  char name[16];

  strcpy(name, pps[0]->value);
  strcat(name, ".obj");

  FILE *f = fopen(name, "r");

  if (!f) {
    exit(EXIT_FAILURE);
  }

  struct obj *o = malloc(sizeof(struct obj));

  if (!o) {
    fclose(f);
    exit(EXIT_FAILURE);
  }

  if (obj_new(o, (struct obj_pps) {.with_seg = false, .with_qud = false, .with_hxd = false})) {
    free(o);
    fclose(f);
    exit(EXIT_FAILURE);
  }

  if (obj_get(o, f)) {
    free(o);
    fclose(f);
    exit(EXIT_FAILURE);
  }

  fclose(f);

  return o;
}

static void test_fde_tdn(void* c) {
  free(c);
}

static MunitResult test_fde(const MunitParameter pps[], void* ctx) { 
  (void)pps;

  struct obj* o = (struct obj*)ctx;
  struct vec x;

  if (vec_new(&x, o->v.len))
    return MUNIT_FAIL;

  if ()

  return MUNIT_OK; 
}

static char* cases[] = { "case1", NULL };
static MunitParameterEnum pms[] = {{ "case", cases }, { NULL, NULL }};

static MunitTest tests[] = {
  { "/fde", test_fde, test_fde_sup, test_fde_tdn, MUNIT_TEST_OPTION_NONE, pms},
  { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }
};

static const MunitSuite suite = {
  "/pde/sse", tests, NULL, 1, MUNIT_SUITE_OPTION_NONE,
};

int main(int argc, char** argv) {
  return munit_suite_main(&suite, NULL, argc, argv);
}

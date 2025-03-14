#include <numx/pde/sse.h>

// clang-format off

vfun dat[] = {0};

static void test_cor_sem(void) {
  struct obj obj;

  obj_new(&obj);

  vcut_add(&obj.vtx, (struct vtx){.x = 10, .y = 0, .z = 0, .ctx = 0});
  vcut_add(&obj.vtx, (struct vtx){.x = 20, .y = 0, .z = 0, .ctx = 0});
  vcut_add(&obj.vtx, (struct vtx){.x = 10, .y = 8, .z = 0, .ctx = 0});
  vcut_add(&obj.vtx, (struct vtx){.x = 20, .y = 8, .z = 0, .ctx = 0});
  vcut_add(&obj.vtx, (struct vtx){.x = 10, .y = 0, .z = 1, .ctx = 0});
  vcut_add(&obj.vtx, (struct vtx){.x = 20, .y = 0, .z = 1, .ctx = 0});
  vcut_add(&obj.vtx, (struct vtx){.x = 10, .y = 8, .z = 1, .ctx = 0});
  vcut_add(&obj.vtx, (struct vtx){.x = 20, .y = 8, .z = 1, .ctx = 0});

  qcut_add(&obj.qud, (struct qud){.vtx = {0, 2, 4, 6}, .ctx = &((struct fctx){.cnd = {.type = NEU, .pps.neu.tta = dat[0]}})});
  qcut_add(&obj.qud, (struct qud){.vtx = {2, 3, 6, 7}, .ctx = &((struct fctx){.cnd = {.type = NEU, .pps.neu.tta = dat[0]}})});
  qcut_add(&obj.qud, (struct qud){.vtx = {4, 5, 6, 7}, .ctx = &((struct fctx){.cnd = {.type = NEU, .pps.neu.tta = dat[0]}})});

  qcut_add(&obj.qud, (struct qud){.vtx = {1, 3, 5, 7}, .ctx = &((struct fctx){.cnd = {.type = ROB, .pps.neu.tta = dat[0]}})});
  qcut_add(&obj.qud, (struct qud){.vtx = {2, 3, 6, 7}, .ctx = &((struct fctx){.cnd = {.type = NEU, .pps.neu.tta = dat[0]}})});
  qcut_add(&obj.qud, (struct qud){.vtx = {4, 5, 6, 7}, .ctx = &((struct fctx){.cnd = {.type = NEU, .pps.neu.tta = dat[0]}})});
}

static void test_cor_mem(void) {}

static void test_cor_irr(void) {}

static void test_cor_cef(void) {}

int main(void) {
  test_cor_sem();
  test_cor_mem();
  test_cor_irr();
  test_cor_cef();

  return 0;
}

#include "lin.h"

double fem_std_lin_apx(struct apx_fun_ctx *ctx, vtx_ptr vtx)
{
    switch (ctx->sim->msh->type) {
        case MSH_C2D:
            return fem_std_lin_c2d_apx(ctx, vtx);
        case MSH_C3D:
            return fem_std_lin_c3d_apx(ctx, vtx);
    }

    return 0;
}

int fem_std_lin_asm(struct sim *sim, struct fem_std_ctx *ctx)
{
    switch (sim->msh->type) {
        case MSH_C2D:
            return fem_std_lin_c2d_asm(sim, ctx);
        case MSH_C3D:
            return fem_std_lin_c3d_asm(sim, ctx);
    }

    return -1;
}
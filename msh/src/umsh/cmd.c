#include <numx/msh/umsh.h>
#include <stdio.h>

int main(int argc, char **argv)
{
    (void)argc;

    char *i = argv[1];
    char *o = argv[2];

    struct umsh msh;

    if (umsh_new(&msh)) {
        fprintf(stderr, "failed to initialize mesh\n");
        return 1;
    }

    char *d = i;
    char *p = strrchr(i, '/') + 1;
    char *t = strrchr(i, '.') + 1;

    *(p - 1) = 0;
    *(t - 1) = 0;

    if (!strcmp(t, "tel")) {
        if (umsh_imp_tel(&msh, d, p)) {
            fprintf(stderr, "failed to import mesh\n");
            return 1;
        }
    } else if (!strcmp(t, "grd")) {
        if (umsh_imp_grd(&msh, d, p)) {
            fprintf(stderr, "failed to import mesh\n");
            return 1;
        }
    } else {
        fprintf(stderr, "failed to import mesh: unknown format\n");
        return 1;
    }

    d = o;
    p = strrchr(o, '/') + 1;
    t = strrchr(o, '.') + 1;

    *(p - 1) = 0;
    *(t - 1) = 0;

    if (!strcmp(t, "cgns")) {
        if (umsh_exp_cgns(&msh, d, p)) {
            fprintf(stderr, "failed to export mesh\n");
            return 1;
        }
    } else {
        fprintf(stderr, "failed to export mesh: unknown format\n");
        return 1;
    }

    umsh_cls(&msh);

    return 0;
}
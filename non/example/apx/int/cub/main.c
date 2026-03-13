#include <stdio.h>

#include <numx/non/apx.h>

int main(int argc, char **argv)
{
    (void)argc;
    (void)argv;

    FILE *in = fopen("non/example/apx/int/cub/points.txt", "r");

    if (!in) {
        return EIO;
    }

    int r = 0;
    int n = 0;

    if (fscanf(in, "%d", &n) != 1) {
        fclose(in);
        return EIO;
    }

    struct vec  x;
    struct imtx k;

    if ((r = vec_new(&x, n))) {
        fclose(in);
        return -r;
    }

    if ((r = imtx_new(&k, (struct imtx_pps){.r = 4, .c = n}))) {
        fclose(in);
        return -r;
    }

    for (int i = 0; i < n; ++i) {
        if (fscanf(in, "%lf %lf", &x.dat[i], &k.dat[0][i]) != 2) {
            fclose(in);
            return EIO;
        }
    }

    fclose(in);

    if ((r = apx_int_cub(&x, &k))) {
        return -r;
    }

    FILE *out = fopen("non/example/apx/int/cub/out/spline.txt", "w+");

    if (!out) {
        return EIO;
    }

    for (int i = 0; i < n - 1; ++i) {
        // xi ai bi ci di
        fprintf(out, "%lf %lf %lf %lf %lf\n", x.dat[i], k.dat[0][i], k.dat[1][i], k.dat[2][i],
            k.dat[3][i]);
    }

    fprintf(out, "%lf %lf\n", x.dat[n - 1], k.dat[0][n - 1]);
    fclose(out);

    return 0;
}

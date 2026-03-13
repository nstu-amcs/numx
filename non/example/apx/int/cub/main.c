#include <stdio.h>
#include <errno.h>

#include <non/apx.h>

int main(int argc, char** argv) {
    FILE* in = fopen("non/example/apx/int/cub/points.txt", "r");

    if (!in) {
        return EIO;
    }

    int r = 0;
    int n = 0;

    if (fscanf(in, "%d", &n) != 1) {
        fclose(in);
        return EIO;
    }

    struct vec x;
    struct imtx k;

    if ((r = vec_new(&x, n))) {
        fclose(in);
        return -r;
    }

    if ((r = imtx_new(&k, (struct imtx_pps){.n = 4, .m = n}))) {
        fclose(in);
        return -r;
    }

    for (int i = 0; i < n; ++i) {
        if (fscanf(in, "%lf %lf", &x.dat[i], &k.dat[0][i]) != 2) {
            fclose(in);
            return EIO;
        }
    }

    return 0;
}
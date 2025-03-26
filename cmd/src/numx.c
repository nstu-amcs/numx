#include <stdio.h>
#include <string.h>

#include <numx/cmd.h>

static void man()
{
    printf("usage: numx {cmd}\n\n");
    printf("COMMAND\n");
    printf("\tpde\t\tSolve partial differential equation;\n");
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        man();
        return 1;
    }

    if (!strcmp("pde", argv[1]))
        return pde(argc, argv);

    man();

    return 0;
}

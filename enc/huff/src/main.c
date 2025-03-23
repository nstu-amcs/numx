#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <huff.h>

static char fname[100];

int man();
int parse_opt(char* opt);
int enc(char* name);
int dec(char* name);

int main(int argc, char* argv[argc]) {
  if (argc == 1)
    return man();
  else if (strstr(argv[1], "--") != NULL)
    return parse_opt(argv[1]);
  else if (strstr(argv[1], ".huff") != NULL)
    return dec(argv[1]);
  else
    return enc(argv[1]);
}

int man() {
  printf("huff [options] filename\n");
  printf("\tFilename: name of the file to compress/decompress.\n");
  printf("\tOptions:\n");
  printf("\t\t--help\tprint this help message to the output stream\n");

  return 0;
}

int parse_opt(char* opt) {
  if (strcmp(opt, "--help"))
    return man();
  else {
    printf("Unrecognized option: %s\n", opt);
    return -1;
  }
}

int enc(char* name) {
  memcpy(fname, name, strlen(name));
  strcat(fname, ".huff");

  FILE* src = fopen(name, "rb");
  FILE* dst = fopen(fname, "wb");

  if (src == NULL || dst == NULL) {
    printf("Error. Could not find or load %s\n", src == NULL ? name : fname);
    return -1;
  }

  if (encode(src, dst) != 0) {
    printf("Error. Could not compress %s\n", name);
    return -1;
  }

  return 0;
}

int dec(char* name) {
  memcpy(fname, name, strlen(name) - strlen(strrchr(name, '.')));

  FILE* src = fopen(name, "rb");
  FILE* dst = fopen(fname, "wb");

  if (src == NULL || dst == NULL) {
    printf("Error. Could not find or load %s\n", src == NULL ? name : fname);
    return -1;
  }

  if (decode(src, dst) != 0) {
    printf("Error. Could not decompress %s\n", name);
    return -1;
  }

  return 0;
}

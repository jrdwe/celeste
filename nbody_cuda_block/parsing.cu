#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>

#include "parsing.cuh"

void parseParameters(char** argv, int argc, int* flag, int* b, double* u, double* l, int* t, char** fstring) {
  int ch;
  while (optind < argc) {
    if ((ch = getopt(argc, argv, "b:f:l:u:t")) != -1 )
      switch (ch) {
        case 'b':
          *b = atoi(optarg);
          break;
        case 'f':
          *fstring = optarg;
          *flag = 1;
          break;
        case 'u':
          *u = atof(optarg); 
          break;
        case 'l':
          *l = atof(optarg); 
          break;
        case 't':
          *t = 1;
          break;
        case '?':
          printf("usage: ./nbody <iters> <dt> ((-b <bodies> -l <lowerbound> -u <upperbound>) | -f <filename>) <testmode> \n");

          // exit program on failure from invalid parameters
          exit(EXIT_FAILURE);
        default:
          continue;
      }
  }
}

int countLines(char* filename) {
  FILE* fptr;
  size_t alloclen = 0; 
  char* line      = NULL;

  fptr = fopen(filename, "r");
  if (!fptr) {
    printf("File doesn't exist!\n");
    exit(EXIT_FAILURE);
  }

  int ctr = 0;
  while (getline(&line, &alloclen, fptr) != -1) 
    ctr++;
  
  free(line); 
  line = NULL;

  fclose(fptr);
  return ctr;
}

void readDataFromFile(char* filename, double* bodies) {
  FILE* fptr;
  size_t alloclen = 0;
  char* line      = NULL;
  const char d[2] = ",";
  int ctr         = 0;

  fptr = fopen(filename, "r");
  if (!fptr) {
    printf("File doesn't exist!\n");
    exit(EXIT_FAILURE);
  }

  while (getline(&line, &alloclen, fptr) != -1) {
    char* token;
    token = strtok(line, d);

    while (token) {
      bodies[ctr++] = atof(token);
      token = strtok(NULL, d);
    }
  }

  free(line);
  line = NULL;

  fclose(fptr);
}


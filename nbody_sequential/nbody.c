#define  _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <time.h>

#include "nbody.h"
#include "parsing.h"

#define NUM_STRUCT_FIELDS 7
#define MIN_ARGUMENTS 5
#define STRUCT_FIELD_MASS_INDEX 6
#define ACCEPTABLE_VARIANCE 1e-1
#define DEFAULT_SIZE 50
#define DEFAULT_LOWER_BOUND -50
#define DEFAULT_UPPER_BOUND 50

void randomiseCelestialBodies(double* body, int size, double lower, double upper) {
  // generate a new seed
  time_t timenow;
  srand((unsigned) time(&timenow));

  // generate random doubles within the provided range
  for (int i = 0; i < size; ++i) {
    double s = (double) rand() / RAND_MAX;

    if (i % NUM_STRUCT_FIELDS == STRUCT_FIELD_MASS_INDEX) {
      // ensures mass isn't negative
      body[i] = upper * s;
    } else {
      body[i] = (upper - lower) * s + lower;
    }
  }
}

void step(struct body* bo, int num, float dt) {
  for (int i = 0; i < num; ++i) {
    double ax = 0.0, ay = 0.0, az = 0.0;

    for (int j = 0; j < num; ++j) {
      if (i == j) { continue; }

      double dx   = bo[j].x - bo[i].x;
      double dy   = bo[j].y - bo[i].y;
      double dz   = bo[j].z - bo[i].z;
      double dist = sqrt((double) ((dx * dx) + (dy * dy) + (dz * dz)));
      double magn = GCONST * ((bo[i].mass * bo[j].mass)/(dist * dist));

      ax          += (((dx / dist) * magn) / bo[i].mass);
      ay          += (((dy / dist) * magn) / bo[i].mass);
      az          += (((dz / dist) * magn) / bo[i].mass);
    }

    bo[i].velocity_x += dt * ax;
    bo[i].velocity_y += dt * ay;
    bo[i].velocity_z += dt * az;
  }
}

double energy(struct body* bo, int num) {
  double e = 0.0f;
  for (int i = 0; i < num; ++i) {
    double matsqd = sqrt((double) (bo[i].velocity_x * bo[i].velocity_x)
      + (bo[i].velocity_y * bo[i].velocity_y) 
      + (bo[i].velocity_z * bo[i].velocity_z));
    double first  = (bo[i].mass * matsqd) / 2;
    double second = 0.0f;

    for (int j = i + 1; j < num; ++j) {
      double dx     = bo[j].x - bo[i].x;
      double dy     = bo[j].y - bo[i].y;
      double dz     = bo[j].z - bo[i].z;
      double dist   = sqrtf(((dx * dx) + (dy*dy) + (dz*dz)));
      second       += (GCONST * bo[i].mass * bo[j].mass) / dist;
    }

    e += first - second;
  }
  
  return e;
}

int main(int argc, char** argv) {

  // check required parameters are provided
  if (argc < MIN_ARGUMENTS) {
    printf("usage: ./nbody <iters> <dt> ((-b <bodies> -l <lowerbound> -u <upperbound>) | -f <filename>) <testmode> \n");
    exit(EXIT_FAILURE);
  }

  // variable declaration
  struct body* bodies;
  double ebefore;
  double eafter;
  clock_t tic;
  clock_t toc;
  char* fstring       = "";
  int f               = 0;
  int testmode        = 0;
  int size            = DEFAULT_SIZE;
  double upper        = DEFAULT_UPPER_BOUND;
  double lower        = DEFAULT_LOWER_BOUND;
  const int iteration = atoi(argv[1]);
  const float dt      = atof(argv[2]);

  // parse all flag information
  parseParameters(argv, argc, &f, &size, &upper, &lower, &testmode, &fstring);

  // get size by reading file
  if (f) 
    size = countLines(fstring);
  
  // malloc required size
  bodies = malloc(sizeof(struct body) * size);
  if (!bodies) { 
    fprintf(stderr, "Unable to malloc requested bodies"); 
    exit(EXIT_FAILURE); 
  }

  if (f) {
    // read in data from provided file
    readDataFromFile(fstring, (double*) bodies);
  } else {
    // randomise all bodies
    randomiseCelestialBodies((double*) bodies, NUM_STRUCT_FIELDS * size, lower, upper);
  }

  if (testmode) {
    // calculate the energy within the system prior
    ebefore = energy(bodies, size);
    
    // record the starting time
    toc = clock();
  }
  
  for (int k = 0; k < iteration; ++k) {
    // execute step function 
    step(bodies, size, dt);

    // update the x & y positions for each data point
    for (int i = 0; i < size; ++i) {
      bodies[i].x += bodies[i].velocity_x * dt;
      bodies[i].y += bodies[i].velocity_y * dt;
      bodies[i].z += bodies[i].velocity_z * dt;
    }
  }
 
  if (testmode) {
    // calculate the energy within the system after
    eafter = energy(bodies, size);

    // calculate the finish time
    tic = clock();
    double time_taken = (double)(tic - toc) / CLOCKS_PER_SEC;
   
    // output time taken to execute function
    printf("Time taken for %d iterations with %d bodies was %f\n", iteration, size, time_taken);

    // check if the energy within the system is acceptable
    if (fabs(ebefore - eafter) > ACCEPTABLE_VARIANCE) 
      printf("Energy did not match %f vs %f\n", ebefore, eafter);
  }
  
  // free allocated memory
  free(bodies);
  return EXIT_SUCCESS;
}

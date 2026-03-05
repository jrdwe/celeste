#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <string.h>

#include "nbody.cuh"

int load_from_dat(const char* path, float4** pos_out, float4** vel_out) {
  FILE* file = fopen(path, "r");
  if (!file) {
    printf("[error]: cannot open '%s'\n", path);
    return 0;
  }

  int count = 0;
  char line[256];
  while (fgets(line, sizeof(line), file) != NULL) {
    count++;
  }

  if (count == 0) {
    printf("[error]: no bodies found in '%s'\n", path);
    fclose(file);
    return 0;
  }

  float4* pos = new float4[count];
  float4* vel = new float4[count];

  rewind(file);

  int i = 0;
  while (fgets(line, sizeof(line), file)) {
    double x, y, z, vx, vy, vz, m;
    if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &x, &y, &z, &vx, &vy, &vz, &m) != 7) {
      printf("[error]: unable to parse line.");
      return 0;
    }

    pos[i] = {(float)x, (float)y, (float)z, (float)m};
    vel[i] = {(float)vx, (float)vy, (float)vz, (float)m};

    i++;
  }

  fclose(file);

  *pos_out = pos;
  *vel_out = vel;

  return count;
}

__global__
void step(float4* d_pos, float4* d_vel, float dt, int n) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i >= n)
    return;

  float4 pos_i = d_pos[i];
  float4 force = {0.0f, 0.0f, 0.0f, 0.0f};

  __shared__ float4 shared_pos[BLOCK_SIZE];

  for (int tile = 0; tile < gridDim.x; tile++) {
    int idx = tile * blockDim.x + threadIdx.x;
    shared_pos[threadIdx.x] = (idx < n) ? d_pos[idx] : make_float4(0.0f, 0.0f, 0.0f, 0.0f);

    __syncthreads();

    for (int j = 0; j < blockDim.x; j++) {
      float3 r;
      r.x = shared_pos[j].x - pos_i.x;
      r.y = shared_pos[j].y - pos_i.y;
      r.z = shared_pos[j].z - pos_i.z;

      double dist_sqr      = (r.x*r.x + r.y*r.y + r.z*r.z) + (DAMPEN * DAMPEN);
      double dist_inv      = rsqrtf(dist_sqr);
      double dist_inv_cube = dist_inv * dist_inv * dist_inv;
      double scale         = shared_pos[j].w * dist_inv_cube;

      force.x += r.x * scale;
      force.y += r.y * scale;
      force.z += r.z * scale;
    }
    __syncthreads();
  }

  d_vel[i].x += dt * force.x;
  d_vel[i].y += dt * force.y;
  d_vel[i].z += dt * force.z;
}

void render_file(unsigned char* image, float4* pos, int size, int count, int max_range) {
  memset(image, 0, WIDTH * HEIGHT * 3);

  const float scale = (WIDTH / 2.0f) / (max_range * 1.2f); // 20% buffer
  for (int i = 0; i < size; i++) {
    int x = (int)(WIDTH  / 2.0f + pos[i].x * scale);
    int y = (int)(HEIGHT / 2.0f + pos[i].y * scale);

    float px = x - WIDTH  / 2.0f;
    float py = y - HEIGHT / 2.0f;

    // distance from center point
    float dist = sqrtf(px*px + py*py) / (WIDTH / 2.0f);

    // exp scale as we move away from the middle
    float brightness = 255.0f * expf(-1.0f * pow(dist, 0.8f));
    unsigned char b = (unsigned char) fmaxf(0.0f, brightness);

    if (x > DOTSIZE && x < WIDTH - DOTSIZE && y > DOTSIZE && y < HEIGHT - DOTSIZE) {
      for (int dx = -DOTSIZE / 2; dx < DOTSIZE / 2; dx++) {
        for (int dy = -DOTSIZE / 2; dy < DOTSIZE / 2; dy++) {
          int idx = 3 * ((y + dy) * WIDTH + (x + dx));

          image[idx]     = b; // R
          image[idx + 1] = b; // G 
          image[idx + 2] = b; // B
        }
      }
    }
  }

  char path[256];
  sprintf(path, "frames/frame_%07d.ppm", count);

  FILE* fp = fopen(path, "wb");
  if (!fp) {
    printf("Error: could not open %s for writing\n", path);
    return;
  }

  fprintf(fp, "P6\n%d %d\n255\n", WIDTH, HEIGHT);
  fwrite(image, 3, WIDTH * HEIGHT, fp);

  fclose(fp);
}

int main(int argc, char** argv) {
  if (argc < 2) {
    printf("[usage]: ./celeste <file>\n");
    exit(EXIT_FAILURE);
  }

  time_t timenow;
  srand((unsigned) time(&timenow));

  float4* pos = nullptr;
  float4* vel = nullptr;

  // load file from the provided dat
  int size = load_from_dat(argv[1], &pos, &vel);
  if (size == 0) {
    printf("[error]: failure reading from dat file.\n");
    return EXIT_FAILURE;
  }

  // used to ensure appropriate scale of display
  float max_range = 0.0f;
  for (int i = 0; i < size; i++) {
      float e = fmaxf(fabsf(pos[i].x), fabsf(pos[i].y));
      if (e > max_range)
        max_range = e;
  }

  // create new folder for ppm files
  if (mkdir("frames", 0777) != 0) {
    printf("[error]: could not create frames directory.\n");
    return EXIT_FAILURE;
  }

  // allocate memory in gpu
  float4 *d_pos, *d_vel;
  if (cudaMalloc(&d_pos, size * sizeof(float4)) != cudaSuccess) {
    printf("[error]: cudaMalloc d_pos failed: %s\n", cudaGetErrorString(cudaGetLastError()));
    return EXIT_FAILURE;
  }

  if (cudaMalloc(&d_vel, size * sizeof(float4)) != cudaSuccess) {
    printf("[error]: cudaMalloc d_vel failed: %s\n", cudaGetErrorString(cudaGetLastError()));
    return EXIT_FAILURE;
  }

  // allocate buffer for ppm frame
  unsigned char* image = (unsigned char*) malloc(WIDTH * HEIGHT * 3);
  if (!image) {
    printf("[error]: failed to allocate image buffer\n");
    return EXIT_FAILURE;
  }

  int nblocks = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;

  clock_t toc = clock();
  for (int k = 0; k < ITERATIONS; ++k) {
    cudaMemcpy(d_pos, pos, size * sizeof(float4), cudaMemcpyHostToDevice);
    cudaMemcpy(d_vel, vel, size * sizeof(float4), cudaMemcpyHostToDevice);

    step<<<nblocks, BLOCK_SIZE>>>(d_pos, d_vel, DT, size);

    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
      printf("[error]: CUDA kernel error on frame %d: %s\n", k, cudaGetErrorString(error));
      break;
    }

    error = cudaDeviceSynchronize();
    if (error != cudaSuccess) {
      printf("[error]: CUDA sync error on frame %d: %s\n", k, cudaGetErrorString(error));
      break;
    }

    cudaMemcpy(pos, d_pos, size * sizeof(float4), cudaMemcpyDeviceToHost);
    cudaMemcpy(vel, d_vel, size * sizeof(float4), cudaMemcpyDeviceToHost);

    for (int i = 1; i < size; i++) {
      pos[i].x += vel[i].x * DT;
      pos[i].y += vel[i].y * DT;
      pos[i].z += vel[i].z * DT;
    }

    if (k % 20 == 0) // save every 20th frame
      render_file(image, pos, size, k / 20, max_range);

    if (k % 200 == 0)
      printf("[info]: processed %d frames.\n", k);
  }

  printf("[info]: time taken for %d iterations with %d bodies was %f\n",
    ITERATIONS, size, (double)(clock() - toc) / CLOCKS_PER_SEC);

  cudaFree(d_pos);
  cudaFree(d_vel);

  free(image);

  delete[] pos;
  delete[] vel;

  return EXIT_SUCCESS;
}

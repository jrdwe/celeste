#ifndef NBODY_H
#define NBODY_H

#define ITERATIONS    100000
#define DT            0.02
#define DAMPEN        0.05
#define DOTSIZE       4
#define WIDTH         1024
#define HEIGHT        1024
#define BLOCK_SIZE    256

/* func: parse dat file of the format: x,y,z,vx,vy,vz,m
 *
 * return: int
 */
int load_from_dat(const char* path, float4** pos_out, float4** vel_out);

/* func: create the ppm file representing a single frame 
 *
 * return: void
 */
void render_file(unsigned char* image, float4* pos, int size, int count, int max_range);

/* func: advance the simulation with time difference (dt)
 *
 * return: void
 */
__global__ void step(float4* d_pos, float4* d_vel, float dt, int n);

#endif

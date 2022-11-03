#ifndef NBODY_H
#define NBODY_H

#define PI (3.141592653589793)
#define SOLARMASS (4 * PI * PI)
#define NDAYS (365.25)
#define GCONST (6.67e-11)

struct body {
	double x;
	double y;
	double z;
	double velocity_x;
	double velocity_y;
	double velocity_z;
	double mass;
};

/* func: generate random double values 
 *       for the celestial bodies
 *
 * return: void
 */
void randomiseCelestialBodies(double* body, int size, double lower, double upper);

/* func: advance the simulation with 
 *       time difference (dt)
 *
 * return: void
 */
__global__
void step(struct body* bodies, int num, float dt);

/* func: calculate the total energy of 
 *       the system
 *
 * return: double
 */
double energy(struct body* bodies, int num);

#endif

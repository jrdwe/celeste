
# n-body physics model

An n-body simulation is an astrophysics model of orbiting bodies within our solar system and galaxy.

This code-base consists of three different implementations of the n-body simulation. A sequential version, an unoptimised CUDA version, and an optimised CUDA version.

Each implementation is designed to run with either an input to generate X random points between two values or receive a file specifying the points. 

The file must be structured with each line in the following form: <x>, <y>, <z>, <velocity-x>, <velocity-y>, <velocity-z>, <mass>

Command line parameters are defined as: ./nbody <iters> <dt> ((-b <bodies> -l <lowerbound> -u <upperbound>) | -f <filename>) <testmode>

# Example executions

Run the sequential approach with 100 iterations, 0.8 dt, and 1000 bodies. Testing is enabled.
./sequential 100 0.8 -b 1000 -t

Run the sequential approach with 100 iterations, 0.8 dt, and 1000 bodies between 100 & -100
./sequential 100 0.8 -b 1000 -u 100 -l -100

Run the cache-block approach with 500 iterations, 0.4 dt, and 100,000 bodies between 2000 & -2000
./cuda_block 500 0.4 -b 100000 -u 2000 -l -2000

Run the naive cuda approach with 250 iterations, 0.2 dt, and a test dataset. Testing is enabled.
./cuda_naive 250 0.2 -f test_datasets/test_dataset1000.txt -t



CC=gcc
GPU=nvcc
CFLAGS=-g -std=c11 -Wall -Werror
TARGET=program
.PHONY: clean
all: $(TARGET)
	
program: 
	make cuda_naive
	make cuda_block
	make sequential

cuda_naive: nbody_cuda_naive/nbody.cu nbody_cuda_naive/parsing.cu
	$(GPU) -g $^ -o $@ -lm

cuda_block: nbody_cuda_block/nbody.cu nbody_cuda_block/parsing.cu
	$(GPU) -g $^ -o $@ -lm

sequential: nbody_sequential/nbody.c nbody_sequential/parsing.c
	$(CC) $(CFLAGS) $^ -o $@ -lpthread -lm

clean:
	rm -f sequential
	rm -f cuda_naive
	rm -f cuda_block

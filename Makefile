GPU=nvcc
TARGET=program
.PHONY: clean
all: $(TARGET)
	
program: 
	make celeste

celeste: src/nbody.cu 
	$(GPU) -g -arch=sm_75 $^ -o $@ -lm

clean:
	rm -f celeste && rm -rf frames

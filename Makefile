all:
	gcc -o solver LU_decomposition.c -fopenmp

clean:
	rm -f solver

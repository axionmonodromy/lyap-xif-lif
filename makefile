clv: lyapunovvectors.c utilities.c
	gcc -Wall -pedantic -o clv lyapunovvectors.c utilities.c -lgsl -lgslcblas -lm 
init: init.c
	gcc -Wall -pedantic -o init init.c utilities.c -lgsl -lgslcblas -lm -fopenmp 


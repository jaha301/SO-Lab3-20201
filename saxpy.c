/**
 * @defgroup   SAXPY saxpy
 *
 * @brief      This file implements an iterative saxpy operation
 *
 * @param[in] <-p> {vector size}
 * @param[in] <-s> {seed}
 * @param[in] <-n> {number of threads to create}
 * @param[in] <-i> {maximum itertions}
 *
 * @author     Danny Munera
 * @date       2020
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

#include <pthread.h>
#include <semaphore.h>


void* compute (void *);

typedef struct _param{
	int ini;
	int end;
	int *n_iter_s;
	int max_iters_s;
	int n_threads_s;
	int its;
	int ps;
	double* Xs;
	double as;
	double* Ys;
	double* Y_avgs_s;

	
} param_t;


//semáforo
sem_t mutex;


int main(int argc, char* argv[]){
	// Variables to obtain command line parameters
	unsigned int seed = 1;
  	int p = 10000000;
  	int n_threads = 2;
  	int max_iters = 1000;
  	// Variables to perform SAXPY operation
	double* X;
	double a;
	double* Y;
	double* Y_avgs;
	int i;
	// Variables to get execution time
	struct timeval t_start, t_end;
	double exec_time;

	// Getting input values
	int opt;
	while((opt = getopt(argc, argv, ":p:s:n:i:")) != -1){
		switch(opt){
			case 'p':
			printf("vector size: %s\n", optarg);
			p = strtol(optarg, NULL, 10);
			assert(p > 0 && p <= 2147483647);
			break;
			case 's':
			printf("seed: %s\n", optarg);
			seed = strtol(optarg, NULL, 10);
			break;
			case 'n':
			printf("threads number: %s\n", optarg);
			n_threads = strtol(optarg, NULL, 10);
			break;
			case 'i':
			printf("max. iterations: %s\n", optarg);
			max_iters = strtol(optarg, NULL, 10);
			break;
			case ':':
			printf("option -%c needs a value\n", optopt);
			break;
			case '?':
			fprintf(stderr, "Usage: %s [-p <vector size>] [-s <seed>] [-n <threads number>]\n", argv[0]);
			exit(EXIT_FAILURE);
		}
	}
	srand(seed);

	printf("p = %d, seed = %d, n_threads = %d, max_iters = %d\n", \
	 p, seed, n_threads, max_iters);

	// initializing data
	X = (double*) malloc(sizeof(double) * p);
	Y = (double*) malloc(sizeof(double) * p);
	Y_avgs = (double*) malloc(sizeof(double) * max_iters);

	for(i = 0; i < p; i++){
		X[i] = (double)rand() / RAND_MAX;
		Y[i] = (double)rand() / RAND_MAX;
	}
	for(i = 0; i < max_iters; i++){
		Y_avgs[i] = 0.0;
	}
	a = (double)rand() / RAND_MAX;

#ifdef DEBUG
	printf("vector X= [ ");
	for(i = 0; i < p-1; i++){
		printf("%f, ",X[i]);
	}
	printf("%f ]\n",X[p-1]);

	printf("vector Y= [ ");
	for(i = 0; i < p-1; i++){
		printf("%f, ", Y[i]);
	}
	printf("%f ]\n", Y[p-1]);

	printf("a= %f \n", a);
#endif

	/*
	 *	Function to parallelize
	 */
	gettimeofday(&t_start, NULL);
	//SAXPY iterative SAXPY mfunction
	
	void *status;
	int auxiliar = 0;
	int pos = 0;
	pthread_t threads[n_threads]; 
	param_t parm[n_threads]; 

	sem_init(&mutex,0,1); 
	while (pos < n_threads)
	{
	parm[pos].ini = (p/n_threads)*pos; 
        	parm[pos].end = (p/n_threads)*(pos+1);
		
	if(pos==n_threads-1){parm[pos].end =p;} //si p no es divisible
		
	parm[pos].Xs = X;
	parm[pos].as = a;
	parm[pos].Ys = Y;
		parm[pos].ps = p;
	//	parm[pos].n_iter_s = 0;
		parm[pos].max_iters_s = max_iters;
		parm[pos].Y_avgs_s = Y_avgs;
		
	pthread_create(&threads[pos], NULL, &compute, &parm[pos]); //creador de los hilos
	auxiliar++;
	pos++;}

	pos=0;
	while (pos < n_threads)
	{
		pthread_join(threads[pos], &status); //join de los hilos
		pos++;
	}
	gettimeofday(&t_end, NULL);

#ifdef DEBUG
	printf("RES: final vector Y= [ ");
	for(i = 0; i < p-1; i++){
		printf("%f, ", Y[i]);
	}
	printf("%f ]\n", Y[p-1]);
#endif

	// Computing execution time
	exec_time = (t_end.tv_sec - t_start.tv_sec) * 1000.0;  // sec to ms
	exec_time += (t_end.tv_usec - t_start.tv_usec) / 1000.0; // us to ms
	printf("Execution time: %f ms \n", exec_time);
	printf("Last 3 values of Y: %f, %f, %f \n", Y[p-3], Y[p-2], Y[p-1]);
	printf("Last 3 values of Y_avgs: %f, %f, %f \n", Y_avgs[max_iters-3], Y_avgs[max_iters-2], Y_avgs[max_iters-1]);
	return 0;
}


void* compute (void *arg){ 
	param_t* par =  (param_t *) arg;
	int ini = par->ini;
	int end = par->end;
	double* Xl = par->Xs;
	double al = par->as;
	double* Yl = par->Ys;
	double* Y_avgs_l = par->Y_avgs_s;
	int itl = par->its;
	int pl = par->ps;
	double Yacc =0.0;	
	int i;
	int max_iters = par->max_iters_s;
	itl = 0;

	while(itl < max_iters){
	Yacc = 0.0; 
		for(i = ini; i < end; i++){
			Yl[i] = Yl[i] + al * Xl[i];
			Yacc += Yl[i];
		}
		
		sem_wait(&mutex); 
		Y_avgs_l[itl] += Yacc / pl; 
		sem_post(&mutex);
		itl++;
		}
	return NULL;
}






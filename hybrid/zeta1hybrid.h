#ifndef ZETA1_HYBRID_H
#define ZETA1_HYBRID_H

double zeta1_hybrid(int n_threads_omp, int iterations, int size, int rank);
void compute(int n, double *values, int n_threads_omp);
void scatter(double *local_values, double *values, int local_count);
double calculate_my_sum(double *local_values, int count);

double zeta1_hybrid(int n_threads_omp, int iterations, int size, int rank) {
	// find local number of iterations
    int local_count = iterations / size + 1; // accomodate for rounding errors
    // setup arrays
    double *values;
    if(rank == 0){
        // allocate enough space for the scatter operation (May this produce error if we do not set values to 0?)
        values = calloc((local_count * size) + size, sizeof(double));
        compute(iterations, values, n_threads_omp);
    }
    double *local_values = calloc(local_count, sizeof(double));
    // partition and send the data to other processors
    scatter(local_values, values, local_count);
    // calculate local sum and reduce them to rank 0
    double my_sum = calculate_my_sum(local_values, local_count);
    double sum = 0;
    MPI_Reduce(&my_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // free memory
    if(rank == 0){
        free (values);
    }
    free(local_values);
    // calculate final result and return
    double result = sqrt(sum*6);
    return result;
}

void compute(int n, double *values, int n_threads_omp) {
    // calculate values using OpenMP
    #pragma omp parallel for num_threads(n_threads_omp)
    for (int i=1; i<=n; i++) {
        values[i-1] = 1.0/((double)i*(double)i);
    }
}

void scatter(double *local_values, double *values, int local_count) {
    // no need for explicit synchronisation because of implicit barrier
    MPI_Scatter(values, local_count, MPI_DOUBLE, local_values, local_count, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
} 

double calculate_my_sum(double *local_values, int count) {
	double sum = 0;
    for (int i=0; i<count; i++) {
        sum += local_values[i];
    }
    return sum;
}

#endif

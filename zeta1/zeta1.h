#ifndef ZETA0_H
#define ZETA0_H

double mpi_zeta1(int iterations, int size, int rank);
void compute(int n, double *values);
void scatter(double *local_values, double *values, int local_count);
double calculate_my_sum(double *local_values, int count);

double mpi_zeta1(int iterations, int size, int rank) {
	// find local number of iterations
    int local_count = iterations / size;
    // setup arrays
    double *values;
    if(rank == 0){
        values = calloc(iterations, sizeof(double));
        compute(iterations, values); // &values ?
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

void compute(int n, double *values) {
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

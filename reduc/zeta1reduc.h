#ifndef ZETA1_REDUC_H
#define ZETA1_REDUC_H

double zeta1reduc(int iterations, int size, int rank);
void compute(int n, double *values);
void scatter(double *local_values, double *values, int local_count);
double calculate_my_sum(double *local_values, int count);
double reduc(int size, int rank, double my_sum);

double zeta1reduc(int iterations, int size, int rank) {
	// find local number of iterations
    int local_count = iterations / size + 1; // accomodate for rounding errors
    // setup arrays
    double *values;
    if(rank == 0){
        // allocate enough space for the scatter operation (May this produce error if we do not set values to 0?)
        values = calloc((local_count * size) + size, sizeof(double));
        compute(iterations, values); // &values ?
    }
    double *local_values = calloc(local_count, sizeof(double));
    // partition and send the data to other processors
    scatter(local_values, values, local_count);
    // calculate local sum and do all reduce
    double my_sum = calculate_my_sum(local_values, local_count);
    double global_sum = reduc(size, rank, my_sum);
    // free memory
    if(rank == 0){
        free (values);
    }
    free(local_values);
    // calculate final result and return
    double result = sqrt(global_sum*6);
    return result;
}

double reduc(int size, int rank, double my_sum) {
    double sum = 0;
    MPI_Allreduce(&my_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sum;
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

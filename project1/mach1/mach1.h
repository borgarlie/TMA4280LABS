#ifndef MACH1_H
#define MACH1_H

double mpi_mach1(int iterations, int size, int rank);
void compute(int n, double *values);
void scatter(double *local_values, double *values, int local_count);
double calculate_my_sum(double *local_values, int count);

double mpi_mach1(int iterations, int size, int rank) {
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
    // calculate local sum and reduce them to rank 0
    double my_sum = calculate_my_sum(local_values, local_count);
    double sum = 0;
    MPI_Reduce(&my_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // free memory
    if(rank == 0){
        free (values);
    }
    free(local_values);
    // return the sum
    return sum;
}

void compute(int n, double *values) {
    for (int i=0; i<n; i++) {
        double c = 1.0;
        if (i % 2 != 0) c = -1.0;
        // calculate first part of the equation
        double x = 0.2;
        double x1 = pow(x, 2*i+1);
        double x2 = 2*i+1;
        double arctanPart1 = (c*x1)/x2;
        // calculate the second part
        double x_2 = (double)1/(double)239;
        double x1_2 = pow(x_2, 2*i+1);
        double x2_2 = 2*i+1;
        double arctanPart2 = (c*x1_2)/x2_2;
        // subtract, scale and put in vector
        values[i] = 16*arctanPart1 - 4*arctanPart2;
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

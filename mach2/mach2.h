#ifndef MACH2_H
#define MACH2_H

double arctan ( int n_threads, double x, int n ) {
    double sum = 0;

    #pragma omp parallel for num_threads(n_threads) reduction(+:sum)
    for (int i = 0; i <= n; i++) {
        double x1 = pow(x, 2*i+1);
        double x2 = 2*i+1;
        double c = 1.0;
        if (i % 2 != 0) c = -1.0;
        sum += (c*x1)/x2;
    }

    return sum;
}

double mach2 ( int n_threads, int n ) {
    double arctan1 = arctan(n_threads, 0.2, n);
    double x = (double)1/(double)239;
    double arctan2 = arctan(n_threads, x, n);
    return 16*arctan1 - 4*arctan2;
}

#endif
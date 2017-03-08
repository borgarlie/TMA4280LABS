#ifndef ZETA2_H
#define ZETA2_H

double zeta2 (int nthreads, int n) {

    double sum = 0;

    #pragma omp parallel for num_threads(nthreads) reduction(+:sum)
    for (int i=1; i<=n; i++) {
        sum += 1.0/((double)i*(double)i);
    }

    double result = sqrt(sum*6);

    return result;
}

#endif
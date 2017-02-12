#ifndef ZETA0_H
#define ZETA0_H

double zeta0 (int n) {

    double sum = 0;

    for (double i=1; i<=n; i++) {
        sum += 1/(i*i);
    }

    double result = sqrt(sum*6);

    return result;
}

#endif
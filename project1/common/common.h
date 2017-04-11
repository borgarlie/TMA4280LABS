#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

#ifndef COMMON_H
#define COMMON_H

int isPowTwo(int x) {
	while (((x % 2) == 0) && x > 1) x /= 2;
 	return (x == 1);
}

#endif
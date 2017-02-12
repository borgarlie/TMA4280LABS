#ifndef COMMON_H
#define COMMON_H

int isPowTwo(int x) {
	while (((x % 2) == 0) && x > 1) x /= 2;
 	return (x == 1);
}

#endif
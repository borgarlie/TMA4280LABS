#ifndef PI
#define PI (3.14159265358979323846)
#endif

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#ifndef COMMON_H
#define COMMON_H

int isPowTwo(int x) {
	while (((x % 2) == 0) && x > 1) x /= 2;
 	return (x == 1);
}

#endif
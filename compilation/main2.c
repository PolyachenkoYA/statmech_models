//#include <stdio.h>

#define DEBUG 1

int main()
{
	#if DEBUG
		int x = 1;
	#else
		int x = 2;
	#endif
	return 0;
}

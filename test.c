#include <stdio.h>

int main()
{
	int test = 0;
	if (test >= 2)
		printf("test greater than 2\n");
	else if (test > 0 && test < 2)
		printf("test between 0 and 2\n");

	return 0;
}

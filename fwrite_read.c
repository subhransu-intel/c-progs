#include<stdio.h>

int main() {
	int buffer;
	int store = 0x1234;
	/* Creating a file and storing an int value */
	FILE * stream;
	int count;

	stream = fopen("filewb.bin", "wb");
	fwrite(&store, sizeof(int), 1, stream);
	fclose(stream);

	// Reading value from file
	stream = fopen("filewb.bin", "rb");
	count = fread(&buffer, sizeof(char), 8, stream);

	printf("count: %d\n",count);
	fclose(stream);
	return(0);
}

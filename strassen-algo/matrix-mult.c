/*
 * Reference (https://en.wikipedia.org/wiki/Strassen_algorithm)
 *
 * | C00 C01 |   | A00 A01| |B00 B01|
 * |         | = |        | |       |
 * | C10 C11 |   | A10 A11| |B10 B11|
 *
 *
 *      | C00 C01 |    | M1 + M4 - M5 + M7      M3 + M5           |
 * C =  |         | =  |                                          |
 *      | C10 C11 |    | M2 + M4	        M1 - M2 + M3 + M6 |
 *
 * where 
 *	M1 = (A00 + A11)(B00 + B11)
 *	M2 = (A10 + A11)B00
 * 	M3 = A00(B01 - B11)
 *	M4 = A11(B10 - B00)
 *	M5 = (A00 + A01)B11
 *	M6 = (A10 - A00)(B00 + B01)
 *	M7 = (A01 - A11)(B10 + B11)
 *
 *
 * Multiplying two 4 x 4 matrix:
 * 	= multiplying 8 2 x 2 matrices
 *
 * Similarly multiplying two N x N matrix:
 * 	= multiplying 8 N/2 x N/2 matrices.
 *
 */

/*
 * Some assumption to keep the program simple:
 *	a. Mathematical operations are limited to int type only.
 *	b. Returns error for any result overflow.
 *	c. Only +ve value for matrix elements assumed.
 *	d. Assumes the matrix is entered in n x n format only. Matrix is
 *	   entered in file a.txt and b.txt.
 *
 * Program can be modified to enter any number of r x c matrix.
 * The matrix then can be padded with 0s to make it n x n matrix which then
 * can use strassen algo.
 *
 * The two matrices to be multiplied can be generated internally or entered
 * through files a.txt and b.txt. matrix A is read from a.txt and B from
 * b.txt
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>

#define DEBUG 0

#if(DEBUG)
#define print_debug(format, ...)		\
do {						\
	printf(format, ##__VA_ARGS__);			\
}while(0)
#else
#define print_debug(format, ...)		\
{						\
	if (0)					\
		printf(format, ##__VA_ARGS__);		\
}
#endif


#define NUM_ELEMS 16

struct matrix {
	int m[NUM_ELEMS][NUM_ELEMS];
	int i;
	int j;
};

void check_overflow(int a, int b, bool add, bool mult)
{
	int s;
	if (add) {
		s = a + b;
		if ((a > 0 && b > 0) && (s < 0)) {
			printf("Addition overflow for a = %d b = %d\n", a, b);
			exit(EXIT_FAILURE);
		} else if ((a < 0 && b < 0) && (s > 0)) {
			printf("Addition overflow for a = %d b = %d\n", a, b);
			exit(EXIT_FAILURE);
		}
	}

	if (mult) {
		s = a * b;

		if (a != 0 && (s/a != b)) {
			printf("multiplication overflow for a = %d b = %d\n", a, b);
			exit(EXIT_FAILURE);
		} else if (b != 0 && (s/b != a)) {
			printf("multiplication overflow for a = %d b = %d\n", a, b);
			exit(EXIT_FAILURE);
		}
	}
}

/* Copy matrix element from m2 to m1 */
void copy_elems_to_quad(struct matrix *m1, struct matrix *m2, int num_elems)
{
	int r, c;
	for (r = 0; r < num_elems; r ++)
		for (c = 0; c < num_elems; c++)
			m1->m[m1->i + r][m1->j + c] = m2->m[m1->i + r][m1->j + c];
}

struct matrix add(struct matrix a, struct matrix b, int n)
{
	struct matrix m;
	int r, c;

	m.i = a.i;
	m.j = a.j;

	print_debug("In add: i= %d j = %d\n", m.i, m.j);
	for (r = 0; r < n; r++) {
		for(c = 0; c < n; c++) {
			check_overflow(a.m[a.i + r][a.j + c], b.m[b.i + r][b.j + c], true, false);
			m.m[m.i + r][m.j + c] = a.m[a.i + r][a.j + c] + b.m[b.i + r][b.j + c];
			print_debug("%d ", m.m[m.i + r][m.j + c]);
		}
		print_debug("\n");
	}
	print_debug("\n");

	return m;
}

struct matrix sub(struct matrix a, struct matrix b, int n)
{
	struct matrix m;
	int r, c;

	m.i = a.i;
	m.j = a.j;

	print_debug("In sub\n");
	for (r = 0; r < n; r++) {
		for(c = 0; c < n; c++) {
			check_overflow(a.m[a.i + r][a.j + c],  -(b.m[b.i + r][b.j + c]), true, false);
			m.m[m.i+r][m.j+c] = a.m[a.i+r][a.j+c] - b.m[b.i+r][b.j+c];
			print_debug("%d ", m.m[m.i + r][m.j + c]);
		}
		print_debug("\n");
	}
	print_debug("\n");

	return m;
}

/**
 * strassen_matrix_multiply: strassen's algo for matrix multiplication.
 * @m: structure holding a,b and c matrix where c = a x b
 * @n: number of row/column for each matrix
 */
struct matrix strassen_matrix_multiply(struct matrix a, struct matrix b, int n)
{
	struct matrix A00, A01, A10, A11; /* Four quadrant of matrix a */
	struct matrix B00, B01, B10, B11; /* Four quadrant of matrix b */
	struct matrix M1, M2, M3, M4, M5, M6, M7;
	struct matrix Q1, Q2, Q3, Q4;
	struct matrix res;
	int r, c, i, j; 

	if (n == 2) {
		int m1, m2, m3, m4, m5, m6, m7;
		struct matrix c;

#if DEBUG
		print_debug("Input for 2 x 2 matrix multiplication:\n");
		for (i = 0; i < 2; i++) {
			for(j = 0; j < 2; j++)
				print_debug("%d ", a.m[a.i + i][a.j + j]);
			print_debug("\n");
		}

		for (i = 0; i < 2; i++) {
			for(j = 0; j < 2; j++)
				print_debug("%d ", b.m[b.i + i][b.j + j]);
			print_debug("\n");
		}
#endif

		/*
		 * Can be optimized by removing repeat addition/multiplications.
		 * Kept it for easier read.
		 */

		/* Check overflow for expressions in m1 */
		check_overflow(a.m[a.i][a.j], a.m[a.i+1][a.j+1], true, false);
		check_overflow(b.m[b.i][b.j], b.m[b.i+1][b.j+1], true, false);
		check_overflow((a.m[a.i][a.j] + a.m[a.i+1][a.j+1]),
				(b.m[b.i][b.j] + b.m[b.i+1][b.j+1]), false, true);

		/* Check overflow for expressions in m2 */
		check_overflow(a.m[a.i+1][a.j], a.m[a.i+1][a.j+1], true, false);
		check_overflow((a.m[a.i+1][a.j] + a.m[a.i+1][a.j+1]), b.m[b.i][b.j],
								false, true);

		/* Check overflow for expressions in m3 */
		check_overflow(b.m[b.i][b.j+1], -(b.m[b.i+1][b.j+1]), true, false);
		check_overflow(a.m[a.i][a.j],
			(b.m[b.i][b.j+1] - b.m[b.i+1][b.j+1]), false, true);

		/* Check overflow for expressions in m4 */
		check_overflow(b.m[b.i+1][b.j], -(b.m[b.i][b.j]), true, false);
		check_overflow(a.m[a.i+1][a.j+1], (b.m[b.i+1][b.j] - b.m[b.i][b.j]),
							false, true);

		/* Check overflow for expressions in m5 */
		check_overflow(a.m[a.i][a.j], a.m[a.i][a.j+1], true, false);
		check_overflow((a.m[a.i][a.j] + a.m[a.i][a.j+1]),
				b.m[b.i+1][b.j+1], false, true);

		/* Check overflow for expressions in m6 */
		check_overflow(a.m[a.i+1][a.j], -(a.m[a.i][a.j]), true, false);
		check_overflow(b.m[b.i][b.j], b.m[b.i][b.j+1], true, false);
		check_overflow((a.m[a.i+1][a.j] - a.m[a.i][a.j]),
				(b.m[b.i][b.j] + b.m[b.i][b.j+1]), false, true);

		/* Check overflow for expressions in m7 */
		check_overflow(a.m[a.i][a.j+1], -(a.m[a.i+1][a.j+1]), true, false);
		check_overflow(b.m[b.i+1][b.j], b.m[b.i+1][b.j+1], true, false);
		check_overflow((a.m[a.i][a.j+1] - a.m[a.i+1][a.j+1]),
				(b.m[b.i+1][b.j] + b.m[b.i+1][b.j+1]),
				false, true);

		/* Strassen's multiplications for a 2 x 2 matrix */
		m1 = (a.m[a.i][a.j] + a.m[a.i+1][a.j+1]) *
				(b.m[b.i][b.j] + b.m[b.i+1][b.j+1]);
		m2 = (a.m[a.i+1][a.j] + a.m[a.i+1][a.j+1]) *
					b.m[b.i][b.j];
		m3 = a.m[a.i][a.j] *
			(b.m[b.i][b.j+1] - b.m[b.i+1][b.j+1]);
		m4 = a.m[a.i+1][a.j+1] *
			(b.m[b.i+1][b.j] - b.m[b.i][b.j]);
		m5 = (a.m[a.i][a.j] + a.m[a.i][a.j+1]) *
			b.m[b.i+1][b.j+1];
		m6 = (a.m[a.i+1][a.j] - a.m[a.i][a.j]) *
			(b.m[b.i][b.j] + b.m[b.i][b.j+1]);
		m7 = (a.m[a.i][a.j+1] - a.m[a.i+1][a.j+1]) *
			(b.m[b.i+1][b.j] + b.m[b.i+1][b.j+1]);

		c.i = a.i;
		c.j = a.j;
	
		/* Check overflow for expressions in c.m[c.i][c.j] */
		check_overflow(m1, m4, true, false);
		check_overflow((m1 + m4), -(m5), true, false);
		check_overflow((m1 + m4 -m5), m7, true, false);

		/* Check overflow for expressions in c.m[c.i][c.j + 1] */
		check_overflow(m3, m5, true, false);

		/* Check overflow for expressions in c.m[c.i + 1][c.j] */
		check_overflow(m2, m4, true, false);

		/* Check overflow for expressions in c.m[c.i + 1][c.j + 1] */
		check_overflow(m1, -(m2), true, false);
		check_overflow((m1 - m2), m3, true, false);
		check_overflow((m1 - m2 + m3), m6, true, false);

		c.m[c.i][c.j] = m1 + m4 - m5 + m7;
		c.m[c.i][c.j+1] = m3 + m5;
		c.m[c.i+1][c.j] = m2 + m4;
		c.m[c.i+1][c.j+1] = m1 - m2 + m3 + m6;

		print_debug("Result 2 x 2 matrix with r = %d c = %d\n", c.i, c.j);
		for (i = 0; i < 2; i++) {
			for(j = 0; j < 2; j++)
				print_debug("%d ", c.m[c.i + i][c.j + j]);
			print_debug("\n");
		}

		return c;
	}

	A00.i = a.i;		A00.j = a.j;
	copy_elems_to_quad(&A00, &a, n/2);

	A01.i = a.i;		A01.j = a.j + (n/2);
	copy_elems_to_quad(&A01, &a, n/2);

	A10.i = a.i + (n/2);	A10.j = a.j;
	copy_elems_to_quad(&A10, &a, n/2);

	A11.i = a.i + (n/2);	A11.j = a.j + (n/2);
	copy_elems_to_quad(&A11, &a, n/2);

	B00.i = b.i;		B00.j = b.j;
	copy_elems_to_quad(&B00, &b, n/2);

	B01.i = b.i; 		B01.j = b.j + (n/2);
	copy_elems_to_quad(&B01, &b, n/2);

	B10.i = b.i + (n/2);	B10.j = b.j;
	copy_elems_to_quad(&B10, &b, n/2);

	B11.i = b.i + (n/2);	B11.j = b.j + (n/2);
	copy_elems_to_quad(&B11, &b, n/2);

	print_debug("\nCalculate M1\n");
	M1 = strassen_matrix_multiply(add(A00, A11, n/2), add(B00, B11, n/2), n/2);

	print_debug("\nCalculate M2\n");
	M2 = strassen_matrix_multiply(add(A10, A11, n/2), B00, n/2);

	print_debug("\nCalculate M3\n");
	M3 = strassen_matrix_multiply(A00, sub(B01, B11, n/2), n/2);

	print_debug("\nCalculate M4\n");
	M4 = strassen_matrix_multiply(A11, sub(B10, B00, n/2), n/2);

	print_debug("\nCalculate M5\n");
	M5 = strassen_matrix_multiply(add(A00, A01, n/2), B11, n/2);

	print_debug("\nCalculate M6\n");
	M6 = strassen_matrix_multiply(sub(A10, A00, n/2), add(B00, B01, n/2), n/2);

	print_debug("\nCalculate M7\n");
	M7 = strassen_matrix_multiply(sub(A01, A11, n/2), add(B10, B11, n/2), n/2);

	Q1 = add(sub(add(M1, M4, n/2), M5, n/2), M7, n/2);
	Q2 = add(M3, M5, n/2);
	Q3 = add(M2, M4, n/2);
	Q4 = add(add(sub(M1, M2, n/2), M3, n/2), M6, n/2);

	for (r = 0, i = 0; r < n/2; r++, i++)
		for (c = 0, j = 0; c < n/2; c++, j++)
			res.m[r][c] = Q1.m[Q1.i + i][Q1.j + j];

	for (r = 0, i = 0; r < n/2; r++, i++)
		for (c = n/2, j = 0; c < n; c++, j++)
			res.m[r][c] = Q2.m[Q2.i + i][Q2.j + j];

	for (r = n/2, i = 0; r < n; r++, i++)
		for (c = 0, j = 0; c < n/2; c++, j++)
			res.m[r][c] = Q3.m[Q3.i + i][Q3.j + j];

	for (r = n/2, i = 0; r < n; r++, i++)
		for (c = n/2, j = 0; c < n; c++, j++)
			res.m[r][c] = Q4.m[Q4.i + i][Q4.j + j];

	return res;
}

void read_from_file(struct matrix *m1, struct matrix *m2, int n)
{
	int i, j;
	FILE *fp;
	int num_line = 0;
	char line[1000];
	char *token;

	fp = fopen("a.txt", "r");
	if (fp == NULL) {
		printf("a.txt open error\n");
		exit(EXIT_FAILURE);
	}

	/* Parse a.txt to read matrix A */
	printf("Elements for matrix A\n");
	i = 0;
	m1->i = m1->j = 0;
	while (fgets(line, 1000, fp) != NULL) {
		j = 0;
		token = strtok(line, " ");

		while(token) {
			m1->m[m1->i + i][m1->j + j] = atoi(token);
			printf("%d ", m1->m[m1->i + i][m1->j + j]);
			if (m1->m[m1->i + i][m1->j + j] < 0)
				exit(EXIT_FAILURE);
			token = strtok(NULL, " ");
			if (++j == n)
				break;
		}
		printf("\n");
		if (++i == n)
			break;
	}
	fclose(fp);

	fp = fopen("b.txt", "r");
	if (fp == NULL) {
		printf("fp open error\n");
		exit(EXIT_FAILURE);
	}

	/* Parse b.txt to read matrix B */
	printf("Elements for matrix B\n");
	i = 0;
	m2->i = m2->j = 0;
	while (fgets(line, 1000, fp) != NULL) {
		j = 0;
		token = strtok(line, " ");

		while(token) {
			m2->m[m2->i + i][m2->j + j] = atoi(token);
			printf("%d ", m2->m[m2->i + i][m2->j + j]);
			if (m2->m[m2->i + i][m2->j + j] < 0)
				exit(EXIT_FAILURE);
			token = strtok(NULL, " ");
			if (++j == n)
				break;
		}
		printf("\n");
		if (++i == n)
			break;
	}

	fclose(fp);
}

void generate_random(struct matrix *m1, struct matrix *m2, int n)
{
	time_t t;
	int i, j;

	srand((unsigned)time(&t));

	m1->i = m1->j = m2->i = m2->j = 0;
	printf("Elements for matrix A\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			m1->m[i][j] = rand()%100;
			printf("%4d ", m1->m[i][j]);
		}
		printf("\n");
	}

	printf("Elements for matrix B\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			m2->m[i][j] = rand()%101;
			printf("%4d ", m2->m[i][j]);
		}
		printf("\n");
	}

}

void print_help()
{
	printf("This program uses strassen's algorithm to multiply two matrices\n");
	printf("Usage: ./a.out <option>\n");
	printf("Options:\n");
	printf("\t-f: 			Read matrix A and B from files a.txt and b.txt respectively\n");
	printf("\t-r: 			Generate matrix A and B internally using rand()\n");
	printf("\t-n <num_row_col>:	Number of row/col\n");
}

int main(int argc, char *argv[])
{
	struct matrix m1, m2, m3;
	int ret = 0;
	int i, j, k, n;
	int input, help, from_file = 0, random = 0;

	if (argc < 4) {
		print_help();
		exit(EXIT_SUCCESS);
	}

	for (i = 0; i < NUM_ELEMS; i++)
		for (j = 0; j < NUM_ELEMS; j++)
			m1.m[i][j] = m2.m[i][j] = m3.m[i][j] = 0;

	while((input = getopt(argc, argv, "frn:")) != -1) {
		switch(input) {
		case 'f':
			from_file = 1;
			break;
		case 'r':
			random = 1;
			break;
		case 'n':
			n = atoi(optarg);
			if (n > NUM_ELEMS) {
				printf("Input is greater than max array row/col elem size %d\n", NUM_ELEMS);
				exit(EXIT_FAILURE);
			}

			if (n > 16) {
				printf("lets restrict to 16 only for testing.\n");
				exit(EXIT_FAILURE);
			}

			break;
		default:
			printf("Invalid option\n");
			help++;
			break;
		}
	}

	if (help || (optind < argc)) {
		print_help();
		exit(EXIT_SUCCESS);
	}

	if (from_file) {
		read_from_file(&m1, &m2, n);
	} else if(random) {
		generate_random(&m1, &m2, n);
	} else {
		print_help();
		exit(EXIT_SUCCESS);
	}

	m3 = strassen_matrix_multiply(m1, m2, n);

	printf("Result with strassen algo: \n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			printf("%d\t", m3.m[m3.i + i][m3.j + j]);
		printf("\n");
	}

	printf("Result with standard multiplication: \n");
	for (i = 0; i < n ; i++) {
		for (j = 0; j < n ; j++) {
			m3.m[i][j] = 0;
			for (k = 0; k < n; k++) {
				m3.m[i][j] += m1.m[i][k] * m2.m[k][j];
			}
			printf("%d\t", m3.m[i][j]);
		}
		printf("\n");
	}


	return 0;
}

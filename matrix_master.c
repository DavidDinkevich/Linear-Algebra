
// Matrix Master
// David Dinkevich
// Created: 3 July 2020

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdbool.h>

// GENERAL UTILS

double smoothenZero(double n) {
	return fabs(n) < (1e-6) ? 0.0 : n;
}

void printInputPrefix() {
	printf(">> ");
}

int readInt() {
	printInputPrefix();
	int n;
	scanf("%d", &n);
	printf("Read value: %d\n", n);
	return n;
}

// VECTORS

double* createVector(int n) {
	return (double *) malloc(sizeof(double) * n);
}

void setVector(double *v1, double *v2, int n) {
	for (int i = 0; i < n; i++) {
		v1[i] = v2[i];
	}
}

double* copyVector(double *v, int n) {
	double *copy = createVector(n);
	setVector(copy, v, n);
	return copy;
}

void zeroVector(double *v, int n) {
	for (int i = 0; i < n; i++) {
		v[i] = 0;
	}
}

double dot(double *vector1, double *vector2, int n) {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += vector1[i] * vector2[i];
	}
	return sum;
}

double magnitude(double *v, int n) {
	return sqrt(dot(v, v, n));
}

void normalizeVector(double *v, int n) {
	double mag = magnitude(v, n);
	for (int i = 0; i < n; i++) {
		v[i] /= mag;
	}
}


void multVector(double *vector, int n, double scalar) {
	for (int i = 0; i < n; i++) {
		vector[i] *= scalar;
	}
}

/*
	Adds the second vector times the scalar to the first vector
	(first vector is changed)
*/
void addVectors(double *v1, double *v2, int n, double scalar) {
	for (int i = 0; i < n; i++) {
		v1[i] += v2[i] * scalar;
	}
}


// MATRICES


double** createMatrix(int n, int m) {
	double **mat = (double **) malloc(sizeof(double*) * n);
	for (int r = 0; r < n; r++) {
		mat[r] = createVector(m);
	}
	return mat;
}

double** copyMatrix(double **mat, int n, int m) {
	double **copy = createMatrix(n, m);
	for (int r = 0; r < n; r++) {
		setVector(copy[r], mat[r], m);
	}
	return copy;
}

void swapRows(double **matrix, int m, int row1, int row2) {
	for (int i = 0; i < m; i++) {
		double temp = matrix[row2][i];
		matrix[row2][i] = matrix[row1][i];
		matrix[row1][i] = temp;
	}
}

int ref(double **matrix, int n, int m) {
	int numRowSwaps = 0;
	for (int c = 0; c < m; c++) {
		for (int r = c + 1; r < n; r++) {
			if (matrix[c][c] == 0 && matrix[r][c] != 0) {
				swapRows(matrix, m, r, c);
				++numRowSwaps;
			}
			else if (matrix[r][c] != 0) {
				addVectors(matrix[r], matrix[c], m, - matrix[r][c] / matrix[c][c]);
			}
		}
	}
	return numRowSwaps;
}

void rref(double **matrix, int n, int m) {
	ref(matrix, n, m);

	for (int c = fmin(n, m) - 1; c >= 0; c--) {
		if (smoothenZero(matrix[c][c]) == 0) {
			continue;
		}
		multVector(matrix[c], m, 1 / matrix[c][c]);
		for (int r = c - 1; r >= 0; r--) {
			double scalar = - matrix[r][c] / matrix[c][c];
			addVectors(matrix[r], matrix[c], m, scalar);
		}
	}
}

void freeMatrix(double **matrix, int n) {
	for (int r = 0; r < n; r++) {
		free(matrix[r]);
	}
	free(matrix);
}

double diagonalProduct(double **mat, int n) {
	double prod = 1;
	for (int r = 0; r < n; r++) {
		prod *= mat[r][r];
	}
	return prod;
}

double** transpose(double **mat, int n, int m) {
	double **trans = createMatrix(m, n);

	for (int r = 0; r < m; r++) {
		for (int c = 0; c < n; c++) {
			trans[r][c] = mat[c][r];
		}
	}
	return trans;
}

int rank(double **mat, int n, int m) {
	double **copy = copyMatrix(mat, n, m);
	ref(copy, n, m);
	int numPivots = 0;
	for (int c = 0; c < fmin(n, m); c++) {
		if (smoothenZero(copy[c][c]) != 0)
			numPivots++;
	}
	freeMatrix(copy, n);
	return numPivots;
}

double** createMinor(double **mat, int n, int m, int i, int j) {
	double **minor = malloc(sizeof(double*) * (n-1));
	for (int r = 0, offsetR = 0; r < n-1; r++) {
		if (r == i)
			++offsetR;

		minor[r] = createVector(m-1);
		for (int c = 0, offsetC = 0; c < n-1; c++) {
			if (c == j) {
				++offsetC;
			}
			minor[r][c] = mat[r + offsetR][c + offsetC];
		}
	}
	return minor;
}

double rowReductionDet(double **mat, int n) {
	double **copy = copyMatrix(mat, n, n);
	int numSwaps = ref(copy, n, n);
	int sign = numSwaps % 2 == 0 ? 1 : -1;
	double det = sign * diagonalProduct(copy, n);
	freeMatrix(copy, n);
	return det;
}

/*
	Params:
	- vector: vector to project
	- orthoBasis: matrix containing the basis vectors of the space
				  onto which the vector will be projected. Rows
				  of matrix are the basis vectors.
	- n: number of basis vectors (dimension of the space)
	- m: length of vectors
*/
double* proj(double *vector, double **orthoBasis, int n, int m) {
	double *result = createVector(m);
	zeroVector(result, m);

	for (int i = 0; i < n; i++) {
		double scalar = dot(vector, orthoBasis[i], m) / dot(orthoBasis[i], orthoBasis[i], m);
		double *tempProj = copyVector(orthoBasis[i], m);
		addVectors(result, tempProj, m, scalar);
		free(tempProj);
	}
	return result;
}

double** gramSchmidt(double **linIndSet, int n, int m) {
	double **copy = copyMatrix(linIndSet, n, m);
	for (int i = 1; i < n; i++) {
		double *prj = proj(linIndSet[i], copy, i, m);
		addVectors(copy[i], prj, m, -1);
		free(prj);
	}
	return copy;
}


double** getCofactorMatrix(double **mat, int n) {
	double **cof = createMatrix(n, n);
	for (int r = 0; r < n; r++) {
		for (int c = 0; c < n; c++) {
			int sign = (r + c) % 2 == 0 ? 1 : -1;
			double **ijMinor = createMinor(mat, n, n, r, c);
			double det = rowReductionDet(ijMinor, n-1);
			freeMatrix(ijMinor, n-1);
			cof[r][c] = sign * det;
		}
	}
	return cof;
}

double** getAdjoint(double **mat, int n) {
	double **cof = getCofactorMatrix(mat, n);
	double **adj = transpose(cof, n, n);
	freeMatrix(cof, n);
	return adj;
}

double** getNullSpace(double **mat, int n, int m, int *dimNullSpc) {
	int rnk = rank(mat, n, m);
	// If full rank, nullspace is 0
	if (rnk == fmin(n, m)) {
		double **zeroSpace = createMatrix(*dimNullSpc = 1, m);
		zeroVector(zeroSpace[0], m);
		return zeroSpace;
	}

	// Row reduce matrix. All the pivotless columns are located on the right
	// of the matrix
	rref(mat, n, m);
	// Get dimension of null space
	*dimNullSpc = m - rnk;
	// Create matrix that will hold null space basis
	double **nullSpcBasis = createMatrix(*dimNullSpc, m);
	// Transpose the matrix for easier access to pivotless-columns
	// (which will now be rows at the bottom of the matrix)
	double **transp = transpose(mat, n, m);

	// Fill in nullSpcBasis with these pivotless cols, negate them,
	// and put 1 where the pivot should have been
	for (int r = 0; r < *dimNullSpc; r++) {
		// Copy the vector into the null space basis matrix
		setVector(nullSpcBasis[r], transp[r + rnk], fmin(n, m));
		// Zero the remainder of the row in the null space basis matrix
		// (if n < m)
		zeroVector(nullSpcBasis[r] + n, m - n);
		// Negate the vector
		multVector(nullSpcBasis[r], m, -1);
		// Place the missing pivot
		nullSpcBasis[r][r + rnk] = 1;
	}

	freeMatrix(transp, m); // m is num rows of the transpose

	return nullSpcBasis;

}

/*
	PRINTING
*/

void printVector(double *vector, int n, bool withParenthesis) {
	if (withParenthesis)
		printf("(");
	for (int i = 0; i < n; i++) {
		printf("%0.3f", smoothenZero(vector[i]));
		if (i < n-1)
			printf(" ");
	}
	if (withParenthesis)
		printf(")");
}

void printMatrix(double **matrix, int n, int m) {
	for (int i = 0; i < n; i++) {
		printVector(matrix[i], m, false);
		printf("\n");
	}
}

void printRowsAsSet(double **matrix, int n, int m) {
	printf("{\n");
	for (int i = 0; i < n; i++) {
		printf("\t");
		printVector(matrix[i], m, true);
		// if (i < n - 1)
		// 	printf(", ");
		printf("\n");
	}
	printf("}");
}

/*
	DIALOGS
*/

int inputVector(double *v, int n) {
	//printInputPrefix();
	for (int c = 0; c < n; c++) {
		if (1 != scanf("%lf", &v[c])) {
			return false;
		}
	}
	return true;
}

bool inputMatrixDims(int *n, int *m) {
	int res = scanf("%d%d", n, m);
	if (res != 2 || !(*n >= 1 && *m >= 1)) {
		printf("Invalid matrix dimensions.\n");
		return false;	
	}
	return true;
}

int inputMatrix(double **mat, int n, int m) {
	for (int r = 0; r < n; r++) {
		if (!inputVector(mat[r], m))
			return false;
	}
	return true;
}

double** runMatrixInputDialog(int *n, int *m) {
	printf("Input matrix dimensions:\n");
	printInputPrefix();
	if (!inputMatrixDims(n, m))
		return NULL;
	printf("Matrix dimensions: %dx%d\n", *n, *m);

	double **matrix = createMatrix(*n, *m);
	printf("Input matrix:\n");
	if (!inputMatrix(matrix, *n, *m))
		return NULL;
	printf("Read matrix:\n");
	printMatrix(matrix, *n, *m);

	return matrix;
}

void runREFDialog() {
	int n, m;
	double **matrix = runMatrixInputDialog(&n, &m);
	if (matrix == NULL) {
		return;
	}
	printf("ref():\n");
	ref(matrix, n, m);
	printMatrix(matrix, n, m);
	freeMatrix(matrix, n);
}

void runRREFDialog() {
	int n, m;
	double **matrix = runMatrixInputDialog(&n, &m);
	if (matrix == NULL) {
		return;
	}
	printf("rref():\n");
	rref(matrix, n, m);
	printMatrix(matrix, n, m);
	freeMatrix(matrix, n);
}

void runDetDialog() {
	int n, m;
	double **matrix = runMatrixInputDialog(&n, &m);
	if (matrix == NULL) {
		return;
	}
	if (n != m) {
		printf("Must input square matrix.\n");
	} else {
		printf("det(): %lf\n", smoothenZero(rowReductionDet(matrix, n)));
	}
	freeMatrix(matrix, n);
}

void runAdjointDialog() {
	int n, m;
	double **matrix = runMatrixInputDialog(&n, &m);
	if (matrix == NULL) {
		return;
	}
	if (n != m) {
		printf("Must input square matrix.\n");
	} else {
		printf("Adjoint:\n");
		double **adj = getAdjoint(matrix, n);
		printMatrix(adj, n, n);
		freeMatrix(adj, n);
	}
	freeMatrix(matrix, n);
}

void runOrthoProjectionDialog() {
	printf("Input vector length:\n");
	int vecLen = readInt();

	printf("Input the dimension of space to project on:\n");
	int spaceDim = readInt();

	printf("Input an orthogonal basis of the space:\n");
	double **basis = createMatrix(spaceDim, vecLen);
	for (int i = 0; i < spaceDim; i++) {
		inputVector(basis[i], vecLen);
	}

	printf("Read value: \n");
	printRowsAsSet(basis, spaceDim, vecLen);
	printf("\n");

	double *toProject = createVector(vecLen);
	printf("Input vector to project:\n");
	inputVector(toProject, vecLen);
	printf("Read value: ");
	printVector(toProject, vecLen, true);
	printf("\n");

	printf("Projection:\n");
	double *prj = proj(toProject, basis, spaceDim, vecLen);
	printVector(prj, vecLen, true);
	printf("\n");
	
	printf("Orthogonal complement:\n");
	addVectors(toProject, prj, vecLen, -1);
	printVector(toProject, vecLen, true);
	printf("\n");

	free(toProject);
	freeMatrix(basis, spaceDim);
	free(prj);
}

void runGramSchmidtDialog() {
	printf("Input vector length:\n");
	int vecLen = readInt();

	printf("Input the size of the set:\n");
	int setLen = readInt();

	printf("Input a linearly independent set of vectors:\n");
	double **set = createMatrix(setLen, vecLen);
	for (int i = 0; i < setLen; i++) {
		inputVector(set[i], vecLen);
	}

	printf("Read value: \n");
	printRowsAsSet(set, setLen, vecLen);
	printf("\n");

	double **orthogonalized = gramSchmidt(set, setLen, vecLen);
	printf("Orthogonalized:\n");
	printRowsAsSet(orthogonalized, setLen, vecLen);
	printf("\n");

	printf("Normalize vectors? (y/n): ");
	char answer;
	scanf(" %c", &answer);
	
	if (answer == 'y') {
		for (int r = 0; r < setLen; r++) {
			normalizeVector(orthogonalized[r], vecLen);
		}
		printf("Normalized:\n");
		printRowsAsSet(orthogonalized, setLen, vecLen);
		printf("\n");
	}

	freeMatrix(set, setLen);
	freeMatrix(orthogonalized, setLen);
}

void runNullSpaceDialog() {
	int n, m;
	double **mat = runMatrixInputDialog(&n, &m);
	if (mat == NULL)
		return;

	int dimNullSpc;
	double **nullSpace = getNullSpace(mat, n, m, &dimNullSpc);
	printf("The null space is the span of:\n");
	printRowsAsSet(nullSpace, dimNullSpc, m);

	freeMatrix(nullSpace, dimNullSpc);	
	freeMatrix(mat, n);
}

int main() {
	printf("Welcome to Matrix Master!\n");

	int input;

	do {
		printf("\nMAIN MENU\n");
		printf("1. Reduce to row echelon form (REF)\n");
		printf("2. Reduce to reduced row echelon form (RREF)\n");
		printf("3. Determinant\n");
		printf("4. Null Space\n");
		printf("5. Adjoint\n");
		printf("6. Orthogonal Projection\n");
		printf("7. Gram-Schmidt\n");
		printf("0: Quit\n");
		printf("Input:\n");
		input = readInt();

		switch (input) {
			case 1: runREFDialog(); break;
			case 2: runRREFDialog(); break;
			case 3: runDetDialog(); break;
			case 4: runNullSpaceDialog(); break;
			case 5: runAdjointDialog(); break;
			case 6: runOrthoProjectionDialog(); break;
			case 7: runGramSchmidtDialog(); break;
			case 0: break;
			default:
				printf("Invalid input. Try again:\n");
		}

	} while (input != 0);

	printf("Process ended, terminating program.");

	return 0;
}

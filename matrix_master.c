
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

/*
	Get whether all entries in the given vector equal the given
	constant "c"
*/
bool allEntriesEqual(double c, double *vec, int n) {
	for (int i = 0; i < n; i++) {
		if (vec[i] != c)
			return false;
	}
	return true;
}


// MATRICES


double** createMatrix(int n, int m) {
	double **mat = (double **) malloc(sizeof(double*) * n);
	for (int r = 0; r < n; r++) {
		mat[r] = createVector(m);
	}
	return mat;
}

void freeMatrix(double **matrix, int n) {
	for (int r = 0; r < n; r++) {
		free(matrix[r]);
	}
	free(matrix);
}

double **createIdentityMatrix(int n) {
	double **identity = createMatrix(n, n);
	for (int r = 0; r < n; r++) {
		for (int c = 0; c < n; c++) {
			identity[r][c] = (r == c) ? 1 : 0;
		}
	}
	return identity;
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

/*
	Represents the row reduction of an augmented matrix.
	"mat1" is of size (n1 x m1)
	"mat2" is of size (n2 x m2)
	If "mat2" is NULL, only "mat1" will be row reduced.
	PRECONDITION: n1 <= n2
*/
int augmentedREF(double **mat1, int n1, int m1, double **mat2, int n2, int m2) {
	if (mat2 != NULL && n1 > n2) {
		printf("ERROR: invalid dimensions for augmented matrix.\n");
		return 0;
	}
	int numRowSwaps = 0;

	int pivotR = 0, pivotC = 0;
	while (pivotR < n1 && pivotC < m1) {		
		// NO PIVOT AT (pivotR, pivotC)
		if (smoothenZero(mat1[pivotR][pivotC]) == 0) {
			// Try to find row with coef in same col and swap
			bool foundAltPivot = false;
			for (int r = pivotR + 1; r < n1; r++) {
				if (mat1[r][pivotC] != 0) {
					foundAltPivot = true;
					// mat1
					swapRows(mat1, m1, pivotR, r);
					// mat2
					if (mat2 != NULL)
						swapRows(mat2, m2, pivotR, r);
					++numRowSwaps;
					break;
				}
			}
			// If we couldn't find an alternative pivot
			if (!foundAltPivot)
				pivotC++; // Move right one column
		}
		// PIVOT EXISTS AT (pivotR, pivotC)
		else {
			// Zero coefs of other rows in this column
			for (int r = pivotR + 1; r < n1; r++) {
				if (mat1[r][pivotC] != 0) {
					double coef = - mat1[r][pivotC] / mat1[pivotR][pivotC];
					// mat1
					addVectors(mat1[r], mat1[pivotR], m1, coef);
					// mat2
					if (mat2 != NULL)
						addVectors(mat2[r], mat2[pivotR], m2, coef);
				}
			}
			// Move to next default pivot location
			pivotC++;
			pivotR++;
		}
	}
	return numRowSwaps;
}

int ref(double **matrix, int n, int m) {
	return augmentedREF(matrix, n, m, NULL, 0, 0); // matrix is not augmented
}

int rank(double **mat, int n, int m) {
	double **copy = copyMatrix(mat, n, m);
	ref(copy, n, m);
	int numPivots = 0;
	int pivotR = 0, pivotC = 0;
	while (pivotR < n && pivotC < m) {
		if (smoothenZero(copy[pivotR][pivotC]) != 0) {
			++numPivots;
			++pivotR;
		}
		++pivotC;
	}
	freeMatrix(copy, n);
	return numPivots;
}

/*
	Represents the row reduction of an augmented matrix.
	"mat1" is of size (n1 x m1)
	"mat2" is of size (n2 x m2)
	If "mat2" is NULL, only "mat1" will be row reduced.
	PRECONDITION: n1 <= n2
*/
void augmentedRREF(double **mat1, int n1, int m1, double **mat2, int n2, int m2) {
	if (mat2 != NULL && n1 > n2) {
		printf("ERROR: invalid dimensions for augmented matrix.\n");
		return;
	}
	augmentedREF(mat1, n1, m1, mat2, n2, m2); // Reduce
	int rnk = rank(mat1, n1, m1);
	
	// Start at first element of last non-zero row
	int pivotR = rnk-1, pivotC = 0;

	while (pivotR >= 0) {
		// NO PIVOT AT (pivotR, pivotC)
		if (smoothenZero(mat1[pivotR][pivotC]) == 0) {
			pivotC = (pivotC + 1) % m1; // Move c to the right
		}
		// PIVOT EXISTS
		else {
			// Divide row by pivot
			double pivot = mat1[pivotR][pivotC];
			// mat1
			multVector(mat1[pivotR], m1, 1 / pivot);
			// mat2
			if (mat2 != NULL)
				multVector(mat2[pivotR], m2, 1 / pivot);

			for (int r = pivotR - 1; r >= 0; r--) {
				double coef = - mat1[r][pivotC];
				// mat1
				addVectors(mat1[r], mat1[pivotR], m1, coef);
				// mat2
				if (mat2 != NULL)
					addVectors(mat2[r], mat2[pivotR], m2, coef);
			}
			// Move to next starting pivot location
			pivotR--;
			pivotC = 0;
		}
	}
}

void rref(double **matrix, int n, int m) {
	augmentedREF(matrix, n, m, NULL, 0, 0); // matrix is not augmented
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

bool isUpperTriangular(double **mat, int n) {
	for (int r = 1; r < n; r++) {
		for (int c = 0; c < r; c++) {
			if (mat[r][c] != 0)
				return false;
		}
	}
	return true;
}

bool isLowerTriangular(double **mat, int n) {
	for (int c = 1; c < n; c++) {
		for (int r = 0; r < c; r++) {
			if (mat[r][c] != 0)
				return false;
		}
	}
	return true;
}

bool isTriangular(double **mat, int n) {
	return isUpperTriangular(mat, n) || isLowerTriangular(mat, n);
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
	Multiplies matrices of dimensions: (n x p)(p x r)
*/
double** multMatrices(double **m1, double **m2, int n, int p, int r) {
	double **res = createMatrix(n, r);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < r; j++) {
			double entry = 0;
			for (int k = 0; k < p; k++) {
				entry += m1[i][k] * m2[k][j];
			}
			res[i][j] = entry;
		}
	}
	return res;
}

double **getInverse(double **mat, int n) {
	double **copy = copyMatrix(mat, n, n);
	double **identity = createIdentityMatrix(n);
	augmentedRREF(copy, n, n, identity, n, n);
	freeMatrix(copy, n);
	return identity;
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

double** gramSchmidt(double **linIndSet, int n, int m, bool normVectors) {
	double **copy = copyMatrix(linIndSet, n, m);
	for (int i = 1; i < n; i++) {
		double *prj = proj(linIndSet[i], copy, i, m);
		addVectors(copy[i], prj, m, -1);
		free(prj);
	}
	if (normVectors) {
		for (int r = 0; r < n; r++) {
			normalizeVector(copy[r], m);
		}
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

	double **transp = transpose(mat, n, m);
	double **identity = createIdentityMatrix(fmax(n, m));
	augmentedRREF(transp, m, n, identity, fmax(n, m), fmax(n, m));

	// Get dimension of null space
	*dimNullSpc = m - rnk;
	double **nullSpcBasis = createMatrix(*dimNullSpc, m);
	for (int r = 0; r < *dimNullSpc; r++) {
		setVector(nullSpcBasis[r], identity[r + rnk], m);
	}

	freeMatrix(transp, m);
	freeMatrix(identity, fmax(n, m));

	return nullSpcBasis;
}

void getQRDecomp(double **mat, int n, double ***Q, double ***R) {
	double **transp = transpose(mat, n, n);
	double **orthogonalized = gramSchmidt(transp, n, n, true);
	*Q = transpose(orthogonalized, n, n);
	// Calculate R:
	*R = createMatrix(n, n);
	for (int r = 0; r < n; r++) {
		for (int c = 0; c < n; c++) {
			if (c < r)
				(*R)[r][c] = 0;
			else {
				(*R)[r][c] = dot(orthogonalized[r], transp[c], n);
			}
		}
	}
	freeMatrix(transp, n);
	freeMatrix(orthogonalized, n);
}

/*
	Precondition: all eigenvalues of matrix are real
*/
double *getEigenvalues(double **mat, int n) {
	double **finalTriangularMat = mat;

	// If the matrix is NOT triangular, execute QR algorithm. Otherwise
	// eigenvalues are the diagonal entries.
	if (!isTriangular(mat, n)) {
		double **Ai = mat;

		for (int k = 0; k < 20; k++) {
			double ***Qi = malloc(sizeof(double **));
			double ***Ri = malloc(sizeof(double **));
			getQRDecomp(Ai, n, Qi, Ri);
			// Toss R because we don't need it
			freeMatrix(*Ri, n); free(Ri);

			double **transpQi = transpose(*Qi, n, n);
			double **tempAi = Ai;
			// Multiply Ai from the left by Qi*
			Ai = multMatrices(transpQi, Ai, n, n, n);
			if (k > 0)
				freeMatrix(tempAi, n);
			freeMatrix(transpQi, n);
			tempAi = Ai;
			// Multiply Ai from the right by Qi
			Ai = multMatrices(Ai, *Qi, n, n, n);
			freeMatrix(tempAi, n);
			freeMatrix(*Qi, n);
			free(Qi);
		}
		// Update
		finalTriangularMat = Ai;
	}

	// Generate spectrum vector (diagonal of finalTriangularMat)
	double *spectrum = createVector(n);
	for (int i = 0; i < n; i++) {
		spectrum[i] = finalTriangularMat[i][i];
	}
	// If we executed the QR algorithm, we need to free the newly
	// generated triangular matrix
	if (finalTriangularMat != mat)
		freeMatrix(finalTriangularMat, n);
	return spectrum;
}

double **getEigenspace(double eigval, double **mat, int n, int *geometricMultiplicity) {
	double **copy = copyMatrix(mat, n, n);
	// Subtract eigenvalue from main diagonal
	for (int i = 0; i < n; i++) {
		copy[i][i] -= eigval;
	}
	double **eigspace = getNullSpace(copy, n, n, geometricMultiplicity);
	freeMatrix(copy, n);
	return eigspace;
}


/*
	PRINTING
*/

void printVector(double *vector, int n, bool withParenthesis) {
	if (withParenthesis)
		printf("(");
	for (int i = 0; i < n; i++) {
		double smoothened = smoothenZero(vector[i]);
		printf(smoothened >= 0 ? " %0.3f" : "%0.3f", smoothened);
		if (i < n-1)
			printf(" ");
	}
	if (withParenthesis)
		printf((n > 0 && vector[0] >= 0) ? " )" : " )");;
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

bool inputMatrix(double **mat, int n, int m) {
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

void runInvertMatrixDialog() {
	int n, m;
	double **matrix = runMatrixInputDialog(&n, &m);
	if (matrix == NULL) {
		return;
	}
	if (n != m) {
		printf("Must input square matrix.\n");
	} else if (smoothenZero(rowReductionDet(matrix, n)) == 0) {
		printf("The given matrix is not invertible.\n");
	} else {
		double **inv = getInverse(matrix, n);
		printf("Inverse:\n");
		printMatrix(inv, n, n);
		freeMatrix(inv, n);
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

	printf("Normalize vectors? (y/n): ");
	char answer;
	scanf(" %c", &answer);

	printf("Input a linearly independent set of vectors:\n");
	double **set = createMatrix(setLen, vecLen);
	for (int i = 0; i < setLen; i++) {
		inputVector(set[i], vecLen);
	}

	printf("Read value: \n");
	printRowsAsSet(set, setLen, vecLen);
	printf("\n");

	// Verify that vectors are linearly independent
	if (rank(set, setLen, vecLen) != fmin(setLen, vecLen)) {
		freeMatrix(set, setLen);
		printf("ERROR: Vectors are not linearly independent.\n");
		return;
	}

	double **orthogonalized = gramSchmidt(set, setLen, vecLen, answer == 'y');
	printf("Orthogonalized:\n");
	printRowsAsSet(orthogonalized, setLen, vecLen);
	printf("\n");
	
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

void runMultMatricesDialog() {
	// Matrices: (nxp)(pxr)
	int n, p, r;
	// Input first matrix
	double **mat1 = runMatrixInputDialog(&n, &p);
	if (mat1 == NULL)
		return;
	// Number of rows in second matrix: must = p
	int temp;
	bool error = false; // Any input errors
	// Input second matrix dimensions
	printf("Input second matrix dimensions:\n");
	printInputPrefix();
	if (!inputMatrixDims(&temp, &r)) {
		error = true;
	}
	else if (temp != p) { // If dimensions are invalid
		printf("Matrices of these sizes cannot be multiplied.\n");
		error = true;
	}
	// Exit if any errors
	if (error)
		return;
	printf("Matrix dimensions: %dx%d\n", p, r);
	// Input second matrix
	double **mat2 = createMatrix(p, r);
	printf("Input second matrix:\n");
	if (!inputMatrix(mat2, p, r)) {
		freeMatrix(mat1, n);
		return;
	}
	printf("Read matrix:\n");
	printMatrix(mat2, p, r);
	// Compute and print product
	double **prod = multMatrices(mat1, mat2, n, p, r);
	printf("Product:\n");
	printMatrix(prod, n, r);
	// Release memory
	freeMatrix(mat1, n);
	freeMatrix(mat2, p);
	freeMatrix(prod, n);
}

void runQRDecompDialog() {
	int n, m;
	double **matrix = runMatrixInputDialog(&n, &m);
	if (matrix == NULL) {
		return;
	}
	if (n != m) {
		printf("Must input square matrix.\n");
	} else {
		double ***Q = malloc(sizeof(double **));
		double ***R = malloc(sizeof(double **));
		getQRDecomp(matrix, n, Q, R);
		printf("Q:\n");
		printMatrix(*Q, n, n);
		printf("R:\n");
		printMatrix(*R, n, n);
		freeMatrix(*Q, n);
		freeMatrix(*R, n);
		free(Q);
		free(R);
	}
	freeMatrix(matrix, n);
}

void runEigenvaluesDialog() {
	int n, m;
	double **matrix = runMatrixInputDialog(&n, &m);
	if (matrix == NULL) {
		return;
	}
	if (n != m) {
		printf("Must input square matrix.\n");
	} else {
		double *spectrum = getEigenvalues(matrix, n);
		printf("\nEigenvalues:\n");

		int sumAlgMults = 0; // Sum of the algebraic multiplicities (in the end = n)

		while (sumAlgMults < n) {
			// GET ALGEBRAIC MULTIPLICITY
			int algebraicMultiplcity = 0;
			for (int j = 0; j < n; j++) {
				// spectrum[sumAlgMults] because we want to skip repititions
				// of the same eigenvalues
				if (spectrum[sumAlgMults] == spectrum[j]) {
					algebraicMultiplcity++;
				}
			}

			// GET EIGENSPACE
			int geometricMultiplicity;
			double **eigenspace = getEigenspace(spectrum[sumAlgMults],
											matrix, n, &geometricMultiplicity);

			// PRINT EIGENVALUE, MULTIPLICITIES, AND EIGENSPACE
			printf("\t%0.3f\n", smoothenZero(spectrum[sumAlgMults]));
			printf("\t\tAlgebraic Multiplicity: %d\n", algebraicMultiplcity);
			printf("\t\tGeometric Multiplicity: %d\n", geometricMultiplicity);
			// Print eigenspace basis
			printf("\t\t%0.3f - eigenspace:\n", spectrum[sumAlgMults]);
				printf("\t\tspan {\n");
			for (int i = 0; i < geometricMultiplicity; i++) {
				printf("\t\t\t");
				printVector(eigenspace[i], m, true);
				printf("\n");
			}
			printf("\t\t}\n\n");

			freeMatrix(eigenspace, geometricMultiplicity);

			sumAlgMults += algebraicMultiplcity;
		}

		free(spectrum);
	}
	freeMatrix(matrix, n);
}

int main() {
	printf("Welcome to Matrix Master!\n");

	int input;

	do {
		printf("\nMAIN MENU\n");
		printf("1. Reduce to row echelon form (REF)\n");
		printf("2. Reduce to reduced row echelon form (RREF)\n");
		printf("3. Determinant\n");
		printf("4. Invert Matrix\n");
		printf("5. Null Space\n");
		printf("6. Adjoint\n");
		printf("7. Orthogonal Projection\n");
		printf("8. Gram-Schmidt\n");
		printf("9. Mult Matrices\n");
		printf("10. QR Decomposition\n");
		printf("11. Eigenvalues\n");
		printf("0: Quit\n");
		printf("Input:\n");
		input = readInt();

		switch (input) {
			case 1: runREFDialog(); break;
			case 2: runRREFDialog(); break;
			case 3: runDetDialog(); break;
			case 4: runInvertMatrixDialog(); break;
			case 5: runNullSpaceDialog(); break;
			case 6: runAdjointDialog(); break;
			case 7: runOrthoProjectionDialog(); break;
			case 8: runGramSchmidtDialog(); break;
			case 9: runMultMatricesDialog(); break;
			case 10: runQRDecompDialog(); break;
			case 11: runEigenvaluesDialog(); break;
			case 0: break;
			default:
				printf("Invalid input. Try again:\n");
		}

	} while (input != 0);

	printf("Process ended, terminating program.");

	return 0;
}

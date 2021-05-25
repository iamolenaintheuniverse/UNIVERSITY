#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <omp.h>

// Function for formatted matrix output
void PrintMatrix(double* pMatrix, int RowCount, int ColCount)
{
	int i, j; // Loop variables
	for (i = 0; i < RowCount; i++)
	{
		for (j = 0; j < ColCount; j++)
			printf("%7.4f ", pMatrix[i * ColCount + j]);
		printf("\n");
	}
}

void DummyDataInitialization(double* pMatrix, int Size)
{
	int i, j; // Loop variables
	for (i = 0; i < Size; i++)
	{
		for (j = 0; j < Size; j++)
		{
			pMatrix[i * Size + j] = i;
		}
	}
}

void RandomDataInitialization(double* pMatrixA, double* pMatrixB, int Size)
{
	int i, j; // Loop variables
	srand(time(NULL));
	for (i = 0; i < Size; i++)
	{
		for (j = 0; j < Size; j++)
		{
			pMatrixA[i * Size + j] = rand() / double(1000);
			pMatrixB[i * Size + j] = rand() / double(1000);
		}
	}
}

void MatrixInitialization(double*& pMatrixA, double*& pMatrixB, double*& pResultC, int& initialMtrSize, int& t)
{
	printf("\nEnter number of threads: ");
	fflush(stdout);
	scanf_s("%d", &t);
	printf("\nEnter size of the initial objects: ");
	fflush(stdout);
	scanf_s("%d", &initialMtrSize);

	pMatrixA = new double[initialMtrSize * initialMtrSize];
	pMatrixB = new double[initialMtrSize * initialMtrSize];
	pResultC = new double[initialMtrSize * initialMtrSize];
	//DummyDataInitialization(pMatrixA, initialMtrSize);
	//DummyDataInitialization(pMatrixB, initialMtrSize);
	RandomDataInitialization(pMatrixA, pMatrixB, initialMtrSize);

	// printf("Matrix A:\n");
	// PrintMatrix(pMatrixA, initialMtrSize, initialMtrSize);
	// printf("Matrix B:\n");
	// PrintMatrix(pMatrixB, initialMtrSize, initialMtrSize);
}

void ParallelResultCalculation(double* pMatrixA, double* pMatrixB, double* pResultC, int initialMtrSize, int t) 
{
	
#pragma omp parallel num_threads(t)
	{
		int i, j, k, t_number;
		t_number = omp_get_thread_num();
		// printf("Thread %d starting matrix multiply.\n", t_number);
#pragma omp for
		for (i = 0; i < initialMtrSize; i++) 
		{
			// printf("Thread %d calculate row = %d.\n", t_number, i);
			for (j = 0; j < initialMtrSize; j++) 
			{
				double cell = 0;
				for (k = 0; k < initialMtrSize; k++) 
				{
					cell += pMatrixA[i * initialMtrSize + k] * pMatrixB[k * initialMtrSize + j];
				}
				pResultC[i * initialMtrSize + j] = cell;
			}
		}
	}

}

void SerialResultCalculation(double* pMatrixA, double* pMatrixB, double* pResultC, int initialMtrSize)
{
	int i, j, k; // Loop variables
	for (i = 0; i < initialMtrSize; i++)
		for (j = 0; j < initialMtrSize; j++)
		{
			double sum = 0;
			for (k = 0; k < initialMtrSize; k++)
				sum = sum + pMatrixA[i * initialMtrSize + k] * pMatrixB[k * initialMtrSize + j];
			pResultC[i * initialMtrSize + j] = sum;
		}
}

void TestResult(double* pMatrixA, double* pMatrixB, double* pResultC, int initialMtrSize)
{
	// Buffer for storing the result of serial matrix multiplication
	double* pSerialResult;
	int equal = 0; // Flag, that shows wheather the vectors are identical
	int i; // Loop variable
	pSerialResult = new double[initialMtrSize * initialMtrSize];
	SerialResultCalculation(pMatrixA, pMatrixB, pSerialResult, initialMtrSize);
	for (i = 0; i < initialMtrSize; i++)
	{
		if (int(pResultC[i]) != int(pSerialResult[i]))
			equal = 1;
	}
	printf("The result of serial matrix multiplication: \n");
	PrintMatrix(pSerialResult, initialMtrSize, initialMtrSize);
	if (equal == 1)
	{
		printf("The result of parallel matrix multiplication: \n");
		PrintMatrix(pResultC, initialMtrSize, initialMtrSize);
		printf("The results of serial and parallel algorithms are NOT identical. Check your code.");
	}
	else
	{
		printf("The result of parallel matrix multiplication: \n");
		PrintMatrix(pResultC, initialMtrSize, initialMtrSize);
		printf("The results of serial and parallel algorithms are identical.");
	}
	delete[] pSerialResult;
}

void ProcessTermination(double* pMatrixA, double* pMatrixB, double* pResultC)
{
	delete[] pMatrixA;
	delete[] pMatrixB;
	delete[] pResultC;
}

int main()
{
	double* pMatrixA; // The first argument - initial matrix A
	double* pMatrixB; // The second argument - initial matrix B
	double* pResultC; // Result matrix for matrix multiplication

	int initialMtrSize; // Sizes of initial matrices

	int t; // number of threads

	time_t start, finish;
	double duration;

	MatrixInitialization(pMatrixA, pMatrixB, pResultC, initialMtrSize, t);

	start = clock();
	ParallelResultCalculation(pMatrixA, pMatrixB, pResultC, initialMtrSize, t);
	finish = clock();
	
	duration = (finish - start) / double(CLOCKS_PER_SEC);
	printf("\n Time of execution: %f\n", duration);

	//TestResult(pMatrixA, pMatrixB, pResultC, initialMtrSize);

	ProcessTermination(pMatrixA, pMatrixB, pResultC);

	return 0;
}
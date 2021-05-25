#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>

void TestDistribution(double* pProcBlockA, double* pProcBlockB, int BlockSize);

int ProcNum; // Number of available processes
int ProcRank; // Rank of current process

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

void ProcessInitialization(double*& pMatrixA, double*& pMatrixB, double*& pResultC, int& initialMtrSize, 
                           double*& pProcBlockA, double*& pProcBlockB, double*& pProcResult, double*& pTemporaryBlock,
                           int& BlockSize, int& GridSize)
{
    if (ProcRank == 0)
    {
        do
        {
            printf("\nEnter size of the initial objects: ");
            fflush(stdout);
            scanf_s("%d", &initialMtrSize);

            if (initialMtrSize < ProcNum)
            {
                printf("Size of the objects must be greater than "
                    "number of processes! \n ");
            }
            if (initialMtrSize % (int)sqrt((double)ProcNum) != 0)
            {
                printf("Size of objects must be divisible by "
                    "number of processes! \n");
            }
        } while ((initialMtrSize < ProcNum) || (initialMtrSize % (int)sqrt((double)ProcNum) != 0));
    }
    MPI_Bcast(&initialMtrSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    GridSize = (int)sqrt((double)ProcNum);
    if (ProcNum != GridSize * GridSize)
        if (ProcRank == 0)
        {
            printf("Number of processes should be perfect square");
            exit(0);
        }
    BlockSize = initialMtrSize / GridSize;

    // Memory allocation
    pProcBlockA = new double[BlockSize * BlockSize];
    pProcBlockB = new double[BlockSize * BlockSize];
    pProcResult = new double[BlockSize * BlockSize];

    pTemporaryBlock = new double[BlockSize * BlockSize];//???

    for (int i = 0; i < BlockSize * BlockSize; i++)
        pProcResult[i] = 0;

    // Obtain the values of initial objects’ elements
    if (ProcRank == 0)
    {
        pMatrixA = new double[initialMtrSize * initialMtrSize];
        pMatrixB = new double[initialMtrSize * initialMtrSize];
        pResultC = new double[initialMtrSize * initialMtrSize];

        //DummyDataInitialization(pMatrixA, initialMtrSize);
        //DummyDataInitialization(pMatrixB, initialMtrSize);
        RandomDataInitialization(pMatrixA, pMatrixB, initialMtrSize);

        printf("Matrix A:\n");
        PrintMatrix(pMatrixA, initialMtrSize, initialMtrSize);
        printf("Matrix B:\n");
        PrintMatrix(pMatrixB, initialMtrSize, initialMtrSize);
    }
}

void CreateGridComm(int GridSize, int* GridCoords, MPI_Comm& GridComm, MPI_Comm& RowComm, MPI_Comm& ColComm)
{
    int DimSize[2] = { GridSize, GridSize };
    int Periodic[2] = { 0,0 };
    int SubDims[2];
    MPI_Cart_create(MPI_COMM_WORLD, 2, DimSize, Periodic, 0, &GridComm);
    MPI_Cart_coords(GridComm, ProcRank, 2, GridCoords);
    SubDims[0] = 0;
    SubDims[1] = 1;
    MPI_Cart_sub(GridComm, SubDims, &RowComm);
    SubDims[0] = 1;
    SubDims[1] = 0;
    MPI_Cart_sub(GridComm, SubDims, &ColComm);
}

void DataDistribution(double* pMatrixA, double* pMatrixB, 
                      double* pProcBlockA, double* pProcBlockB, 
                      int initialMtrSize, int BlockSize,
                      MPI_Datatype& MPI_BLOCK, MPI_Comm GridComm,
                      MPI_Request *request, MPI_Status *status)
{
    MPI_Type_vector(BlockSize, BlockSize, initialMtrSize, MPI_DOUBLE, &MPI_BLOCK);
    MPI_Type_commit(&MPI_BLOCK);
    if (ProcRank == 0)
    {
        for (int r = 0; r < ProcNum; r++)
        {
            int c[2];
            MPI_Cart_coords(GridComm, r, 2, c);
            MPI_Isend(pMatrixA + c[0] * initialMtrSize * BlockSize + c[1] * BlockSize, 1,
                MPI_BLOCK, r, 111, GridComm, request);
            MPI_Isend(pMatrixB + c[0] * initialMtrSize * BlockSize + c[1] * BlockSize, 1,
                MPI_BLOCK, r, 222, GridComm, request);
        }
    }

    MPI_Recv(pProcBlockA, BlockSize * BlockSize, MPI_DOUBLE, 0, 111, GridComm, status);
    MPI_Recv(pProcBlockB, BlockSize * BlockSize, MPI_DOUBLE, 0, 222, GridComm, status);
}

void pProcBlockA_SendRow(int iter, double* pProcBlockA, double* pTemporaryBlock, 
                         int GridSize, int BlockSize, int* GridCoords,
                         MPI_Comm RowComm)
{
    int p = (GridCoords[0] + iter) % GridSize;
    if (GridCoords[1] == p)
        for (int i = 0; i < BlockSize * BlockSize; i++)
            pTemporaryBlock[i] = pProcBlockA[i];
    MPI_Bcast(pTemporaryBlock, BlockSize * BlockSize, MPI_DOUBLE, p, RowComm);
}

void pBlockMult(double* pTemporaryBlock, double* pProcBlockB, double* pProcResult, int BlockSize)
{
    for (int i = 0; i < BlockSize; i++)
        for (int j = 0; j < BlockSize; j++)
        {
            double t = 0;
            for (int k = 0; k < BlockSize; k++)
                t += pTemporaryBlock[i * BlockSize + k] * pProcBlockB[k * BlockSize + j];
            pProcResult[i * BlockSize + j] += t;
        }
}

void pProcBlockB_SendCol(double* pProcBlockA, double* pProcBlockB, int BlockSize, int GridSize, int* GridCoords, 
                         MPI_Comm ColComm, MPI_Request* req1, MPI_Request* req2, MPI_Status* status)
{
    MPI_Status s;

    int NextProc = GridCoords[0] + 1;
    if (GridCoords[0] == GridSize - 1)
        NextProc = 0;

    int PrevProc = GridCoords[0] - 1;
    if (GridCoords[0] == 0)
        PrevProc = GridSize - 1;

    MPI_Sendrecv_replace(pProcBlockB, BlockSize * BlockSize, MPI_DOUBLE, PrevProc, 0, NextProc, 0, ColComm, &s);
    //MPI_Irecv(pProcBlockB, BlockSize * BlockSize, MPI_DOUBLE, NextProc, 333, ColComm, req1);
    //MPI_Isend(pProcBlockB, BlockSize * BlockSize, MPI_DOUBLE, PrevProc, 333, ColComm, req2);

    //MPI_Wait(req1, status);
    //MPI_Wait(req2, status);

    //TestDistribution(pProcBlockA, pProcBlockB, BlockSize);
}

// Process rows and vector multiplication
void ParallelResultCalculation(double* pProcBlockA, double* pProcBlockB, double* pProcResult, double* pTemporaryBlock,
                               int BlockSize, int GridSize, int* GridCoords,
                               MPI_Comm RowCom, MPI_Comm ColCom, MPI_Status* status)
{
    MPI_Request req1, req2;
    for (int i = 0; i < GridSize; i++)
    {
        pProcBlockA_SendRow(i, pProcBlockA, pTemporaryBlock, GridSize, BlockSize, GridCoords, RowCom);
        pBlockMult(pTemporaryBlock, pProcBlockB, pProcResult, BlockSize);
        pProcBlockB_SendCol(pProcBlockA, pProcBlockB, BlockSize, GridSize, GridCoords, ColCom, &req1, &req2, status);
    }
}

// Result matrix replication
void ResultReplication(double* pProcResult, double* pResultC, int initialMtrSize, int BlockSize, 
                       MPI_Comm GridComm, MPI_Datatype MPI_BLOCK,
                       MPI_Request* request, MPI_Status* status)
{
    //PrintMatrix(pProcResult, BlockSize, BlockSize);
    int request_complete = 0;
    MPI_Isend(pProcResult, BlockSize * BlockSize, MPI_DOUBLE, 0, 444, GridComm, request);

    if (ProcRank == 0)
    {
        for (int r = 0; r < ProcNum; r++)
        {
            int c[2];
            //MPI_Probe(r, 444, GridComm, status);
            MPI_Cart_coords(GridComm, r, 2, c);
            MPI_Recv(pResultC + c[0] * initialMtrSize * BlockSize + c[1] * BlockSize, 1,
                    MPI_BLOCK, r, 444, GridComm, status);
        }
        //PrintMatrix(pResultC, initialMtrSize, initialMtrSize);
    }
}

void TestDistribution(double* pProcBlockA, double* pProcBlockB, int BlockSize)
{
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < ProcNum; i++)
    {
        if (ProcRank == i)
        {
            printf("\nProcRank = %d \n", ProcRank);
            printf(" Matrix A Stripe:\n");
            PrintMatrix(pProcBlockA, BlockSize, BlockSize);
            printf(" Matrix B Stripe:\n");
            PrintMatrix(pProcBlockB, BlockSize, BlockSize);
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void TestPartialResults(double* pProcResult, int BlockSize) 
{
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < ProcNum; i++) 
    {
        if (ProcRank == i) 
        {
            printf("ProcRank = %d \n", ProcRank);
            printf("Part of result matrix: \n");
            PrintMatrix(pProcResult, BlockSize, BlockSize);
        }
        MPI_Barrier(MPI_COMM_WORLD);
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
    MPI_Barrier(MPI_COMM_WORLD);
    if (ProcRank == 0)
    {
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
}

void ProcessTermination(double* pMatrixA, double* pMatrixB, double* pResultC, 
                        double* pTemporaryBlock, double* pProcBlockA, double* pProcBlockB, double* pProcResult)
{
    if (ProcRank == 0)
    {
        delete[] pMatrixA;
        delete[] pMatrixB;
        delete[] pResultC;
    }
    delete[] pProcBlockA;
    delete[] pProcBlockB;
    delete[] pProcResult;
    delete[] pTemporaryBlock;
}

int main(int argc, char* argv[])
{
    MPI_Comm GridComm, ColComm, RowComm;

    MPI_Datatype MPI_BLOCK;

    double* pMatrixA; // The first argument - initial matrix A
    double* pMatrixB; // The second argument - initial matrix B
    double* pResultC; // Result matrix for matrix multiplication
    int initialMtrSize; // Sizes of initial matrices

    double* pProcBlockA; // Block of the matrix A on current process
    double* pProcBlockB; // Block of the matrix B on current process
    double* pProcResult; // Block of result matrix on current process
    double* pTemporaryBlock;
    
    int BlockSize;
    int GridSize;

    int GridCoords[2];

    MPI_Request request;
    MPI_Status status;

    double Start, Finish, Duration;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if (ProcRank == 0)
        printf("Parallel matrix multiplication program\n");

    ProcessInitialization(pMatrixA, pMatrixB, pResultC, initialMtrSize, pProcBlockA, pProcBlockB, pProcResult, pTemporaryBlock, BlockSize, GridSize);

    CreateGridComm(GridSize, GridCoords, GridComm, RowComm, ColComm);

    Start = MPI_Wtime();
    DataDistribution(pMatrixA, pMatrixB, pProcBlockA, pProcBlockB, initialMtrSize, BlockSize, MPI_BLOCK, GridComm, &request, &status);
    //TestDistribution(pProcBlockA, pProcBlockB, BlockSize);

    ParallelResultCalculation(pProcBlockA, pProcBlockB, pProcResult, pTemporaryBlock, BlockSize, GridSize, GridCoords, RowComm, ColComm, &status);
    //TestPartialResults(pProcResult, BlockSize);

    ResultReplication(pProcResult, pResultC, initialMtrSize, BlockSize, GridComm, MPI_BLOCK, &request, &status);
    TestResult(pMatrixA, pMatrixB, pResultC, initialMtrSize);
    Finish = MPI_Wtime();
    
    Duration = Finish - Start;
    if (ProcRank == 0) 
    {
        printf("\nTime of execution = % f\n", Duration);
    }

    //TestResult(pMatrixA, pMatrixB, pResultC, initialMtrSize);
    MPI_Wait(&request, &status);
    
    ProcessTermination(pMatrixA, pMatrixB, pResultC, pTemporaryBlock, pProcBlockA, pProcBlockB, pProcResult);

    MPI_Finalize();
}
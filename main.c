#include <stdio.h>
#include <malloc.h>
#include "mpi.h"

#define N1 4
#define N2 4
#define N3 4
#define t 0.00001f
#define eps 0.00001f

#define RANK_ROOT 0
#define CONTINUE 1
#define EXIT 0


void fill_matrix(double *A, const int height, const int width) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (i == j)
                A[i * width + j] = 2.0f;
            else
                A[i * width + j] = 1.0f;
        }
    }
}

double *create_matrix() {
    unsigned long long matr_size = (unsigned long long) N1 * N2;
    double *A = calloc(matr_size, sizeof(double));
    fill_matrix(A, N1, N2);
    return A;
}


int *set_chunk_sizes(int size, int matrix_height) {
    int *chunk_sizes = calloc(size, sizeof(int));
    int quot = matrix_height / size;
    int remain = matrix_height % size;

    for (int i = 0; i < size; i++) {
        if (i < size - remain) {
            chunk_sizes[i] = quot;
        } else {
            chunk_sizes[i] = quot + 1;
        }
    }
    return chunk_sizes;
}


int *set_offset(int size, const int *chunk_size_arr) {
    int *offset_arr = calloc(size, sizeof(int));
    offset_arr[0] = 0;
    for (int i = 1; i < size; i++) {
        offset_arr[i] = offset_arr[i - 1] + chunk_size_arr[i - 1];
    }
    return offset_arr;
}


int *set_matr_chunk_sizes(const int *chunk_size_arr, int size, int param) {
    int *matrix_chunk_size_arr = calloc(size, sizeof(int));
    for (int i = 0; i < size; i++) {
        matrix_chunk_size_arr[i] = chunk_size_arr[i] * param;
    }
    return matrix_chunk_size_arr;
}


int *set_matr_offset(const int *offset_arr, int size, int param) {
    int *matrix_offset_arr = calloc(size, sizeof(int));
    for (int i = 0; i < size; i++) {
        matrix_offset_arr[i] = offset_arr[i] * param;
    }
    return matrix_offset_arr;
}


double *send_matrix_from_root(int comm_rank, int comm_size, double *A, int *chunk_size_arr, int *offset_arr) {
    int *matrix_chunk_size_arr = set_matr_chunk_sizes(chunk_size_arr, comm_size, N1);
    int *matrix_offset_arr = set_matr_offset(offset_arr, comm_size, N1);

    double *matrix_chunk = calloc(matrix_chunk_size_arr[comm_rank], sizeof(double));

    MPI_Scatterv(A, matrix_chunk_size_arr, matrix_offset_arr, MPI_DOUBLE, matrix_chunk,
                 matrix_chunk_size_arr[comm_rank],
                 MPI_DOUBLE, RANK_ROOT, MPI_COMM_WORLD);
    return matrix_chunk;
}

int run(int comm_size, int comm_rank) {

    int *chunk_size_arr = set_chunk_sizes(comm_size, N1);
    int *offset_arr = set_offset(comm_size, chunk_size_arr);

    int comm_chunk_size = chunk_size_arr[comm_rank];

    double *A = NULL;

    double *matrix_chunk = NULL;

    if (comm_rank == RANK_ROOT) {
        A = create_matrix();
    }

    matrix_chunk = send_matrix_from_root(comm_rank, comm_size, A, chunk_size_arr, offset_arr);
}

MPI_Comm create_new_comm(int size, int rank, int *dims) {
    int periods[2] = {0, 0}, coords[2], reorder = 1;

    int sizey, sizex, ranky, rankx;
    int prevy, prevx, nexty, nextx;
    MPI_Comm comm2d;

    MPI_Dims_create(size, 2, dims);
    sizex = dims[0];
    sizey = dims[1];
    if (rank == RANK_ROOT)
        printf("Dimensions (%d, %d)\n", sizex, sizey);

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &comm2d);
    MPI_Comm_rank(comm2d, &rank);
    MPI_Cart_get(comm2d, 2, dims, periods, coords);
    printf("Process %d has coordinates (%d, %d)\n", rank, coords[0], coords[1]);

//    MPI_Cart_shift(comm2d, 0, 1, &prevy, &nexty);
//    MPI_Cart_shift(comm2d, 1, 1, &prevx, &nextx);
    return comm2d;
}

//void create_new_groups(MPI_Comm comm2d, MPI_Comm *groups, int *dims) {
//    int periods[2] = {0, 0}, coords[2], reorder = 1;
//    MPI_Cart_get(comm2d, 2, dims, periods, coords);
//
//    MPI_Comm_split(comm2d, coords[0], 1, &groups[coords[0]]);
//    MPI_Comm_split(comm2d, coords[1] + dims[0], 1, &groups[coords[1] + dims[0]]);
//}

int main(int argc, char **argv) {

//    double *A = calloc(N1 * N2, sizeof(double));
//    fill_matrix(A, N1, N2);
//    double *B = calloc(N2 * N3, sizeof(double));
//    fill_matrix(B, N2, N3);

    MPI_Init(&argc, &argv);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    double *arr;
    double *arr_recv_row;
    double *arr_recv_col;
    if (world_rank == RANK_ROOT) {
        arr = calloc(world_size, sizeof(double));
        for (int i = 0; i < world_size; i++) {
            arr[i] = i;
        }
        for (int i = 0; i < world_size; i++) {
            printf("%lf ", arr[i]);
        }
        puts("\n");
    }
    int dims[2] = {0, 0}; //x y

    MPI_Comm comm2d = create_new_comm(world_size, world_rank, dims);
    int periods[2] = {0, 0}, coords[2];
    MPI_Cart_get(comm2d, 2, dims, periods, coords);

    MPI_Comm row_comm, col_comm;
    MPI_Comm_split(comm2d, coords[0], coords[1], &col_comm);
    MPI_Comm_split(comm2d, coords[1], coords[0], &row_comm);
    MPI_Barrier(comm2d);
    int row_rank, row_size;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_size(row_comm, &row_size);
    int col_rank, col_size;
    MPI_Comm_rank(col_comm, &col_rank);
    MPI_Comm_size(col_comm, &col_size);

    arr_recv_row = calloc(row_size, sizeof(double));
    arr_recv_col = calloc(col_size, sizeof(double));
    // printf("SIZES: %d %d ",row_size,col_size);
    if (coords[0] == 0) {

        MPI_Scatter(arr, row_size, MPI_DOUBLE, arr_recv_row, row_size, MPI_DOUBLE, 0, col_comm);
    }
    if (coords[1] == 0) {
        MPI_Scatter(arr, col_size, MPI_DOUBLE, arr_recv_col, col_size, MPI_DOUBLE, 0, row_comm);

    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(arr_recv_row, row_size, MPI_DOUBLE, 0, row_comm);
    MPI_Bcast(arr_recv_col, col_size, MPI_DOUBLE, 0, col_comm);

    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < col_size; i++) {
        printf("%lf ", arr_recv_col[i]);
    }
    puts("\n");
    for (int i = 0; i < row_size; i++) {
        printf("%lf ", arr_recv_row[i]);
    }
    puts("\n\n");

    //   create_new_groups(comm2d, groups, dims);

//    int new_rank, new_size;
//    MPI_Comm_rank(row_comm, &new_rank);
//    MPI_Comm_size(row_comm, &new_size);
//
//    if (rank == 0) {
//        printf("Number of processes in MPI_COMM_WORLD: %d\n", size);
//        printf("Number of processes in new_comm: %d\n", new_size);
//    }
//
//    for (int i = 0; i < new_size; i++) {
//        MPI_Barrier(row_comm); // синхронизация процессов в новом коммуникаторе
//        if (new_rank == i) {
//            printf("Process %d in row_comm: Hello from (%d, %d) in Cart_comm !\n", i, coords[0], coords[1]);
//        }
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    MPI_Comm_rank(col_comm, &new_rank);
//    MPI_Comm_size(col_comm, &new_size);
//
//    if (rank == 0) {
//        printf("Number of processes in MPI_COMM_WORLD: %d\n", size);
//        printf("Number of processes in new_comm: %d\n", new_size);
//    }
//
//    for (int i = 0; i < new_size; i++) {
//        MPI_Barrier(col_comm); // синхронизация процессов в новом коммуникаторе
//        if (new_rank == i) {
//            printf("Process %d in row_comm: Hello from (%d, %d) in Cart_comm !\n", i, coords[0], coords[1]);
//        }
//    }

//
//     MPI_Comm_rank(row_comm, &rank);
//     printf("Comm row rank: %d\n", rank);
//     MPI_Comm_rank(col_comm, &rank);
//     printf("Comm col rank: %d\n", rank);
//    // MPI_Comm_rank(groups[coords[1] + dims[0]], &rank);
//   // printf("Comm rawf: %d\n", rank);
//    puts("\n");

    // run(size, rank);

    //   MPI_Comm_size(MPI_COMM_WORLD, &size);
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    if (rank == RANK_ROOT) {
//
    //printf("Comm size: %d\n", size);
//    }

    MPI_Finalize();
    return 0;
}


/* 2 2 2   1 | 1 | 1
 * _____
 * 3 3 3   2 | 2 | 2   ==
 * _____
 * 2 2 2   1 | 1 | 1*/
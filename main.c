#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "mpi.h"

#define N1 4
#define N2 4
#define N3 4

#define RANK_ROOT 0
#define CONTINUE 1
#define EXIT 0

typedef struct {
    double *data;
    int height;
    int width;
} Matrix;


void fill_matrix(Matrix *A) {
    for (int i = 0; i < A->height; i++) {
        for (int j = 0; j < A->width; j++) {
            if (i == j)
                A->data[i * A->width + j] = 2.0f;
            else
                A->data[i * A->width + j] = 1.0f;
        }
    }
}

void old_fill_matrix(double *A, int height, int width) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (i == j)
                A[i * width + j] = 2.0f;
            else
                A[i * width + j] = 1.0f;
        }
    }
}

void print_matrix(Matrix *A) {
    for (int i = 0; i < A->height; i++) {
        for (int j = 0; j < A->width; j++) {
            printf("%lf ", A->data[i * A->width + j]);
        }
        puts("\n");
    }
}

void old_print_matrix(double *A, int height, int width) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            printf("%lf ", A[i * width + j]);
        }
        puts("\n");
    }
}
//
//double *create_matrix() {
//    unsigned long long matr_size = (unsigned long long) N1 * N2;
//    double *A = calloc(matr_size, sizeof(double));
//    fill_matrix(A);
//    return A;
//}



void matrix_multiply(double *A, int A_height, int A_width, double*B, int B_height, int B_width){
    if(A_width!= B_height)
        return;




}


MPI_Comm create_new_comm(int size, int rank, int *dims) {
    int periods[2] = {0, 0}, coords[2], reorder = 1;

    int sizey, sizex, ranky, rankx;
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
    MPI_Init(&argc, &argv);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    double *arr;
    double *arr_recv_row;
    double *arr_recv_col;
    Matrix *A;
    Matrix *B;
    Matrix *A_recv;
    Matrix *B_recv;
    double *A1;
    double *B1;
    if (world_rank == RANK_ROOT) {
        A1 = calloc(N1 * N2, sizeof(double));
        old_fill_matrix(A1, N1, N2);
        A = (Matrix *) calloc(1, sizeof(Matrix));
        A->height = N1;
        A->width = N2;
        A->data = calloc(N1 * N2, sizeof(double));
        fill_matrix(A);
        B1 = calloc(N1 * N2, sizeof(double));
        old_fill_matrix(B1, N1, N2);



//        B = calloc(1, sizeof(Matrix));
//        B->height = N2;
//        B->width = N3;
//        B->data = calloc(N2 * N3, sizeof(double));
//        fill_matrix(B);
//        print_matrix(A);
//        arr = calloc(world_size, sizeof(double));
//        for (int i = 0; i < world_size; i++) {
//            arr[i] = i;
//        }
//        for (int i = 0; i < world_size; i++) {
//            printf("%lf ", arr[i]);
//        }
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

    if (row_size != dims[0] || col_size != dims[1]) {
        perror("Bad size of topology");
        exit(1);
    }
    if (N1 % col_size != 0 || N3 % row_size != 0) {
        perror("Bad size of matrix");
        exit(1);
    }

    int A_split_size = N1 / col_size;
    int B_split_size = N3 / row_size;
    //  printf("Slices: A:%d B:%d\n", A_split_size, B_split_size);
    A_recv = (Matrix *) calloc(1, sizeof(Matrix));
    A_recv->height = A_split_size;
    A_recv->width = N2;
    double *A_recv1 = calloc(A_split_size * N2, sizeof(double));
    double *B_recv1 = calloc(B_split_size * N2, sizeof(double));
    A_recv->data = calloc(A_split_size * N2, sizeof(double));
    // print_matrix(A_recv);
//    B_recv = malloc(sizeof (Matrix));
//    A_recv->height = A_split_size;
//    A_recv->width = N2;
//    B_recv->data = calloc(col_size * B_split_size, sizeof(double));
    arr_recv_row = calloc(row_size, sizeof(double));
    arr_recv_col = calloc(col_size, sizeof(double));
    //  printf("ALO:%d ", N2 * A_split_size);
    // printf("SIZES: %d %d ", row_size, col_size);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Datatype column;
    MPI_Type_vector(N2, B_split_size, N3, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);
    MPI_Datatype column_resized;
    MPI_Type_create_resized(column, 0, (int) (B_split_size * sizeof(double)), &column_resized);
    MPI_Type_commit(&column_resized);

    if (coords[0] == 0) {
//        if (world_rank ==RANK_ROOT) {
//            printf("%lf\n", A->data[A_split_size * N2-1]);
//            printf("%lf\n", A->data[A_split_size * N2 * 2-1]);
//        }
        MPI_Scatter(A1, N2 * A_split_size, MPI_DOUBLE, A_recv1, N2 * A_split_size, MPI_DOUBLE, 0, col_comm);
        //    MPI_Scatter(A->data, A_split_size * N2, MPI_DOUBLE, A_recv->data, A_split_size * N2, MPI_DOUBLE, 0, col_comm);
    }
    if (coords[1] == 0) {
        MPI_Scatter(B1, 1, column_resized, B_recv1, N2 * B_split_size, MPI_DOUBLE, 0, row_comm);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(A_recv1, N2 * A_split_size, MPI_DOUBLE, 0, row_comm);
    MPI_Bcast(B_recv1, N2 * B_split_size, MPI_DOUBLE, 0, col_comm);

    MPI_Barrier(MPI_COMM_WORLD);
    sleep(1 + world_rank);
//    for (int i = 0; i < col_size; i++) {
//        printf("%lf ", arr_recv_col[i]);
//    }
    puts("\n");
    sleep(1 + world_rank);
    printf("RANK (%d,%d): ", coords[0], coords[1]);
    old_print_matrix(A_recv1, A_split_size, N2);
    puts("\n\n");
    old_print_matrix(B_recv1, N2, B_split_size);

    matrix_multiply(A_recv1, A_split_size, N2, B_recv1, N2, B_split_size);
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
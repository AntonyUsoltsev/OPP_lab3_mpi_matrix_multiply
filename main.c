#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
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

//void old_fill_matrix(double *A, int height, int width) {
//    for (int i = 0; i < height; i++) {
//        for (int j = 0; j < width; j++) {
//            if (i == j)
//                A[i * width + j] = 2.0f;
//            else
//                A[i * width + j] = 1.0f;
//        }
//    }
//}
void old_fill_matrix(double *A, int height, int width) {
    for (int i = 0; i < height * width; i++) {
        A[i] = i + 1;
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

void
matrix_multiply(const double *A, int A_height, int A_width, const double *B, int B_height, int B_width, double *C) {
    if (A_width != B_height)
        return;

    for (int i = 0; i < A_height; i++) {
        for (int j = 0; j < B_width; j++) {
            double sum = 0;
            for (int k = 0; k < A_width; k++) {
                sum += A[i * A_width + k] * B[k * B_width + j];
            }
            C[i * B_width + j] = sum;
        }
    }
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


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    double *A1;
    double *B1;
    double *C;
    if (world_rank == RANK_ROOT) {
        A1 = calloc(N1 * N2, sizeof(double));
        old_fill_matrix(A1, N1, N2);
        B1 = calloc(N2 * N3, sizeof(double));
        old_fill_matrix(B1, N2, N3);
        C = calloc(N1 * N3, sizeof(double));

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
    int B_split_size = N3 / row_size;;
    double *A_recv1 = calloc(A_split_size * N2, sizeof(double));
    double *B_recv1 = calloc(B_split_size * N2, sizeof(double));

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
    // sleep(1 + world_rank);
//    for (int i = 0; i < col_size; i++) {
//        printf("%lf ", arr_recv_col[i]);
//    }
    //puts("\n");
    // sleep(1 + world_rank);
//     printf("RANK (%d,%d): ", coords[0], coords[1]);
//     old_print_matrix(A_recv1, A_split_size, N2);
//    puts("\n\n");
//    old_print_matrix(B_recv1, N2, B_split_size);
//     puts("\n\n");
    double *C_part = calloc(A_split_size * B_split_size, sizeof(double));
    matrix_multiply(A_recv1, A_split_size, N2, B_recv1, N2, B_split_size, C_part);
    //sleep(1 + world_rank);
    printf("World_rank %d ,RANK (%d,%d) , Displs: %d\n", world_rank, coords[0], coords[1],
           coords[0] * B_split_size + A_split_size * coords[1] * N3);

    for (int i = 0; i < A_split_size * B_split_size ; i++) {
        printf("%lf ", C_part[i]);
    }
    //  old_print_matrix(C_part, A_split_size, B_split_size);
    puts("\n");


    int *block_lengths = calloc(A_split_size, sizeof(int));
    long *matrix_displace = calloc(A_split_size, sizeof(int));
    MPI_Datatype *types = calloc(A_split_size, sizeof(MPI_Datatype));
    for (int i = 0; i < A_split_size; i++) {
        block_lengths[i] = B_split_size;
        matrix_displace[i] = (long) (i * N3 * sizeof(double));
        types[i] = MPI_DOUBLE;
    }

    MPI_Datatype structure;
    MPI_Type_create_struct(A_split_size, block_lengths, matrix_displace, types, &structure);
    MPI_Type_commit(&structure);

    MPI_Datatype matrix_back, matrix_back_resized;
    MPI_Type_vector(A_split_size, B_split_size, N3, MPI_DOUBLE, &matrix_back);
    MPI_Type_commit(&matrix_back);
    MPI_Type_create_resized(matrix_back, 0, (int) (B_split_size * sizeof(double)), &matrix_back_resized);
    MPI_Type_commit(&matrix_back_resized);


    // A_split_size = N1/dims[1];
    // B_split_size = N3/dims[0];
    int *displs;
    int *sizes;
    if (world_rank == RANK_ROOT) {
        displs = calloc(world_size, sizeof(int));
        sizes = calloc(world_size, sizeof(int));
        for (int i = 0; i < dims[1]; i++) {
            for (int j = 0; j < dims[0]; j++) {
                displs[i + j * dims[1]] = (j * B_split_size + A_split_size * i * N3);
            }
            //  displs[i] = (coords[1] * N3 + coords[0]) * B_split_size;
            //  displs[i] = coords[1] * B_split_size + A_split_size * coords[0] * N3;
        }
        for (int i = 0; i < world_size; i++) {
            sizes[i] = A_split_size * B_split_size;
        }
        for (int i = 0; i < world_size; i++) {
            printf("\nDispls %d, Sizes %d\n", displs[i], sizes[i]);
        }
    }
    MPI_Gatherv(C_part, 1, structure, C, sizes, displs, MPI_DOUBLE, RANK_ROOT, MPI_COMM_WORLD);

    if (world_rank == RANK_ROOT) {
        old_print_matrix(C, N1, N3);
    }
    MPI_Type_free(&matrix_back);
    MPI_Type_free(&matrix_back_resized);
    MPI_Type_free(&column);
    MPI_Type_free(&column_resized);
    MPI_Type_free(&structure);

    MPI_Finalize();
    return 0;
}



/* 2 2 2   1 | 1 | 1
 * _____
 * 3 3 3   2 | 2 | 2   ==
 * _____
 * 2 2 2   1 | 1 | 1*/
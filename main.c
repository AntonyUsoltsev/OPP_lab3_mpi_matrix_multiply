#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N1 3360
#define N2 3360
#define N3 3360

#define RANK_ROOT 0
#define X 0
#define Y 1

void fill_matrix(double *A, int height, int width) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (i == j)
                A[i * width + j] = 2.0f;
            else
                A[i * width + j] = 1.0f;
        }
    }
}

//void fill_matrix(double *A, int height, int width) {
//    for (int i = 0; i < height * width; i++) {
//        A[i] = i + 1;
//    }
//}


void file_print_matrix(double *A, int height, int width, char *file_name) {
    FILE *file = fopen(file_name, "w");
    fprintf(file, "[");
    for (int i = 0; i < height; i++) {
        fprintf(file, "[");
        for (int j = 0; j < width; j++) {
            fprintf(file, "%lf, ", A[i * width + j]);
        }
        fprintf(file, "],");
        fputs("\n", file);
    }
    fprintf(file, "]");
    fclose(file);
}

void
matrix_multiply(const double *A, int A_height, int A_width, const double *B, int B_height, int B_width, double *C) {
    if (A_width != B_height)
        return;
    for (int i = 0; i < A_height; i++) {
        for (int k = 0; k < A_width; k++) {
            for (int j = 0; j < B_width; j++) {
                C[i * B_width + j] += A[i * A_width + k] * B[k * B_width + j];
            }
        }
    }
}

MPI_Comm create_new_comm(int size, int rank, int *dims) {
    int periods[2] = {0, 0}, coords[2], reorder = 1;
    int sizey, sizex;
    MPI_Comm comm2d;

    MPI_Dims_create(size, 2, dims);
    sizex = dims[X];
    sizey = dims[Y];
    if (rank == RANK_ROOT)
        printf("Dimensions (%d, %d)\n", sizex, sizey);

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &comm2d);
    MPI_Comm_rank(comm2d, &rank);
    MPI_Cart_get(comm2d, 2, dims, periods, coords);
    //  printf("Process %d has coordinates (%d, %d)\n", rank, coords[0], coords[1]);
    return comm2d;
}

void
send_matrix(double *A_recv, double *B_recv, int A_split_size, int B_split_size, double *A, double *B, const int *coords,
            MPI_Comm row_comm, MPI_Comm col_comm) {

    MPI_Datatype column;
    MPI_Type_vector(N2, B_split_size, N3, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);
    MPI_Datatype column_resized;
    MPI_Type_create_resized(column, 0, (int) (B_split_size * sizeof(double)), &column_resized);
    MPI_Type_commit(&column_resized);

    if (coords[X] == 0) {
        MPI_Scatter(A, N2 * A_split_size, MPI_DOUBLE, A_recv, N2 * A_split_size, MPI_DOUBLE, 0, col_comm);
    }
    if (coords[Y] == 0) {
        MPI_Scatter(B, 1, column_resized, B_recv, N2 * B_split_size, MPI_DOUBLE, 0, row_comm);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(A_recv, N2 * A_split_size, MPI_DOUBLE, RANK_ROOT, row_comm);
    MPI_Bcast(B_recv, N2 * B_split_size, MPI_DOUBLE, RANK_ROOT, col_comm);

    MPI_Type_free(&column);
    MPI_Type_free(&column_resized);

}

void receive_matrix(int A_split_size, int B_split_size, int world_rank, int world_size, const int *dims, double *C,
                    const double *C_part, MPI_Comm comm2d) {
    MPI_Datatype send_matrix_part, send_matrix_part_resized;
    MPI_Type_vector(A_split_size, B_split_size, N3, MPI_DOUBLE, &send_matrix_part);
    MPI_Type_commit(&send_matrix_part);
    MPI_Type_create_resized(send_matrix_part, 0, (int) (B_split_size * sizeof(double)), &send_matrix_part_resized);
    MPI_Type_commit(&send_matrix_part_resized);

    int *displs = NULL, *sizes = NULL;
    if (world_rank == RANK_ROOT) {
        displs = calloc(world_size, sizeof(int));
        sizes = calloc(world_size, sizeof(int));
        for (int i = 0; i < dims[Y]; i++) {
            for (int j = 0; j < dims[X]; j++) {
                displs[i + j * dims[Y]] = (j * B_split_size + A_split_size * i * N3) / B_split_size;
            }
        }
        for (int i = 0; i < world_size; i++) {
            sizes[i] = 1;
        }
    }
    MPI_Gatherv(C_part, B_split_size * A_split_size, MPI_DOUBLE, C, sizes, displs, send_matrix_part_resized, RANK_ROOT,
                comm2d);
    MPI_Type_free(&send_matrix_part);
    MPI_Type_free(&send_matrix_part_resized);
    if (world_rank == RANK_ROOT) {
        free(displs);
        free(sizes);
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    double *A = NULL, *B = NULL, *C = NULL;
    if (world_rank == RANK_ROOT) {
        A = calloc(N1 * N2, sizeof(double));
        fill_matrix(A, N1, N2);
        //    file_print_matrix(A, N1, N2,"./matrix_A.txt");
        B = calloc(N2 * N3, sizeof(double));
        fill_matrix(B, N2, N3);
        //  file_print_matrix(B, N2, N3,"./matrix_B.txt");
        C = calloc(N1 * N3, sizeof(double));
    }

    double start = MPI_Wtime();

    int dims[2] = {0, 0}, periods[2] = {0, 0}, coords[2]; //x y
    MPI_Comm comm2d = create_new_comm(world_size, world_rank, dims);

    MPI_Cart_get(comm2d, 2, dims, periods, coords);
    MPI_Comm row_comm = NULL, col_comm = NULL; //x y
    MPI_Comm_split(comm2d, coords[0], coords[1], &col_comm);
    MPI_Comm_split(comm2d, coords[1], coords[0], &row_comm);

    int row_rank, row_size;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_size(row_comm, &row_size);
    int col_rank, col_size;
    MPI_Comm_rank(col_comm, &col_rank);
    MPI_Comm_size(col_comm, &col_size);

    if (row_size != dims[X] || col_size != dims[Y]) {
        perror("Bad size of topology");
        exit(1);
    }
    if (N1 % col_size != 0 || N3 % row_size != 0) {
        perror("Bad size of matrix");
        exit(1);
    }

    int A_split_size = N1 / dims[Y];
    int B_split_size = N3 / dims[X];
    double *A_recv = calloc(A_split_size * N2, sizeof(double));
    double *B_recv = calloc(B_split_size * N2, sizeof(double));

    send_matrix(A_recv, B_recv, A_split_size, B_split_size, A, B, coords, row_comm, col_comm);

    double *C_part = calloc(A_split_size * B_split_size, sizeof(double));
    matrix_multiply(A_recv, A_split_size, N2, B_recv, N2, B_split_size, C_part);

    receive_matrix(A_split_size, B_split_size, world_rank, world_size, dims, C, C_part, comm2d);

    double end = MPI_Wtime();

    if (world_rank == RANK_ROOT) {
//        file_print_matrix(C, N1, N3,"./matrix_C.txt");
//        char script[100];
//        sprintf(script, "/mnt/c/'Python 3.8.2'/python.exe ./check_matrix.py");
//        if (system(script) != 0) {
//            perror("Script didn't run");
//        }
        printf("%lf sec\n", end - start);
        free(A);
        free(B);
        free(C);
    }
    free(A_recv);
    free(B_recv);
    MPI_Finalize();
    return 0;
}
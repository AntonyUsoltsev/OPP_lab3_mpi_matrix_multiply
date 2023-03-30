# OPP_lab3_mpi_matrix_multiply
Compile command:

    mpicc main.c -o main -O3 --std=c99

Run command:

    mpirun -oversubscribe -np $proc_count ./main

import sys

import numpy


def fill_matrix_2(height, width, matrix):
    k = 1.0
    for i in range(0, height):
        tmp_matr = []
        for j in range(0, width):
            tmp_matr.append(k)
            k += 1.0
        matrix[i] = (tmp_matr)


def fill_matrix_1(height, width, matrix):
    for i in range(0, height):
        tmp_matr = []
        for j in range(0, width):
            if i == j:
                tmp_matr.append(2.0)
            else:
                tmp_matr.append(1.0)
        matrix[i] = (tmp_matr)


argv = sys.argv
n1, n2, n3, type = map(int, argv[1::])

matrix1 = [[]] * n1
matrix2 = [[]] * n2

if type == 1:
    fill_matrix_1(n1, n2, matrix1)
    fill_matrix_1(n1, n2, matrix2)
if type == 2:
    fill_matrix_2(n1, n2, matrix1)
    fill_matrix_2(n1, n2, matrix2)

matr_mul = (numpy.matmul(matrix1, matrix2))
file = open("matrix.txt", "r")
str = file.read()
c_matrix = eval(str)
c_matrix = numpy.array(c_matrix)

print((c_matrix == matr_mul).all())

import numpy

file = open("matrix_A.txt", "r")
str = file.read()
A_matrix = eval(str)
A_matrix = numpy.array(A_matrix)

file = open("matrix_B.txt", "r")
str = file.read()
B_matrix = eval(str)
B_matrix = numpy.array(B_matrix)


matr_mul = (numpy.matmul(A_matrix, B_matrix))
file = open("matrix_C.txt", "r")
str = file.read()
C_matrix = eval(str)
C_matrix = numpy.array(C_matrix)

print((C_matrix == matr_mul).all())

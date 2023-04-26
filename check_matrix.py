import numpy

file = open("matrix_A.txt", "r")
data = file.read()
A_matrix = eval(data)
A_matrix = numpy.array(A_matrix)

file = open("matrix_B.txt", "r")
data = file.read()
B_matrix = eval(data)
B_matrix = numpy.array(B_matrix)

matrix_mul = (numpy.matmul(A_matrix, B_matrix))
file = open("matrix_C.txt", "r")
data = file.read()
C_matrix = eval(data)
C_matrix = numpy.array(C_matrix)

print((C_matrix == matrix_mul).all())

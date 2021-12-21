# Алгоритм взят с сайта https://www.kazedu.kz/referat/200522
# 1) Читаем матрицу  В и А из файлов и проверяем последнюю на симметричность и на ненулевой определитель
# 2) По специальным формулам вычисляем вспомогательную треугольную матрицу Т (T * T^T = А)
# 3) Вводим вспомогательный вектор Y, той же размерности, что и матрица А 
# 4) T^T*Y = В, T*X = Y, где T^T - транспонированная матрица T, В - матрица В
# 5) По специальным формулам находим Y, а потом X, последний записываем в файл
# 6) Генерируем матрицы Гильберта для n = 2, 3, 4, 5, 6, 7, 8, для каждой находим число обусловленности
# 7) Число обусловленности матрицы я взял как произведение  ее нормы и нормы матрицы, обратной к ней.
# 8) Норму матрицы я определил как максимальный элемент из списка, который образован
# путем построчного сложения модулей элементов матрицы.
# Важное дополнение! Если А[0][0] = 0, то при поиске матрицы T возникает деление на 0.
# Чтобы решить проблему, ищется ненулевой элемент на главной диагонали и он ставится на место A[0][0]
# путем перестановки 2-х строк и 2-х столбцов.
# Если на главной диагонали нет ненулевых чисел, в таком случае этот алгоритм не применим.
import copy
import math
import cmath
import sys
import numpy as np
from functools import reduce

def cut_matrix (matrix, index_column, index_string): # возвращает матрицу без index_column-го столбца и index_string-й строки
    new_matrix = copy.deepcopy(matrix) # создаем копию матрицы
    k = len(new_matrix[0]) # определяем размерность матрицы
    for i in range(k): # цикл от о до размерности - 1, который удаляет колонку
        del new_matrix[i][index_column] 
    del new_matrix[index_string] # удаляем строку
    return new_matrix


def make_gilbert_matrix (n): # Метод возвращает матрицу Гильберта размерности n
    matrix = []
    for i in range(n):
        matrix.append([])
        for j in range(n):
            matrix[i].append(1 / (i + j + 1))
    return matrix
    

def calculate_determinant (matrix): # вычисляет определитель
    if len(matrix) == 1:
        return matrix[0][0]
    if len(matrix[0]) == 2:  # если размерность  матрицы = 2, то вычислить просто:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    else:
        minor = 0;
        for i in range(len(matrix[0])): # иначе раскладываем по первой строке
            minor = minor + matrix[0][i] * calculate_determinant(cut_matrix(matrix, i, 0)) * ((-1) ** i)
        return minor


def make_inverse_matrix (matrix): # возвращает обратную матрицу
    inverse_matrix = np.zeros((len(matrix), len(matrix)))
    for i in range(len(matrix[0])):
        for j in range(len(matrix[0])):
            inverse_matrix[i][j] = calculate_determinant(cut_matrix(matrix, j, i)) * ((-1) ** (i + j))
    # осталось траснспонировать и умножить на 1/Determinant
    temp_matrix = copy.deepcopy(inverse_matrix)
    determinant = calculate_determinant(matrix)
    for i in range(len(matrix[0])):
        for j in range(len(matrix[0])):
            inverse_matrix[i][j] = temp_matrix[j][i] / determinant
    return inverse_matrix


def get_column(matrix, index):
    result = []
    for i in range(len(matrix)):
        result.append(matrix[i][index])
    return result


def print_matrix_in_file(matrix, file):
    string = ""
    for i in range(len(matrix)): #  выводим матрицу
        for j in range(len(matrix)):
            if type(matrix[i][j]) == np.complex128:
                if abs(matrix[i][j].imag) <=  sys.float_info.epsilon:
                    string = str(round(matrix[i][j].real, 5))
                if abs(matrix[i][j].real) <=  sys.float_info.epsilon:
                    string = str(round(matrix[i][j].imag, 5)) + "j"           
                if abs(matrix[i][j].imag) <=  sys.float_info.epsilon \
                       and abs(matrix[i][j].real) <=  sys.float_info.epsilon:
                    string = "0"
            else:
                string = str(round(matrix[i][j], 5))
            file.write(string)
            for k in range (len(string), 15):
                file.write(" ")
        file.write("\n")


def is_symmetrical (matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if matrix[i][j] != matrix[j][i]:
                return 0
    return 1


def delete_column(matrix, index):
    result = copy.deepcopy(matrix)
    for i in range(len(matrix)):
        del result[i][index]
    return result

def multiply_matrix (matrix1, matrix2):
    rez = np.zeros((len(matrix1),len(matrix1)), dtype=np.complex)  
    for i in range(len(matrix1)):
        for j in range(len(matrix1)):
            rez[i][j] = reduce(lambda X, Y: X + Y, [x*y for x,y in zip(matrix1[i], get_column(matrix2, j))]) 
    return rez

def make_triangular_matrix(matrix_a):
    triangular_matrix = np.zeros((len(matrix_a),len(matrix_a)), dtype = np.complex)
        # вычисляем треугольную матрицу, по тем формулам
    for i in range(len(triangular_matrix)):
        for j in range (len(triangular_matrix)):
            if i == 0 and j == 0:
                triangular_matrix[i][j] = cmath.sqrt(matrix_a[i][j])             
            if i == 0 and j > 0:
                triangular_matrix[0][j] = matrix_a[0][j] / triangular_matrix[0][0]
            if i > 0 and i == j:
                summ = 0
                for k in range(i):
                   summ += (triangular_matrix[k][i] ** 2)
                triangular_matrix[i][i] = cmath.sqrt(matrix_a[i][i] - summ) 
            if i < j and i > 0:
                summ = 0
                for k in range(i):
                    summ += (triangular_matrix[k][i] * triangular_matrix[k][j])
                triangular_matrix[i][j] = (matrix_a[i][j] - summ) / triangular_matrix[i][i]
            if i > j:
                triangular_matrix[i][j] = 0
    return triangular_matrix


if __name__ == "__main__":
    matrix_a = []
    matrix_b = []
    output_file_name = "output.txt"
    with open("input.txt") as file_a: #  вводим матрицу А из файла
        for line in file_a:
            matrix_a.append([float(x) for x in line.split()])

    matrix_b = get_column(matrix_a, len(matrix_a))
    matrix_a = delete_column(matrix_a, len(matrix_a[0]) - 1)

    with open(output_file_name,'w') as f:
        f.write("Введенная матрица:\n")
        print_matrix_in_file(matrix_a, f)
    
    is_correct = 1
    if is_symmetrical(matrix_a) == 0:
        with open(output_file_name,'a') as f:
            f.write("Матрица не симметричная :( \n")
        is_correct = 0
    
    if calculate_determinant(matrix_a) == 0:
        with open(output_file_name,'a') as f:
            f.write("У матрицы А нулевой определитель :( \n")
        is_correct = 0

    backup_a = copy.deepcopy(matrix_a)
    is_swapped = 0           #  флаг перестановки 
    index_for_swapped = -1    #  индекс  строки и столбца с теми, с кем свапаем 

    if matrix_a[0][0] == 0:                  #  проверяем первый элемент на главной диагонали  на ноль
        for i in range(1, len(matrix_a)):    # если да, находим  ненулевой элемент на главной диагонали
            if matrix_a[i][i] != 0:          #  нашли
                matrix_a[0], matrix_a[i] = matrix_a[i], matrix_a[0] # свапаем строки
                matrix_b[0], matrix_b[i] = matrix_b[i], matrix_b[0] # свапаем элементы в векторе В
                for j in range(len(matrix_a)):                   #  свапаем столбцы
                    matrix_a[j][0], matrix_a[j][i] = matrix_a[j][i], matrix_a[j][0]
                is_swapped = 1               #  указываем, что перестановка была выполнена
                index_for_swapped = i        #  и индекс
                break

    if is_swapped == 1:
        with open(output_file_name,'a') as f:
            f.write("A[0][0] = 0, поэтому переставляем строки и столбцы, " + \
                    "чтобы исправить это, сохраняя симметричность:\n")
            print_matrix_in_file(matrix_a, f)

    if matrix_a[0][0] == 0:
        with open(output_file_name,'a') as f:
            f.write("Эту систему нельзя решить таким методом :( \n")
        is_correct = 0
    if is_correct == 1:
        triangular_matrix = make_triangular_matrix(matrix_a)
        # вычисляем треугольную матрицу, по тем формулам
        matrix_y = []
        for i in range(len(triangular_matrix)):
            if i == 0:
                matrix_y.append(matrix_b[0] / triangular_matrix[0][0])
            else:
                summ = 0
                for k in range(i):
                    summ += triangular_matrix[k][i] * matrix_y[k]
                matrix_y.append((matrix_b[i] - summ) / triangular_matrix[i][i])
            #  ищем решение
        matrix_x = [0 for _ in range(len(triangular_matrix))] # заполняме вектор Х нулями
        for i in range(len(triangular_matrix)- 1, -1, -1): # вычисляем Х
            if i == len(triangular_matrix)- 1:
                matrix_x[i] = matrix_y[i] / triangular_matrix[i][i]
            else:
                summ = 0;
                for k in range(i + 1, len(triangular_matrix)):
                    summ += triangular_matrix[i][k] * matrix_x[k]
                matrix_x[i] = (matrix_y[i] - summ) / triangular_matrix[i][i]

        # решение найдено, выводим результаты:
        with open(output_file_name,'a') as f:
            f.write("Треугольная матрица: \n")
            print_matrix_in_file(triangular_matrix, f)
            f.write("Транспонированная треугольная матрица: \n")
            print_matrix_in_file(np.transpose(triangular_matrix),  f)
            f.write("Произведение транспонированной и обычной: \n")
            print_matrix_in_file(multiply_matrix(np.transpose(triangular_matrix), triangular_matrix), f)
            f.write("Вспомогательный вектор Y: \n")
            for i in range(len(matrix_y)):
                f.write("Y" + str(i + 1) + " = " + str(matrix_y[i]) + "\n")
            f.write('Решение системы уравнений: \n')
            if is_swapped == 1:
                matrix_x[0], matrix_x[index_for_swapped] = matrix_x[index_for_swapped], matrix_x[0]
                matrix_b[0], matrix_b[index_for_swapped] = matrix_b[index_for_swapped], matrix_b[0]
                matrix_a = copy.deepcopy(backup_a)
            for i in range(len(matrix_x)):
                f.write("X"+ str(i + 1) +  " = " + str(round(matrix_x[i].real, 5)) + '\n')

            # сделаю проверку:

            f.write("Сделаем проверку найденного решения: \n")
            for i in range(len(matrix_a)):
                string = ""
                for j in range(len(matrix_a)):
                    string += "( " + str(round(matrix_a[i][j].real, 5)) + " * " \
                                   + str(round(matrix_x[j].real, 5)) + ") + "
                string = string[:-2]
                Summ = reduce(lambda X, Y: X + Y, [x*y for x,y in zip(matrix_a[i],matrix_x)])
                string += " = " + str(Summ.real) + " | "  + str(matrix_b[i]) + "\n"
                f.write(string)              

            # делаем дополнительное задание, результат вычислений - в файле. 
            f.write("Исследование зависимости матрицы Гильберта и числа обусловленности:\n")
            for n in range (2, 8): #  цикл по размерностым матрицы
                f.write("n = " + str(n)  + ". Матрица: \n")
                gilbert_matrix = make_gilbert_matrix(n) #  делаем матрицу
                print_matrix_in_file(gilbert_matrix, f)
                inversed_gilbert_matrix = make_inverse_matrix(gilbert_matrix) #  делаем обратную
                summ_str_in_gilbert_matrix = []
                summ_str_in_inversed_gilbert_matrix = []
                for i in range(n): # ищем число обусловленности
                    summ_str_in_inversed_gilbert_matrix.append(reduce(lambda x, y: \
                        abs(x) + abs(y), inversed_gilbert_matrix[i]))

                    summ_str_in_gilbert_matrix.append(reduce(lambda x, y: abs(x) + abs(y), gilbert_matrix[i]))
                f.write("Число: " + str(max(summ_str_in_gilbert_matrix) *\
                        max(summ_str_in_inversed_gilbert_matrix)) + "\n") # выводим его
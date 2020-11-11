import sys

from newtonraphson import Newtonraphson
from desacoplado import Newtondesacoplado
from desacopladomod import Desacopladomod
from fastnewton import Fastnewton
from gaussseidel import Gaussseidel
from gaussjacobi import Gaussjacobi


class Solution:
    def __init__(self, y_matrix, g_matrix, b_matrix, data_bus, line_data):
        self.y_matrix = y_matrix
        self.g_matrix = g_matrix
        self.b_matrix = b_matrix
        self.bus_data = data_bus
        self.line_data = line_data

    def get_solution_method(self):
        print("Métodos de solução do fluxo de carga:")
        print("[0] Newton Raphson")
        print("[1] Newton Desacoplado")
        print("[2] Newton Desacoplado Modificado (divide a eq do fluxo de potência por V[k])")
        print("[3] Newton Desacoplado Rápido BB")
        print("[4] Newton Desacoplado Rápido XB - mais usado")
        print("[5] Newton Desacoplado Rápido BX")
        print("[6] Newton Desacoplado Rápido XX")
        print("[7] Gauss-Seidel")
        print("[8] Gauss-Jacobi")
        print("[9] Fluxo de carga linearizado")

        method = int(input("Como deseja solucionar o fluxo de carga? "))

        return method

    def get_solution(self, method):

        if method == 0:
            solution = Newtonraphson(self.y_matrix, self.g_matrix, self.b_matrix, self.bus_data)
            solution.get_newtonraphson_solution()

        elif method == 1:
            solution = Newtondesacoplado(self.y_matrix, self.g_matrix, self.b_matrix, self.bus_data)
            solution.get_newtondesacoplado_solution()

        elif method == 2:
            solution = Desacopladomod(self.y_matrix, self.g_matrix, self.b_matrix, self.bus_data)
            solution.get_desacopladomod_solution()

        elif method == 3:
            solution = Fastnewton(self.y_matrix, self.g_matrix, self.b_matrix, self.bus_data,
                                  self.line_data, method)
            solution.get_fastnewton_solution()

        elif method == 4:
            solution = Fastnewton(self.y_matrix, self.g_matrix, self.b_matrix, self.bus_data,
                                  self.line_data, method)
            solution.get_fastnewton_solution()

        elif method == 5:
            solution = Fastnewton(self.y_matrix, self.g_matrix, self.b_matrix, self.bus_data,
                                  self.line_data, method)
            solution.get_fastnewton_solution()

        elif method == 6:
            solution = Fastnewton(self.y_matrix, self.g_matrix, self.b_matrix, self.bus_data,
                                  self.line_data, method)
            solution.get_fastnewton_solution()

        elif method == 7:
            solution = Gaussseidel(self.y_matrix, self.bus_data)
            solution.get_gaussseidel_solution()

        elif method == 8:
            solution = Gaussjacobi(self.y_matrix, self.bus_data)
            solution.get_gaussjacobi_solution()

        elif method == 9:
            pass

        else:
            print("Método de solução não existente. Encerrando processo.")
            sys.exit()

    def get_analise_option(self):
        print("Opções de análise:")
        print("[0] Aplicar valores de Base")
        print("[1] Perdas nas linhas")
        print("[2] Correntes nas linhas")
        option = int(input("O que deseja analisar? "))

        return option

    def get_analise(self, option):
        pass


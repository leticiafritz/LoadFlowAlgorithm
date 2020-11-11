import cmath
import math

import numpy as np


class Gaussseidel:

    def __init__(self, y_matrix, bus_data):
        self.y_matrix = y_matrix
        self.bus_data = bus_data

    def get_init_data(self):
        v = np.zeros(max(self.bus_data.iloc[:, 0]), dtype=complex)
        p = np.zeros(max(self.bus_data.iloc[:, 0]))
        q = np.zeros(max(self.bus_data.iloc[:, 0]))
        code = np.zeros(max(self.bus_data.iloc[:, 0]))  # 0-> slack, 1-> PV, 2->PQ
        print("Dados de Entrada")
        for n in range(len(self.bus_data)):
            if self.bus_data.iloc[n, 1] == "slack":
                v[n], code[n] = cmath.rect(self.bus_data.iloc[n, 2], self.bus_data.iloc[n, 3]), 0
                print("Barra ", n + 1, " (", self.bus_data.iloc[n, 1], "): ",
                      "Ê = ", v[n], " pu // Incógnitas: P e Q")
                v_init = v[n]
            elif self.bus_data.iloc[n, 1] == "PQ":
                p[n] = float(self.bus_data.iloc[n, 4]) - float(self.bus_data.iloc[n, 6])
                q[n] = float(self.bus_data.iloc[n, 5]) - float(self.bus_data.iloc[n, 7])
                code[n] = 2
                print("Barra ", n + 1, " (", self.bus_data.iloc[n, 1], "): ",
                      "P = ", "%.4f" % p[n], " pu, Q = ", "%.4f" % q[n], " pu, // Incógnitas: V e teta")

        return v, v_init, p, q, code

    def get_incognitas_vector(self, code, v_init):
        variable_vector = []
        print("Valores iniciais: ")
        i = 0
        for n in range(len(code)):
            if code[n] == 2:
                variable_vector.append(v_init)
                print("Ê", self.bus_data.iloc[n, 0], ": ", variable_vector[i])
                i += 1

        return variable_vector

    def get_gaussseidel_solution(self):
        print("\n")
        print("Iniciando o método de Gaus-Jacobi ...")

        # Organizando os dados
        v, v_init, p, q, code = self.get_init_data()

        # Vetor de incognitas
        variable_vector = self.get_incognitas_vector(code, v_init)
        i = 0
        for k in range(len(v)):
            if v[k] != 0:
                v[k] = v[k]
            else:
                v[k] = variable_vector[i]
                i += 1

        # Processo iterativo
        tolerance = float(input("Qual a tolerância? "))
        iteration = 0
        control = tolerance + 1
        v_old = v

        while control > tolerance:
            v_current = np.zeros(len(v), dtype=complex)
            for k in range(len(code)):
                if k == 0:
                    v_current[k] = v_init
            print("--------------------------------------------------")
            print("Iteração: ", iteration)
            for k in range(len(code)):
                sum_ye_old = 0
                sum_ye_current = 0
                if code[k] == 2:
                    for m in range(len(code)):
                        if m < k:
                            sum_ye_current = sum_ye_current + self.y_matrix[k, m] * v_current[m]
                        elif m > k:
                            sum_ye_old = sum_ye_old + self.y_matrix[k, m] * v_old[m]
                    v_current[k] = (1 / self.y_matrix[k, k]) * \
                                   (((p[k] - q[k] * 1j) / (v_old[k].real - v_old[k].imag * 1j)) -
                                    sum_ye_old - sum_ye_current)
                    print("Barra", k + 1, ": Ê = ", np.abs(v_current[k]), ",", math.degrees(np.angle(v_current[k])),
                          "º")
            # Incluindo a tensão da barra slack no vetor de tensão
            for k in range(len(code)):
                if k == 0:
                    v_current[k] = v_init
            # Descobrindo corrente da barra slack
            current = 0
            for k in range(len(code)):
                if code[k] == 0:
                    for m in range(len(code)):
                        current = current + (self.y_matrix[k, m] * v_current[m])
            print("Corrente da Barra slack: ", current, "pu")
            # Descobrindo potência aparente da barra slack
            for k in range(len(code)):
                if code[k] == 0:
                    s = v_current[k] * (current.real - current.imag * 1j) - self.bus_data.iloc[k, 6] - \
                        self.bus_data.iloc[k, 7] * 1j
            print("Potência aparente gerada da Barra slack: ", s, "pu")
            # Testando controle
            v_control = v_current - v_old
            control = max(abs(v_control))
            v_old = v_current
            iteration += 1
        # Resultado:
        print("--------------------------------------------------")
        print("Processo convergiu com", iteration - 1, " iterações, pelo método de Gauss Seidel.")
        print("Solução: ")
        print("Variáveis de estado:")
        for k in range(len(code)):
            print("Barra", k + 1, "(", self.bus_data.iloc[k, 1], ") Ê =", np.abs(v_current[k]),
                  "[pu],", math.degrees(np.angle(v_current[k])), "[º]")
        print("Potência líquida da barra slack: Sg = ", s, "[pu]")
        print("--------------------------------------------------")
        print("\n")

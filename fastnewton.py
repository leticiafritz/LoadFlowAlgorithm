import math

import pandas as pd
import numpy as np

from systemmatrixsolve import Systemmatrixsolve


class Fastnewton:
    def __init__(self, y_matrix, g_matrix, b_matrix, data_bus, data_line, type):
        self.y_matrix = y_matrix
        self.g_matrix = g_matrix
        self.b_matrix = b_matrix
        self.data_bus = pd.DataFrame(data_bus)
        self.data_line = pd.DataFrame(data_line)
        self.type = type

    def get_init_data(self):
        v = np.zeros(max(self.data_bus.iloc[:, 0]))
        teta = np.zeros(max(self.data_bus.iloc[:, 0]))
        p = np.zeros(max(self.data_bus.iloc[:, 0]))
        q = np.zeros(max(self.data_bus.iloc[:, 0]))
        code = np.zeros(max(self.data_bus.iloc[:, 0]))  # 0-> slack, 1-> PV, 2->PQ
        print("Dados de Entrada")
        for n in range(len(self.data_bus)):
            if self.data_bus.iloc[n, 1] == "slack":
                v[n], teta[n], code[n] = self.data_bus.iloc[n, 2], math.radians(self.data_bus.iloc[n, 3]), 0
                print("Barra ", n + 1, " (", self.data_bus.iloc[n, 1], "): ",
                      "V = ", "%.4f" % v[n], " pu, teta = ", "%.4f" % teta[n], " rad, // Incógnitas: P e Q")
                teta_init = teta[n]
                v_init = v[n]
            elif self.data_bus.iloc[n, 1] == "PV":
                v[n] = self.data_bus.iloc[n, 2]
                p[n] = float(self.data_bus.iloc[n, 4]) - float(self.data_bus.iloc[n, 6])
                code[n] = 1
                print("Barra ", n + 1, " (", self.data_bus.iloc[n, 1], "): ",
                      "V = ", "%.4f" % v[n], " pu, P = ", "%.4f" % p[n], " pu, // Incógnitas: Q e teta")
            elif self.data_bus.iloc[n, 1] == "PQ":
                p[n] = float(self.data_bus.iloc[n, 4]) - float(self.data_bus.iloc[n, 6])
                q[n] = float(self.data_bus.iloc[n, 5]) - float(self.data_bus.iloc[n, 7])
                code[n] = 2
                print("Barra ", n + 1, " (", self.data_bus.iloc[n, 1], "): ",
                      "P = ", "%.4f" % p[n], " pu, Q = ", "%.4f" % q[n], " pu, // Incógnitas: V e teta")

        return v, v_init, teta, teta_init, p, q, code

    def get_incognitas_vector(self, code, teta_init, v_init):
        variable_teta_vector = []
        variable_v_vector = []
        p_calc_vector = []
        q_calc_vector = []
        mismatches_teta_vector = []
        mismatches_v_vector = []
        mismatches_p_vector = []
        mismatches_q_vector = []
        print("Valores iniciais: ")
        i = 0
        for n in range(len(code)):
            if code[n] == 1 or code[n] == 2:
                variable_teta_vector.append(teta_init)
                p_calc_vector.append(0)
                print("teta", self.data_bus.iloc[n, 0], ": ", "%.4f" % variable_teta_vector[i])
                mismatches_teta_vector.append(teta_init)
                mismatches_p_vector.append(0)
                i += 1
        i = 0
        for n in range(len(code)):
            if code[n] == 2:
                variable_v_vector.append(v_init)
                q_calc_vector.append(0)
                print("V", self.data_bus.iloc[n, 0], ": ", "%.4f" % variable_v_vector[i])
                mismatches_v_vector.append(v_init)
                mismatches_q_vector.append(0)
                i += 1

        return variable_teta_vector, variable_v_vector, mismatches_teta_vector, mismatches_v_vector, \
               mismatches_p_vector, mismatches_q_vector

    def get_power_balance(self, code, teta, v, variable_vector, p, q, mismatches_power_vector, ref):
        # Modificando v e teta para incluir o chute inicial
        i = 0
        for k in range(len(code)):
            if code[k] == 1 or code[k] == 2:
                teta[k] = variable_vector[i]
                i += 1
        for k in range(len(code)):
            if code[k] == 2:
                v[k] = variable_vector[i]
                i += 1

        # Balanço de Potência
        p_calc_vector = np.zeros(len(code))
        q_calc_vector = np.zeros(len(code))
        mismatches_power_p = np.zeros(len(code))
        mismatches_power_q = np.zeros(len(code))
        for k in range(len(code)):
            p_sum = 0
            q_sum = 0
            for m in range(len(code)):
                p_sum = p_sum + v[m] * (self.g_matrix[k, m] * math.cos(teta[k] - teta[m]) +
                                        self.b_matrix[k, m] * math.sin(teta[k] - teta[m]))
                q_sum = q_sum + v[m] * (self.g_matrix[k, m] * math.sin(teta[k] - teta[m]) -
                                        self.b_matrix[k, m] * math.cos(teta[k] - teta[m]))
            p_calc_vector[k] = v[k] * p_sum
            q_calc_vector[k] = v[k] * q_sum
            mismatches_power_p[k] = p[k] - p_calc_vector[k]
            mismatches_power_q[k] = q[k] - q_calc_vector[k]

        # Organizando vetor de mismatches
        i_mismatche = 0
        for k in range(len(code)):
            if code[k] == 1 or code[k] == 2:
                if ref == 0:
                    mismatches_power_vector[i_mismatche] = mismatches_power_p[k]
                    print("Barra", k + 1, " Pesp = ", "%.4f" % p[k], "Pcal = ", p_calc_vector[k],
                          ", mismatche Pesp-Pcalc = ", "%.4f" % mismatches_power_vector[i_mismatche])
                i_mismatche += 1
        for k in range(len(code)):
            if code[k] == 2:
                if ref == 1:
                    mismatches_power_vector[i_mismatche] = mismatches_power_q[k]
                    print("Barra", k + 1, " Qesp = ", "%.4f" % q[k], " Qesp = ", "%.4f" % q_calc_vector[k],
                            ", mismatche Qesp-Qcalc = ", "%.4f" % mismatches_power_vector[i_mismatche])
                i_mismatche += 1

        return teta, v, p_calc_vector, q_calc_vector, mismatches_power_vector

    def get_b_linha_matrix(self, code):
        npq = list(code).count(2)
        npv = list(code).count(1)
        b_linha = np.zeros((npq + npv, npq + npv))

        row = 0
        col = 0
        if self.type == 3 or self.type == 5:
            for k in range(len(code)):
                for m in range(len(code)):
                    if k == m:
                        if code[k] == 1 or code[k] == 2:
                            b_linha[row, col] = - self.b_matrix[k, k]
                            if col < (b_linha.shape[1] - 1):
                                col += 1
                            else:
                                col = 0
                                row += 1
                    else:
                        if code[k] == 1 or code[k] == 2:
                            if code[m] == 1 or code[m] == 2:
                                b_linha[row, col] = - self.b_matrix[k, m]
                                if col < (b_linha.shape[1] - 1):
                                    col += 1
                                else:
                                    col = 0
                                    row += 1
        elif self.type == 4 or self.type == 6:
            for k in range(len(code)):
                for m in range(len(code)):
                    if k == m:
                        if code[k] == 1 or code[k] == 2:
                            b_linha[row, col] = 0
                            for m_aux in range(len(code)):
                                if (self.data_line.iloc[m_aux, 0] == k) or \
                                        (self.data_line.iloc[m_aux, 1] == k):
                                    b_linha[row, col] = b_linha[row, col] + \
                                                        (self.data_line.iloc[k, m_aux] ** (-1))
                            if col < (b_linha.shape[1] - 1):
                                col += 1
                            else:
                                col = 0
                                row += 1
                    else:
                        if code[k] == 1 or code[k] == 2:
                            if code[m] == 1 or code[m] == 2:
                                b_linha[row, col] = -  (self.data_line.iloc[k, m] ** (-1))
                                if col < (b_linha.shape[1] - 1):
                                    col += 1
                                else:
                                    col = 0
                                    row += 1

        return b_linha

    def get_b_2linha_matrix(self, code):
        npq = list(code).count(2)
        b_2linha = np.zeros((npq, npq))

        row = 0
        col = 0
        if self.type == 3 or self.type == 4:
            for k in range(len(code)):
                for m in range(len(code)):
                    if k == m:
                        if code[k] == 2:
                            b_2linha[row, col] = - self.b_matrix[k, k]
                            if col < (b_2linha.shape[1] - 1):
                                col += 1
                            else:
                                col = 0
                                row += 1
                    else:
                        if code[k] == 2:
                            if code[m] == 2:
                                b_2linha[row, col] = - self.b_matrix[k, m]
                                if col < (b_2linha.shape[1] - 1):
                                    col += 1
                                else:
                                    col = 0
                                    row += 1
        elif self.type == 5 or self.type == 6:
            for k in range(len(code)):
                for m in range(len(code)):
                    if k == m:
                        if code[k] == 2:
                            b_2linha[row, col] = 0
                            for m_aux in range(len(code)):
                                if (self.data_line.iloc[m_aux, 0] == k) or \
                                        (self.data_line.iloc[m_aux, 1] == k):
                                    b_2linha[row, col] = b_2linha[row, col] + \
                                                        (self.data_line.iloc[k, m_aux] ** (-1))
                            if col < (b_2linha.shape[1] - 1):
                                col += 1
                            else:
                                col = 0
                                row += 1
                    else:
                        if code[k] == 2:
                            if code[m] == 2:
                                b_2linha[row, col] = -  (self.data_line.iloc[k, m] ** (-1))
                                if col < (b_2linha.shape[1] - 1):
                                    col += 1
                                else:
                                    col = 0
                                    row += 1
        return b_2linha

    def get_fastnewton_solution(self):
        print("\n")
        print("Iniciando o método de Newton Desacoplado ...")

        # Organizando os dados
        v, v_init, teta, teta_init, p, q, code = self.get_init_data()

        # Vetor de incognitas
        variable_teta_vector, variable_v_vector, mismatches_teta_vector, mismatches_v_vector, \
        mismatches_p_vector, mismatches_q_vector = self.get_incognitas_vector(code, teta_init, v_init)
        variable_vector = variable_teta_vector + variable_v_vector
        mismatches_power_vector = mismatches_p_vector + mismatches_q_vector

        # Processo iterativo
        tolerance = float(input("Qual a tolerância? "))
        iteration_p, iteration_q = 0, 0
        control_p, control_q = tolerance + 1, tolerance + 1

        while control_p > tolerance or control_q > tolerance:
            if control_p > tolerance:
                print("--------------------------------------------------")
                print("Meia Iteração Ativa: ", iteration_p)
                # Balanço de potência
                teta, v, p_calc_vector, q_calc_vector, mismatches_power_vector = \
                    self.get_power_balance(code, teta, v, variable_vector, p, q, mismatches_power_vector, 0)
                # Matriz B'
                b_linha = self.get_b_linha_matrix(code)
                print("Matriz B':")
                print(np.around(b_linha, 4))
                # Resolvendo deltaP  = H*deltaTeta, para achar o vetor mismatche das variaveis de estado
                for n in range(len(mismatches_p_vector)):
                    mismatches_p_vector[n] = mismatches_power_vector[n]
                solvesystem = Systemmatrixsolve(b_linha, mismatches_p_vector)
                mismatches_teta_vector = solvesystem.linearsystem()
                variable_teta_vector = variable_teta_vector + mismatches_teta_vector
                for n in range(len(variable_teta_vector)):
                    variable_vector[n] = variable_teta_vector[n]
                print("Atualização das variáveis de estado:")
                i = 0
                for n in range(len(code)):
                    if code[n] == 1 or code[n] == 2:
                        print("Dteta", n + 1, ":", "%.4f" % mismatches_teta_vector[i], "rad, teta atualizado: ",
                              "%.4f" % variable_teta_vector[i])
                        i += 1
                print("Vetor de mismatches de potência:")
                i = 0
                for n in range(len(code)):
                    if code[n] == 1 or code[n] == 2:
                        print("Barra", n + 1, " deltaP: ", mismatches_p_vector[i], "pu, P calculado: ",
                              p_calc_vector[n])
                        i += 1
                control_p = max([abs(item) for item in mismatches_p_vector])
                iteration_p += 1
            if control_q > tolerance:
                print("--------------------------------------------------")
                print("Meia Iteração Reativa: ", iteration_q)
                # Balanço de potência
                teta, v, p_calc_vector, q_calc_vector, mismatches_power_vector = \
                    self.get_power_balance(code, teta, v, variable_vector, p, q, mismatches_power_vector, 1)
                # Matriz B"
                b_2linha = self.get_b_2linha_matrix(code)
                print("Matriz B'':")
                print(np.around(b_2linha, 4))
                # Resolvendo deltaQ  = L*deltaV, para achar o vetor mismatche das variaveis de estado
                for n in range(len(mismatches_q_vector)):
                    mismatches_q_vector[n] = mismatches_power_vector[n + len(mismatches_p_vector)]
                solvesystem = Systemmatrixsolve(b_2linha, mismatches_q_vector)
                mismatches_v_vector = solvesystem.linearsystem()
                variable_v_vector = variable_v_vector + mismatches_v_vector
                for n in range(len(variable_v_vector)):
                    variable_vector[n + len(mismatches_p_vector)] = variable_v_vector[n]
                print("Atualização das variáveis de estado:")
                i = 0
                for n in range(len(code)):
                    if code[n] == 2:
                        print("DV", n + 1, ":", "%.4f" % mismatches_v_vector[i], "pu, V atualizado: ",
                              "%.4f" % variable_v_vector[i])
                        i += 1
                print("Vetor de mismatches de potência:")
                i = 0
                for n in range(len(code)):
                    if code[n] == 2:
                        print("Barra", n + 1, " deltaQ: ", mismatches_q_vector[i], "pu, Q calculado: ",
                              q_calc_vector[n])
                        i += 1
                control_q = max([abs(item) for item in mismatches_q_vector])
                iteration_q += 1
        # Resultado:
        print("--------------------------------------------------")
        print("Processo convergiu com", iteration_p - 1, " meia iteração ativa e ", iteration_q - 1,
              " meia iteração reativa, pelo método de Newton Desacoplado.")
        print("Solução: ")
        print("Variáveis de estado:")
        for k in range(len(code)):
            print("Barra", k + 1, "(", self.data_bus.iloc[k, 1], ") teta:", "%.4f" % teta[k],
                  "rad, V:", "%.4f" % v[k], "pu")
        print("Potência líquida das barras:")
        for k in range(len(code)):
            print("Barra", k + 1, "(", self.data_bus.iloc[k, 1], ") P = ", "%.4f" % p_calc_vector[k], "[pu] Q = ",
                  "%.4f" % q_calc_vector[k], "[pu]")
        print("--------------------------------------------------")
        print("\n")

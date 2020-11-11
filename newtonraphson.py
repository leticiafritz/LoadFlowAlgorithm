import math

import pandas as pd
import numpy as np

from systemmatrixsolve import Systemmatrixsolve


class Newtonraphson:
    def __init__(self, y_matrix, g_matrix, b_matrix, data_bus):
        self.y_matrix = y_matrix
        self.g_matrix = g_matrix
        self.b_matrix = b_matrix
        self.data = pd.DataFrame(data_bus)

    def get_init_data(self):
        v = np.zeros(max(self.data.iloc[:, 0]))
        teta = np.zeros(max(self.data.iloc[:, 0]))
        p = np.zeros(max(self.data.iloc[:, 0]))
        q = np.zeros(max(self.data.iloc[:, 0]))
        code = np.zeros(max(self.data.iloc[:, 0]))  # 0-> slack, 1-> PV, 2->PQ
        print("Dados de Entrada")
        for n in range(len(self.data)):
            if self.data.iloc[n, 1] == "slack":
                v[n], teta[n], code[n] = self.data.iloc[n, 2], math.radians(self.data.iloc[n, 3]), 0
                print("Barra ", n + 1, " (", self.data.iloc[n, 1], "): ",
                      "V = ", "%.4f" % v[n], " pu, teta = ", "%.4f" % teta[n], " rad, // Incógnitas: P e Q")
                teta_init = teta[n]
                v_init = v[n]
            elif self.data.iloc[n, 1] == "PV":
                v[n] = self.data.iloc[n, 2]
                p[n] = float(self.data.iloc[n, 4]) - float(self.data.iloc[n, 6])
                code[n] = 1
                print("Barra ", n + 1, " (", self.data.iloc[n, 1], "): ",
                      "V = ", "%.4f" % v[n], " pu, P = ", "%.4f" % p[n], " pu, // Incógnitas: Q e teta")
            elif self.data.iloc[n, 1] == "PQ":
                p[n] = float(self.data.iloc[n, 4]) - float(self.data.iloc[n, 6])
                q[n] = float(self.data.iloc[n, 5]) - float(self.data.iloc[n, 7])
                code[n] = 2
                print("Barra ", n + 1, " (", self.data.iloc[n, 1], "): ",
                      "P = ", "%.4f" % p[n], " pu, Q = ", "%.4f" % q[n], " pu, // Incógnitas: V e teta")

        return v, v_init, teta, teta_init, p, q, code

    def get_incognitas_vector(self, code, teta_init, v_init):
        variable_vector = []
        p_calc_vector = []
        q_calc_vector = []
        print("Valores iniciais: ")
        i = 0
        for n in range(len(code)):
            if code[n] == 1 or code[n] == 2:
                variable_vector.append(teta_init)
                p_calc_vector.append(0)
                print("teta", self.data.iloc[n, 0], ": ", "%.4f" % variable_vector[i])
                i += 1
        for n in range(len(code)):
            if code[n] == 2:
                variable_vector.append(v_init)
                q_calc_vector.append(0)
                print("V", self.data.iloc[n, 0], ": ", "%.4f" % variable_vector[i])
                i += 1
        mismatches_state_vector = np.zeros(len(variable_vector))
        mismatches_power_vector = np.zeros(len(variable_vector))

        return variable_vector, mismatches_state_vector, mismatches_power_vector

    def get_power_balance(self, code, teta, v, variable_vector, p, q, mismatches_power_vector):
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
                mismatches_power_vector[i_mismatche] = mismatches_power_p[k]
                print("Barra", k + 1, " Pesp = ", "%.4f" % p[k], "Pcal = ", p_calc_vector[k],
                      ", mismatche Pesp-Pcalc = ", "%.4f" % mismatches_power_vector[i_mismatche])
                i_mismatche += 1
        for k in range(len(code)):
            if code[k] == 2:
                mismatches_power_vector[i_mismatche] = mismatches_power_q[k]
                print("Barra", k + 1, " Qesp = ", "%.4f" % q[k], " Qesp = ", "%.4f" % q_calc_vector[k],
                      ", mismatche Qesp-Qcalc = ", "%.4f" % mismatches_power_vector[i_mismatche])
                i_mismatche += 1

        return teta, v, p_calc_vector, q_calc_vector, mismatches_power_vector

    def get_jacobiana_matrix(self, iteration, code, p_calc_vector, q_calc_vector, v, teta):
        # Estrutura da matriz Jacobiana
        npq = list(code).count(2)
        npv = list(code).count(1)
        h_matrix = np.zeros((npq + npv, npq + npv))
        n_matrix = np.zeros((npq + npv, npq))
        m_matrix = np.zeros((npq, npq + npv))
        l_matrix = np.zeros((npq, npq))
        j_matrix = np.zeros((2 * npq + npv, 2 * npq + npv))

        # Submatriz H
        row = 0
        col = 0
        for k in range(len(code)):
            for m in range(len(code)):
                if k == m:
                    if code[k] == 1 or code[k] == 2:
                        h_matrix[row, col] = -q_calc_vector[k] - (v[k] ** 2) * self.b_matrix[k, k]
                        if col < (h_matrix.shape[1] - 1):
                            col += 1
                        else:
                            col = 0
                            row += 1
                else:
                    if code[k] == 1 or code[k] == 2:
                        if code[m] == 1 or code[m] == 2:
                            h_matrix[row, col] = v[k] * v[m] * \
                                                 (self.g_matrix[k, m] * math.sin(teta[k] - teta[m]) -
                                                  self.b_matrix[k, m] * math.cos(teta[k] - teta[m]))
                            if col < (h_matrix.shape[1] - 1):
                                col += 1
                            else:
                                col = 0
                                row += 1

        # Submatriz N
        row = 0
        col = 0
        for k in range(len(code)):
            for m in range(len(code)):
                if k == m:
                    if code[k] == 2:
                        n_matrix[row, col] = (v[k] ** (-1)) * (p_calc_vector[k] +
                                                               (v[k] ** 2) * self.g_matrix[k, k])
                        if col < (n_matrix.shape[1] - 1):
                            col += 1
                        else:
                            col = 0
                            row += 1
                else:
                    if code[k] == 1 or code[k] == 2:
                        if code[m] == 2:
                            n_matrix[row, col] = v[k] * (self.g_matrix[k, m] * math.cos(teta[k] - teta[m]) +
                                                         self.b_matrix[k, m] * math.sin(teta[k] - teta[m]))
                            if col < (n_matrix.shape[1] - 1):
                                col += 1
                            else:
                                col = 0
                                row += 1

        # Submatriz M
        row = 0
        col = 0
        for k in range(len(code)):
            for m in range(len(code)):
                if k == m:
                    if code[k] == 2:
                        m_matrix[row, col] = p_calc_vector[k] - (v[k] ** 2) * self.g_matrix[k, k]
                        if col < (m_matrix.shape[1] - 1):
                            col += 1
                        else:
                            col = 0
                            row += 1
                else:
                    if code[k] == 2:
                        if code[m] == 1 or code[m] == 2:
                            m_matrix[row, col] = -v[k] * v[m] * \
                                                 (self.g_matrix[k, m] * math.cos(teta[k] - teta[m]) +
                                                  self.b_matrix[k, m] * math.sin(teta[k] - teta[m]))
                            if col < (m_matrix.shape[1] - 1):
                                col += 1
                            else:
                                col = 0
                                row += 1

        # Submatriz L
        row = 0
        col = 0
        for k in range(len(code)):
            for m in range(len(code)):
                if k == m:
                    if code[k] == 2:
                        l_matrix[row, col] = (v[k] ** (-1)) * (q_calc_vector[k] -
                                                               (v[k] ** 2) * self.b_matrix[k, m])
                        if col < (l_matrix.shape[1] - 1):
                            col += 1
                        else:
                            col = 0
                            row += 1
                else:
                    if code[k] == 2:
                        if code[m] == 2:
                            l_matrix[row, col] = v[k] * (self.g_matrix[k, m] * math.sin(teta[k] - teta[m]) -
                                                         self.b_matrix[k, m] * math.cos(teta[k] - teta[m]))
                            if col < (l_matrix.shape[1] - 1):
                                col += 1
                            else:
                                col = 0
                                row += 1

        # Montando a matriz Jacobiana
        j_matrix_aux1 = np.concatenate((h_matrix, n_matrix), axis=1)
        j_matrix_aux2 = np.concatenate((m_matrix, l_matrix), axis=1)
        j_matrix = np.concatenate((j_matrix_aux1, j_matrix_aux2), axis=0)

        print("Matriz Jacobiana da iteração", iteration, ":")
        print(np.around(j_matrix, 4))

        return j_matrix, h_matrix, n_matrix, m_matrix, l_matrix

    def get_newtonraphson_solution(self):
        print("\n")
        print("Iniciando o método de Newton Raphson ...")

        # Organizando os dados
        v, v_init, teta, teta_init, p, q, code = self.get_init_data()

        # Vetor de incognitas
        variable_vector, mismatches_state_vector, mismatches_power_vector = \
            self.get_incognitas_vector(code, teta_init, v_init)

        # Processo iterativo
        tolerance = float(input("Qual a tolerância? "))
        iteration = 0
        print("--------------------------------------------------")
        print("Iteração: ", iteration)
        # Balanço de potência
        teta, v, p_calc_vector, q_calc_vector, mismatches_power_vector = \
            self.get_power_balance(code, teta, v, variable_vector, p, q, mismatches_power_vector)
        control = max(abs(mismatches_power_vector))

        while control > tolerance:
            iteration += 1
            print("--------------------------------------------------")
            print("Iteração: ", iteration)
            # Matriz Jacobiana
            j_matrix, h_matrix, n_matrix, m_matrix, l_matrix = \
                self.get_jacobiana_matrix(iteration, code, p_calc_vector, q_calc_vector, v, teta)
            # Resolvendo deltaP  = J*deltaEstado, para achar o vetor mismatche das variaveis de estado
            solvesystem = Systemmatrixsolve(j_matrix, mismatches_power_vector)
            mismatches_state_vector = solvesystem.linearsystem()
            variable_vector = variable_vector + mismatches_state_vector
            # Balanço de potência
            teta, v, p_calc_vector, q_calc_vector, mismatches_power_vector = \
                self.get_power_balance(code, teta, v, variable_vector, p, q, mismatches_power_vector)
            print("Atualização das variáveis de estado:")
            i = 0
            for n in range(len(code)):
                if code[n] == 1 or code[n] == 2:
                    print("Dteta", n + 1, ":", "%.4f" % mismatches_state_vector[i], "rad, teta atualizado: ",
                          variable_vector[i])
                    i += 1
            for n in range(len(code)):
                if code[n] == 2:
                    print("DV", n + 1, ":", "%.4f" % mismatches_state_vector[i], "pu, V atualizado: ",
                          variable_vector[i])
                    i += 1
            print("Vetor de mismatches de potência:")
            i = 0
            for n in range(len(code)):
                if code[n] == 1 or code[n] == 2:
                    print("Barra", n + 1, " deltaP: ", mismatches_power_vector[i])
                    i += 1
            for n in range(len(code)):
                if code[n] == 2:
                    print("Barra", n + 1, " deltaQ: ", mismatches_power_vector[i])
                    i += 1
            control = max(abs(mismatches_power_vector))

        # Resultado:
        print("--------------------------------------------------")
        print("Processo convergiu com", iteration, "pelo método de Newton Raphson.")
        print("Solução: ")
        print("Variáveis de estado:")
        for k in range(len(code)):
            print("Barra", k + 1, "(", self.data.iloc[k, 1], ") teta:", "%.4f" % teta[k],
                  "rad, V:", "%.4f" % v[k], "pu")
        print("Potência líquida das barras:")
        for k in range(len(code)):
            print("Barra", k + 1, "(", self.data.iloc[k, 1], ") P = ", "%.4f" % p_calc_vector[k], "[pu] Q = ",
                  "%.4f" % q_calc_vector[k], "[pu]")
        print("--------------------------------------------------")
        print("\n")

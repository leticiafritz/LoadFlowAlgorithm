import math

import pandas as pd
import numpy as np

from systemmatrixsolve import Systemmatrixsolve


class Newtondesacoplado:
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
                print("teta", self.data.iloc[n, 0], ": ", "%.4f" % variable_teta_vector[i])
                mismatches_teta_vector.append(teta_init)
                mismatches_p_vector.append(0)
                i += 1
        i = 0
        for n in range(len(code)):
            if code[n] == 2:
                variable_v_vector.append(v_init)
                q_calc_vector.append(0)
                print("V", self.data.iloc[n, 0], ": ", "%.4f" % variable_v_vector[i])
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

    def get_jacobiana_matrix(self, code, p_calc_vector, q_calc_vector, v, teta):
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

        return j_matrix, h_matrix, n_matrix, m_matrix, l_matrix

    def get_newtondesacoplado_solution(self):
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
                # Matriz Jacobiana
                j_matrix, h_matrix, n_matrix, m_matrix, l_matrix = \
                    self.get_jacobiana_matrix(code, p_calc_vector, q_calc_vector, v, teta)
                print("Submatriz H da Jacobiana, iteração", iteration_p, ":")
                print(np.around(h_matrix, 4))
                # Resolvendo deltaP  = H*deltaTeta, para achar o vetor mismatche das variaveis de estado
                for n in range(len(mismatches_p_vector)):
                    mismatches_p_vector[n] = mismatches_power_vector[n]
                solvesystem = Systemmatrixsolve(h_matrix, mismatches_p_vector)
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
                # Matriz Jacobiana
                j_matrix, h_matrix, n_matrix, m_matrix, l_matrix = \
                    self.get_jacobiana_matrix(code, p_calc_vector, q_calc_vector, v, teta)
                print("Submatriz L da Jacobiana, iteração", iteration_q, ":")
                print(np.around(l_matrix, 4))
                # Resolvendo deltaQ  = L*deltaV, para achar o vetor mismatche das variaveis de estado
                for n in range(len(mismatches_q_vector)):
                    mismatches_q_vector[n] = mismatches_power_vector[n + len(mismatches_p_vector)]
                solvesystem = Systemmatrixsolve(l_matrix, mismatches_q_vector)
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
            print("Barra", k + 1, "(", self.data.iloc[k, 1], ") teta:", "%.4f" % teta[k],
                  "rad, V:", "%.4f" % v[k], "pu")
        print("Potência líquida das barras:")
        for k in range(len(code)):
            print("Barra", k + 1, "(", self.data.iloc[k, 1], ") P = ", "%.4f" % p_calc_vector[k], "[pu] Q = ",
                  "%.4f" % q_calc_vector[k], "[pu]")
        print("--------------------------------------------------")
        print("\n")

############################################################
# PROJETO: SOLUÇÃO DO FLUXO DE CARGA
# NOME: LETÍCIA FRITZ HENRIQUE
# E-MAIL: LETICIA.HENRIQUE@UFJF.BR
# VERSÃO: 1.0
############################################################

# BIBLIOTECAS
import os
import sys

import pandas as pd
import numpy as np

from ybus import Ybus
from solution import Solution


# FUNÇÃO PRINCIPAL
def main():
    # Dados de Entrada
    print("ALGORITMO PARA SOLUÇÃO DO FLUXO DE CARGA")
    print("-----")
    line_data = pd.read_excel(os.path.dirname(sys.argv[0]) + '/grid_data.xlsx',
                              index_col=None, sheet_name='line')
    bus_data = pd.read_excel(os.path.dirname(sys.argv[0]) + '/grid_data.xlsx',
                             index_col=None, sheet_name='bus')

    # Calculo da Y barra
    ybus_matrix = Ybus(line_data, bus_data)
    y_matrix, g_matrix, b_matrix = ybus_matrix.get_ybus_matrix()

    # Método de resolução do Fluxo de Carga
    solve = Solution(y_matrix, g_matrix, b_matrix, bus_data, line_data)
    method = solve.get_solution_method()
    solve.get_solution(method)

    # Análise do fluxo de carga
    option = solve.get_analise_option()
    solve.get_analise(option)


# RUN
if __name__ == '__main__':
    main()

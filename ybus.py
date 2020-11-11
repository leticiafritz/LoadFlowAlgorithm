import pandas as pd
import numpy as np
import math


class Ybus:
    def __init__(self, data_line, data_bus):
        self.data = pd.DataFrame(data_line)
        self.b_sh_bus = list(pd.DataFrame(data_bus).iloc[:, 8])

    def get_ybus_matrix(self):
        y_bus = np.zeros((max(self.data.iloc[:, 1]), max(self.data.iloc[:, 1])), dtype=complex)
        # Preenche a matriz admitância fora da diagonal principal
        for k in range(max(self.data.iloc[:, 1])):
            row = self.data.iloc[k, 0] - 1
            col = self.data.iloc[k, 1] - 1
            y_bus[row, col] = -(self.data.iloc[k, 5]) * (math.e ** (0-self.data.iloc[k, 6]*1j)) * \
                              (1 / (complex(self.data.iloc[k, 2], self.data.iloc[k, 3])))
            y_bus[col, row] = -(self.data.iloc[k, 5]) * (math.e ** (0+self.data.iloc[k, 6]*1j)) * \
                              (1 / (complex(self.data.iloc[k, 2], self.data.iloc[k, 3])))

        # Preenche a matriz admitância dentro da diagonal principal
        for k in range(max(self.data.iloc[:, 1])):
            y_sum = 0
            for circuit in range(len(self.data.iloc[:, 1])):
                if (self.data.iloc[circuit, 0] == k + 1) or (self.data.iloc[circuit, 1] == k + 1):
                    y_sum = y_sum + (((self.data.iloc[k, 5] ** 2) *
                                      (1 / complex(self.data.iloc[circuit, 2], self.data.iloc[circuit, 3]))) +
                                     (complex(0, self.data.iloc[circuit, 4]) / 2))
            y_bus[k, k] = y_sum + complex(0, float(self.b_sh_bus[k]))

        # Encontrando matriz G e B
        g_bus = y_bus.real
        b_bus = y_bus.imag

        # Output
        print("Matriz Admitância Nodal:")
        print(np.around(y_bus, 4))
        print("\n")
        print("Matriz G:")
        print(np.around(g_bus, 4))
        print("\n")
        print("Matriz B:")
        print(np.around(b_bus, 4))
        print("\n")

        return y_bus, g_bus, b_bus

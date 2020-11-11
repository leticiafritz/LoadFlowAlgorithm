import sys

import numpy as np
from scipy.linalg import solve
np.set_printoptions(suppress=True)


class Systemmatrixsolve:
    def __init__(self, a_matrix, b_vector):
        self.a_matrix = a_matrix
        self.b_vector = b_vector

    def linearsystem(self):

        return solve(self.a_matrix, self.b_vector)



import numpy as np

def rabi_model(x, tDec, T, A, B, C):
    """
    __ Description __
    Fits Rabi oscillations of the format
    A e^(-x/tDec) cos(2π x / T + B) + C
    """

    return A * np.sin(2 * np.pi * x / T + B) * np.exp(-x / tDec) + C

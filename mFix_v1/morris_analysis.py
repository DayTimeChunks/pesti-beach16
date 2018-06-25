# -*- coding: utf-8 -*-
from SALib.analyze import morris as moa
import numpy as np


def getInputVector(row, sample_matrix):
    """
    :param row: relevant sample row
    :param sample_matrix: numpy sample matrix
    :return: a numpy row with required input parameters
    """
    test_vector = sample_matrix[row]
    return test_vector


def getTestOutput(variable, runs):
    """
    :param variable: name of the variable to extract
    :param runs: number of runs
    :return: the output vector, length = no. of runs
    """
    Y = np.zeros([runs])
    folder = 1
    for i in range(runs):
        path = str(folder) + "/" + variable + ".tss"
        # res = pd.read_table(path, header=None)
        res = np.loadtxt(path, float, skiprows=4, usecols=[1])
        Y[i] = res[-1]
        folder += 1
    return Y


def getSiList(ouput_vars, problem, param_values, grid_jump, p, runs):
    morris_results = []
    for indx in range(len(ouput_vars)):
        var_name = ouput_vars[indx]
        Y = getTestOutput(var_name, runs)
        Si = moa.analyze(problem, param_values, Y, num_resamples=1000, print_to_console=True,
                         grid_jump=grid_jump, num_levels=p)
        morris_results.append(Si)
    return morris_results


def saveInputMatrix(param_values):
    np.savetxt("input_vectors.txt", param_values)


def getInputMatrix(path):
    path += "input_vectors.txt"
    return np.loadtxt(path, float)









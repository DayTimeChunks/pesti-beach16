from SALib.sample import latin
from SALib.test_functions import Ishigami
from SALib.analyze import delta
import numpy as np


# problem = {
#         'num_vars': 3,
#         'names': ['z3_factor', 'cZ0Z1', 'cZ'],
#         'bounds': [[0.75, 0.99], [0.01, 1.0], [0.01, 1.0]]
#     }

problem = {
    'num_vars': 3,
    'names': ['x1', 'x2', 'x3'],
    'bounds': [[-3.14159265359, 3.14159265359],
               [-3.14159265359, 3.14159265359],
               [-3.14159265359, 3.14159265359]]
}

X = latin.sample(problem, 100)
Y = Ishigami.evaluate(X)
Si = delta.analyze(problem, X, Y, print_to_console=True)

# print(latin.sample(problem, 6))


# Parameter delta delta_conf S1 S1_conf
# x1 0.169099 0.021929 0.274358 0.030497
# x2 0.273540 0.027849 0.324121 0.059738
# x3 0.149089 0.014088 0.009655 0.018261
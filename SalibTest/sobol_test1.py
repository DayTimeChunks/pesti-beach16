from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.test_functions import Ishigami
import numpy as np


# Source
# http://salib.readthedocs.io/en/latest/basics.html

# Define the model inputs
problem = {
    'num_vars': 3,
    'names': ['x1', 'x2', 'x3'],
    'bounds': [[-3.14159265359, 3.14159265359],
               [-3.14159265359, 3.14159265359],
               [-3.14159265359, 3.14159265359]]
}

# Generate samples
# Store in a 8000-samples by 3-variables numpy matrix
param_values = saltelli.sample(problem, 1000)
# print(param_values.shape) # -> (8000, 3)
# print(param_values.shape[0]) # -> 8000

# TODO: Is saltelli.sample() also used in Morris?
# ...to generate the matrix of test cases?

# Run model (example)
Y = Ishigami.evaluate(param_values)

# Define a function that extracts the output of the model,
# based on a vector of input parameters.
def evaluate_model(input_vector):
    pass

# We need to provide our output-Y based on our specific model.
my_own_model = False
if my_own_model:
    # Y should also be a numpy matrix (n-samples tested by 2-column)
    # So initiate the numpy matrix with zeroes.
    Y = np.zeros([param_values.shape[0]])

    # Store the value of each model simulation:
    for i, X in enumerate(param_values):
        Y[i] = evaluate_model(X)



# Perform analysis
Si = sobol.analyze(problem, Y, print_to_console=True)

# Print the first-order sensitivity indices
print Si['S1']
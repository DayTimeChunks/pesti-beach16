from SALib.sample import latin


# TODO: Test if this works!

problem = {
        'num_vars': 3,
        'names': ['z3_factor', 'cZ0Z1', 'cZ'],
        'bounds': [[0.75, 0.99], [0.01, 1.0], [0.01, 1.0]]
    }

print(latin.sample(problem, 4))


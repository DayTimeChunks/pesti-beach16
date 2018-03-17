
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd

# Which computer's directory?
PC = True
# What model version?
version = 'v4'

cols = pd.read_table('colOutlet.txt', sep=' ', header=None)

print('1739: ', cols[1737][0])
print('1739: ', cols[1738][0])
print('1739: ', cols[1739][0])


print(len(cols))
zeros = 0
ones = 0
col_nums = 0
for col in cols:
    col_nums += 1
    if cols[col][0] == 0:
        zeros += 1
    elif cols[col][0] == 1:
        ones += 1
        outlet = col
        print("outlet col: ", outlet)

# print("values: ", zeros, ones, outlet)


test = 'test'
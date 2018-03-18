
import csv
import json
path = 'initial.csv'
with open(path, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    data_list = list()
    my_dict = {}
    for row in reader:
        my_dict[row[0].strip()] = float(row[1])
        # data_list.append(row)
# data = [dict(zip(data_list[0], row)) for row in data_list]
# data.pop(0)
# s = json.dumps(data)
print (my_dict)

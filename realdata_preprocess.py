from tqdm import tqdm, trange
import random


data = []

with open('flights.csv') as f:
    for line in tqdm(f):
        line = line.strip().split(',')
        if line[3].strip() and line[9].strip() and line[11].strip() and line[14].strip() and line[17].strip() and line[22].strip():
            data.append([line[3].strip(), line[9].strip(), line[11].strip(), line[14].strip(), line[17].strip(), line[22].strip()])

for i in trange(1, len(data)):
    if int(data[i][0]) <= 5:
        data[i][0] = 1
    else:
        data[i][0] = 0
    
    if int(data[i][1]) >= 700 and int(data[i][1]) <= 1800:
        data[i][1] = 1
    else:
        data[i][1] = 0
        
    if int(data[i][2]) <= 0:
        data[i][2] = 1
    else:
        data[i][2] = 0
    
    data[i][3] = round((int(data[i][3]) - 18) / 700, 4)
    data[i][4] = round((int(data[i][4]) - 31) / (4983 - 31), 4)
    
    tmp = int(data[i][5])
    if tmp <= 0:
        data[i][5] = 1
    elif tmp <= 5:
        data[i][5] = 2
    elif tmp <= 15:
        data[i][5] = 3
    else:
        data[i][5] = 4

tmp = data[1:]
random.shuffle(tmp)
import csv
with open('flights_clean.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerows([data[0]] + tmp)

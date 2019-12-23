import numpy as np
import pylab as plt
import xlrd
import pandas as pd
import xlsxwriter
from sklearn import preprocessing

def wgn(x, snr):
    snr = pow(10,(snr/10.0))
    x_p = np.sum(pow(x,2)/len(x))
    n_p = x_p/snr
    return np.random.randn(len(x)) * np.sqrt(n_p)

def excel_to_matrix(path):
    table = xlrd.open_workbook(path).sheets()[0]
    row = table.nrows
    col = table.ncols
    datamatrix = np.zeros((row, col))
    for x in range(col):
        cols = np.matrix(table.col_values(x))  
        datamatrix[:, x] = cols #

    min_max_scaler = preprocessing.MinMaxScaler()
    datamatrix  = min_max_scaler.fit_transform(datamatrix)
    return datamatrix


d = excel_to_matrix("THP1_data.xlsx")


array_d = []

for i in range(len(d)):
    for j in range(len(d[0])):
        array_d.append(d[i,j])

v_d = np.array(array_d)

#print(v.T)


d_noise = wgn(v_d.T,6)
d_with_noise = v_d.T + d_noise

noise_data = np.zeros((len(d),len(d[0])))
idx = 0

for k in range(len(d)):
    for l in range(len(d[0])):
        noise_data[k][l] = d_with_noise[idx]
        idx += 1
        

print(noise_data)
print(d_with_noise[-2])
print(len(noise_data[0]))


workbook = xlsxwriter.Workbook('THP1_noise.xlsx')
worksheet = workbook.add_worksheet()



for row in range(len(noise_data)):
    for col in range(len(noise_data[0])):
        worksheet.write_number(row, col, noise_data[row][col])

workbook.close()







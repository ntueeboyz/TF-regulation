#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 21:02:39 2019

@author: Kuo
"""

import numpy as np
import xlsxwriter
from wgn import wgn
from excel_to_matrix.py import excel_to_matrix


def add_noise(data,snr):
    
    # Input array
    inputdata = data + '.xlsx'
    
    d = excel_to_matrix(inputdata)
    array_d = []

    for i in range(len(d)):
        for j in range(len(d[0])):
            array_d.append(d[i,j])

    v_d = np.array(array_d)


    # Add noise with SNR
    d_noise = wgn(v_d.T,snr)
    d_with_noise = v_d.T + d_noise

    noise_data = np.zeros((len(d),len(d[0])))
    idx = 0

    for k in range(len(d)):
        for l in range(len(d[0])):
            noise_data[k][l] = d_with_noise[idx]
            idx += 1
    
    # Output
    outputdata = data + '_noise_snr' + str(snr) + '.xlsx'

    workbook = xlsxwriter.Workbook(outputdata)
    worksheet = workbook.add_worksheet()



    for row in range(len(noise_data)):
        for col in range(len(noise_data[0])):
            worksheet.write_number(row, col, noise_data[row][col])

    workbook.close()
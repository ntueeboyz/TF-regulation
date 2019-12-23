#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 21:16:37 2019

@author: Kuo
"""
import xlrd
import numpy as np
from sklearn import preprocessing

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
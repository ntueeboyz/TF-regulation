#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 21:14:43 2019

@author: Kuo
"""

import numpy as np

def wgn(x, snr):
    snr = pow(10,(snr/10.0))
    x_p = np.sum(pow(x,2)/len(x))
    n_p = x_p/snr
    return np.random.randn(len(x)) * np.sqrt(n_p)
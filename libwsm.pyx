# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 19:21:39 2020

@author: gionuno
"""

import cython;
from libcpp.vector cimport vector;
import numpy as np;

cdef extern from "libwsm.hpp":
    cdef vector[vector[double]] wsm_cpp(vector[vector[double]] &,vector[double] &,int);

def weighted_sort_means(X,w,K):
    return np.array(wsm_cpp(X,w,K));
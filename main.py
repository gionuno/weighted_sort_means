# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 19:30:30 2020

@author: gionuno
"""

import numpy as np;
import numpy.random as rd;

from libwsm import weighted_sort_means;

import matplotlib.pyplot as plt;
import matplotlib.image  as img;

im = img.imread("zebras.jpg");

H = np.zeros((256,256,256));
for x in range(im.shape[0]):
    for y in range(im.shape[1]):
        H[im[x,y,0],im[x,y,1],im[x,y,2]] += 1.0;

X = [];
w = [];
for i in range(256):
    for j in range(256):
        for k in range(256):
            if H[i,j,k] > 0:
                X.append(np.array([1.0*i,1.0*j,1.0*k]));
                w.append(H[i,j,k]);

K = 16;
C = weighted_sort_means(np.array(X),np.array(w),K);

jm = np.zeros(im.shape);
for x in range(im.shape[0]):
    for y in range(im.shape[1]):
        k = 0;
        min_dist = 1e15;
        for l in range(C.shape[0]):
            aux_dist = np.linalg.norm(im[x,y,:]-C[l,:]);
            if min_dist > aux_dist:
                min_dist = aux_dist;
                k = l;
        jm[x,y,:] = C[k];


fig, ax = plt.subplots(2,2)

ax[0,0].imshow(im/255.0);
ax[0,0].set_axis_off();
ax[0,1].imshow(jm/255.0);
ax[0,1].set_axis_off();
ax[1,0].imshow(np.abs(im/255.0-jm/255.0));
ax[1,0].set_axis_off();
ax[1,1].imshow(C.reshape((4,4,3))/255.0);
ax[1,1].set_axis_off();

plt.show()
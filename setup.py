# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 19:26:59 2020

@author: gionuno
"""

from distutils.core import setup;
from distutils.extension import Extension;
from Cython.Build import cythonize;

e_modules = cythonize([Extension("libwsm",["libwsm.pyx"],extra_compile_args=['-std=c++11'],libraries=[],language="c++")]);
setup(name="libwsm",ext_modules = e_modules,script_args=['build'])
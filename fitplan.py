"""
Code to calculate the ODR fit of the data
"""
import numpy as np
nmax=20
nmax_one=400
nmax_two=2000
dim_zero=3
dim_one=4
dim=20
dim_two=21
dim_three=420
scale=22.0
x = np.ndarray(shape=(dim_one, nmax_one), dtype=float, order='F')
mat = np.ndarray(shape=(dim, dim_two), dtype=float, order='F') # Values of dim by dim_two
vec = np.ndarray(shape=dim, dtype=float, order='F') # Size of dim
#  coef(DIM1,DIM1)
coef = np.ndarray(shape=(dim_one, dim_one), dtype=float, order='F')
# covmat(DIM0,DIM,DIM2), scoef(DIM1,DIM1)
covmat = np.ndarray(shape=(dim_zero, dim, dim_two), dtype=float, order='F')
scoef = np.ndarray(shape=(dim_one, dim_one), dtype=float, order='F')
# var(DIM1)
var = np.ndarray(shape=dim_one, dtype=float, order='F')
# work(DIM3), res(NMAX1)
work = np.ndarray(shape=(dim_three), dtype=float, order='F')
res = np.ndarray(shape=(nmax_one), dtype=float, order='F')
# real cond
cond = 0.0

# integer i, errcode, norder, info
i = 0
errcode = 0
norder = 0
info = 0
# integer nclus, nlda, ntotal, nfit
nclus = 0
nlda = 0
ntotal = 0
nfit = 0
# logical vflag
vlfag = False

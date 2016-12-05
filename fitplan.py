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
x = [][] # Values of dim_one by nmax_one
mat = [][] # Values of dim by dim_two
vec = [] # Size of dim
#  coef(DIM1,DIM1)
coef = [][]
# covmat(DIM0,DIM,DIM2), scoef(DIM1,DIM1)
covmat = [][][]
scoef = [][]
# var(DIM1)
var = []
# work(DIM3), res(NMAX1)
work = []
res = []
# real cond
cond = 0.0

# integer i, errcode, norder, info
# integer nclus, nlda, ntotal, nfit

# logical vflag

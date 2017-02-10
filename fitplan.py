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


########################### fitplan.f converted code #########################
"""c Get the dimension  and flag for output. Read data from STDIN
      read(*,*) norder, vflag
      if(norder .eq. 2 .or. norder .eq. 3) then
        do 30 i=1,NMAX1,1
           read(*,*,err=940,end=40) (x(j,i), j=1,norder)
c Find the column with expected mu in it and transform
      if(norder .eq. 3 .and. i .eq. 1) then
        if( x(1,i) .gt. x(2,i) .and. x(1,i) .gt. x(3,i) ) then
           ncol=1
        else if( x(2,i) .gt. x(1,i) .and. x(2,i) .gt. x(3,i) ) then
           ncol=2
        else if( x(3,i) .gt. x(1,i) .and. x(3,i) .gt. x(2,i) ) then
           ncol=3
        endif
      endif
      if(norder .eq. 3) then
        x(ncol,i) = x(ncol,i)-SCALE
      endif
  30    continue
      else
        write(*,*) "Wrong dimension used, must be 2 or 3"
        stop
      endif

  40    continue

      if(vflag) then
        write(*,*) "==> End of reading <=="
      endif
      write(*,*) "ncol=",ncol

      ntotal = i-1
"""
# Get dimensions and flags for output


# Initialize x(norder+1,j) = 1
for j in range(1, ntotal):
    x[norder+1, j] = 1.0

# initialize coefficients
for i in range(1, norder):
    for j in range(1, norder):
        coef[i,j] = 0.0
    coef[i, norder+1] = -1.0

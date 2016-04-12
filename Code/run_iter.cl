# 2016apr12 Various fits showing how to use iterfitlin
# Real runs would have nboot=100
#
# Requirements and how to run this: 

# Compile random.f using the make file. Put the executable random in your bin directory
# Edit random.cl to point to that executable instead of to my bin directory
# Start iraf
# Then in the directory with all the files, declare the two cl-scripts as tasks
# cl>    task iterfitlin=iterfitlin.cl
# cl>task random=random.cl

# Declare this driver script which as no parameters that is why the declaration looks different
# cl> task $run_iter=run_iter.cl
# Then run it from the IRAF prompt
# cl> run_iter


# The data files:
# comafit.fits  Coma cluster sample
# rxj0152allfit.fits  RXJ0152 sample
# rxj1226allfit.fits  RXJ1226 sample

delete input ver-
printf("comafit.fits\n", > "input")
iterfitlin input y_col=lreJB_kpc x1_col=lsigma_cor x2_col=lIeJB_cor nboot=10 seed=14

delete input ver-
printf("comafit.fits\n", > "input")
iterfitlin input y_col=lML_JB x1_col=lMass x2_col="" nboot=10 seed=14

delete input ver-
printf("rxj0152allfit.fits\nrxj1226allfit.fits\n", > "input")
iterfitlin input y_col=lreJB_kpc_DEV x1_col=lsigma_cor x2_col=lIeJB_DEV nboot=10 seed=14

delete input ver-
printf("rxj0152allfit.fits\nrxj1226allfit.fits\n", > "input")
iterfitlin input y_col=lML_JB_DEV x1_col=lMass_DEV x2_col="" nboot=10 seed=14


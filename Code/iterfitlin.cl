# Fit plane or line iteratively 
#   if more than one cluster, each cluster in separate STSDAS table
#
# Residuals are calculated according to resalgo
#      per:  perpendicular to the fitted relation
#            delta = (y - a*x1 - b*x2 - c)/sqrt(1+a^2+b^2)
#      y  :  in y     delta = (y - a*x1 - b*x2 - c)
#      x1 :  in x1
#      x2 :  in x2
#
# Minization algorithms:  (minalgo)
#    quartile : distance between upper (0.75) quartile point 
#               and lower (0.25) quartile point
#    delta100 : SUM abs(res)
#    delta60  : SUM abs(res) of the central 60 percent of distribution
#    rms100   : SUM res^2
#    rms60    : SUM res^2 of the central 60 percent of distribution
#
# Zero points : median (use with delta and quartile)
#               mean   (use with rms)
# 
# Bootstrap uncertainties on coefficients if nboot>0
# Seed for random is in each step   seed-(nboot-step)
#
# WARNING: This is a rather slow algorithm. 
# FP for 200 galaxies in 10 clusters takes 5min on Sparc 10 
# if restart=yes and nboot=0 (no uncertainties)
# 
# Version: 21.06.95
#          16.10.95  Solaris roeskva
#          14.02.96  changed resalgo=y initial values, faster conv.
#          10.11.97  use rename instead of trename in 2.11
#          21.01.98  bootstrap changed to the proper way of doing it
#          06.03.98  fixed bug in bootstrap for 2 parameter fits
#                    improved random number generation
#          27.07.98  stop if itera>3*maxiter || niterb>3*maxiter
#                    no matter the other parameter's niter
#          30.03.04  handle INDEF automatically, works under IRAF 2.12.2 Redhat9.0
# Inger Jorgensen, Gemini Observatory
# e-mail: inger@gemini.edu

procedure iterfitlin(tablist)

char tablist {prompt="List of input STSDAS tables"}
char minalgo {"delta100",min="quartile|delta100|delta60|rms100|rms60",
 prompt="Distance to minimize (delta100,delta60,rms100,rms60,quartile)"}
char zeropoint {"median",min="mean|median",
 prompt="Zero points (median,mean)"}
char resalgo {"per",min="per|y|x1|x2",
 prompt="Residual to minimize (per,y,x1,x2)"}
char y_col  {"lre_GR_sc",prompt="Column name for y"}
char x1_col {"lsig_re",prompt="Column name for x1"}
char x2_col  {"lIe_GR_sccor",prompt="Column name for x2 (optional)"}
char gal_col {"galaxy",prompt="Column name for galaxy"}
char grp_col {"group",prompt="Column name for group"}
real facta   {0.05,prompt="Starting factor for changes in a"}
real factb   {0.02,prompt="Starting factor for changes in b"}
int  maxiter {20,prompt="Maximum number of iterations"}
bool restart {yes,prompt="Restart iteration with smaller factors"}
int  nboot   {0,prompt="Number of estimates for bootstrap"}
int  seed    {1,prompt="Seed for random used in bootstrap"}
struct *scanfile {"",prompt="List directed struct"}

begin

char n_tablist, n_recol, n_Iecol, n_sigcol, n_galcol, n_grpcol
char n_minalgo, n_zeropoint, n_resalgo
real n_factain, n_factbin
int  n_maxiter, n_nboot, n_seed
bool n_restart

char tmpall, tmpsel, tmpout, tmpexp, n_table, n_expression, n_allcol
char tmpran, n_taball, tmpboo
int  n_clus, n_totgal, n_i, n_member, n_rich, n_totboot
real n_a, n_b, n_zero, n_norm, n_signa, n_signb, n_facta, n_factb
real n_ma, n_mb, n_factas, n_factbs
real n_rms, n_lrerms
real n_nlow, n_nhigh, n_vlow, n_vhigh, n_delta, n_tabval
real n_vlowout, n_vhighout, n_deltaout, n_aout, n_bout
real n_ea, n_eb, n_ssa, n_sa, n_ssb, n_sb
int  n_itera, n_iterb
bool n_fla, n_flb, n_flend, n_flfirst, n_twopar, n_flprint, n_flboot

n_tablist=tablist ; n_recol=y_col ; n_Iecol=x2_col ; n_sigcol=x1_col
n_galcol=gal_col ; n_grpcol=grp_col
n_allcol=n_galcol//","//n_grpcol//","//n_recol//","//n_sigcol//","//n_Iecol
n_factain=facta ; n_factbin=factb
n_factas=n_factain ; n_factbs=n_factbin   # save values
n_minalgo=minalgo ; n_maxiter=maxiter ; n_zeropoint=zeropoint
n_restart=restart ; n_resalgo=resalgo ; n_nboot=nboot ; n_seed=seed
# make sure seed is positive and larger than nboot
if(n_seed<=0 || n_seed<=n_nboot)
   n_seed=max(n_seed*n_seed,(n_seed+1)*(n_nboot+1))

n_flprint=yes  # printing on
if(n_resalgo=="x2" && (n_Iecol=="" || n_Iecol==" ")) {
  print("Coordinate chosen for minimization is not present. Stopped")
  bye
}
tmpall = mktemp("tmpall")
tmpsel = mktemp("tmpsel")
tmpout = mktemp("tmpout")
tmpexp = mktemp("tmpexp")

# check if two parameter fit
n_twopar=no
if(n_Iecol=="" || n_Iecol==" ") {
  n_twopar=yes ; n_b=0 ; n_factbin=0
}

cache("tinfo","tstat")

# put everything in one table with nclus as the cluster number
scanfile=n_tablist
n_clus=0
while(fscan(scanfile,n_table)!=EOF) {
 n_clus+=1
 tproject(n_table,tmpsel,columns=n_allcol)
 tcalc(tmpsel,"nclus",n_clus,colfmt="i3")
 if(!access(tmpall//".tab") ) 
   rename(tmpsel//".tab",tmpall//".tab")
 else {
   tmerge(tmpall//","//tmpsel,tmpout,"append")
   tdelete(tmpall//","//tmpsel,verify=no)
   rename(tmpout//".tab",tmpall//".tab")
 }
}
tinfo(tmpall,ttout=no)
n_totgal = tinfo.nrows
n_taball=tmpall          # always as used for deletion

# Initialize bootstrap parameters, 
#    set n_flboot=no until bootstrap starts
n_flboot=no
if(n_nboot>0) {    
  tmpran = mktemp("tmpran")
  tmpboo = mktemp("tmpboo")
  n_totboot=n_nboot    
  n_ssa=0 ; n_sa=0 ; n_ssb=0 ; n_sb=0
  tcalc(tmpall,"row","ROWNUM",colfmt="i4")
}

print("")
printf("Fitting technique : iterative, %s %s minimized, %s zero points\n",
  n_resalgo,n_minalgo,n_zeropoint)
printf("Number of clusters: %4d\n",n_clus)
printf("Number of galaxies: %4d\n",n_totgal)
printf(" (n_facta,n_bfact)=(%7.4f,%7.4f)\n",n_factain,n_factbin)
print("Columns           : "//n_allcol)
print("")

scanfile=""

bootstrap:   # Jump point for bootstrap

# fit the cluster with most data to get an idea of where to 
# iterate from
for(n_i=1 ; n_i<=n_clus ; n_i+=1) {
 tstat(tmpall,"nclus",lowlim=n_i,highlim=n_i, >> "/dev/null")
 if(n_i==1) {
    n_member=tstat.nrows
    n_rich=n_i
 }
 else  {
    if(tstat.nrows>n_member) {
      n_member=tstat.nrows
      n_rich=n_i
    }
 }
}
if(n_flprint) {
printf("Cluster number with most data         :  %3d\n",n_rich)
printf("Number of galaxies in this cluster    :  %3d\n",n_member)
 }

# fitting this cluster - get rid of any INDEF points
if(n_twopar)
  tselect(tmpall,tmpsel,"nclus=="//n_rich//" && "//n_recol//"<99999. && "//n_sigcol//"<99999.")
else
  tselect(tmpall,tmpsel,"nclus=="//n_rich//" && "//n_recol//"<99999. && "//n_sigcol//"<99999. && "//n_Iecol//"<99999.")
tinfo(tmpsel,ttout-)
if(n_flprint)
  printf("Number of galaxies fit in this cluster:  %3d\n",tinfo.nrows)
tfitlin(tmpsel,n_recol,n_sigcol,n_Iecol,rows="-",verbose=no, > tmpout)
# get the RMA coefficients
if(n_resalgo=="y") {
head(tmpout,nlines=6) | fields("STDIN","3-4",lines="6") | \
scan(n_a,n_b) 
} else {
tail(tmpout,nlines=3) | fields("STDIN","2-3",lines="1") | \
scan(n_a,n_b)
}
if(n_twopar)
  n_b=0.0
if(n_flprint) {
printf("Initial values               (a,b)=(%7.4f,%7.4f)\n",n_a,n_b)
print("")
}
delete(tmpout,verify=no)
tdelete(tmpsel,verify=no)

# if two parameter fit - make face zero column
if(n_twopar) {
  tcalc(tmpall,"tmp","0.",colfmt="f2.0")
  n_Iecol="tmp"
}

restart:    # Jump point for restarting with smaller steps

# Iteration to minimize some measurement of the residuals
# residuals defined according to resalgo

if(n_flprint) {
if(n_minalgo=="quartile") {
printf(" ----a----   ----b----\n"
printf("  i  fact     i  fact     a       b       low q.    high q.    delta\n")
}
if(n_minalgo=="delta100") {
printf(" ----a----   ----b----\n"
printf("  i  fact     i  fact     a       b       delta       rms\n")
}
if(n_minalgo=="delta60") {
printf(" ----a----   ----b----\n"
printf("  i  fact     i  fact     a       b       delta       rms     N(delta)\n")
}
if(n_minalgo=="rms100") {
printf(" ----a----   ----b----\n"
printf("  i  fact     i  fact     a       b         rms\n")
}
if(n_minalgo=="rms60") {
printf(" ----a----   ----b----\n"
printf("  i  fact     i  fact     a       b         rms      N(rms)\n")
}
}  # end of n_flprint

#bootstrap:   # Jump point for bootstrap

n_flend=no
n_itera=1 ; n_iterb=1 
n_fla=yes ; n_flb=no ; n_facta=n_factain; n_factb=n_factbin
n_signa=1. ; n_signb=1 ; n_flfirst=yes
nextres: 

# derive the zero points and calculate the residuals
# the zero points are median or mean as defined by n_zeropoint
n_norm=1.                         # min in y
if(n_resalgo=="per")
  n_norm=sqrt(1.+n_a**2+n_b**2)   # min perpendicular
if(n_resalgo=="x1")
  n_norm=-1.*n_a                  # min in x1
if(n_resalgo=="x2")
  n_norm=-1.*n_b                  # min in x2

tcalc(tmpall,"res","0.",colfmt="f6.3")

###################### ZEROPOINT TRUMP ############################################
for(n_i=1 ; n_i<=n_clus ; n_i+=1) {
# delta y
 n_expression = "("//n_recol//"-"//n_a//"*"//n_sigcol//"-"//n_b//"*"//n_Iecol//")"
 print(n_expression, > tmpexp)
 tcalc(tmpall,"z"//n_i,"@"//tmpexp,colfmt="f6.3")
 delete(tmpexp,verify=no)
 n_expression="z"//n_i//"*(nclus=="//n_i//")+1000.*(nclus!="//n_i//")"
 print(n_expression, >tmpexp)
 tcalc(tmpall,"z"//n_i,"@"//tmpexp,colfmt="f6.3")
 delete(tmpexp,verify=no)
 tstat(tmpall,"z"//n_i,lowlim=INDEF,highlim=100., >> "/dev/null")
 if(n_zeropoint=="median")
   n_zero=tstat.median
 else
   n_zero=tstat.mean
# printf("Zero point for cluster  %-3d : %8.5f\n",n_i,n_zero)
# residuals normalized
 n_expression="((z"//n_i//"-"//n_zero//")*(nclus=="//n_i//"))/"//n_norm//"+1000.*(nclus!="//n_i//")"
 print(n_expression, > tmpexp)
 tcalc(tmpall,"r"//n_i,"@"//tmpexp,colfmt="f6.3")
 delete(tmpexp,verify=no)
 n_expression="res+((z"//n_i//"-"//n_zero//")*(nclus=="//n_i//"))/"//n_norm
 print(n_expression, > tmpexp)
 tcalc(tmpall,"res","@"//tmpexp,colfmt="f6.3"
 delete(tmpexp,verify=no)
}

# ----------- quartile minimized -----------------
if(n_minalgo=="quartile") {
# get lower and upper quatile point (mean of two closest points)
tsort(tmpall,"res")
n_nlow = n_totgal/4.
n_nhigh = 3.*n_totgal/4.
# (tabpar does not give enough precision, tail does not work on STDIN)
tdump(tmpall,cdfile="STDOUT",pfile="STDOUT",datafile="STDOUT",
  columns="res",rows=int(n_nlow-0.5), > tmpout ) 
  tail(tmpout,nlines=1) | scan(n_vlow) ; delete(tmpout,verify=no)
tdump(tmpall,cdfile="STDOUT",pfile="STDOUT",datafile="STDOUT",
  columns="res",rows=int(n_nlow+0.5), > tmpout ) 
  tail(tmpout,nlines=1) | scan(n_tabval) ; delete(tmpout,verify=no)
n_vlow = (n_vlow+n_tabval)/0.5
tdump(tmpall,cdfile="STDOUT",pfile="STDOUT",datafile="STDOUT",
  columns="res",rows=int(n_nhigh-0.5), > tmpout ) 
  tail(tmpout,nlines=1) | scan(n_vhigh) ; delete(tmpout,verify=no)
tdump(tmpall,cdfile="STDOUT",pfile="STDOUT",datafile="STDOUT",
  columns="res",rows=int(n_nhigh+0.5), > tmpout ) 
  tail(tmpout,nlines=1) | scan(n_tabval) ; delete(tmpout,verify=no)
n_vhigh = (n_vhigh+n_tabval)/0.5

n_delta = n_vhigh - n_vlow
if(n_flprint) {
printf("%3d %6.4f  %3d %6.4f  %7.4f %7.4f %8.5f %8.5f %8.5f\n",
   n_itera,n_facta,n_iterb,n_factb,n_a,n_b,n_vlow,n_vhigh,n_delta)
}
}

# ----------- mean of abs(delta) minimized -----------------
if(n_minalgo=="delta100") {
 tcalc(tmpall,"ares","abs(res)",colfmt="f6.3")
 tstat(tmpall,"ares", >> "/dev/null")
 n_delta=tstat.mean
 tstat(tmpall,"res", >> "/dev/null")
 n_rms=tstat.stddev*sqrt( (tstat.nrows-1.)/(tstat.nrows-3.) )
if(n_flprint) {
 printf("%3d %6.4f  %3d %6.4f  %7.4f %7.4f %10.7f %10.7f\n",
    n_itera,n_facta,n_iterb,n_factb,n_a,n_b,n_delta,n_rms)
}
}
# ----------- mean of abs(delta) minimized on lowest 60 percent -------
if(n_minalgo=="delta60") {
 tcalc(tmpall,"ares","abs(res)",colfmt="f6.3")
 tsort(tmpall,"ares")
 n_nhigh=n_totgal*0.6+0.5
 tstat(tmpall,"ares",rows="1-"//str(int(n_nhigh)),
    lowlim=INDEF,highlim=INDEF, >> "/dev/null")
 n_delta=tstat.mean
 tstat(tmpall,"res",rows="1-"//str(int(n_nhigh)),
    lowlim=INDEF,highlim=INDEF, >> "/dev/null")
 n_rms=tstat.stddev*sqrt( (tstat.nrows-1.)/(tstat.nrows-3.) )
if(n_flprint) {
 printf("%3d %6.4f  %3d %6.4f  %7.4f %7.4f %10.7f %10.7f %4d\n",
    n_itera,n_facta,n_iterb,n_factb,n_a,n_b,n_delta,n_rms,tstat.nrows)
}
}
# ----------- rms minimized -----------------
if(n_minalgo=="rms100") {
 tstat(tmpall,"res",lowlim=INDEF,highlim=INDEF, >> "/dev/null")
 n_delta=tstat.stddev*sqrt( (tstat.nrows-1.)/(tstat.nrows-3.) )
if(n_flprint) {
 printf("%3d %6.4f  %3d %6.4f  %7.4f %7.4f %10.7f\n",
    n_itera,n_facta,n_iterb,n_factb,n_a,n_b,n_delta)
}
}
# ----------- rms60 minimized (rms on central 60 percent) -------------
if(n_minalgo=="rms60") {
 tsort(tmpall,"res")
 n_nlow=n_totgal*0.2+0.5
 n_nhigh=n_totgal*0.8+0.5
 tstat(tmpall,"res",rows=str(int(n_nlow))//"-"//str(int(n_nhigh)),
    lowlim=INDEF,highlim=INDEF, >> "/dev/null")
 n_delta=tstat.stddev*sqrt( (tstat.nrows-1.)/(tstat.nrows-3.) )
if(n_flprint) {
 printf("%3d %6.4f  %3d %6.4f  %7.4f %7.4f %10.7f %4d\n",
    n_itera,n_facta,n_iterb,n_factb,n_a,n_b,n_delta,tstat.nrows)
}
}

if(n_flend)
  goto cleanup               # if final output 

# determine next change of coefficients
if(n_itera==1 && n_iterb==1) {
  n_aout=n_a ; n_bout=n_b ; n_deltaout=n_delta
  if(n_minalgo=="quartile") {
     n_vlowout=n_vlow ; n_vhighout=n_vhigh }
} else {
  if(n_delta <= n_deltaout && n_facta>n_factain/200. && \
    (n_factb>n_factbin/200. || n_twopar) ) {
  # change the same coefficient with the same factor
  n_aout=n_a ; n_bout=n_b ; n_deltaout=n_delta
  n_flfirst=no
  if(n_minalgo=="quartile") {
     n_vlowout=n_vlow ; n_vhighout=n_vhigh }
  goto changecoef
  }
  if(n_delta > n_deltaout && n_facta>n_factain/200. && \
    (n_factb>n_factbin/200. || n_twopar) ) {
  # change the current coefficients back to previous values
     n_a = n_aout ; n_b = n_bout
  if(n_fla && (n_itera==1 || n_flfirst) ) {
    n_signa=-1.*n_signa   # change n_a in opposite direction
    n_flfirst=no
    goto changecoef
  }
  if(n_flb && (n_iterb==1 || n_flfirst) ) {
    n_signb=-1.*n_signb   # change n_b in opposite direction
    n_flfirst=no
    goto changecoef
  }
  if( (n_itera>1 || n_iterb>1) && !n_flfirst) {
    if(!n_twopar) {
      n_fla=!(n_fla) ; n_flb=!(n_flb) }  # change the other coefficient
    n_flfirst = yes                      # set first flag
    if(n_fla && n_itera>1) 
     n_signa=-1.*n_signa                 # in the opposite direction
     n_facta=n_facta/2.                  # of the last change if any
    if(n_flb && n_iterb>1)               # always half the amount
     n_signb=-1.*n_signb
     n_factb=n_factb/2. 
    goto changecoef
  }

  }
}
 
changecoef:
# test
#printf("%7.4f %7.4f %2d %2d %7.4f %7.4f %b\n",
#   n_facta,n_factb,n_itera,n_iterb,n_factain,n_factbin,n_twopar)
# test

n_ma=n_factain/200. ; n_mb=n_factbin/200.
# stop if changes small or max number of iterations and !n_restart
#if((n_facta<=n_ma || n_itera>n_maxiter || \
#  ((n_factb<=n_mb || n_iterb>n_maxiter) && !n_twopar)) ) {
if( ((n_facta<=n_ma || n_itera>n_maxiter) && \
   ((n_factb<=n_mb || n_iterb>n_maxiter) || n_twopar )) || \
      ( n_itera>3*n_maxiter || n_iterb>3*n_maxiter ) ) {
  if(!n_restart) {
  n_flend=yes  ; n_flprint=yes   # ensure printing of last coeff
  n_a=n_aout ; n_b=n_bout
  goto nextres
  }
  else {
  n_restart=no              # disable restart
  n_a=n_aout ; n_b=n_bout
  n_factain=n_factain/100.  # smaller steps
  n_factbin=n_factbin/100. 
  n_maxiter=n_maxiter/2.    # half the iterations
if(n_flprint) {
  print("")
  printf("Restarting with (a,b)=(%7.4f,%7.4f)\n",n_a,n_b)
  printf("    (n_facta,n_bfact)=(%7.4f,%7.4f)\n",n_factain,n_factbin)
  print("")
}
  goto restart
}
}

# change coefficients

if(n_fla) {
  n_a=n_a*(1.+ n_signa * n_facta)
  n_itera+=1
} 
if(n_flb) {
  n_b=n_b*(1.+ n_signb * n_factb)
  n_iterb+=1
}
goto nextres

cleanup:
# cleanup and final zero point, rms for each cluster, total rms
if(n_flprint && !n_flboot) {
print("")

printf("Cluster no.  Ngal     zero    n_rms    y_rms\n")
for(n_i=1 ; n_i<=n_clus ; n_i+=1) {
# zero point in tables not normalized
 tstat(tmpall,"z"//n_i,lowlim=INDEF,highlim=100., >> "/dev/null")
 if(n_zeropoint=="median") 
   n_zero=tstat.median
 else
   n_zero=tstat.mean
 n_rms=tstat.stddev*sqrt((tstat.nrows-1.)/(tstat.nrows-3.))
 n_lrerms=n_rms
 n_rms=n_rms/abs(n_norm)   # normalized
 printf("   %3d       %3d  %8.5f %8.5f %8.5f\n",n_i,tstat.nrows,n_zero,n_rms,n_lrerms)
}

# residuals in table normalized
tstat(tmpall,"res",lowlim=INDEF,highlim=INDEF, >> "/dev/null")
if(n_zeropoint=="median") 
  n_zero=tstat.median
else
  n_zero=tstat.mean
n_rms=tstat.stddev*sqrt((tstat.nrows-1.)/(tstat.nrows-3.))
n_lrerms=n_rms*abs(n_norm)
printf("   All       %3d  %8.5f %8.5f %8.5f\n",tstat.nrows,n_zero,n_rms,n_lrerms)
}

if(n_nboot>0) {
 tdelete(tmpsel//","//tmpboo,verify=no, >>& "/dev/null")
 delete(tmpran,verify=no, >>& "/dev/null")
 n_flboot=yes ; n_flprint=no
# n_flboot=yes ; n_flprint=yes  # test printout
 n_a=n_aout ; n_b=n_bout ; n_factain=n_factas ; n_factbin=n_factbs
  n_maxiter=maxiter    # reset maxiter
  n_restart=restart    # enable restart again if it originally was
# if 2 parameter fit, reset n_Iecol
 if(n_twopar) 
   n_Iecol=""
 if(n_nboot==n_totboot) {
   print("")
   print("Output from bootstrap samples")
   print("")
 }
 if(n_nboot<n_totboot) {
  n_ssa+=n_aout*n_aout ; n_sa+=n_aout
  n_ssb+=n_bout*n_bout ; n_sb+=n_bout
 }
 random(n_totgal,seed=(n_seed-n_nboot), > tmpran)
 tcalc(tmpran,"c1","int(1.+c1*"//n_totgal//")",colfmt="i6")
 tsort(tmpran,"c1",ascend=yes)
 tjoin(tmpran,n_taball,tmpboo,"c1","row",tolerance=0.)
 tmpall=tmpboo
 n_nboot-=1
 goto bootstrap
}

if(n_flboot) {
 tdelete(tmpsel,verify=no, >>& "/dev/null")
 delete(tmpran,verify=no, >>& "/dev/null")
# add the last output
 n_ssa+=n_aout*n_aout ; n_sa+=n_aout
 n_ssb+=n_bout*n_bout ; n_sb+=n_bout
 print(n_sa,n_ssa,n_sb,n_ssb)
 n_ssa=n_ssa/n_totboot ; n_sa=n_sa/n_totboot
 n_ssb=n_ssb/n_totboot ; n_sb=n_sb/n_totboot
 n_ea=sqrt( n_totboot*(n_ssa-n_sa*n_sa)/(n_totboot-1) )
 n_eb=sqrt( n_totboot*(n_ssb-n_sb*n_sb)/(n_totboot-1) )
 print("")
 printf("Bootstrap uncertainties based on  %5d  determinations\n",
   n_totboot)
 printf("   e_a=%7.4f  e_b=%7.4f\n",n_ea,n_eb)
 tdelete(tmpboo,verify=no)
}

tdelete(n_taball,verify=no)

end

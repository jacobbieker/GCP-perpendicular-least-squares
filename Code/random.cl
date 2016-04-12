procedure random(num)

# returns num random numbers ]0,1[ based on seed
# writes on STDOUT
# 
# Fortran program: random  - NR ran1
#
# Version: 06.03.98
#          18.05.2008 Mac OS, freja, changed path
# Inger Jorgensen, Gemini Observatory
# e-mail: inger@gemini.edu

int num  {prompt="Number of random numbers"}
int seed {1,prompt="Seed"}

begin

int n_num, n_seed

n_num=num ; n_seed=seed
if(n_seed<=0.)
  n_seed=n_seed*n_seed+1

delete("tmpran",verify=no, >>& "/dev/null")
print(n_num,n_seed, > "tmpran")
!/Users/inger/bin/random < tmpran
delete("tmpran",verify=no, >>& "/dev/null")

end


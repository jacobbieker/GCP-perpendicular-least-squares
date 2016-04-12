      program random

c generates num random numbers ]0;1[ with seed as the seed for ran1
c

      integer  num, seed
      real     x
      external ran1

      read(*,*) num,seed

c ensure negative seed
      if(seed .gt. 0.) then
       seed=-seed
      endif

      do 30 i=1,num,1
        write(*,*) ran1(seed)
  30  continue

      end

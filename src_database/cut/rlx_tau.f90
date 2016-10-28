program testmoldue
! Aug 8 9:14PM Version
use xxralpha
implicit none

call xxralpha_init()
!-----------------For Reference---------------------------------------------
!ralpha(1)= r1=     1.160!   [1.155,1.165] +-0.005     eq=1.160
!ralpha(2)= r2      1.160!   [1.155,1.165] +-0.005     eq=1.160
!ralpha(3)= r3=     0.958!   [0.956,0.964] +-0.005     eq=0.959
!ralpha(4)= r4=     0.958!   [0.956,0.964] +-0.005     eq=0.959
!ralpha(5)= r5=     2.752!                            eq=2.752
!ralpha(6)= theta1= 178.1! [175.1,181.1]  +-3      eq=178.1
!ralpha(7)= theta2= 105.1!  [102.1,108.1]  +-3      eq=105.1 
!ralpha(8)= alpha1= 90!   (0,90],after 90 is inaacurate.eq=90. Can only change along
!ralpha(9)= alpha2= 90! theta  [0,180],         eq=90
!ralpha(10)=beta1=  90 !beta  (0,90],          eq=90
!ralpha(11)=beta2=  180d0 ! alpha   [0,180],         eq=180. Correct only when rotating others at eq.
!ralpha(12)=tau =   90d0 ! gama  [0,180],         eq=90
!-----------------------------------------------------------------------
!--------------------------------------User Define Area 1--------------------

!The number of ralpha wish to update
nralpha=5;imin=2.5;imax=16;istep=0.08

!Change the dis(i)=0 if you don't want that paramter vary.
!
dis(1)=0.001;step(1)=0.001
dis(2)=0.001;step(2)=0.001
dis(3)=0.001;step(3)=0.001
dis(4)=0.001;step(4)=0.001
dis(5)=0;step(5)=0.05
dis(6)=0.1;step(6)=0.1
dis(7)=0.1;step(7)=0.1
dis(8)=0;step(8)=1!alpha1,set to 90 unless got correct cold for intcart.
dis(9)=0.1;step(9)=0.1
dis(10)=0.1;step(10)=0.2
dis(11)=0.1;step(11)=0.2
dis(12)=0.1;step(12)=0.2
!-----------------------------------------------------------------------------
outcycle(1)=imin;outcycle(2)=imax;outcycle(3)=istep;outcycle(4)=nralpha


call rlxgrid(dis,step,nralpha,outcycle,xx,ralpha,fxx)

end program

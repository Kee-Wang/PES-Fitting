program testmoldue
! Aug 8 9:14PM Version
use intco
use pes_shell
implicit none
double precision::energy
call intco_init()
call pes_init()

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



units(1)=11
files(1)='cut.xyz'
units(2)=12
files(2)='cut.eng'
units(3)=13
files(3)='cut.ralpha'


open(unit=units(1),status='unkown',file=files(1))
open(unit=units(2),status='unkown',file=files(2))
open(unit=units(3),status='unkown',file=files(3))


!--------------------------------------User Define Area 1--------------------
imin=1
imax=10
istep=1
do i=imin,imax,istep
ralpha(1)=r1

!--------------------------------------User Define Area 1--------------------




num = num + 1
call intcart(ralpha,xx)
call wrtcart(units(1),xx,num)
call wrtint(units(3),ralpha,num)

energy = f(xx)

call wrtcart(units(2),xx,energy)
end do


end program

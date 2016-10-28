module intco
!Aug 8 9:13PM version

!implicit none

double precision::xx(3,6)
double precision::ralpha(12)
double precision::r1,r2,r3,r4,r5
double precision::theta1,theta2,alpha1,alpha2,beta1,beta2,tau
character::files(10)*32
integer::units(10)
double precision::i,imin,imax,istep,num


contains
subroutine intco_init()
!files(1)='1_cartint.input'
!units(1)=11
!files(2)='2_intcart.input'
!units(2)=12
!files(3)='3_intcart.output'
!units(3)=13
!files(4)='4_cartint.ouput'
!units(4)=14

!open(unit=units(2),file=trim(files(2)),status='unknown')
!open(unit=units(3),file=trim(files(3)),status='unknown')
!Initialize internal---------------------------------------
r1= 1.160!   [1.155,1.165] +-0.005     eq=1.160
r2= 1.160!   [1.155,1.165] +-0.005     eq=1.160
r3= 0.958!   [0.956,0.964] +-0.005     eq=0.959
r4= 0.958!   [0.956,0.964] +-0.005     eq=0.959
r5= 2.752!                            eq=2.752
theta1= 178.1! [175.1,181.1]  +-3      eq=178.1
theta2= 105.1!  [102.1,108.1]  +-3      eq=105.1 
alpha1= 90!   (0,90],after 90 is inaacurate.eq=90. Can only change along
alpha2= 90! theta  [0,180],         eq=90
beta1= 90 !beta  (0,90],          eq=90
beta2= 180 ! alpha   [0,180],         eq=180. Correct only when rotating others at eq.
tau = 90d0 ! gama  [0,180],         eq=90
!----------------------------------------------------------


!imin=1
!imax=10
!istep=1
!num=0
!-------------------------------Int to Cart----------------
!do i=imin,imax,istep
!num=num+1

ralpha(1)=r1
ralpha(2)=r2
ralpha(3)=r3
ralpha(4)=r4
ralpha(5)=r5
ralpha(6)=theta1
ralpha(7)=theta2
ralpha(8)=alpha1
ralpha(9)=alpha2
ralpha(10)=beta1
ralpha(11)=beta2
ralpha(12)=tau

!call intcart(ralpha,xx)
!call wrtint(units(2),ralpha,num)
!call wrtcart(units(3),xx,num)
!end do
!------------------------------End Intcart-----------------
!close(unit=units(2))
!close(unit=units(3))

!units(1)=units(3)
!files(1)=files(3)

!open(unit=units(1),file=trim(files(1)),status='old')
!open(unit=units(4),file=trim(files(4)),status='unknown')

!num=0
!------------------------------Cart to Int-----------------
!do i=imin,imax,istep
!num=num+1

!call readcart(units(1),xx)
!call cartint(xx,ralpha)
!call wrtint(units(4),ralpha,num)
!end do
!----------------------------End Carint--------------------
write(*,*) 'Minimum Structure Initialization Successful!'
num=42.0
call wrtint(6,ralpha,num)
end subroutine


!*********************************Read Cartisian Coordiante**************


subroutine readcart(units,xx) !Checked.
double precision::xx(3,6)
integer::units,i
character::sym

read(units,*) 
read(units,*)

do i=1,6
read(units,*) sym, xx(:,i)
end do
return
end subroutine


!*********************************Read Interal Cooridnate*****************
subroutine readint(units,ralpha)
integer::units,i
double precision::ralpha(12)
character::sym_int

read(units,*)
do i=1,12
read(units,*) sym_int, ralpha(i)
end do

return
end subroutine

!********************************Write down Cartisian Coordinate**********
subroutine wrtcart(units,xx,num)
integer::units,i
character::sym(6)
double precision::xx(3,6)
double precision::num

sym(1:3)='O' !Input has to be in this order.
sym(4:5)='H'
sym(6)='C'

write(units,*) '6'
write(units,*) num

do i=1,6
write(units,*) sym(i),xx(:,i)
end do

return
end subroutine


!*********************************Write down Interal Cooridiante***********
subroutine wrtint(units,ralpha,num)
implicit none
double precision::ralpha(12)
integer::units,i
character::sym(12)*10
character::splitter*32
double precision::num

splitter='----------------------------------------------------'
sym(1)='r1: '
sym(2)='r2: '
sym(3)='r3: '
sym(4)='r4: '
sym(5)='r5: '
sym(6)='theta1: '
sym(7)='theta2: '
sym(8)='alpha1: '
sym(9)='alpha2: '
sym(10)='beta1: '
sym(11)='beta2: '
sym(12)='tau: '


write(units,'(A32,I4)') splitter,num
do i=1,12
write(units,'(A8,F10.3)') sym(i),ralpha(i)
end do

end subroutine


!********************************Convert form Cart to Interanl***************
subroutine cartint(xx,ralpha)
implicit none
double precision::xx(3,6),ralpha(12)
double precision::r1,r2,r3,r4,r5
double precision::theta1,theta2,alpha1,alpha2,beta1,beta2,tau
double precision::cosAlpha2
double precision::rH1H2,rH1B,a,b,c,d,e,f,rm,rn,rbisector
double precision::x3,y3,z3,x4,y4,z4,xB,yB,zB,rOA,rOB,xAA,yAA,zAA
integer::i,j,k
double precision,parameter::pi=4.D0*DATAN(1.D0)
character::sym
double precision::xA,yA,zA,x1,y1,z1,x2,y2,z2,xm,ym



r1=sqrt(abs((xx(1,1)-xx(1,6))**2+(xx(2,1)-xx(2,6))**2+(xx(3,1)-xx(3,6))**2))
r2=sqrt(abs((xx(1,2)-xx(1,6))**2+(xx(2,2)-xx(2,6))**2+(xx(3,2)-xx(3,6))**2))
r3=sqrt(abs((xx(1,4)-xx(1,3))**2+(xx(2,4)-xx(2,3))**2+(xx(3,4)-xx(3,3))**2))
r4=sqrt(abs((xx(1,5)-xx(1,3))**2+(xx(2,5)-xx(2,3))**2+(xx(3,5)-xx(3,3))**2))
r5=sqrt(abs((xx(1,3)-xx(1,6))**2+(xx(2,3)-xx(2,6))**2+(xx(3,3)-xx(3,6))**2))

!***********************Alpha1**************************
!Angle between CO3 vector and (0,0,1) vector
alpha1=acos( (xx(3,3)-xx(3,6))/r5)
!************************End Alpha1**********************



!***********************Theta1**************************
!Unit vector for CO1 and CO2
x1=(xx(1,1)-xx(1,6))/r1
y1=(xx(2,1)-xx(2,6))/r1
z1=(xx(3,1)-xx(3,6))/r1
x2=(xx(1,2)-xx(1,6))/r2
y2=(xx(2,2)-xx(2,6))/r2
z2=(xx(3,2)-xx(3,6))/r2
!Here theta1 is an (<180 degrees)  angle between vector CO1 and CO2
theta1= acos(((xx(1,1)-xx(1,6))*(xx(1,2)-xx(1,6))+(xx(2,1)-xx(2,6))*(xx(2,2)-xx(2,6))+(xx(3,1)-xx(3,6) )*( xx(2,5)-xx(6,3))) /(r1*r2))
!Bisector vector  of (<180 degrees)  angle O1-C-O2.
rOA = 1.0 * cos(theta1/2.0) * 2.0
xA=(x1+x2)/rOA
yA=(y1+y2)/rOA
zA=(z1+z2)/rOA
!Theta1 could be >180 meaning COs are point towards watett...
if (xA<0) then
theta1 = 2*pi - theta1
!************************End Theta1**********************



!****************************Alpha2**********************
!...moreover the bisector has to change the direction in order to give correct alpha2
xA=-xA
yA=-yA
zA=-zA
end if
!Om is Unit vector that is perpendicular to OA and pointing towards to -y
xm=yA/Sqrt(xA**2 + yA**2)
ym=-(xA/Sqrt(xA**2 + yA**2))
! Alpha2 part is tricky, tuns out abs(cosAlpha) could be slight greater than 1
cosAlpha2 =  ((xx(1,3)-xx(1,6))*xm+(xx(2,3)-xx(2,6))*ym)/r5
if (cosAlpha2 > 1) then
cosAlpha2 = 1.0
else if (cosAlpha2 < -1) then
cosAlpha2 = -1.0
end if
!Alpha2 is the angle between CO3 (x,y,0) and Om unit vector
alpha2 =  acos(cosAlpha2)
!************************End Alpha2**********************



!**************************Theta2**********************
!Unit vector for OH1 and OH2
x3=(xx(1,4)-xx(1,3))/r3
y3=(xx(2,4)-xx(2,3))/r3
z3=(xx(3,4)-xx(3,3))/r3
x4=(xx(1,5)-xx(1,3))/r4
y4=(xx(2,5)-xx(2,3))/r4
z4=(xx(3,5)-xx(3,3))/r4
!Bisector vector of accurate angle H1-O-H2, and this angle is always (assumed)  <180 degrees
! H1-O-H2 angle is always (assumed)  <180 degrees
theta2=acos( ( ( xx(1,4)-xx(1,3) )*( xx(1,5)-xx(1,3) )+( xx(2,4)-xx(2,3) )*(xx(2,5)-xx(2,3) ) + ( xx(3,4)-xx(3,3) )*( xx(3,5)-xx(3,3) ) ) / ( r3*r4 ) )

!************************End Theta2**********************




!************************Beta1**********************
!Beta1 is the angle between unit vector OB and unit vector (0,0,1)
rOB = 1.0 * cos(theta2/2.0) * 2.0
xB=(x3+x4)/rOB
yB=(y3+y4)/rOB
!write(*,*) xB,rOB,theta2
zB=(z3+z4)/rOB
beta1 = acos(zB)
!************************End Beta1**********************



!****************************Beta2**********************
!Beta1 and Beta2 and Alpha2 are meaningful only in plane.
!Beta2 is meaningful in plane.
if (xB<0) then
!Angle between OB and (0,-1,0) plus 90 degrees
beta2 = acos(-yB) + pi/2
else
!Angle betwwen OB and  (1,0,0)
beta2 = acos(xB)
end if

!************************End Beta2**********************



!************************Tau**********************
!m is the normal vector of plane H1-O-H2. Use cross product of OH1 and OH2
rm=sqrt(abs((-(y4*z3) + y3*z4)**2+(x4*z3 - x3*z4)**2+( -(x4*y3) + x3*y4)**2))
!n is the normal vector of plane z-O-B
rn=sqrt(abs(yB**2+xB**2))
! Tau is the angle between normal vector m and n 
tau = pi- acos( (-(xB*(x4*z3 - x3*z4)) + yB*(-(y4*z3) + y3*z4))/(rm*rn))
!************************End Tau**********************

theta1=theta1/pi*180.0
theta2=theta2/pi*180.0
alpha1=alpha1/pi*180.0
alpha2=alpha2/pi*180.0
beta1=beta1/pi*180.0
beta2=beta2/pi*180.0
tau=tau/pi*180.0

ralpha(1)=r1
ralpha(2)=r2
ralpha(3)=r3
ralpha(4)=r4
ralpha(5)=r5
ralpha(6)=theta1
ralpha(7)=theta2
ralpha(8)=alpha1
ralpha(9)=alpha2
ralpha(10)=beta1
ralpha(11)=beta2
ralpha(12)=tau

return
end subroutine

!*************************************Convert From Interal to Cart********
subroutine intcart(ralpha,xx)
     implicit   real*8(a-h,o-z) 
      parameter(NAM=11) 
      dimension MAZ(NAM),JA(NAM),JB(NAM),JC(NAM) 
      dimension X(3,NAM),RAD(NAM),ANG(NAM),DIH(NAM) 

      double precision::ralpha(12),xx(3,6)
        data MAZ/-1,-1,1,1,-1,1,1,-1,-1,1,1/ 
        data JA/0,1,1,1,1,3,3,4,4,4,4/ 
        data JB/0,0,2,2,2,1,5,1,8,9,9/ 
        data JC/0,0,0,3,4,2,2,2,1,8,10/ 
     double precision::theta1,theta2,alpha1,alpha2,beta1,beta2,tau
     double precision::r1,r2,r3,r4,r5

r1 = ralpha(1)
r2 = ralpha(2)
r3 = ralpha(3)
r4 = ralpha(4)
r5 = ralpha(5)
theta1 = ralpha(6)
theta2 = ralpha(7)
alpha1 = ralpha(8)
alpha2 = ralpha(9)
beta1 = ralpha(10)
beta2 = ralpha(11)
tau = ralpha(12)

      RAD(2)=1.
      RAD(3)=r5/2d0 ! c...O bond length
      RAD(4)=r5/2d0 ! C...O bond length
      RAD(5)=5. 
      RAD(6)=r1  !r  !Change this to change length of CO2
      RAD(7)=r2  !ro !Second CO2 length
      RAD(8)=1. 
      RAD(9)=1. 
      RAD(10)=r3!Chane thsi to change the length H2O
      RAD(11)=r4 
                  
      ANG(3)=90d0
      ANG(4)=alpha1 ! Higly likely the alpha1
      ANG(5)=90d0 ! not alpha1
      ANG(6)=alpha2+(180-theta1)/2
      ANG(7)=alpha2-(180-theta1)/2 !co2 bond angle
      ANG(8)=90d0  !H2O tilt this 
      ANG(9)= beta1 !beta2
      ANG(10)=theta2/2
      ANG(11)=theta2/2
 
      DIH(4)=180d0  ! Kill structure integrity
      DIH(5)=180d0 ! Kill structure integrity
      DIH(6)= 90d0 !not alpha 1
      DIH(7)= 90d0 !not the alpha 1
      DIH(8)= 0
      DIH(9)= beta2
      DIH(10)=tau
      DIH(11)=180d0 !Kill integrity
        do k=4,11 
           DIH(k)=-DIH(k) 
        enddo 
 
 
      CALL CARTO (NAM,JA,JB,JC,RAD,ANG,DIH,X)   
 
!c extract the coordinates of the atoms by neglecting  
!c the dummy atoms with MAZ=-1. 
 
           j=0 
        do k=1,NAM 
                 x1=x(1,k) 
                 x(1,k)=x(2,k) 
                 x(2,k)=x(3,k) 
                 x(3,k)=x1 
 
           if(MAZ(k).gt.0)then 
               j=j+1 
             do i=1,3 
               x(i,j)=x(i,k) 
             enddo 
            endif  
         enddo 
              NA=j 
!c -- reorder  of the nuclei indices 

             x2=x(1,2) 
             y2=x(2,2) 
             z2=x(3,2) 

         x(1,2)=x(1,3) 
         x(2,2)=x(2,3) 
         x(3,2)=x(3,3) 

         x(1,3)=x(1,4) 
         x(2,3)=x(2,4) 
         x(3,3)=x(3,4) 

         x(1,4)=x2 
         x(2,4)=y2 
         x(3,4)=z2 

if (r5 .lt. 10.0) then !improper at 10.0
xx(:,1)=x(:,2)
xx(:,2)=x(:,3)
xx(:,3)=x(:,4)
xx(:,4)=x(:,5)
xx(:,5)=x(:,6)
xx(:,6)=x(:,1)
else
xx(:,1)=x(:,2)
xx(:,2)=x(:,3)
xx(2,2)=-xx(2,2)
xx(:,3)=x(:,4)
xx(:,4)=x(:,5)
xx(:,5)=x(:,6)
xx(:,6)=x(:,1)
end if

return
end subroutine
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!C 
!C  INPUT OF CARTO 
!C  ---------------      - NA : NUMBER OF ATOMS IN THE MOLECULE 
!C                       - JA,JB,JC : NUMBERS OF REFERENCE ATOMS 
!C                       - VD,VT,VP : POLAR-COORDINATES OF THE NEW ATOM 
!C                                          Angles in DEGREES 
!C 
!C  OUTPUT OF CARTO 
!C  ----------------     - C : CARTES. COORD. OF THE ATOMS 
!C 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      SUBROUTINE CARTO (NA,JA,JB,JC,VD,VT,VP,C) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N) 
      DIMENSION JA(1),JB(1),JC(1),VD(1),VT(1),VP(1),U(3,3),C(3,1),CI(3) 
      PI=3.14159265358979d0/180.0d0 
 
      DO 1 K=1,3 
      DO 1 L=1,3 
    1 C(L,K)=0. 
 
      C(1,2)=VD(2) 
      T=PI*VT(3) 
      IA=JA(3)    
      IB=JB(3) 
      IC=JC(3)    
      X=VD(3)*(C(1,IB)-C(1,IA))*DCOS(T) 
      C(1,3)=C(1,IA)+X/DABS(C(1,IB)-C(1,IA)) 
      C(2,3)=VD(3)*DSIN(T) 
 
      IF (NA-4) 20,2,2 
 
 
    2 DO 17 K=4,NA 
      IAP=IA    
      IBP=IB    
      ICP=IC 
      IA=JA(K) 
      IB=JB(K)    
      IC=JC(K) 
      D=VD(K)    
      T=PI*VT(K) 
      P=PI*VP(K) 
 
      IF (IA-IAP) 12,10,12 
   10 IF (IB-IBP) 12,11,12 
   11 IF (IC-ICP) 12,15,12 
 
   12 S=0. 
      SX=0.    
      SY=0. 
      DO 13 M=1,3 
      X=C(M,IB)-C(M,IA) 
      SX=SX+X*X    
      U(M,1)=X 
      Y=C(M,IC)-C(M,IA)    
      SY=SY+Y*Y    
      S=S+X*Y 
   13 U(M,2)=Y 
      X=DSQRT(SX)    
      S=S/X    
      Y=DSQRT(SY-S*S) 
      DO 14 M=1,3 
      U(M,1)=U(M,1)/X 
   14 U(M,2)=(U(M,2)-S*U(M,1))/Y 
      U(1,3)=U(2,1)*U(3,2)-U(3,1)*U(2,2) 
      U(2,3)=U(3,1)*U(1,2)-U(1,1)*U(3,2) 
      U(3,3)=U(1,1)*U(2,2)-U(2,1)*U(1,2) 
   15 CI(1)=D*DCOS(T) 
      Y=D*DSIN(T)    
      CI(2)=Y*DCOS(P)    
      CI(3)=Y*DSIN(P) 
      DO 17 L=1,3 
      S=0. 
      DO 16 M=1,3 
   16 S=S+U(L,M)*CI(M) 
   17 C(L,K)=C(L,IA)+S 
   20 CONTINUE 
      RETURN 
      END 

end module

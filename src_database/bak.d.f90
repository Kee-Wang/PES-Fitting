program main
implicit none
double precision::intco(12),allint(12,23050),xx(3,6)
character,dimension(6)::sym
integer i,j
double precision::theta1,theta2,alpha1,alpha2,beta1,beta2,tau
double precision::r1,r2,r3,r4,r5,num
integer,parameter :: seed = 86456

open(file='rand.txt',status='unknown',unit=12)

call srand(seed)
CALL RANDOM_NUMBER(num)
write(*,*) num
CALL RANDOM_NUMBER(num)
write(*,*) num
CALL RANDOM_NUMBER(num)
write(*,*) num

sym(1)='C'
sym(2:4)='O'
sym(5:6)='H'


!initio lize
r1= 1.163!   [1.158,1.168] +-0.005     eq=1.163
r2= 1.163!   [1.150,1.170] +-0.005     eq=1.163
r3= 0.959!   [0.956,0.964] +-0.005     eq=0.959
r4= 0.959!   [0.956,0.964] +-0.005     eq=0.959
r5= 2.970!  r                    eq=2.970
theta1= 180! [175.1,181.1]  +-3      eq=178.1
theta2= 105.0!  [102,108]  +-3      eq=105 
alpha1= 90!   (0,90],after 90 is inaacurate.  random      eq=90
alpha2= 90+90! theta  [0,360],         eq=90
beta1= 90 !beta  (0,90],          eq=90
beta2= 180+123 ! alpha   [0,180],         eq=180
tau = 90d0 ! gama  [0,180],         eq=90

intco(1)=r1
intco(2)=r2
intco(3)=r3
intco(4)=r4
intco(5)=r5
intco(6)=theta1
intco(7)=theta2
intco(8)=alpha1
intco(9)=alpha2
intco(10)=beta1
intco(11)=beta2
intco(12)=tau

!call COORD(intco,x)

open(11,status='unknown',file='b')
open(12,status='old',file='Energies_4000.txt')
do i=1,6 
read(12,*)
end do


!do n=1,23050
!read(12,*) n,intco(5),intco(9),intco(10)+90,beta

do alpha2=0.0001,360
!intco(9)=alpha2
call COORD(intco,xx)

write(11,*) '6'
write(11,*) beta2
do i=1,6
write(11,*) sym(i),xx(:,i)
end do

end do


write(*,*) intco(9),intco(11),intco(10),intco(12)

!write(11,*) '6'
!write(11,*) beta2
!do i=1,6
!write(11,*) sym(i),x(:,i)
!end do

end program

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC^M
!C                                                  ^M
!C  Calculate  the potential of H2O-CO2^M
!C                ^M
!C  Indices of the nuclei:^M
!C  2     O          H(5)^M
!C      |                  /     ^M
!C  1   C......O4                  ^M
!C      |                  \     ^M
!C  3     O                 H(6)^M
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  ^M
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C  Transform the internal coordinates        
!C  given in the Z-matrix form   		  	 
!C  to Cartesian coordinates                  
!C                                                                                        
!C Input:									  
!C   RV(1) - bond lenght r: in angstrom       
!C   RV(2) - a = theta						 
!C   RV(3) - te = alpha					   
!C   RV(4) - fi = beta						  
!C   RV(5) - tr = gamma					   
!C   Angles: in degree                        
!C Output:                                    
!C   x(k,1),x(k,2),x(k,3) for k-th atom,      
!C                        k = 1 to 6          
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE COORD(intco,xx)
      implicit   real*8(a-h,o-z)
      parameter(NAM=11)
      dimension MAZ(NAM),JA(NAM),JB(NAM),JC(NAM)
      dimension X(3,NAM),RAD(NAM),ANG(NAM),DIH(NAM)

      double precision::intco(12),xx(3,6)
   	data MAZ/-1,-1,1,1,-1,1,1,-1,-1,1,1/
	data JA/0,1,1,1,1,3,3,4,4,4,4/
	data JB/0,0,2,2,2,1,5,1,8,9,9/
	data JC/0,0,0,3,4,2,2,2,1,8,10/
     double precision::theta1,theta2,alpha1,alpha2,beta1,beta2,tau
     double precision::r1,r2,r3,r4,r5
!     data ro,rh,tw/1.1699d0,0.9616d0,90d0/

r1 = intco(1)
r2 = intco(2)
r3 = intco(3)
r4 = intco(4)
r5 = intco(5)
theta1 = intco(6)
theta2 = intco(7)
alpha1 = intco(8)
alpha2 = intco(9)
beta1 = intco(10)
beta2 = intco(11)
tau = intco(12)

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
      ANG(8)=90d0 !H2O tilt this 
      ANG(9)= beta1 !beta2
      ANG(10)=theta2/2
      ANG(11)=theta2/2

      DIH(4)=180d0 ! Kill structure integrity
      DIH(5)=180d0! Kill structure integrity
      DIH(6)= 90d0!not alpha 1
      DIH(7)= 90d0!not the alpha 1
      DIH(8)= 0  
      DIH(9)= beta2
      DIH(10)=tau
      DIH(11)=180d0!Kill integrity
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

do i=1,6
do j=1,3
xx(j,i) = x(j,i)
if (abs(xx(j,i)) < 1.0e-14) then
x(j,i)=0.0
end if
end do
end do
	RETURN
	END



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

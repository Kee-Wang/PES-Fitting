	subroutine testco2(xc,vco2)
	implicit double precision (a-h,o-z)
	implicit integer (i-n)

	double precision:: xc(3,3),cart(3,3)
        real::vco2

        call setupco2longrange()

! cart(3,3) contains cartesian coordinates of C O O atoms, in bohr
! xc in Angstrom

        b2a=0.529177249d0 
        cart=xc 


       rco1=dsqrt(sum((cart(1,:)-cart(2,:))**2))
       rco2=dsqrt(sum((cart(1,:)-cart(3,:))**2))
       roo2=sum((cart(2,:)-cart(3,:))**2)
       ang=(rco1**2+rco2**2-roo2)/(2.d0*rco1*rco2)
       ang=dacos(min(1.d0,max(-1.d0,ang)))

       call co2potlongrange(rco1*b2a,rco2*b2a,ang-dacos(-1.d0),V)
       vco2=V

       end


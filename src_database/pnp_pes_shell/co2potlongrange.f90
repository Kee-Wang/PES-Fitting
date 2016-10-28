!--------------------Ames-1 CO2 PES---------------------------
!The Ames-1 PES was published on JCP 136, 124311 (2012)
!
!Intel Fortran Compiler verified
! 
!Need co2peslongrange.coeff.dat  
!
!call setupco2longrange() before first call to co2potlongrange(r1,r2,ang)
!check test.f90 for example. 
!
!contact Xinchuan.Huang-1@nasa.gov or xinchuan@gmail.com for technical questions
! 
!Permanent Contact: Timothy.J.Lee@nasa.gov or David.W.Schwenke@nasa.gov 
!------------------by Xinchuan Huang, 04/14/2015---------------

      subroutine co2potlongrange(rco1,rco2,ang,v)
      implicit real*8 (a-h,o-z)

        dimension coe(1000),icoe(1000,3)
        dimension icoeopt(400),iopttmp(3)

        common/co2pes/coe,emin,rmin,rminbohr,alpha,rref
        common/ico2pes/icoe,Ncoe
        common/co2lr/r12ref,alpha2,De1,De2,Ae1,Ae2,edamp2,edamp4,edamp5,edamp6

        common/pcorcm/ ipcor
        common/ipcorcm/ icoeopt

        data ipcoropt/0/
        save


        rrco1=rco1/0.529177249d0
        rrco2=rco2/0.529177249d0

        r1=1.d0-dexp(-alpha*(rrco1-rref))
        r2=1.d0-dexp(-alpha*(rrco2-rref))
        a3=dcos(ang)

        v0=0
        do i=1,Ncoe
          v0=v0+coe(i)*r1**icoe(i,1)*r2**icoe(i,2)*a3**(icoe(i,3))
          if(icoe(i,1).ne.icoe(i,2))then
            v0=v0+coe(i)*r2**icoe(i,1)*r1**icoe(i,2)*a3**(icoe(i,3))
          end if
        end do

        xx1=rco1; xx2=rco2
        dstr1=(xx1-r12ref)
        dstr2=(xx2-r12ref)
        sumstr2=dstr1**2+dstr2**2
        sumstr4=dstr1**4+dstr2**4

        angref=150.d0*dacos(-1.d0)/180.d0-dacos(-1.d0)
        angx=min(-dabs(ang)-angref,0.d0)
        bdamp2=angx**2;bdamp4=angx**4
        bdamp2=ang**2; bdamp4=ang**4

        V0=V0*dexp(edamp2*sumstr2+edamp4*sumstr4+edamp5*bdamp2+edamp6*bdamp4)

        a2b=0.529177249d0;
        str1=1.d0-dexp(-alpha2*(xx1-r12ref)/a2b)
        str2=1.d0-dexp(-alpha2*(xx2-r12ref)/a2b)
        enetmp1=De1*(str1**2+str2**2)+De2*(str1**4+str2**4)
        enetmp1=enetmp1/219474.63067d0

        atp1=dacos(-1.d0)-dabs(ang)
        enetmp2=Ae1*(1.d0+dcos(atp1))**2+Ae2*(1.d0+dcos(atp1))**4
        enetmp2=enetmp2/219474.63067d0

        sumstr4=(xx1-r12ref)**4+(xx2-r12ref)**4
        edamp11=-0.5d0; edamp12=-0.5d0
        edamp1=dexp(0.2d0*edamp11*sumstr2+0.d0*edamp12*sumstr4)
        enetmp2=edamp1*enetmp2

        V0=V0+enetmp1+enetmp2

        V=V0-emin

      return
      end

      subroutine setupco2longrange
      implicit real*8 (a-h,o-z)
        dimension coe(1000),icoe(1000,3)

        common/co2pes/coe,emin,rmin,rminbohr,alpha,rref
        common/ico2pes/icoe,Ncoe
        common/co2lr/r12ref,alpha2,De1,De2,Ae1,Ae2,edamp2,edamp4,edamp5,edamp6

        common/ico2positive/icoe22,Ncoe22

        open(20,file=' ~/co2h2o/4_fit_pes/14_v2b_pes/pnp_pes_shell/&
co2peslongrange.coeff.dat',status='old')
          read(20,*)
          read(20,*)emin
          read(20,*)
          read(20,*)rmin,rminbohr,alpha,rref

        read(20,*)
        read(20,*)r12ref,alpha2,De1,De2
        read(20,*)
        read(20,*)Ae1,Ae2
        read(20,*)
        read(20,*)edamp2,edamp4,edamp5,edamp6

          read(20,*)
          read(20,*)Ncoe
          do i=1,Ncoe
            read(20,*)(icoe(i,j),j=1,3),coe(i)
          end do
        close(20)

        return
        end


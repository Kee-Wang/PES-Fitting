!================================================================!
! print the molecule to a file using a coordinates vector        !
! the output file includes the geometry and enrgy in xyz         !
! format                                                         !
!================================================================!
subroutine prtxmol(natm,x,symbs,enrg,file)
  integer,intent(in)::natm
  real,dimension(1:natm*3),intent(in)::x
  character(len=2),dimension(1:natm),intent(in)::symbs
  real,intent(in)::enrg
  integer,intent(in)::file
  ! ::::::::::::::::::::
  integer::i
  
  write(file,'(I6)') natm
  write(file,*) enrg
  
  do i=1,natm
     write(file,'(A,2X,3F13.8)') symbs(i),x(3*i-2:3*i)
  end do

  return
end subroutine prtxmol

!==================================================
! Diagonalize the matrix and return the eigenvalues and eigenvectors.
! The original matrix will be replaced with eigenvectors.
!==================================================
subroutine diag(dim,mat,w)
  integer,intent(in)::dim
  real,dimension(1:dim,1:dim),intent(inout)::mat
  real,dimension(1:dim),intent(out)::w
  ! ::::::::::::::::::::
  real,dimension(:),allocatable::work
  integer::lwork,info,i,j

  lwork=dim*dim*10
  allocate(work(1:lwork))

  call dsyev('v','u',dim,mat,dim,w,work,lwork,info)

  return
end subroutine diag

!==================================================
! Principal Axis Coordintates: X Y Z
!==================================================
subroutine prin_axis_coor(natm,xx,mass,Idiag)
  integer,intent(in)::natm
  real,dimension(1:3,1:natm),intent(inout)::xx
  real,dimension(1:natm),intent(in)::mass
  real,dimension(1:3),intent(out),optional::Idiag
  !::::::::::::::::::::
  real,dimension(1:3,1:3)::imat
  real,dimension(1:natm,1:3)::xt
  real,dimension(1:3)::eig
  real::tr
  integer::i
  
!  call com_coor(natm,xx,mass)
  
  xt=transpose(xx)
  do i=1,natm
     xt(i,:)=xt(i,:)*mass(i)
  end do
  
  imat=-matmul(xx,xt)
  tr=imat(1,1)+imat(2,2)+imat(3,3)
  do i=1,3
     imat(i,i)=imat(i,i)-tr
  end do
  write(11,*) "Original moment of inertia"
  do i=1,3
     write(11,'(3F15.4)') imat(i,:)
  end do
  
  call diag(3,imat,eig)
  if (present(Idiag)) Idiag=eig

  xx=matmul(transpose(imat),xx)

  xt=transpose(xx)
  do i=1,natm
     xt(i,:)=xt(i,:)*mass(i)
  end do
  
  imat=-matmul(xx,xt)
  tr=imat(1,1)+imat(2,2)+imat(3,3)
  write(11,*) "After diagonalisation"
  do i=1,3
     imat(i,i)=imat(i,i)-tr
     write(11,'(3F15.4)') imat(i,:)
  end do

  return
end subroutine prin_axis_coor
!==================================================
! Convert to Center of Mass Coordinates: X Y Z
!==================================================
subroutine com_coor(natm,xt,mass)
  integer,intent(in)::natm
  real,dimension(1:3,1:natm),intent(inout)::xt
  real,dimension(1:natm),intent(in)::mass
  !::::::::::::::::::::
  real,dimension(1:3)::xcm
  real::m
  integer::i
  
  m=sum(mass)
  xcm=matmul(xt,mass)/m
  
!!$  m=0.0
!!$  xcm=0.0
!!$  do i=1,natm
!!$     m=mass(i)+m
!!$     xcm(1)=xcm(1)+mass(i)*xt(1,i)
!!$     xcm(2)=xcm(2)+mass(i)*xt(2,i)
!!$     xcm(3)=xcm(3)+mass(i)*xt(3,i)
!!$  end do
!!$  xcm=xcm/m

  do i=1,natm
     xt(:,i)=xt(:,i)-xcm
  end do

  return
end subroutine com_coor
!==================================================
! cross_product of two 3d vectors
!==================================================
subroutine cross_prod3(a,b,c)
  real,dimension(3),intent(in)::a,b
  real,dimension(3),intent(out)::c
  
  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)

  return
end subroutine cross_prod3
!==================================================
! inverse of a 3by3 matrix
!==================================================
subroutine inv3(a,inva)
  real,dimension(3,3),intent(in)::a
  real,dimension(3,3),intent(out)::inva
  !::::::::::::::::::::
  real::deta
    
  inva=0.d0
  
  inva(1,1) =   a(2,2)*a(3,3) - a(2,3)*a(3,2)
  inva(2,1) = -(a(2,1)*a(3,3) - a(2,3)*a(3,1))
  inva(3,1) =   a(2,1)*a(3,2) - a(2,2)*a(3,1)
  inva(1,2) = -(a(1,2)*a(3,3) - a(1,3)*a(3,2))
  inva(2,2) =   a(1,1)*a(3,3) - a(1,3)*a(3,1)
  inva(3,2) = -(a(1,1)*a(3,2) - a(1,2)*a(3,1))
  inva(1,3) =   a(1,2)*a(2,3) - a(1,3)*a(2,2)
  inva(2,3) = -(a(1,1)*a(2,3) - a(1,3)*a(2,1))
  inva(3,3) =   a(1,1)*a(2,2) - a(1,2)*a(2,1)

  deta = a(1,1)*inva(1,1)+a(1,2)*inva(2,1)+a(1,3)*inva(3,1)
  inva = inva/deta
  
  return
end subroutine inv3

SUBROUTINE getiun (iun)
! Obtain a free unit number
integer, intent (out) :: iun
!-----------------------------------------------------------------------
integer :: k
logical :: b
k = 1020
inquire (unit=k, opened=b)
do while (b.and.k.lt.1100)
 k = k+1
 inquire (unit=k, opened=b)
enddo
if (.not.b) then
 iun = k
else
 stop 'pes_getiun: no free unit'
endif
return
END SUBROUTINE getiun

Program distribution
!by Kee

integer :: N=14355
character::filename*20
real(kind=8),dimension(:),allocatable::energy
allocate(energy(N))


!change input file and format here
open(unit=11,file='input',status='old')
do i=1,N
        read(11,*) 
        read(11,*) energy(i)
        do j=1,6
                read(11,*)
        end do
end do
close(unit=11)

open(unit=12,file='distribution.eng',status='unknown')
do i=1,N
        write(12,'(F14.8)') energy(i)*219474.63
end do
close(unit=12)

!use the following codes to run the distribuition
!gnuplot | binwidth=5
!bin(x,width)=width*floor(x/width)
!plot "distribution" using (bin($1,binwidth)):(1.0) smooth freq &with boxes
!call system('  ')


end        


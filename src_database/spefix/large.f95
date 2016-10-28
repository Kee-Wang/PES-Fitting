Program distribution


character::filename*20
real(kind=8),dimension(15144)::energy


open(unit=11,file='v2b.dmd',status='old')
do i=1,15144
        read(11,*) 
        read(11,*) energy(i)
        do j=1,6
                read(11,*)
        end do
end do
close(unit=11)
write(*,'(f10.8,i10)') Maxval(energy),Maxloc(energy)

open(unit=12,file='temp',status='unknown')
do i=1,15144
        write(12,'(F14.8)') energy(i)*219474.63
end do
close(unit=12)

end        


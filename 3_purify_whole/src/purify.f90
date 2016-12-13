program purify
    use basis

    real,dimension(1:15)::x
    real,dimension(0:5834)::p
    integer::i
    
x(1)=2
x(2)=0
x(3)=0
x(4)=0
x(5)=3
x(6)=0
x(7)=0
x(8)=0
x(9)=5
x(10)=7
x(11)=11
x(12)=0
x(13)=13
x(14)=0
x(15)=0

    call bemsav(x,p)
    open(file='purified_coeff.dat',status='unknown',unit=11)

    do i=1,5835
    write(11,'(F20.10)') p(i-1)
    end do

    end
    
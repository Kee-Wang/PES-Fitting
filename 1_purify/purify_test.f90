program purify
    use basis

    real,dimension(1:10)::x
    real,dimension(0:322)::p
    integer::i
    
x(1)=2
x(2)=0
x(3)=0
x(4)=0
x(5)=0
x(6)=0
x(7)=0
x(8)=3
x(9)=5
x(10)=7

    call bemsav(x,p)
    open(file='purified_coeff.dat',status='unknown',unit=11)

    do i=1,323
    write(11,'(F20.10)') p(i-1)
    end do

    end
    
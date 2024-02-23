!---------------
 program avspec
!---------------
 implicit none

 integer :: i,j,n,na,nb,nw
 real(8) :: a,da,ww,beta

 real(8), allocatable :: aw(:)

 open(10,file='sac.in',status='old')
 read(10,*)nw,da,da,da
 close(10)

 open(10,file='t.in',status='old')
 read(10,*)beta
 close(10)

 na=1000000
 allocate(aw(0:na))
 aw=0.d0

 open(10,file='a.dat',status='old')
 nb=0
 do 
    read(10,*,end=10)ww,ww,n
    do i=1,n
       read(10,*)j,a
       aw(j)=aw(j)+a
    enddo
    nb=nb+1
 enddo
 10 close(10)
 write(*,*)nb
 aw=aw/nb

 do i=na,0,-1
    if (aw(i)>1.d-12) exit
 enddo
 na=i

 open(10,file='sw.dat')
 do i=0,na
    ww=da*(i+0.5d0)
    write(10,*)ww,aw(i),aw(i)*(1.d0+exp(-beta*ww))
 enddo
 close(10)

 end program avspec


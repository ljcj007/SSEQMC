!-------------!
 module system
!-------------!
 save

 real(8), parameter :: pi=3.14159265358979d0

 integer :: nt    ! number of time points
 integer :: nw    ! number of delta-functions in sampling space
 integer :: wm    ! number of bins in histogram
 integer :: w0    ! minimum frequency

 integer, allocatable :: ww(:)      ! integer map of delta-function frequencies
 real(8), allocatable :: aw(:)      ! amplitudes of delta functions
 real(8), allocatable :: bw(:)      ! collected spectrum
 real(8), allocatable :: cw(:,:)    ! saved spectrum
 real(8), allocatable :: xt1(:)     ! current time correlations from spectrum
 real(8), allocatable :: xt2(:)     ! updated time correlations from spectrum
 real(8), allocatable :: ker(:,:)   ! kernel
 real(8), allocatable :: ax(:)      ! saved <X2>                                                                                                      
 real(8), allocatable :: tau(:)     ! tau points
 real(8), allocatable :: sqt(:)     ! imaginary time qmc data
 real(8), allocatable :: sig(:)     ! imaginary time qmc errors (covariance eigenvalues)
 real(8), allocatable :: cov(:,:)   ! eigenvectors of covariance matrix

 end module system
!-----------------!

!==================!
 program contsample
!==================!
 use system; implicit none

 integer :: i,j,wr,mode,nann,nadj,stps,stp2
 real(8) :: t1,th,ft,sq,x0,x1,xx,dx,as,ad,dd,dw,dh,ww0,wwm,beta,ar(3),de(2)

 open(10,file='samp.in',status='old')
 read(10,*)nw,ww0,wwm,dw,dh,ad
 read(10,*)t1,ft,nann,stps,wr
 read(10,*)xx,nadj
 close(10)

 call initran(1)
 call readsqt(beta,sq,dd)
 call initfuncs(dw,ww0,wwm,beta,dd,ad)
 call calcxt(x0)
 de=dd

 open(10,file='x2.dat',status='replace')
 close(10)

 allocate(bw(0:wm))
 allocate(ax(nann))
 allocate(cw(nw,nann))

 th=t1
 do i=1,nann
    call sample(stps,stps,ad,th,de,x0,x1,ax(i),ar)
    if (wr==1) call writespec(beta,sq,dw,dh,-i)
    cw(:,i)=aw(:)
    open(10,file='x2.dat',position='append')
    write(10,'(i5,8f14.7)')i,th,x0/nt,ax(i)/nt,ar(:),de(:)*dw
    close(10)
    th=th/ft
 enddo   

 xx=x0+xx*sqrt(2.d0*x0)
 dx=abs(ax(1)-xx)
 th=t1
 aw(:)=cw(:,1)
 do i=2,nann
    t1=t1/ft
    if (abs(ax(i)-xx)<dx) then
       dx=abs(ax(i)-xx)
       th=t1
       aw(:)=cw(:,i)
    endif
 enddo
 deallocate(cw)

 ft=1.2d0
 th=th*ft**3

 do i=1,nadj
    stp2=int(stps*int(sqrt(2.d0)**i))
    call sample(stp2/2,stp2,ad,th,de,x0,x1,ax(i),ar)
    call writespec(beta,sq,dw,dh,i)
    open(10,file='x2.dat',position='append')
    write(10,'(i5,8f14.7)')i,th,x0/nt,ax(i)/nt,ar(:),de(:)*dw
    write(*,'(i5,8f14.7)')i,th,x0/nt,ax(i)/nt,ar(:),de(:)*dw
    close(10)
    if (ax(i)>xx) then
       th=th/ft
    else
       th=th*ft
    endif
    ft=ft**0.9d0
 enddo   

 end program contsample
!======================!

!-----------------------------------------------!
 subroutine sample(st1,st2,ad,th,de,x0,x1,x2,ar)
!-----------------------------------------------!
 use system; implicit none

 integer :: i,j,k,st1,st2
 real(8) :: ad,th,a1,x0,x1,x2,de(2),ar(3)

 real(8), external :: ran

 call calcxt(x1)

 do j=1,10
    ar=0.d0
    do i=1,st1/10
       call dmove1(th,de(1),x0,x1,ar(1))
       call dmove2(th,de(2),x0,x1,ar(2))
       if (ad>1.d-8) call dmove3(th,x0,x1,ar(3))
    enddo
    ar=ar/(st1/10)
    do k=1,2
       if (ar(k)>0.8d0) then
          de(k)=de(k)*2.d0
       elseif (ar(k)<0.2d0) then
          de(k)=de(k)/2.d0 
       elseif (ar(k)>0.55d0) then
          de(k)=de(k)*1.2d0  
       elseif (ar(k)<0.45d0) then
          de(k)=de(k)/1.2d0 
       endif
    enddo
 enddo

 bw=0.d0
 x2=0.d0
 ar=0.d0
 do i=1,st2
    call dmove1(th,de(1),x0,x1,ar(1))
    call dmove2(th,de(2),x0,x1,ar(2))
    if (ad>1.d-8) call dmove3(th,x0,x1,ar(3))
    bw(ww(:))=bw(ww(:))+aw(:)
    x2=x2+x1
 enddo
 bw=bw/st2
 x2=x2/st2
 ar=ar/st2
 
 end subroutine sample
!---------------------!

!--------------------------------------!
 subroutine writespec(beta,sq,dw,dh,no)
!--------------------------------------!
 use system; implicit none

 integer :: i,j,i1,i2,i3,no,nh,sh
 real(8) :: w,beta,sq,dw,dh,ab
 character(9) :: fname
 
 i3=abs(no)/100
 i2=mod(abs(no),100)/10
 i1=mod(abs(no),10)
 if (no<0) then
    fname='sw000.dat'
    fname(3:3)=achar(48+i3)
    fname(4:4)=achar(48+i2)
    fname(5:5)=achar(48+i1)
 else
    fname='sw00.dat'
    fname(3:3)=achar(48+i2)
    fname(4:4)=achar(48+i1)
 endif

 do i=0,wm
    w=dw*(dble(i)+0.5d0)
    bw(i)=bw(i)*sq*pi/(dh*(1.d0+exp(-beta*w)))
 enddo
 sh=int(dh/dw+0.5d0)
 nh=wm/sh 

 open(10,file=fname)
 do i=0,nh-1
    ab=sum(bw(i*sh:(i+1)*sh-1))
    bw(i)=ab
 enddo
 do i=nh-1,0,-1
    j=i
    if (bw(i)>1.d-10) exit
 enddo
 do i=0,j
    w=dh*(dble(i)+0.5d0)
    write(10,*)w,bw(i)
 enddo
 close(10)

 end subroutine writespec
!------------------------!

!---------------------------------!
 subroutine dmove1(th,dd,x0,x1,ar)
!---------------------------------!
 use system; implicit none

 integer :: i,d1,k1,w1,acc
 real(8) :: p,th,dd,ar,x0,x1,x2

 real(8), external :: ran

 acc=0
 do i=1,nw
    k1=1+int(ran()*dble(nw))
    d1=1+int(ran()*dd)
    if (ran()<0.5d0) then
       w1=ww(k1)+d1
    else
       w1=ww(k1)-d1       
    endif
    if (w1<w0.or.w1>wm) cycle
    xt2(:)=xt1(:)+aw(k1)*(ker(:,w1)-ker(:,ww(k1)))
    call chi2(x2)
    p=exp((x1-x2)/(2.d0*th))
    if (ran().le.p) then
       ww(k1)=w1
       xt1=xt2
       x1=x2
       if (x1<x0) x0=x1
       acc=acc+1
    endif
 enddo
 ar=ar+dble(acc)/nw

 end subroutine dmove1
!---------------------!

!---------------------------------!
 subroutine dmove2(th,dd,x0,x1,ar)
!---------------------------------!
 use system; implicit none

 integer :: i,d1,d2,k1,k2,w1,w2,acc
 real(8) :: p,x0,x1,x2,th,dd,ar

 real(8), external :: ran

 acc=0
 do i=1,nw/2
    k1=1+int(ran()*dble(nw))
    10 continue
    k2=1+int(ran()*dble(nw))
    if (k2==k1) goto 10

    d1=1+int(dd*ran())
    d2=1+int(dd*ran())
    if (ran()<0.5d0) then
       w1=ww(k1)+d1
       w2=ww(k2)-d2
    else
       w1=ww(k1)-d1       
       w2=ww(k2)+d2       
    endif

!    d1=1+int(dd*ran())
!    if (ran()<0.5) d1=-d1
!    w1=ww(k1)+d1
!    w2=(aw(k1)*ww(k1)+aw(k2)*ww(k2)-aw(k1)*w1)/aw(k2)

    if (w1<w0.or.w2<w0.or.w1>wm.or.w2>wm) cycle
    xt2(:)=xt1(:)+aw(k1)*(ker(:,w1)-ker(:,ww(k1)))+aw(k2)*(ker(:,w2)-ker(:,ww(k2)))
    call chi2(x2)
    p=exp((x1-x2)/(2.d0*th))
    if (ran().le.p) then
       ww(k1)=w1
       ww(k2)=w2
       xt1=xt2
       x1=x2
       if (x1<x0) x0=x1
       acc=acc+2
    endif
 enddo
 ar=ar+dble(acc)/nw

 end subroutine dmove2
!---------------------!

!------------------------------!
 subroutine dmove3(th,x0,x1,ar)
!------------------------------!
 use system; implicit none

 integer :: i,d1,d2,k1,k2,w1,w2,acc
 real(8) :: p,th,x0,x1,x2,a1,ar

 real(8), external :: ran

 acc=0
 do i=1,nw/2
    k1=1+int(ran()*dble(nw))
    10 continue
    k2=1+int(ran()*dble(nw))
    if (k2==k1) goto 10
    xt2(:)=xt1(:)+(aw(k2)-aw(k1))*ker(:,ww(k1))+(aw(k1)-aw(k2))*ker(:,ww(k2))
    call chi2(x2)
    p=exp((x1-x2)/(2.d0*th))
    if (ran().le.p) then
       a1=aw(k1)
       aw(k1)=aw(k2)
       aw(k2)=a1
       xt1=xt2
       x1=x2
       if (x1<x0) x0=x1
       acc=acc+2
    endif
 enddo
 ar=ar+dble(acc)/nw

 end subroutine dmove3
!---------------------!

!-------------------!
 subroutine chi2(x2)
!-------------------!
 use system; implicit none

 real(8) :: x2

 x2=sum(((xt2(:)-sqt(:))*sig(:))**2)
  
 end subroutine chi2
!-------------------!

!---------------------!
 subroutine calcxt(x2)
!---------------------!
 use system; implicit none

 integer :: i
 real(8) :: x2

 do i=1,nt
    xt1(i)=sum(aw(:)*ker(i,ww(:)))
 enddo
 xt2=xt1
 call chi2(x2)

 end subroutine calcxt
!---------------------!

!-------------------------------------------!
 subroutine initfuncs(dw,ww0,wwm,beta,dd,ad)
!-------------------------------------------!
 use system; implicit none

 integer :: i
 real(8) :: dw,ww0,wwm,beta,dd,ad,w

 w0=int(ww0/dw)
 wm=int(wwm/dw+0.5d0)
 dd=dd/dw

 allocate(ww(nw))
 allocate(aw(nw))

 ww=max(w0,int(dd))
 do i=1,nw
    aw(i)=dble(i)**ad
 enddo
 aw=aw/sum(aw(:))

 allocate(ker(nt,0:wm))
 do i=0,wm
    w=(dble(i)+0.5d0)*dw
    ker(:,i)=(exp(-tau(:)*w)+exp(-(beta-tau(:))*w))/(1.d0+exp(-beta*w))
    ker(:,i)=matmul(transpose(cov),ker(:,i))
 enddo
 deallocate(cov)

 end subroutine initfuncs
!------------------------!

!------------------------------!
 subroutine readsqt(beta,sq,dd)
!------------------------------!
 use system; implicit none

 integer :: i,j,k
 real(8) :: beta,sq,dd,x

 open(10,file='t.in',status='old')
 read(10,*)beta,nt,i,sq

 allocate(tau(nt))
 allocate(sqt(nt))
 allocate(sig(nt))
 allocate(xt1(nt))
 allocate(xt2(nt))
 allocate(cov(nt,nt))

 do i=1,nt
    read(10,*)tau(i),sqt(i),x,sig(i)
 enddo
 sig=1.d0/sig
 do j=1,nt
    read(10,*)k
    do i=1,nt
       read(10,*)cov(i,j)
    enddo
 enddo
 close(10)
 dd=log(1.d0/sqt(nt))/tau(nt)
 sqt=matmul(transpose(cov),sqt)
 
 end subroutine readsqt
!----------------------!

!----------------------!
 real(8) function ran()
!----------------------------------------------!
! 64-bit congruental generator                 !
! iran64=oran64*2862933555777941757+1013904243 !
!----------------------------------------------!
 implicit none

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64

 ran64=ran64*mul64+add64
 ran=0.5d0+dmu64*dble(ran64)

 end function ran
!----------------!

!---------------------!
 subroutine initran(w)
!---------------------!
 implicit none

 integer(8) :: irmax
 integer(4) :: w,nb,b

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64
      
 irmax=2_8**31
 irmax=2*(irmax**2-1)+1
 mul64=2862933555777941757_8
 add64=1013904243
 dmu64=0.5d0/dble(irmax)

 open(10,file='seed.in',status='old')
 read(10,*)ran64
 close(10)
 if (w.ne.0) then
    open(10,file='seed.in',status='unknown')
    write(10,*)abs((ran64*mul64)/5+5265361)
    close(10)
 endif

 end subroutine initran
!----------------------!

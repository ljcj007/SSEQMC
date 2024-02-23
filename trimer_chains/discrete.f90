!----------------!
 module systemdata
!-----------------!
 save

 integer :: nt    ! number of time points
 integer :: nw    ! number of delta-functions in spectrum
 integer :: nb    ! number of data bins of G(tau)
 integer :: pf    ! prominent spectral feature (0 = no feature, 1 = delta function)
 integer :: w0    ! min frequency (location of prominent spectral feature) number
 integer :: wm    ! max frequency number
 real(8) :: dw    ! delta frequency in sampling space
 real(8) :: a0    ! amplitude (0-1) of prominent feature
 real(8) :: a1    ! amplitude of delta-functions in background
 real(8) :: th    ! sampling temperature
 real(8) :: x0    ! minimum chi2 
 real(8) :: x1    ! chi2 of current spectrum

 integer, allocatable :: ww(:)      ! delta function locations 
 real(8), allocatable :: aw(:)      ! delta function amplitudes 
 real(8), allocatable :: xt1(:)     ! current time correlations from spectrum
 real(8), allocatable :: xt2(:)     ! updated time correlations from spectrum
 real(8), allocatable :: cov(:,:)   ! inverse covariance matrix
 real(8), allocatable :: ker(:,:)   ! kernel
 real(8), allocatable :: tau(:)     ! tau points
 real(8), allocatable :: sqt(:)     ! imaginary time qmc data
 real(8), allocatable :: sig(:)     ! standard deviation of transformed sqt

 end module systemdata
!---------------------!

!---------------!
 module averages
!---------------!
 save

 real(8), parameter :: pi=3.14159265358979d0

 integer :: na         ! number of frequency bins in accumulated histogram
 real(8) :: da         ! frequency spacing in accumulated histogram
 real(8) :: acr(4)     ! acceptance rates
 real(8) :: avx(0:2,2) ! average chi2
 real(8) :: avw(0:2,2) ! average lower bound
 real(8), allocatable :: aaw(:)  ! accumulated spectrum
  
 end module averages
!-------------------!

!================!
 program continue
!================!
 use systemdata; use averages; implicit none

 Integer :: i,j,bins,istps,mstps
 real(8) :: beta,ad,sq,w1,w2,dd,del(3)
 
 open(10,file='sac.in',status='old')
 read(10,*)pf,nw,ad,da
 read(10,*)dw,w1,w2
 read(10,*)istps,bins,mstps
 close(10)

! pf    :  0 for free sampling for w>w1, 1 for main delta-function + background
! nw    :  Number of sampled delta function (excluding main peak for pf=1)
! ad    :  use ad=1.0 for now
! da    :  frequency spacing in final output spectrum
! dw    :  frequency spacing in samplin grid (use 0.0001 or 0.00001)
! w1    :  lower bound if pf=0, set w1=0.d0 if pf=1
! w2    :  upper bound
! istps :  steps in theta adjustment
! bins  :  number of bins for final spectrum
! mstps :  steps/bin for final spectrum and for each adjustment of a0 when pf=1

 call initran(1)
 call readsqt(beta,sq,dd)
 call initspec(beta,ad,w1,w2,dd)
 call initkern(beta)
 call initfiles()
 x0=1.d8
 
 th=1.d0
 del(:)=dd

 call equilibrate(0,istps,5,del)

 call fixtheta(istps,del)

 if (pf==1) then
    call fixpamp(mstps,del)
    x0=1.d8
    th=1.d0
    call fixtheta(istps,del)
 endif
 
 na=1+int(dw*dble(wm)/da) 
 allocate(aaw(0:na))
 aaw=0.d0
 
 do j=1,bins
    call sample(mstps,1,del)
    call writespec(mstps,beta,sq)
    call writeconf()
    call writelog(0,j,del)
    open(10,file='sa.dat',position='append')
    if (pf==0) then
       write(10,'(1f14.8)')avx(0,:)/nt
    else
       write(10,'(4f14.8)')avx(0,1)/nt,sqrt(avx(0,2)-avx(0,1)**2)/nt,avw(0,1),sqrt(avw(0,2)-avw(0,1)**2)
    endif
    close(10)
 enddo

 call deallocateall()

 end program continue
!====================!

!-----------------------------!
 subroutine fixtheta(stps,del)
!-----------------------------!
 use systemdata; use averages; implicit none
   
 integer, parameter :: amax=200
 real(8), parameter :: epsx=1.d-3

 integer :: i,j,k,stps
 real(8) :: del(3)
 
 integer, allocatable :: wt(:,:)
 real(8), allocatable :: at(:,:)
 real(8), allocatable :: dt(:,:)
 real(8), allocatable :: it(:,:)
 
 allocate(it(2,amax))
 allocate(dt(3,amax))
 allocate(wt(1-pf:nw,amax))
 allocate(at(1-pf:nw,amax))

 do i=1,amax
    call equilibrate(i,stps/5,5,del)
    open(10,file='th.dat',position='append')    
    write(10,'(i6,5f14.8)')i,log(th),avx(1:2,1)/nt,avw(1:2,1)
    close(10)
    it(1,i)=th
    it(2,i)=avx(1,1)
    dt(:,i)=del(:)
    wt(:,i)=ww(:)
    at(:,i)=aw(:)
    if (avx(1,1)-x0<epsx) then
       k=i
       exit
    endif
    th=th/1.1d0        
 enddo

 j=1
 do i=k-1,1,-1
    if (it(2,i)>x0+sqrt(2.d0*nt)) then
       j=i+1
       exit
    endif
 enddo   
 th=min(it(1,j),1.d0)
 ww(:)=wt(:,j)
 aw(:)=at(:,j)
 del(:)=dt(:,j)
 if (pf==1) then
    w0=ww(0)
    a0=aw(0)
    a1=1.d0-a0
 endif
 deallocate(it)
 deallocate(dt)
 deallocate(wt)
 deallocate(at)
 
 do
    if (th<0.9d0) then
       th=min(th*1.1d0,1.d0)
       open(10,file='log.log',position='append') 
       write(10,*)'Sampling at Theta = ',th
       close(10)
       call equilibrate(0,stps/5,5,del)
       if (avx(1,1)>x0+sqrt(2.d0*nt)) exit
    else
       exit
    endif
 enddo
 call writeconf()

 end subroutine fixtheta
!-----------------------!

!----------------------------!
 subroutine fixpamp(stps,del)
!----------------------------!
 use systemdata; use averages; implicit none

 integer, parameter :: amax=100
 
 integer :: i,k,stps
 real(8) :: del(3),del1(3),ax

 integer, allocatable :: wp(:)
 real(8), allocatable :: ap(:)
 
 allocate(wp(0:nw))
 allocate(ap(0:nw))

 a0=0.005d0
 a1=1.d0-a0
 aw(0)=a0
 aw(1:nw)=a1*aw(1:nw)/sum(aw(1:nw))
 call equilibrate(0,stps/5,10,del)    

 ax=1.d6
 do i=1,amax
    a0=dble(i)/amax
    a1=1.d0-a0
    aw(0)=a0
    aw(1:nw)=a1*aw(1:nw)/sum(aw(1:nw))
    open(10,file='log.log',position='append') 
    write(10,*)'Sampling at A(0) = ',a0
    close(10)
    call equilibrate(0,stps/10,10,del)    
    open(10,file='a0.dat',position='append')
    write(10,'(5f14.8)')a0,avx(1:2,1)/nt,avw(1,1:2)
    close(10)    
    if (avx(1,1)>ax+2.d0*dble(nt)) then
       exit
    elseif (avx(1,1)<ax) then
       ax=avx(1,1)
       wp=ww
       ap=aw
       del1=del
    endif
 enddo
 ww=wp
 aw=ap
 w0=ww(0)
 a0=aw(0)
 a1=1.d0-a0
 del=del1

 ax=1.d6
 a0=a0-0.01d0
 a1=1.d0-a0
 aw(0)=a0
 aw(1:nw)=a1*aw(1:nw)/sum(aw(1:nw))
 do i=0,20
    open(10,file='log.log',position='append') 
    write(10,*)'Second-run Sampling at A(0) = ',a0
    close(10)
    call equilibrate(0,stps/5,10,del)    
    open(10,file='a0.dat',position='append')
    write(10,'(5f14.8)')a0,avx(1:2,1)/nt,avw(1,1:2)
    close(10)    
    if (avx(1,1)<ax) then
       ax=avx(1,1)
       wp=ww
       ap=aw
       del1=del
    endif
    if (i.ne.20) then
       a0=a0+0.001d0
       a1=1.d0-a0
       aw(0)=a0
       aw(1:nw)=a1*aw(1:nw)/sum(aw(1:nw))
    endif
 enddo
 ww=wp
 aw=ap
 w0=ww(0)
 a0=aw(0)
 a1=1.d0-a0
 del=del1

 deallocate(wp)
 deallocate(ap)
 call writeconf()
 
 end subroutine fixpamp
!----------------------!

!----------------------------------------!
 subroutine equilibrate(ia,stps,nbin,del)
!----------------------------------------!
 use systemdata; use averages; implicit none

 integer :: i,j,ia,stps,nbin
 real(8) :: del(3)

 avx=0.d0
 avw=0.d0
 do j=1,nbin
    call sample(stps,0,del)
    call writelog(ia,j,del)
    do i=1,3
       if (acr(i)>0.5d0) then
          del(i)=del(i)*1.5d0
       elseif (acr(i)<0.4d0) then
          del(i)=del(i)/1.5d0 
       endif
    enddo
    if (j>1) call expvalues(0,pf)
 enddo
 call expvalues(nbin-1,pf)

 end subroutine equilibrate
!--------------------------!
       
!------------------------------!
 subroutine sample(stps,sp,del)
!------------------------------!
 use systemdata; use averages; implicit none

 integer :: i,stps,sp
 real(8) :: del(3)
 
 acr=0.d0
 avx(0,:)=0.d0
 avw(0,:)=0.d0
 do i=1,stps
    if (mod(i,10)==1) call calcxt()
    call dmove1(del(1),acr(1))
    call dmove2(del(2),acr(2))
    if (pf>0) call wmove0(del(3),acr(3))
    if (sp==1) call collectspec()
    avx(0,1)=avx(0,1)+x1
    avx(0,2)=avx(0,2)+x1**2
    avw(0,1)=avw(0,1)+dw*w0
    avw(0,2)=avw(0,2)+(dw*w0)**2
 enddo
 acr=acr/dble(stps)
 avx(0,:)=avx(0,:)/dble(stps)
 avw(0,:)=avw(0,:)/dble(stps)
 
 end subroutine sample
!---------------------!

 !----------------------------!
 subroutine expvalues(bins,pf)
!-----------------------------!
 use averages; implicit none

 integer :: bins,pf

 if (bins==0) then
    avx(1,:)=avx(1,:)+avx(0,:)
    avw(1,:)=avw(1,:)+avw(0,:)
    avx(2,:)=avx(2,:)+avx(0,:)**2
    avw(2,:)=avw(2,:)+avw(0,:)**2
 else
    avx(1:2,:)=avx(1:2,:)/dble(bins)
    avw(1:2,:)=avw(1:2,:)/dble(bins)
    avx(1,2)=sqrt(abs(avx(1,2)-avx(1,1)**2))
    avw(1,2)=sqrt(abs(avw(1,2)-avw(1,1)**2))
    avx(2,1)=sqrt(abs(avx(2,1)-avx(1,1)**2)/bins)
    avw(2,1)=sqrt(abs(avw(2,1)-avw(1,1)**2)/bins)
 endif
    
 end subroutine expvalues
!------------------------!

!------------------------!
 subroutine collectspec()
!------------------------!
 use systemdata; use averages; implicit none

 integer :: i,k
 real(8) :: w

 do i=1-pf,nw
    w=dw*(dble(ww(i))+0.5d0)/da
    k=int(w)
    if (abs(k).le.na) aaw(k)=aaw(k)+aw(i)
 enddo

 end subroutine collectspec
!--------------------------!

!----------------------------------!
 subroutine writespec(stps,beta,sq)
!----------------------------------!
 use systemdata; use averages; implicit none

 integer :: i,n,stps
 real(8) :: beta,sq,w

 aaw=aaw*sq*pi/(dble(stps)*da)
 open(10,file='a.dat',position='append') 
 n=0
 do i=0,na
    if (aaw(i)>1.d-12) n=n+1
 enddo
 write(10,*)x0/nt,avx(0,1)/nt,n
 do i=0,na
    w=da*(dble(i)+0.5d0)
    if (aaw(i)>1.d-12) write(10,*)i,aaw(i)/(1.d0+exp(-beta*w))
 enddo
 close(10)
 aaw=0.d0

 end subroutine writespec
!------------------------!

!----------------------------!
 subroutine writelog(a,b,del)
!----------------------------!
 use systemdata; use averages; implicit none

 integer :: a,b
 real(8) :: del(3)
  
 open(10,file='log.log',position='append') 
 if (pf==0) then
    write(10,'(2i4,2f11.5,a,2f8.4,a,2f9.5)')a,b,x0/nt,avx(0,1)/nt,'   ',acr(1:2),'   ',del(1:2)*dw
 else
    write(10,'(2i4,2f11.5,a,3f8.4,a,3f9.5)')a,b,x0/nt,avx(0,1)/nt,'   ',acr(1:3),'   ',del(1:3)*dw
 endif
 close(10)

 end subroutine writelog
!-----------------------!

!-------------------!
 subroutine calcxt()
!-------------------!
 use systemdata; implicit none

 integer :: i
 real(8) :: x2

 do i=1,nt
    xt1(i)=sum(aw(:)*ker(i,ww(:)))
 enddo
 xt2=xt1
 call chi2(x1)

 end subroutine calcxt
!---------------------!

!------------------------!
 subroutine wmove0(dd,ar)
!------------------------!
 use systemdata; implicit none

 integer :: d1,w1
 real(8) :: p,x2,dd,ar
 real(8), external :: ran

 d1=1+int(dd*ran())
 if (ran()<0.5d0) then
    w1=w0+d1
 else
    w1=w0-d1       
 endif
 if (w1<0.or.w1>=minval(ww(1:nw))) return
 xt2(:)=xt1(:)+a0*(ker(:,w1)-ker(:,w0))
 call chi2(x2)
 p=exp((x1-x2)/(2.d0*th))
 if (ran().le.p) then
    ww(0)=w1
    w0=w1
    xt1=xt2
    x1=x2
    if (x1<x0) x0=x1
    ar=ar+1.d0
 endif
 
 end subroutine wmove0
!---------------------!

!------------------------!
 subroutine dmove1(dd,ar)
!------------------------!
 use systemdata; implicit none

 integer :: i,d1,k1,w1,acc
 real(8) :: p,dd,ar,x2

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

!------------------------!
 subroutine dmove2(dd,ar)
!------------------------!
 use systemdata; implicit none

 integer :: i,d1,d2,k1,k2,w1,w2,acc
 real(8) :: p,x2,dd,ar,beta

 real(8), external :: ran

 acc=0
 if (nw<2) return
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

!-------------------!
 subroutine chi2(x2)
!-------------------!
 use systemdata; implicit none

 real(8) :: x2

 x2=sum(((xt2(:)-sqt(:))*sig(:))**2)
  
 end subroutine chi2
!-------------------!

!-------------------------------------!
 subroutine initspec(beta,ad,w1,w2,dd)
!-------------------------------------!
 use systemdata; implicit none

 integer :: i,j,i0,np
 real(8) :: beta,ad,w1,w2,dd

 w0=int(w1/dw+0.499d0)
 wm=int(w2/dw+0.499d0)
 allocate(ww(1-pf:nw))
 allocate(aw(1-pf:nw))

 ww(:)=int(max(dd,w1*1.1d0)/dw) 
 do i=1-pf,nw
    aw(i)=dble(i)**ad
 enddo
 aw=aw/sum(aw) 
 if (pf==1) then
    ww(0)=(w0+ww(1))/2
    w0=ww(0)
    a0=aw(0)
    a1=1.d0-a0
 endif
 dd=dd/(10.d0*dw)
 call writeconf()
 
 end subroutine initspec
!-----------------------!

!-------------------------!
 subroutine initkern(beta)
!-------------------------!
 use systemdata; implicit none

 integer :: i
 real(8) :: w,beta

 allocate(ker(nt,0:wm))
 do i=0,wm
    w=(dble(i)+0.5d0)*dw
    ker(:,i)=(exp(-tau(:)*w)+exp(-(beta-tau(:))*w))/(1.d0+exp(-beta*w))
    ker(:,i)=matmul(transpose(cov),ker(:,i))
 enddo
 deallocate(cov)

 end subroutine initkern
!-----------------------!

!----------------------!
 subroutine writeconf()
!----------------------!
 use systemdata; implicit none

 integer :: i

 open(10,file='conf',status='replace')
 write(10,*)nw
 do i=1-pf,nw
    write(10,'(i9,3f14.8)')ww(i),dble(ww(i))*dw,aw(i)
 enddo
 close(10)

 end subroutine writeconf
!------------------------!

!------------------------------!
 subroutine readsqt(beta,sq,dd)
!------------------------------!
 use systemdata; implicit none

 integer :: i,j,k
 real(8) :: beta,sq,dd,x

 open(10,file='t.in',status='old')
 read(10,*)beta,nt,nb,sq

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
 subroutine initfiles()
!----------------------!

 open(10,file='log.log',status='replace') 
 close(10)
 open(10,file='th.dat',status='replace') 
 close(10)
 open(10,file='a0.dat',status='replace') 
 close(10)
 open(10,file='sa.dat',status='replace') 
 close(10)
 open(10,file='a.dat',status='replace') 
 close(10)

 end subroutine initfiles
!------------------------!

!--------------------------!
 subroutine deallocateall()
!--------------------------!
 use systemdata; use averages

 deallocate(ww)
 deallocate(aw)
 deallocate(xt1)
 deallocate(xt2)
 deallocate(ker)
 deallocate(tau)
 deallocate(sqt)
 deallocate(sig)
 deallocate(aaw)

 end subroutine deallocateall
!----------------------------!

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

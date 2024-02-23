!--------------!
 module cordata
!--------------!

 integer :: nq,nt,nb,nbt
 integer, allocatable :: ntq(:)
 integer, allocatable :: qlb(:)
 real(8), allocatable :: tau(:),cor(:,:,:),samp(:,:)
 real(8), allocatable :: sqt(:,:),sig(:,:),sq(:)

 end module cordata
!------------------!

!=====================!
 program tcorrelations
!=========================================================!                                                                 
! Do statistics of time correlations stored in 'cor.dat'.
! Write averages and errors to:
! - sq.dat   : static structure factor for each q
! - tq.dat   : time correlations for all times, all q
! - qxxx.dat : normalized time correlations and inverse of
!              covariance matrix for q=xxx 
!----------------------------------------------------------! 
 use cordata; implicit none

 integer :: ns,qq,sk,rb
 real(8) :: beta

 open(10,file='tres.in',status='old')
 read(10,*)nq    !  Number of q-points in file
 read(10,*)beta  !  Beta 
 read(10,*)qq    !  if 0 do covariance for all q, if not use q only
 read(10,*)nb    !  use number of bins (0=all)
 read(10,*)rb    !  Rebinning factor
 read(10,*)sk    !  use every sk time point
 read(10,*)nbt   !  number of bootstrap samples
 close(10)

 open(10,file='tres.log',status='replace'); close(10)

 call initran(1)
 call readindata(sk,rb)
 call computemeans()
 call covariance(qq,beta)

 end program tcorrelations
!=========================!

!-----------------------!
 subroutine computemeans
!-----------------------!
 use cordata; implicit none

 integer :: i,j

 allocate(sq(nq))
 allocate(ntq(nq))
 allocate(sqt(0:nt,nq))
 allocate(sig(0:nt,nq))
 allocate(samp(0:nt,0:nbt))

 open(10,file='sq.dat',status='replace')
 open(20,file='tq.dat',status='replace')
 do j=1,nq
    ntq(j)=nt
    call bootstraps(j,0)
    do i=0,nt       
       sqt(i,j)=samp(i,0)
       sig(i,j)=sqrt(sum((samp(i,1:nbt)-samp(i,0))**2)/dble(nbt))
       if (ntq(j)==nt.and.sig(i,j)/sqt(i,j).gt.0.1d0) ntq(j)=i-1
    enddo
 enddo
 do j=1,nq
    sq(j)=sqt(0,j)
    write(10,'(i5,2f15.9)')qlb(j),sqt(0,j),sig(0,j)
    write(20,*)j
    do i=0,nt
       write(20,'(3f15.9)')tau(i),sqt(i,j),sig(i,j)
    enddo
 enddo     
 close(10)
 close(20)

 end subroutine computemeans
!---------------------------!

!------------------------------!
 subroutine covariance(qq,beta)
!------------------------------!
 use cordata; implicit none

 integer :: i,j,k,n,q,q1,q2,qq,i1,i2,i3,ie
 real(8) :: beta
 character(8) :: fname

 real(8), allocatable :: cov(:,:),ss0(:)

 allocate(cov(nt,nt))
 allocate(ss0(nt))

 fname='q000.dat'
 if (qq==0) then
    q1=1
    q2=nq
 else
    q1=qq
    q2=qq
 endif
 do q=q1,q2
    n=ntq(q)
    i3=qlb(q)/100
    i2=mod(qlb(q),100)/10
    i1=mod(qlb(q),10)
    fname(4:4)=achar(48+i1)
    fname(3:3)=achar(48+i2)
    fname(2:2)=achar(48+i3)
    call bootstraps(q,1)    
    do i=1,n       
       sqt(i,q)=samp(i,0)
       sig(i,q)=sqrt(sum((samp(i,1:nbt)-samp(i,0))**2)/dble(nbt))
    enddo
    call computecov(n,cov(1:n,1:n))
    call diasym(cov(1:n,1:n),ss0(1:n),n)
    open(10,file=fname,status='replace')
    write(10,'(f12.3,2i8,2f16.10)')beta,ntq(q),nbt,sq(q)
    do i=1,n
       write(10,'(f12.9,4f20.15)')tau(i),sqt(i,q),sig(i,q),sqrt(ss0(i)/nbt)
    enddo
    do j=1,n
       write(10,*)j
       do i=1,n
          write(10,*)cov(i,j)
       enddo
    enddo
 enddo

 end subroutine covariance
!-------------------------!

!----------------------------!
 subroutine computecov(n,cov)
!----------------------------!
 use cordata; implicit none

 integer :: i,j,n
 real(8) :: cov(n,n)

 do j=1,n
 do i=1,n
    cov(i,j)=sum((samp(i,1:nbt)-samp(i,0))*(samp(j,1:nbt)-samp(j,0)))
 enddo
 enddo  

 end subroutine computecov
!-------------------------!

!-----------------------------!
 subroutine bootstraps(q,norm)
!-------------------------------------------!
! Generates nb bootstrap samples for given q
! Normalize to 1 at tau=0 if nr==1
!-------------------------------------------!
 use cordata; implicit none

 integer :: i,j,r,q,norm

 real(8), external :: ran

 samp=0.d0
 do i=1,nb
    samp(:,0)=samp(:,0)+cor(:,i,q)
 enddo
 do j=1,nbt
 do i=1,nb
    r=int(ran()*dble(nb))+1
    samp(:,j)=samp(:,j)+cor(:,r,q)
 enddo
 enddo
 samp=samp/dble(nb)
 if (norm==1) samp=samp/sq(q)

 end subroutine bootstraps
!-------------------------!

!---------------------------!
 subroutine readindata(s,rb)
!---------------------------!
 use cordata; implicit none

 integer :: i,j,k,b,q,s,u,nu,rb,kk
 real(8) :: c,beta

 open(10,file='tgrid.dat',status='old')
 read(10,*)nt
 nu=-1
 do i=0,nt
    read(10,*)c
    if (mod(i,s)==0) then
       nu=nu+1
    endif
 enddo
 rewind(10)
 read(10,*)nt
 allocate(tau(0:nu))
 nu=-1
 do i=0,nt
    read(10,*)c
    if (mod(i,s)==0) then
       nu=nu+1
       tau(nu)=c
    endif
 enddo

 allocate(qlb(1:nq))
 open(10,file='cor.dat',status='old') 
 b=0
 do 
    do j=1,nq
       read(10,*,end=10)q
       qlb(j)=q
       do i=0,nt
          read(10,*)c
       enddo
    enddo
    b=b+1
 enddo
 10 rewind(10)
 if (nb==0) then
    nb=b
 else
    nb=min(b,nb)
 endif
 nb=(nb-mod(nb,rb))/rb
 allocate(cor(0:nu,nb,0:nq)); cor=0.d0
 do k=0,rb*nb-1
    kk=1+k/rb
    do j=1,nq
       read(10,*)q
       u=-1
       do i=0,nt
          read(10,*)c
          if (mod(i,s)==0) then
             u=u+1
             cor(u,kk,j)=cor(u,kk,j)+c
          endif
       enddo
    enddo
 enddo
 close(10)
 cor=cor/dble(rb)
 nt=nu

 open(10,file='tres.log',position='append')
 write(10,*)'Number of bins read      : ',b
 write(10,*)'Using number of bins     : ',nb*rb
 write(10,*)'Rebinned to              : ',nb
 write(10,*)'Number of times          : ',nt
 close(10)
 
 end subroutine readindata
!-------------------------!

!---------------------------!
 subroutine diasym(aa,eig,n)
!---------------------------!
 implicit none

 integer :: n,l,inf
 real(8) ::  aa(n,n),eig(n),work(n*(3+n/2))

 l=n*(3+n/2)
 call dsyev('V','U',n,aa,n,eig,work,l,inf)

 end subroutine diasym
!---------------------!

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

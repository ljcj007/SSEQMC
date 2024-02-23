!-------------------------!
!-a program to calclulta--!
!-criticalilty in checker-!
!----by nusen MA----------!
!-------------------------!
!-------------!
 module system
!-------------!
 save

 integer :: ll
 integer :: nn
 integer :: n2
 integer :: nh
 integer :: mm
 integer :: maxcut
 real(8) :: prob(2)
 real(8) :: beta

 integer, allocatable :: spin(:)
 integer, allocatable :: oper(:)
 integer, allocatable :: site(:,:) 
 integer, allocatable :: otyp(:)
 integer, allocatable :: disc(:)

 integer, allocatable :: frst(:)
 integer, allocatable :: last(:)
 integer, allocatable :: vrtx(:)

 end module system
!-----------------!

!------------!
 module mdata
!------------!
 save

 integer, allocatable :: cos0(:)
 real(8), allocatable :: cos1(:,:)
 real(8), allocatable :: sin1(:,:)
 real(8), allocatable :: bres(:)

 real(8) :: nmsr=0.d0
 real(8) :: enrg=0.d0
 real(8) :: usus=0.d0
 real(8) :: usus2=0.d0
 real(8) :: asus=0.d0
 real(8) :: str0=0.d0
 real(8) :: str1(2)=0.d0
 real(8) :: rhosx=0.d0
 real(8) :: rhosy=0.d0
 real(8) :: binder=0.d0
 real(8) :: amag2=0.d0
 real(8) :: amag4=0.d0
!------------------------!
!---calculate slope------!
!------------------------!

 real(8) :: nstr=0.d0
 real(8) :: mn4=0.d0
 real(8) :: mn2=0.d0
 real(8) :: wxn2=0.d0
 real(8) :: wyn2=0.d0
 real(8) :: zn2=0.d0

!-------results---------!
 real(8) :: rxslope
 real(8) :: ryslope
 real(8) :: xslope
 real(8) :: bslope


 end module mdata
!----------------!

!===================!
 program ssechecker
!===================!
 use system; implicit none

 integer :: i,j,n,b,c,ens,nbins,rbins
 integer :: Ne, Nm,order,mulg
 integer :: cx,cy
 real(8) :: gg

 open(10,file='read.in',status='old')
 read(10,*)nbins,Ne,Nm,maxcut !order=0 or order=1 equilibrium or not
 read(10,*)cx,cy
 read(10,*)order
 read(10,*)ll
 read(10,*)beta
 read(10,*)mulg
 close(10)

 call initran(1)
 call initarrays()
 call initconf(cx,cy)
 gg=dble(mulg)/10000.d0

 open(10,file='stop',status='replace')
 write(10,*)0
 close(10)
 call bprobabilities()
 if (order==0)then
    do i=1,Ne
        call mcsweep()
        call adjustcutoff()
    enddo
    open(unit=10,file="spinconf.dat",status='new')
    write(10,"(I5)")spin(:)
    close(10)
    open(unit=10,file="bondconf.dat",status='new')
    write(10,"(2I20)")mm,nh
    write(10,"(I10)")oper(:)
    close(10)
      do j=1,nbins
         call writelog()
          do i=1,Nm
             call mcsweep()
             call measure()
          enddo
        open(unit=10,file="spinconf.dat",status='replace')
        write(10,"(I5)")spin(:)
        close(10)
        open(unit=10,file="bondconf.dat",status='replace')
        write(10,"(2I20)")mm,nh
        write(10,"(I10)")oper(:)
        close(10)
        call writedata(ll,gg)
     enddo
  else
    open(unit=10,file="spinconf.dat",status='old')
    read(10,"(I5)")spin(:)
    close(10)
    open(unit=10,file="bondconf.dat",status='old')
    read(10,"(2I20)")mm,nh
    if (allocated(oper)) then
        deallocate(oper)
        deallocate(vrtx)
    endif
    allocate(oper(0:mm-1))
    allocate(vrtx(0:4*mm-1))
    read(10,"(I10)")oper(:)
    close(10)
    do j=1,nbins
        call writelog()
        do i=1,Nm
            call mcsweep()
            call measure()
        enddo
        open(unit=10,file="spinconf.dat",status='replace')
        write(10,"(I5)")spin(:)
        close(10)
        open(unit=10,file="bondconf.dat",status='replace')
        write(10,"(2I20)")mm,nh
        write(10,"(I10)")oper(:)
        close(10)
        call writedata(ll,gg)
    enddo
  endif
 call deallocateall()

 contains

   !---------------------------!
    subroutine bprobabilities()
   !---------------------------!

       prob(1)=0.5d0*beta*dble(n2)
       prob(2)=gg*prob(1)

    end subroutine bprobabilities
   !-----------------------------!

  !---------------------!
   subroutine writelog()
  !---------------------!

   open(10,file='log.txt',status='replace')
   write(10,*)' Disorder realization : ',c
   write(10,*)' Beta value #         : ',beta
   write(10,*)' Bin                  : ',b
   write(10,*)' String length        : ',nh
   close(10)

   end subroutine writelog
  !-----------------------!

  !----------------------!
   subroutine checkstop()
  !----------------------!
   integer :: s   

   open(10,file='stop',status='old')
   read(10,*)s
   close(10)
   if (s/=0) then
      open(10,file='log.txt',position='append')
      write(10,*)' Stopping by request'
      stop
   endif

   end subroutine checkstop
  !------------------------!

 end program ssechecker
!===========================!


!--------------------!
 subroutine double_beta
!--------------------!
 use system; implicit none
  integer :: i,n

   if (mm>maxcut/2)then
    print *,"can't double beta", beta, mm
    stop
  end if
  n=2*mm
  vrtx(0:mm-1)=oper(:)
  deallocate(oper)
  allocate(oper(0:n-1))
   oper(0:mm-1)=vrtx(0:mm-1)
  do i=mm,n-1
   oper(i)=oper(n-i-1)
  end do
  deallocate(vrtx)
  allocate(vrtx(0:4*n-1))
  mm=mm*2
  nh=nh*2 
  beta=beta*2.d0
end subroutine double_beta
!--------------------!

!--------------------!
 subroutine mcsweep()
!--------------------!
 use system; implicit none

 call diagonalupdate()
 call linkvertices()
 call loopupdate()
 end subroutine mcsweep
!----------------------!

!---------------------------!
 subroutine diagonalupdate()
!---------------------------!
 use system; implicit none

 integer :: i,b
 real(8),external :: ran

 do i=0,mm-1
    if (oper(i)==0) then       
       b=int(ran()*dble(n2))+1
       if (spin(site(1,b))/=spin(site(2,b))) then
          if (ran()*dble(mm-nh)<prob(otyp(b))) then
             oper(i)=2*b
             nh=nh+1 
          endif
       endif
    elseif (mod(oper(i),2)==0) then       
       b=oper(i)/2
       if (dble(mm-nh+1)>ran()*prob(otyp(b))) then
          oper(i)=0
          nh=nh-1
       endif
    else
       b=oper(i)/2
       spin(site(1,b))=-spin(site(1,b))
       spin(site(2,b))=-spin(site(2,b))     
    endif
 enddo

 end subroutine diagonalupdate
!-----------------------------!

!--------------------------!
 subroutine  linkvertices()
!--------------------------!
 use system; implicit none

 integer :: i,b,p,s1,s2,v0,v1,v2
 real(8), external :: ran

 last(:)=-1
 frst(:)=-1

 do p=0,mm-1
    v0=4*p
    b=oper(p)/2
    if (b/=0) then      
       s1=site(1,b)
       s2=site(2,b)
       v1=last(s1)
       v2=last(s2)
       if (v1/=-1) then
          vrtx(v1)=v0
          vrtx(v0)=v1
       else
          frst(s1)=v0
       endif
       if (v2/=-1) then
          vrtx(v2)=v0+1
          vrtx(v0+1)=v2
       else
          frst(s2)=v0+1
       endif
       last(s1)=v0+2
       last(s2)=v0+3       
    else
       vrtx(v0:v0+3)=-1
    endif
 enddo

 do i=1,nn
    if (frst(i)/=-1) then
       vrtx(frst(i))=last(i)
       vrtx(last(i))=frst(i)
    endif
 enddo

 end subroutine linkvertices
!---------------------------!

!-----------------------!
 subroutine loopupdate()
!-----------------------!
 use system; implicit none

 integer :: i,p,v0,v1,v2
 real(8), external :: ran

 do v0=0,4*mm-1,2       
    if (vrtx(v0)<0) cycle
    v1=v0
    if (ran()<0.5d0) then
       do 
          vrtx(v1)=-1
          v2=ieor(v1,1)          
          v1=vrtx(v2)          
          vrtx(v2)=-1
          if (v1==v0) exit
       enddo
    else
       do 
          p=v1/4
          oper(p)=ieor(oper(p),1)
          vrtx(v1)=-2
          v2=ieor(v1,1)
          v1=vrtx(v2)
          vrtx(v2)=-2
          if (v1==v0) exit
       enddo
    endif
 enddo

 do i=1,nn
    if (frst(i)/=-1) then       
       if (vrtx(frst(i))==-2) spin(i)=-spin(i)       
    elseif (ran()<0.5d0) then
       spin(i)=-spin(i)
    endif
 enddo

 end subroutine loopupdate
!-------------------------!

!--------------------!
 subroutine measure()
!--------------------!
 use system; use mdata; implicit none

 integer :: i,j,b,s1,s2,jj(0:1),strb
 real(8) :: am,xx0,ss0,ss4,cm1(2),sm1(2),ss1(2)

 xx0=0.d0
 ss0=0.d0
 ss1(:)=0.d0
 ss4=0.d0
 am=0.d0
 strb=0
 do i=1,nn
    am=am+dble(spin(i)*cos0(i))/2
 enddo
 do i=1,2
    cm1(i)=0.5d0*sum(dble(spin(:))*cos1(i,:))
    sm1(i)=0.5d0*sum(dble(spin(:))*sin1(i,:))
 enddo
 jj(:)=0
 do i=0,mm-1
    b=oper(i)/2
    if(otyp(b)==2) strb=strb+1
    if (mod(oper(i),2)==1) then
       b=oper(i)/2
       s1=site(1,b)
       s2=site(2,b)
       spin(s1)=-spin(s1)
       spin(s2)=-spin(s2)
       am=am+2.d0*dble(spin(s1)*cos0(s1))
       jj((b-1)/nn)=jj((b-1)/nn)+spin(s2)
       cm1(:)=cm1(:)+dble(spin(s1))*cos1(:,s1)+dble(spin(s2))*cos1(:,s2)
       sm1(:)=sm1(:)+dble(spin(s1))*sin1(:,s1)+dble(spin(s2))*sin1(:,s2)
    endif
    if (oper(i)/=0) then
       xx0=xx0+am
       ss0=ss0+am**2
       ss4=ss4+am**4
       ss1(:)=ss1(:)+cm1(:)**2+sm1(:)**2
    endif
 enddo

 nmsr=nmsr+1.d0
 enrg=enrg+dble(nh)
 nstr=nstr+dble(strb)
 usus=usus+0.25d0*dble(sum(spin(:))**2)
 zn2=zn2+0.25d0*dble(sum(spin(:))**2)*dble(strb)
 asus=asus+(xx0**2+ss0)/(dble(nh)*dble(nh+1))
 amag2=amag2+ss0/dble(nh)
 mn2=mn2+(ss0*strb)/dble(nh)
 amag4=amag4+ss4/dble(nh)
 mn4=mn4+(ss4*strb)/dble(nh)
 str0=str0+ss0/dble(nh)
 str1(:)=str1(:)+ss1(:)/dble(nh)
 rhosx=rhosx+dfloat(jj(0))**2
 rhosy=rhosy+dfloat(jj(1))**2
 wxn2=wxn2+(dfloat(jj(0))**2)*strb
 wyn2=wyn2+(dfloat(jj(1))**2)*strb

 end subroutine measure
!----------------------!

!-----------------------------------------!
 subroutine writedata(l,gg)
!-----------------------------------------!
 use system; use mdata; implicit none
    
 integer :: i,bin,nt(2),l
 real(8) :: gg,aa

 aa=1/(l*gg)
 enrg=enrg/nmsr
 str0=str0/nmsr
 str1=str1/nmsr
 asus=asus/nmsr
 usus=usus/nmsr
 rhosx=rhosx/nmsr
 rhosy=rhosy/nmsr
 amag2=amag2/nmsr
 amag4=amag4/nmsr
 
 nstr=nstr/nmsr
 zn2=zn2/nmsr
 mn2=mn2/nmsr
 mn4=mn4/nmsr
 wxn2=wxn2/nmsr
 wyn2=wyn2/nmsr

 open(10,file='binder.dat',status='unknown',position='append')
 write(10,'(f10.5,5f38.10)')gg,amag2,amag4,mn4,mn2,nstr
 close(10)
 open(10,file="stiff.dat",status='unknown',position='append')
 write(10,'(f10.5,5f30.10)')gg,wxn2,wyn2,rhosx,rhosy,nstr
 close(10)

 open(10,file="suscep.dat",status='unknown',position='append')
 write(10,'(f10.5,3f30.10)')gg,zn2,usus,nstr
 close(10)
 !rxslope=(aa*(wxn2-nstr*rhosx))/beta
 !ryslope=(aa*(wyn2-nstr*rhosy))/beta
 !xslope=aa*beta*(zn2-usus*nstr)
 !bslope=(mn4+nstr*amag4-(2*amag4*mn2)/amag2)/(gg*(amag2**2))


 nt(:)=0
 do i=1,n2
    nt(otyp(i))=nt(otyp(i))+1
 enddo
 enrg=enrg/beta
 enrg=enrg-(nt(1)+gg*nt(2))/4.d0

 enrg=enrg/dble(nn)
 str0=str0/dble(nn)
 str1=str1/dble(nn)
 asus=beta*asus/dble(nn)
 usus=beta*usus/dble(nn)
 rhosx=1.5d0*rhosx/(beta*nn)
 rhosy=1.5d0*rhosy/(beta*nn)



 
 bres(1)=enrg
 bres(2)=usus
 bres(3)=asus
 bres(4)=str0
 bres(5)=str1(1)
 bres(6)=str1(2)
 bres(7)=rhosx
 bres(8)=rhosy
 bres(9)=0.5d0*(rhosx+rhosy)





 open(10,file='res.dat',status='unknown',position='append')
 write(10,'(f10.5,9f18.10)')gg,bres(:)
 close(10)

 !open(10,file='slope.dat',status='unknown',position='append')
 !write(10,'(f10.5,4f18.10)')gg,rxslope,ryslope,xslope,bslope
 !close(10)


nmsr=0.d0
enrg=0.d0
usus=0.d0
asus=0.d0
str0=0.d0
str1=0.d0
rhosx=0.d0
rhosy=0.d0
amag2=0.d0
amag4=0.d0
binder=0.d0

mn2=0.d0
mn4=0.d0
nstr=0.d0
zn2=0.d0
wxn2=0.d0
wyn2=0.d0
!xslope=0.d0
!bslope=0.d0
!rxslope=0.d0
!ryslope=0.d0

 end subroutine writedata
!------------------------!

!-----------------------------!
 subroutine average(nr,dat,av)
!-----------------------------!
 implicit none

 integer :: nr
 real(8) :: dat(nr),av(2)

 av(1)=sum(dat(:))/dble(nr)
 av(2)=sum(dat(:)**2)/dble(nr)
 av(2)=sqrt(abs(av(2)-av(1)**2)/dble(nr-1))

 end subroutine average
!----------------------!

!-----------------------------------------!
 subroutine corrlength(ll,nr,dat0,dat1,av)
!-----------------------------------------!
 implicit none

 real(8), parameter :: pi=3.141592653589793d0

 integer :: i,j,k,nr,nb,ll
 real(8) :: dat0(nr),dat1(nr),a0,a1,av(2),nc
 real(8),external :: ran

 nb=nr*10
 a0=sum(dat0(:))/dble(nr)
 a1=sum(dat1(:))/dble(nr)
 av(1)=sqrt((a0/a1-1.d0))/(2.d0*pi/dble(ll))
 nc=av(1)
 av(1)=av(1)/dble(ll)
 av(2)=0.d0
 do i=1,nb
    a0=0.d0
    a1=0.d0
    do j=1,nr
       k=int(ran()*nr)+1
       a0=a0+dat0(k)
       a1=a1+dat1(k)
    enddo
    a0=a0/dble(nr)   
    a1=a1/dble(nr)   
    a1=sqrt(a0/a1-1.d0)/(2.d0*pi/dble(ll))
    av(2)=av(2)+(nc-a1)**2
 enddo
 av(2)=av(2)/dble(nb)
 av(2)=sqrt(abs(av(2)))

 end subroutine corrlength
!-------------------------!

!-------------------------!
 subroutine adjustcutoff()
!-------------------------!
 use system; implicit none

 integer :: m1

 m1=nh+nh/3
 if (m1<=mm) then 
  return
 else if (m1>=maxcut) then
   print*, "m exceeded maximum cutoff"
   stop
 end if
 vrtx(0:mm-1)=oper(:)
 deallocate(oper)
 allocate(oper(0:m1-1))
 oper(0:mm-1)=vrtx(0:mm-1)
 oper(mm:m1-1)=0
 deallocate (vrtx)
 allocate(vrtx(0:4*m1-1))
 mm=m1

 end subroutine adjustcutoff
!---------------------------!

!---------------------------!
 subroutine initconf(cx,cy)
!---------------------------!
 use system; implicit none

 integer :: i,j,x,y,b1,b2,cx,cy
 real(8),external :: ran

 do i=1,nn
    spin(i)=(-1)**(mod(i-1,ll)+(i-1)/ll)
 enddo
 if (allocated(oper)) then
    deallocate(oper)
    deallocate(vrtx)
 endif    
 mm=20
 allocate(oper(0:mm-1))
 allocate(vrtx(0:4*mm-1))
 oper(:)=0
 nh=0

 otyp(:)=1

    do y=0,ll-1
    do x=0,(ll-1)/cx
       i=1+cx*x+y*ll
       otyp(i)=2
    enddo
    enddo
    do x=0,ll-1
    do y=0,(ll-1)/cy
      i=1+x+cy*y*ll
      j=i+nn
      otyp(j)=2
    enddo 
    enddo

 end subroutine initconf
!-----------------------!

!-------------------------------!
 subroutine initarrays()
!-------------------------------!
 use system; use mdata; implicit none

 real(8), parameter :: pi=3.141592653589793d0

 integer :: i,x1,x2,y1,y2,s1,s2,s3,ens

 nn=ll**2
 n2=2*nn

 allocate(site(2,n2))
 do y1=0,ll-1
 do x1=0,ll-1
    x2=mod(x1+1,ll)
    y2=mod(y1+1,ll)
    s1=1+x1+y1*ll
    s2=1+x2+y1*ll
    s3=1+x1+y2*ll
    site(1,s1)=s1
    site(2,s1)=s2
    site(1,s1+nn)=s1
    site(2,s1+nn)=s3
 enddo
 enddo

 allocate(cos0(nn))
 allocate(cos1(2,nn))
 allocate(sin1(2,nn))
 do y1=0,ll-1
 do x1=0,ll-1
    i=1+x1+y1*ll
    cos0(i)=(-1)**(x1+y1)
    cos1(1,i)=cos(dble(x1)*(pi-2.d0*pi/dble(ll)))*(-1.d0)**y1
    sin1(1,i)=sin(dble(x1)*(pi-2.d0*pi/dble(ll)))*(-1.d0)**y1
    cos1(2,i)=cos(dble(y1)*(pi-2.d0*pi/dble(ll)))*(-1.d0)**x1
    sin1(2,i)=sin(dble(y1)*(pi-2.d0*pi/dble(ll)))*(-1.d0)**x1
 enddo
 enddo

 allocate(otyp(n2))
 allocate(spin(nn))
 allocate(frst(nn))
 allocate(last(nn))
 allocate(bres(9))


 end subroutine initarrays
!-------------------------!

!--------------------------!
 subroutine deallocateall()
!--------------------------! 
 use system; use mdata; implicit none

 deallocate(spin)
 deallocate(oper)
 deallocate(site)
 deallocate(otyp)
 deallocate(frst)
 deallocate(last)
 deallocate(vrtx)
 deallocate(cos0)
 deallocate(cos1)
 deallocate(sin1)
 deallocate(bres)

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

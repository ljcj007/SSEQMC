!-------------------------------------------!
! Basic 2D S=1/2 Heisenberg SSE program     !
! Author: Anders Sandvik, Boston University !
! This version dated: 21.01.2009            !
!-------------------------------------------!
 module configuration
!--------------------!
 save

 integer :: lx
 integer :: ly
 integer :: nn
 integer :: nb
 integer :: nh
 integer :: mm

 real(8) :: beta
 real(8) :: aprob
 real(8) :: dprob
 real(8) :: j1

 integer, allocatable :: spin(:)
 integer, allocatable :: bsites(:,:) 
 real(8), allocatable :: jsites(:) 

 integer, allocatable :: opstring(:)

 integer, allocatable :: frstspinop(:)
 integer, allocatable :: lastspinop(:)
 integer, allocatable :: vertexlist(:)

 end module configuration
!------------------------!

!----------------------!
 module measurementdata
!----------------------!
 save

 real(8) :: enrg1=0.d0
 real(8) :: enrg2=0.d0
 real(8) :: amag1=0.d0 
 real(8) :: amag2=0.d0 
 real(8) :: asusc=0.d0 
 real(8) :: stiff=0.d0 
 real(8) :: ususc=0.d0
 real(8) :: data1(7)=0.d0
 real(8) :: data2(7)=0.d0

!--------------------------!
 end module measurementdata
!--------------------------!

!============================!
 program basic_heisenberg_sse
!============================!
 use configuration; implicit none

 integer :: i,j,nbins,msteps,isteps

 open(10,file='sse.in',status='old')
 read(10,*)lx,ly,beta,nbins,msteps,isteps
 read(10,*)j1
 close(10)

 call initran(1)
 call makelattice()
 call initconfig()

 aprob=0.5d0*beta*nb
 dprob=1.d0/(0.5d0*beta*nb)

 do i=1,isteps
    call diagonalupdate()
    call loopupdate()
    call adjustcutoff(i)
 enddo

 open(10,file='results.txt',status='replace')
 write(10,*)'Finished equilibration, M = ',mm
 close(10)

 do j=1,nbins
    do i=1,msteps
       call diagonalupdate()
       call loopupdate()
       call measureobservables()
    enddo
    call writeresults(msteps,j)
 enddo
   
 call deallocateall()

 end program basic_heisenberg_sse
!================================!

!---------------------------!
 subroutine diagonalupdate()
!---------------------------!
 use configuration; implicit none

 integer :: i,b,op
 real(8), external :: ran
 real(8) :: p

 do i=0,mm-1
    op=opstring(i)
    if (op==0) then       
       b=min(int(ran()*nb)+1,nb)
       if (spin(bsites(1,b))/=spin(bsites(2,b))) then
	  p=aprob/jsites(b)
          if (aprob>=dfloat(mm-nh).or.aprob>=ran()*(mm-nh)) then
             opstring(i)=2*b
             nh=nh+1 
          endif
       endif
    elseif (mod(op,2)==0) then        
       p=dprob*(mm-nh+1)*jsites(b)
       if (p>=1.d0.or.p>=ran()) then
          opstring(i)=0
          nh=nh-1
       endif
    else
       b=op/2
       spin(bsites(1,b))=-spin(bsites(1,b))
       spin(bsites(2,b))=-spin(bsites(2,b))
    endif
 enddo

 end subroutine diagonalupdate
!-----------------------------!

!-----------------------!
 subroutine loopupdate()
!-----------------------!
 use configuration; implicit none

 integer :: i,n,l,b,op,s1,s2,v0,v1,v2
 real(8), external :: ran
 frstspinop(:)=-1
 lastspinop(:)=-1

 do v0=0,4*mm-1,4
    op=opstring(v0/4)
    if (op/=0) then
       b=op/2
       s1=bsites(1,b)
       s2=bsites(2,b)
       v1=lastspinop(s1)
       v2=lastspinop(s2)
       if (v1/=-1) then
          vertexlist(v1)=v0
          vertexlist(v0)=v1
       else
          frstspinop(s1)=v0
       endif
       if (v2/=-1) then
          vertexlist(v2)=v0+1
          vertexlist(v0+1)=v2
       else
          frstspinop(s2)=v0+1
       endif
       lastspinop(s1)=v0+2
       lastspinop(s2)=v0+3
    else
       vertexlist(v0:v0+3)=0
    endif
 enddo
 do s1=1,nn
    v1=frstspinop(s1)
    if (v1/=-1) then
        v2=lastspinop(s1)
        vertexlist(v2)=v1
        vertexlist(v1)=v2
    endif
 enddo

 do v0=0,4*mm-1,2       
    if (vertexlist(v0)<1) cycle
    v1=v0
    if (ran()<0.5d0) then
       do 
          opstring(v1/4)=ieor(opstring(v1/4),1)
          vertexlist(v1)=-1
          v2=ieor(v1,1)
          v1=vertexlist(v2)
          vertexlist(v2)=-1
          if (v1==v0) exit
       enddo
    else
       do 
          vertexlist(v1)=0
          v2=ieor(v1,1)
          v1=vertexlist(v2)
          vertexlist(v2)=0
          if (v1==v0) exit
       enddo
    endif
 enddo

 do i=1,nn
    if (frstspinop(i)/=-1) then
       if (vertexlist(frstspinop(i))==-1) spin(i)=-spin(i)
    else
       if (ran()<0.5) spin(i)=-spin(i)
    endif
 enddo

 end subroutine loopupdate
!-------------------------!

!-------------------------------!
 subroutine measureobservables()
!-------------------------------!
  use configuration; use measurementdata; implicit none

 integer :: i,b,op,s1,s2,am,jj(0:1)
 real(8) :: am1,am2,ax1

 am=0
 do i=1,nn
    am=am+spin(i)*(-1)**(mod(i-1,lx)+(i-1)/lx)
 enddo      
 am=am/2
 am1=0.d0
 am2=0.d0
 ax1=0.d0
 jj(:)=0
 do i=0,mm-1
    op=opstring(i)
    if (op==0) then
        cycle
    elseif (mod(op,2)==1) then        
       b=op/2
       s1=bsites(1,b)
       s2=bsites(2,b)
       spin(s1)=-spin(s1)
       spin(s2)=-spin(s2)
       jj((b-1)/nn)=jj((b-1)/nn)+spin(s2)
       am=am+2*spin(s1)*(-1)**(mod(s1-1,lx)+(s1-1)/lx)
    endif
    ax1=ax1+dfloat(am)
    am1=am1+dfloat(abs(am))
    am2=am2+dfloat(am)**2
 enddo
 if (nh/=0) then
    ax1=(ax1**2+am2)/(dfloat(nh)*dfloat(nh+1))
    am1=am1/nh
    am2=am2/nh
 else
    am1=dfloat(abs(am))
    am2=dfloat(am)**2
    ax1=am2
 endif

 enrg1=enrg1+dfloat(nh)
 enrg2=enrg2+dfloat(nh)**2
 amag1=amag1+am1
 amag2=amag2+am2
 asusc=asusc+ax1
 stiff=stiff+0.5d0*(dfloat(jj(0))**2+dfloat(jj(1))**2)
 ususc=ususc+dfloat(sum(spin)/2)**2

 end subroutine measureobservables
!---------------------------------!

!------------------------------------!
 subroutine writeresults(msteps,bins)
!------------------------------------!
 use configuration; use measurementdata; implicit none

 integer :: i,msteps,bins
 real(8) :: wdata1(7),wdata2(7)

 enrg1=enrg1/msteps
 enrg2=enrg2/msteps
 amag1=amag1/msteps
 amag2=amag2/msteps
 asusc=asusc/msteps
 stiff=stiff/msteps
 ususc=ususc/msteps

 enrg2=(enrg2-enrg1*(enrg1+1.d0))/nn
 enrg1=enrg1/(beta*nn)-0.5d0
 amag1=amag1/nn
 amag2=amag2/nn
 asusc=beta*asusc/nn    
 ususc=beta*ususc/nn
 stiff=stiff/(beta*nn)

 data1(1)=data1(1)+enrg1
 data1(2)=data1(2)+enrg2
 data1(3)=data1(3)+amag1
 data1(4)=data1(4)+amag2
 data1(5)=data1(5)+asusc
 data1(6)=data1(6)+stiff
 data1(7)=data1(7)+ususc

 data2(1)=data2(1)+enrg1**2
 data2(2)=data2(2)+enrg2**2
 data2(3)=data2(3)+amag1**2
 data2(4)=data2(4)+amag2**2
 data2(5)=data2(5)+asusc**2
 data2(6)=data2(6)+stiff**2
 data2(7)=data2(7)+ususc**2

 do i=1,7
    wdata1(i)=data1(i)/bins
    wdata2(i)=data2(i)/bins
    wdata2(i)=sqrt(abs(wdata2(i)-wdata1(i)**2)/bins)
 enddo

 open(10,file='results.txt',status='replace')
 write(10,*)' Cut-off L : ',mm
 write(10,*)' Number of bins completed : ',bins
 write(10,*)' ========================================='
 write(10,10)' -E/N       : ',wdata1(1),wdata2(1)
 write(10,10)'  C/N       : ',wdata1(2),wdata2(2)
 write(10,10)'  <|m|>     : ',wdata1(3),wdata2(3)
 write(10,10)'  S(pi,pi)  : ',wdata1(4),wdata2(4)
 write(10,10)'  X(pi,pi)  : ',wdata1(5),wdata2(5)
 write(10,10)'  rho_s     : ',wdata1(6),wdata2(6)
 write(10,10)'  X(0,0)    : ',wdata1(7),wdata2(7)
 write(10,*)' ========================================='
 10 format(1x,a,2f14.8)
 close(10)

 enrg1=0.d0
 enrg2=0.d0
 amag1=0.d0
 amag2=0.d0
 asusc=0.d0
 stiff=0.d0
 ususc=0.d0

 end subroutine writeresults
!---------------------------!

!-----------------------------!
 subroutine adjustcutoff(step)
!-----------------------------!
 use configuration; implicit none

 integer, allocatable :: stringcopy(:)
 integer :: mmnew,step

 mmnew=nh+nh/3
 if (mmnew<=mm) return

 allocate(stringcopy(0:mm-1))
 stringcopy(:)=opstring(:)
 deallocate(opstring)
 allocate(opstring(0:mmnew-1))
 opstring(0:mm-1)=stringcopy(:)
 opstring(mm:mmnew-1)=0
 deallocate(stringcopy)

 mm=mmnew
 deallocate (vertexlist)
 allocate(vertexlist(0:4*mm-1))

 open(unit=10,file='results.txt',status='replace')
 write(10,*)' Step: ',step,'  Cut-off L: ',mm
 close(10)

 end subroutine adjustcutoff
!---------------------------!

!-----------------------!
 subroutine initconfig()
!-----------------------!
 use configuration; implicit none

 integer :: i
 real(8), external :: ran
 allocate(spin(nn))
 do i=1,nn
    spin(i)=2*int(2.*ran())-1
 enddo

 mm=20
 allocate(opstring(0:mm-1))
 opstring(:)=0
 nh=0

 allocate(frstspinop(nn))
 allocate(lastspinop(nn))
 allocate(vertexlist(0:4*mm-1))

 end subroutine initconfig
!-------------------------!
!------------------------!
 subroutine makelattice()
!------------------------!
 use configuration; implicit none

 integer :: s,x1,y1,s1,s2

 nn=lx*ly
 nb=8*nn

 allocate(bsites(2,nb))
 allocate(jsites(nb))

 do y1=0,ly-1
 do x1=0,lx-1
    s=1+x1+y1*lx
    bsites(1,4*s-3)=4*s-3
    bsites(2,4*s-3)=4*s-2
    bsites(1,4*s-2)=4*s-2
    bsites(2,4*s-2)=4*s-1
    bsites(1,4*s-1)=4*s-1
    bsites(2,4*s-1)=4*s
    bsites(1,4*s)=4*s
    bsites(2,4*s)=4*s-3
    jsites(4*s-3)=j1
    jsites(4*s-2)=j1
    jsites(4*s-1)=j1
    jsites(4*s)=j1
    jsites(s+4*nn)=1d0
    jsites(s+5*nn)=1d0
    jsites(s+6*nn)=1d0
    jsites(s+7*nn)=1d0
    s1=1+mod(x1+1,lx)+y1*lx
    s2=1+mod(x1-1+lx,lx)+y1*lx
    bsites(1,s+4*nn)=4*s-2
    bsites(2,s+4*nn)=4*s1
    bsites(1,s+5*nn)=4*s
    bsites(2,s+5*nn)=4*s2-2
    s1=1+x1+mod(y1+1,ly)*lx
    s1=1+x1+mod(y1+ly-1,ly)*lx
    bsites(1,s+6*nn)=4*s-1
    bsites(2,s+6*nn)=4*s1-3   
    bsites(1,s+7*nn)=4*s-3
    bsites(2,s+7*nn)=4*s2-1     
 enddo
 enddo
 nn=nn*4
 end subroutine makelattice
!--------------------------!

!--------------------------!
 subroutine deallocateall()
!--------------------------! 
 use configuration; implicit none

 deallocate (spin)
 deallocate (bsites)
 deallocate (opstring)
 deallocate (frstspinop)
 deallocate (lastspinop)
 deallocate (vertexlist)

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

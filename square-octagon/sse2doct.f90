!--------------------!
 module configuration
!--------------------!
 save

 integer :: lx       ! system length in x direction
 integer :: ly       ! system length in y direction
 integer :: nn       ! number of sites,nn=lx*ly
 integer :: nn1       ! number of sites,nn=lx*ly
 integer :: nb       ! number of bonds ,nb=2*nn
 integer :: ns       ! number of substrings (time slices)
 integer :: ms       ! maximum sub-string length
 integer :: mm       ! total string length (mm=ns*ms)

 real(8) :: beta     ! inverse temperature
 real(8) :: prob(2)

 integer, allocatable :: nh(:)        ! number of operators in each substring
 integer, allocatable :: spin(:)      ! spin state
 integer, allocatable :: bsites(:,:)  ! list of sites bsites(1,b),bsites(2,b) at bond b
 integer, allocatable :: opstring(:)  ! operator string
 integer, allocatable :: otyp(:)

 integer, allocatable :: vertexlist(:) ! list of vertex links
 integer, allocatable :: frstspinop(:) ! location in vertexlist() of first operation on a site 
 integer, allocatable :: lastspinop(:) ! location in vertexlist() of last operation on a site 

 end module configuration
!------------------------!

!----------------------!
 module measurementdata
!----------------------!
 save

 integer :: nqx        
 integer :: nqy          
 integer :: nq           ! number of q-correlations to measure
 integer :: nq1           ! number of q-correlations to measure
 integer :: ntau         ! number of time separations to measure
 integer :: tstp         ! time-slice step size in time averaging

 integer :: nms1=0       ! number of static measurements
 integer :: nms2=0       ! number of time-correlation measurements
 real(8) :: energ=0.d0   ! energy
 real(8) :: amag2=0.d0   ! squared sublattice magnetization
 real(8) :: ususc=0.d0   ! uniform susceptibility
 real(8),allocatable::crr(:)

 integer, allocatable :: qptsx(:)    ! list of q-points in x direction 
 integer, allocatable :: qptsy(:)    ! list of q-points in y direction
 integer, allocatable :: qpts(:)    ! list of q-points 
 integer, allocatable :: tgrd(:)    ! grid of time points
 real(8), allocatable :: tcor(:,:)  ! momentum-time correlation functions
 real(8), allocatable :: ref(:,:)   ! real part of fourier transform of spins
 real(8), allocatable :: imf(:,:)   ! imaginary part of fourier transform of spins
 real(8), allocatable :: phi(:,:,:) ! phase factors for fourier transforms
 real(8), allocatable :: tcor1(:,:)  ! momentum-time correlation functions
 real(8), allocatable :: ref1(:,:)   ! real part of fourier transform of spins
 real(8), allocatable :: imf1(:,:)   ! imaginary part of fourier transform of spins
 real(8), allocatable :: phi1(:,:,:) ! phase factors for fourier transforms
 real(8), allocatable :: tc(:)      ! working array for time correlations

!--------------------------!
 end module measurementdata
!--------------------------!

!==============!
 program hchecker
!==============!
 use configuration; use measurementdata; implicit none

 integer :: i,j,k,nbin,mstps,istps,tmsr,gtype
 real(8) :: gg
 real(8) :: dtau,tmax

 open(10,file='sse.in',status='old')
 read(10,*)lx,ly,beta,dtau           
 read(10,*)nbin,mstps,istps     
 read(10,*)gtype,tmax              
 read(10,*)tstp,tmsr 
 read(10,*)gg          
 close(10)

 nq=lx*ly

 allocate(qpts(nq))
 do k=1,nq
  qpts(k)=k
 enddo
! lx,ly  =  number of spins in x and y direction
! beta   =  inverse temperature
! dtau   =  time-slize widt
! nbin   =  number of bins
! mstps  =  steps per bin for measurements
! istps  =  initial (equlibration) steps
! gtype  =  time-grid type (1,2,3)
! tmax   =  maximum delta-time for correlation functions
! nqx,nqy=  number of q-points in x and y direction for time correlations
! tstp   =  time-step when computing time averages of correlations  
! tmsr   =  measure time correlations after every tmsr MC sweep
! qpts   =  q-values to do (among values 0,....,lx*ly/4, actual q is this times 2*pi/lx/ly)
 call initran(1)
 call makelattice()
 call initconfig(dtau)
 call taugrid(dtau,gtype,tmax)
 call qarrays()

  prob(1)=0.5d0*dtau*dble(nb)
  prob(2)=gg*prob(1)

 do i=1,istps
    call diagonalupdate()
    call linkvertices()
    call loopupdate()
    call adjustcutoff(i)
 enddo

 do j=1,nbin
    do i=1,mstps
       call diagonalupdate()
       call linkvertices()
       call loopupdate()
       call measure1()
       if (mod(i,tmsr)==0) call measure2()
    enddo
    call writedata()
 enddo
   
 call deallocateall()

 end program hchecker
!==================!

!---------------------------!
 subroutine diagonalupdate()
!---------------------------!
 use configuration; implicit none

 integer :: i,j,b,s,op

 real(8), external :: ran

 i=0
 do s=1,ns
 do j=1,ms
    op=opstring(i)
    if (op==0) then       
       b=int(ran()*nb)+1
       if (spin(bsites(1,b))/=spin(bsites(2,b))) then
          if (ran()*(ms-nh(s))<=prob(otyp(b))) then
             opstring(i)=2*b
             nh(s)=nh(s)+1 
          endif
       endif
    elseif (mod(op,2)==0) then        
       if (ran()*prob(otyp(b))<=(ms-nh(s)+1)) then
          opstring(i)=0
          nh(s)=nh(s)-1
       endif
    else
       b=op/2
       spin(bsites(1,b))=-spin(bsites(1,b))
       spin(bsites(2,b))=-spin(bsites(2,b))
    endif
    i=i+1
 enddo
 enddo

 end subroutine diagonalupdate
!-----------------------------!

!-------------------------!
 subroutine linkvertices()
!-------------------------!
 use configuration; implicit none

 integer :: i,j,k,b,op,s1,s2,v0,v1,v2

 frstspinop(:)=-1
 lastspinop(:)=-1

 do i=0,mm-1
    op=opstring(i)
    v0=4*i
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
       vertexlist(v0:v0+3)=-1
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
 
 end subroutine linkvertices
!---------------------------!

!-----------------------!
 subroutine loopupdate()
!-----------------------!
 use configuration; implicit none

 integer :: i,v0,v1,v2,nl
 real(8), external :: ran

 do v0=0,4*mm-1,2       
    if (vertexlist(v0)<0) cycle
    v1=v0
    if (ran()<0.5d0) then
       do 
          opstring(v1/4)=ieor(opstring(v1/4),1)
          vertexlist(v1)=-2
          v2=ieor(v1,1)
          v1=vertexlist(v2)
          vertexlist(v2)=-2
          if (v1==v0) exit
       enddo
    else
       do 
          vertexlist(v1)=-1
          v2=ieor(v1,1)
          v1=vertexlist(v2)
          vertexlist(v2)=-1
          if (v1==v0) exit
       enddo
    endif
 enddo

 do i=1,nn
    if (frstspinop(i)/=-1) then
       if (vertexlist(frstspinop(i))==-2) spin(i)=-spin(i)
    else
       if (ran()<0.5) spin(i)=-spin(i)
    endif
 enddo

 end subroutine loopupdate
!-------------------------!

!---------------------!
 subroutine measure1()
!-----------------------------!
! standard static measurements
!-----------------------------!
 use configuration; use measurementdata; implicit none

 integer :: i,b,op,s1,s2,am,x,y,r
 real(8) :: am2

 am=0
 do i=1,nn
    am=am+spin(i)*(-1)**(mod(i-1,lx)+(i-1)/lx)
 enddo      
 am=am/2
 am2=0.d0
 do i=0,mm-1
    op=opstring(i)
    if (mod(op,2)==1) then        
       b=op/2
       s1=bsites(1,b)
       s2=bsites(2,b)
       spin(s1)=-spin(s1)
       spin(s2)=-spin(s2)
       am=am+2*spin(s1)*(-1)**(mod(s1-1,lx)+(s1-1)/lx)
    endif
    am2=am2+dfloat(am)**2
 enddo
 am2=am2/dble(mm)

 do s1=1, nn ! nn: spin number
   x=mod(s1-1,lx) !x index
   y=int((s1-1)/lx) !y index
   do r=0, lx/2
     s2=mod(x+r,lx)+y*lx+1
     crr(r)=crr(r)+spin(s1)*spin(s2)
   end do
 end do

 energ=energ+dble(sum(nh))
 amag2=amag2+am2
 ususc=ususc+dble(sum(spin)/2)**2
 nms1=nms1+1
 crr=crr+0.25d0*crr

 end subroutine measure1
!-----------------------!

!---------------------!
 subroutine measure2()
!---------------------------------------------!
! time dependent q-space correlation functions
!---------------------------------------------!
 use configuration; use measurementdata; implicit none

 integer :: i,j,s,q,b,op,s1,s2,t1,t2,dt,dd
 integer, allocatable :: spinq(:,:) 
 allocate(spinq(4,nn1))

 i=0
 do s=1,ns
    do j=1,ms
       op=opstring(i)
       if (mod(op,2)==1) then        
          b=op/2
          s1=bsites(1,b)
          s2=bsites(2,b)
          spin(s1)=-spin(s1)
          spin(s2)=-spin(s2)
       endif
       i=i+1
    enddo
    
    do q=1,nn1
        spinq(1,q)=spin(4*q-3)
        spinq(2,q)=spin(4*q-2)
        spinq(3,q)=spin(4*q-1)
        spinq(4,q)=spin(4*q)
    enddo
    
    do q=1,nq
       ref(s,q)=sum(phi(1,:,q)*dble(spin(:)))  ! real part of the ft at slice s
       imf(s,q)=sum(phi(2,:,q)*dble(spin(:)))  ! imag part of the ft at slice s
    enddo
    
    do q=1,nq1
       ref1(s,q)=sum(phi1(1,:,q)*dble(spinq(1,:)))  
       imf1(s,q)=sum(phi1(2,:,q)*dble(spinq(1,:)))  
       ref1(s,q+nq1)=sum(phi1(1,:,q)*dble(spinq(2,:)))  
       imf1(s,q+nq1)=sum(phi1(2,:,q)*dble(spinq(2,:)))  
       ref1(s,q+nq2*2)=sum(phi1(1,:,q)*dble(spinq(3,:)))  
       imf1(s,q+nq2*2)=sum(phi1(2,:,q)*dble(spinq(3,:)))  
       ref1(s,q+nq3*2)=sum(phi1(1,:,q)*dble(spinq(4,:)))  
       imf1(s,q+nq3*2)=sum(phi1(2,:,q)*dble(spinq(4,:)))  
    enddo
 enddo
 do q=1,nq            ! loop over q points
    dd=0
    tc=0.d0
    do t1=1,ns,tstp   ! average time point t1 over time slices
       do dt=0,ntau
          t2=mod(t1+tgrd(dt)-1,ns)+1   ! add time difference in tgrd(dt) to get time t2
          tc(dt)=tc(dt)+ref(t1,q)*ref(t2,q)+imf(t1,q)*imf(t2,q)  ! the real part of the correlation
       enddo
       dd=dd+1        ! count the number of measurements for t-averages
    enddo
    tcor(:,q)=tcor(:,q)+0.75d0*tc(:)/dble(dd)  
 enddo
 do q=1,nq            ! loop over q points
    dd=0
    tc=0.d0
    do t1=1,ns,tstp   ! average time point t1 over time slices
       do dt=0,ntau
          t2=mod(t1+tgrd(dt)-1,ns)+1   
          tc(dt)=tc(dt)+ref1(t1,q)*ref1(t2,q)+imf1(t1,q)*imf1(t2,q)  
       enddo
       dd=dd+1        
    enddo
    tcor1(:,q)=tcor1(:,q)+0.75d0*tc(:)/dble(dd)  
 enddo
 
 nms2=nms2+1
 deallocate (spinq)

 end subroutine measure2
!-----------------------!

!----------------------!
 subroutine writedata()
!--------------------------------------------!
! writes the simple measurements to 'res.dat'
! writes time correlations to 'cor.dat'
!--------------------------------------------!
 use configuration; use measurementdata; implicit none

 integer :: i,q,t,r,nt(2)
 real(8) :: gg

 energ=energ/dble(nms1)
 amag2=amag2/dble(nms1)
 ususc=ususc/dble(nms1)
 crr=crr/dble(nms1)


 nt(:)=0
 do i=1,nb
    nt(otyp(i))=nt(otyp(i))+1
 enddo
 energ=energ/beta
 energ=energ-(nt(1)+gg*nt(2))/4.d0

 energ=energ/dble(nn)
 amag2=3.d0*amag2/dble(nn)**2
 ususc=beta*ususc/dble(nn)
 crr=crr/dble(nn)

 open(10,file='res.dat',position='append')
 write(10,'(3f15.10)')energ,amag2,ususc,r,crr(r)
 close(10)
 energ=0.d0
 amag2=0.d0
 ususc=0.d0
 nms1=0
 crr(:)=0.d0

 tcor=tcor/dble(nms2) 
 open(10,file='cor.dat',position='append')
 do q=1,nq
    write(10,'(2i8)')q,qpts(q)
    do t=0,ntau
       write(10,'(f20.12)')tcor(t,q)
    enddo
 enddo
 close(10)
 tcor=0.d0
 
 tcor1=tcor1/dble(nms2) 
 open(10,file='cor1.dat',position='append')
 do q=1,nq
    write(10,'(2i8)')q,qpts(q)
    do t=0,ntau
       write(10,'(f20.12)')tcor1(t,q)
    enddo
 enddo
 close(10)
 tcor1=0.d0
 
 nms2=0

 end subroutine writedata
!------------------------!

!-----------------------------!
 subroutine adjustcutoff(step)
!-----------------------------!
 use configuration; implicit none

 integer, allocatable :: stringcopy(:)
 integer :: i,nm1,nm2,nm3,step

 nm1=maxval(nh)
 nm2=nm1+nm1/5+1
 nm1=sum(nh)
 nm3=nm1+nm1/2
 if (nm2<=ms.and.nm3<=mm) return
 if (nm2*ns<nm3) nm2=nm3/ns+1

 allocate(stringcopy(0:mm-1))
 stringcopy(:)=opstring(:)
 deallocate(opstring)
 allocate(opstring(0:ns*nm2-1)); opstring=0
 do i=1,ns
    opstring((i-1)*nm2:(i-1)*nm2+ms-1)=stringcopy((i-1)*ms:i*ms-1)
 enddo    
 deallocate(stringcopy)

 ms=nm2
 mm=ns*ms
 deallocate(vertexlist)
 allocate(vertexlist(0:4*mm))

 open(unit=10,file='info.txt',position='append')
 write(10,*)' Step: ',step,'  Cut-off L: ',mm,'  Current nh: ',sum(nh)
 close(10)

 end subroutine adjustcutoff
!---------------------------!

!---------------------------!
 subroutine initconfig(dtau)
!---------------------------!
 use configuration;use measurementdata; implicit none

 integer :: i,j,x,y
 real(8) :: ran,dtau

 external :: ran

 allocate(spin(nn))
 do i=1,nn
    spin(i)=(-1)**(mod(i-1,lx)+(i-1)/lx)
 enddo

 ms=4
 ns=int(beta/dtau+0.1d0)
 mm=ms*ns
 allocate(nh(ns)); nh=0
 allocate(opstring(0:mm-1)); opstring(:)=0
 
 allocate(frstspinop(nn))
 allocate(lastspinop(nn))
 allocate(vertexlist(0:4*mm))
 allocate(crr(0:nn/2))


 end subroutine initconfig
!-------------------------!

!-----------------------------------!
 subroutine taugrid(dtau,gtype,tmax)
!-------------------------------------------------------------------------!
! constructs array tgrd(0,...,ntau) of time separations (in units of dtau)
! gtype = 1 for uniform grid of all separatiins up to tmax
!         2 for quadratic grid t=dtau*n^2/4 up to beta/2
!         3 for uniform grid up to tmax followed by quadratic to beta/2
! actual time points including dtau factor are written to 'tgrid.dat'
!--------------------------------------------------------------------------!
 use configuration; use measurementdata; implicit none

 integer :: i,t1,t2,gtype,nsm
 real(8) :: dtau,tmax

 nsm=int((tmax+1.d-5)/dtau)
 open(unit=10,file='tgrid.dat',status='replace')
 if (gtype==1) then
    ntau=nsm
    write(10,'(i4)')ntau
    allocate(tgrd(0:ntau)) 
    do i=0,ntau
       tgrd(i)=i
       write(10,'(f15.8)')tgrd(i)*dtau
    enddo
 elseif (gtype==2) then    
    ntau=0
    t1=0
    do
       t2=((ntau+1)**2)/4
       if (t2==t1) t2=t1+1
       if (t2.le.nsm) then
          ntau=ntau+1
          t1=t2
       else
          exit
       endif
    enddo
    write(10,'(i4)')ntau
    allocate(tgrd(0:ntau)) 
    tgrd(0)=0
    write(10,'(f15.8)')0.d0
    do i=1,ntau
       tgrd(i)=(i**2)/4
       if (tgrd(i)==tgrd(i-1)) tgrd(i)=tgrd(i-1)+1       
       write(10,'(f15.8)')tgrd(i)*dtau
    enddo
 elseif (gtype==3) then    
    ntau=nsm
    t1=0
    i=0
    do
       t2=((i+1)**2)/4
       if (t2==t1) t2=t1+1
       if (t2.gt.nsm.and.t2.le.ns/2) then
          ntau=ntau+1
          t1=t2
       elseif (t2.gt.ns/2) then
          exit
       endif
       i=i+1
    enddo
    write(10,'(i4)')ntau
    allocate(tgrd(0:ntau)) 
    tgrd(0)=0
    do i=0,nsm
       tgrd(i)=i
       write(10,'(f15.8)')tgrd(i)*dtau
    enddo
    ntau=nsm
    t1=0
    i=0
    do
       t2=((i+1)**2)/4
       if (t2==t1) t2=t1+1
       if (t2.gt.nsm.and.t2.le.ns/2) then
          ntau=ntau+1
          tgrd(ntau)=t2
          write(10,'(f15.8)')tgrd(ntau)*dtau
          t1=t2
       elseif (t2.gt.ns/2) then
          exit
       endif
       i=i+1
    enddo
 endif
 close(10) 
 allocate(tc(0:ntau)) 
 allocate(tcor(0:ntau,nq)) 
 tcor=0.d0
 allocate(tcor1(0:ntau,nq)) 
 tcor1=0.d0

 end subroutine taugrid
!----------------------!

!--------------------! 
 subroutine qarrays()
!-----------------------------------------------------------!
! phase factors for fourier transforms of spin configuration
!-----------------------------------------------------------!
 use configuration; use measurementdata; implicit none

 real(8), parameter :: pi=3.14159265358979323d0

 integer :: q1,q2,rx,ry,r,q,ix,iy
 real(8) :: qx,qy
 nqx=lx/2
 nqy=ly/2
 nq1=nqx*nqy
 
 allocate(phi(2,nn,nq))
 allocate(ref(ns,nq))
 allocate(imf(ns,nq))

 do q=1,nq
    qx=2.d0*pi*dble(qpts(mod(q-1,lx)+1))/dble(2*lx)
    qy=2.d0*pi*dble(qpts((q-1)/lx+1)))/dble(2*ly)
    do r=1,nn1
       rx=mod(r-1,lx)+1
       ry=(r-1)/lx+1
       phi(1,4*r-3,q)=cos(dble(2*rx)*qx+dble(2*ry+0.5)*qy)/sqrt(dble(nn))
       phi(2,4*r-3,q)=sin(dble(2*rx)*qx+dble(2*ry+0.5)*qy)/sqrt(dble(nn))
       phi(1,4*r-2,q)=cos(dble(2*rx+0.5)*qx+dble(2*ry)*qy)/sqrt(dble(nn))
       phi(2,4*r-2,q)=sin(dble(2*rx+0.5)*qx+dble(2*ry)*qy)/sqrt(dble(nn))
       phi(1,4*r-1,q)=cos(dble(2*rx)*qx+dble(2*ry-0.5)*qy)/sqrt(dble(nn))
       phi(2,4*r-1,q)=sin(dble(2*rx)*qx+dble(2*ry-0.5)*qy)/sqrt(dble(nn))
       phi(1,4*r,q)=cos(dble(2*rx-0.5)*qx+dble(2*ry)*qy)/sqrt(dble(nn))
       phi(2,4*r,q)=sin(dble(2*rx-0.5)*qx+dble(2*ry)*qy)/sqrt(dble(nn))
    enddo
 enddo

 allocate(phi1(2,nn1,nq1))
 allocate(ref1(ns,nq1))
 allocate(imf1(ns,nq1))
 do q=1,nq1
    qx=2.d0*pi*dble(qpts(mod(q-1,nqx)+1))/dble(lx)
    qy=2.d0*pi*dble(qpts((q-1)/nqx+1)))/dble(ly)
    do r=1,nn1
       rx=mod(r-1,lx)+1
       ry=(r-1)/lx+1
       phi1(1,r,q)=cos(dble(rx)*qx+dble(ry)*qy)/sqrt(dble(nn1))       
       phi1(2,r,q)=sin(dble(rx)*qx+dble(ry)*qy)/sqrt(dble(nn1)) 
    enddo
 enddo

 end subroutine qarrays
!----------------------!

!------------------------!
 subroutine makelattice()
!------------------------!
 use configuration; implicit none

 integer :: x1,x2,y1,y2,s,s1,s2

 nn1=lx*ly
 nn=nn1
 nb=nn1*8
 allocate(bsites(2,nb))
 allocate(otyp(nb))
 
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
    otyp(4*s-3)=2
    otyp(4*s-2)=2
    otyp(4*s-1)=2
    otyp(4*s)=2
    otyp(s+4*nn)=1
    otyp(s+5*nn)=1
    otyp(s+6*nn)=1
    otyp(s+7*nn)=1
    s1=1+mod(x1+1,lx)+y1*lx
    s2=1+mod(x1-1+lx,lx)+y1*lx
    bsites(1,s+4*nn)=4*s-2
    bsites(2,s+4*nn)=4*s1
    bsites(1,s+5*nn)=4*s
    bsites(2,s+5*nn)=4*s2-2
    s1=1+x1+mod(y1+1,ly)*lx
    s2=1+x1+mod(y1+ly-1,ly)*lx
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
 use configuration; use measurementdata; implicit none

 deallocate (nh)
 deallocate (spin)
 deallocate (bsites)
 deallocate (opstring)
 deallocate (frstspinop)
 deallocate (lastspinop)
 deallocate (vertexlist)
 deallocate (qpts)
 deallocate (qptsx)
 deallocate (qptsy)
 deallocate (tcor)
 deallocate (phi)
 deallocate (ref)
 deallocate (imf)
 deallocate (tcor1)
 deallocate (phi1)
 deallocate (ref1)
 deallocate (imf1)
 deallocate (tc)
 deallocate (crr)
 deallocate (otyp)

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

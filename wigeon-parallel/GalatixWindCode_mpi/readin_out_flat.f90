module readin_out
  use datablock
  use hydro_stat

contains

!	 --------------------------------------------------------------------------
!     Name:      data_io.for
!     Function:  read in some parameters
!     ---------------------------------------------------------------------------

  subroutine setup

    rLbox=15.0   ! kpc
    radc=0.30    ! kpc for Mg=1.e8 Msun
    radd=6.53    !kpc   for Mg=1.e8 Msun
    rhoc=6.0e+7  ! Msun kpc^-3 for Mg=1.e+8 Msun
    cfl=0.32      
    tint=0.      ! Myr
!    tend=.5   ! MyrI  yeh
    tend=100.   ! Myr
    tsupend=50.0 ! Myr
    cs_init=8.0  ! km s^-1
    p_igm=1.0*k_boltzmann  ! cm^-3 K
    rho0=527.42
    ist=4000
    gamma=1.6666666667
    gm1=gamma-1.
    call Phys_Unit
    radc=radc/l_unit
    radd=radd/l_unit
    rhoc=rhoc/(rho_unit*1.48e+31)
    rho0=1.0
    tint=tint/t_unit
    tend=tend/t_unit
    tsupend=tsupend/t_unit
    cs_init=cs_init/v_unit
    g0=(1.e+6*365*86400)**2*g_gravity*rho_unit*t_unit*t_unit
    dedp=1.e+38*t_unit*1.e+6*365.*86400./ene_unit/(l_unit*3.09e+21)**3/405.0
    dmdp=gm1*dedp/1000.0/cs_init**2

! table the output time sequence 

    do i=1,n_out
       t_out(i)=i*tend/n_out
    enddo
    if(myid.eq.0) then
       open (10,file=datadir(1:lblnk(datadir))//'coolrate_1.0',status='old',form='unformatted')
       read(10) coolrate
       close(10)
    endif
    call mpi_bcast(coolrate,n_cool,mpi_double_precision,0,mpi_comm_world,ierr)

    call TVD_RK_parameter

    
    do iz=1,nLz 
       do iy=1,nLy
          do ix=1,nLx
             rix=real(ix+nxst)
             riy=real(iy+nyst)
             riz=real(iz+nzst)
             
             xx=sqrt((rix-real(Lx+1)/2.)**2+(riy-real(Ly+1)/2.)**2+(riz-real(Lz+1)/2.)**2)/radc
             if(xx==0.) then
                phi=4.0*pi*g0*rhoc/3.0
             else    
                phi=4.0*pi*g0*rhoc*(xx-atan(xx))/xx**3
             endif
             frc(1,ix,iy,iz)=phi*(rix-real(Lx+1)/2.)
             frc(2,ix,iy,iz)=phi*(riy-real(Ly+1)/2.)
             frc(3,ix,iy,iz)=phi*(riz-real(Lz+1)/2.)
          enddo
       enddo
    enddo
    open (69,file='frctest.dat')
!    write(69,*)frc(1,1:Lx,1:Ly,1:Lz)
    close(69)
! -----------------------------------------------------------------------------

  end subroutine setup

!     -----------------------------------------------------------------------
!     Name:      initialize.f
!     Function:  set up the grid and initial condition ( u(x,y,z) at t= 0 )
!     -----------------------------------------------------------------------

  subroutine initialize
    include'mpif.h'

    real,allocatable,dimension(:,:,:)::uct
    character (len=60):: datafile
    integer,dimension(3)::ndim
    
    allocate(uct(nLx,nLy,nLz),STAT=istat)

!    write(*,'(''enter the status of the initial condition'')')
!    i_run=integerread('initial at t=0. (0) / restart from a restart file (1) = ')
    i_run=0 
    if (i_run==0) then
       
       nt=0
       tnum = 0.
       i_out=1
       supin=.false.
       write(datafile,'(i4''.0000'')')ist 
       uc(:,:,:,:)=0.0
!  initialize the gas component
!$omp parallel do default(shared) private(i,j,k,rix,riy,riz,rrr,rrr2,zz,fw,fr,f0,ee,eks1,rho_init,vcrr)
       do k=1,nLz
          do j=1,nLy
             do i=1,nLx 
                rix=real(i+nxst)
                riy=real(j+nyst)
                riz=real(k+nzst)

                rrr=sqrt((rix-real(Lx+1)/2.)**2+(riy-real(Ly+1)/2.)**2)
                rrr2=sqrt(rrr**2+(riz-real(Lz+1)/2.)**2)
                zz=(riz-real(Lz+1)/2.)
                fw=feks(rrr2)
                fr=feks(rrr)
                f0=feks(0.0)
                ee=0.93/exp((zz*l_unit/2.5)**2)
                eks1=fw-ee**2*fr-(1.-ee**2)*f0
                rho_init=rho0*exp(-293.0*eks1/64.0)
                if((rho_init>=(1.0/527.42)).and.(rrr<=(radd-4.))) then
                   uc(1,i,j,k)=rho_init
                   vcrr=ee*sqrt(rrr*sqrt(frc(1,i,j,k)**2+frc(2,i,j,k)**2))
                   uc(2,i,j,k)=uc(1,i,j,k)*vcrr*real(-(riy-real(Ly+1)/2.))/rrr
                   uc(3,i,j,k)=uc(1,i,j,k)*vcrr*real(rix-real(Lx+1)/2.)/rrr
                else if(rrr>(radd-4.)) then
                   rho_init=rho_init*exp(-0.8*(rrr+4.-radd))*exp(-1.5*abs(zz)*l_unit/2.5)        
                   uc(1,i,j,k)=max(rho_init,3.0e-5)
                else
                   rho_init=rho_init*exp(-1.5*abs(zz)*l_unit/2.5)     
                   uc(1,i,j,k)=max(rho_init,3.0e-5)
                endif
                
                uc(5,i,j,k)=uc(1,i,j,k)*cs_init*cs_init/gamma/gm1+(uc(2,i,j,k)**2+uc(3,i,j,k)**2+uc(4,i,j,k)**2)/(uc(1,i,j,k)+1.e-8)/2.0
             enddo
          enddo
       enddo
       
!$omp end parallel do 

!   initialize points out of caculate region according to flat boundary
!   condition  

       if(myid.eq.0) call head_logfile(i_run)
!       if(myid.eq.0) call data_output(.true.)
!       call data_output(.true.)
       
    else if (i_run==1) then
       
!       open (2,file=datadir(1:lblnk(datadir))//'ucdump.dat',status='old',form='unformatted')
!       read (2) tnum
!       read (2) nt
!       read (2)  i_out
!       do ic=1,mnp
!          read (2) uct
!          uc(ic,1:Lx,1:Ly,1:nLz)=uct(1:Lx,1:Ly,1:nLz)
!       enddo
!       close(2)
!       if(tnum<tsupend) then
!          supin=.true.
!       else
!          supin=.false.
!       endif
       
!       do i=-2,0
!          uc(:,i,1:Ly,1:Lz)=uc(:,1,1:Ly,1:Lz)
!          uc(:,1:Lx,i,1:Lz)=uc(:,1:Lx,1,1:Lz)
!          uc(:,1:Lx,1:Ly,i)=uc(:,1:Lx,1:Ly,1)
!       enddo
       
!       do i=Ls+1,Ls+3
!          uc(:,i,1:Ly,1:Lz)=uc(:,Lx,1:Ly,1:Lz)
!          uc(:,1:Lx,i,1:Lz)=uc(:,1:Lx,Ly,1:Lz)
!          uc(:,1:Lx,1:Ly,i)=uc(:,1:Lx,1:Ly,Lz)
!       enddo
!    if(myid.eq.0) call head_logfile(i_run)

    endif

    deallocate(uct,STAT=istat)

  end subroutine initialize
  
  subroutine head_logfile(i_run)
    character(len=60) datafile
    
    write(datafile,'(''prun''i4''.log'')')ist
    
    select case (i_run) 
    case(0)
       open (68,file=datadir(1:lblnk(datadir))//datafile(1:lblnk(datafile)),status='unknown')
    case(1)
       open (68,file=datadir(1:lblnk(datadir))//datafile(1:lblnk(datafile)),&
&       status='old',position='append')
       write(68,'(''----------------------------- RESTART ----------------------------------'')')
    end select
    close(68)
  end subroutine head_logfile

! -----------------------------------------------------------------------------

  subroutine data_output(write_disk)
    character fname*6,ucout*20,sts*8
    logical write_disk    
!     -------------------------------------------------------------------------
!     Name:      data_output.f
!     Function:  save the solution
!     -------------------------------------------------------------------------
    if (write_disk) then 
       write(ucout,'(''uc'',i4.4,''.'',i4.4,''.'',i6.6)')ist,i_out-1,myid
       sts='unknown'
       if(myid.eq.0) write(*,'(''save the numerical solution'',''step no='',i4)')nt
    else     
       write(fname,'(i6.6)')myid
       ucout='ucdump.'//fname//'.dat'
       sts='unknown'
    endif
    
    open(10,file=datadir(1:lblnk(datadir))//ucout,form='unformatted')
!    if(myid.eq.0) then
       write(10) tnum
       write(*,*) 'in output tnum=',tnum
       write(10) nt
       write(10) i_out
       write(*,*) "write_log time",tnum,nt,i_out
!    endif

    do ic=1,mnp
       write(10) uc(ic,1:nLx,1:nLy,1:nLz)
    enddo

    close(10)
    return
  end subroutine data_output

! --------------------- write log file and output on screen ------------------------

  subroutine write_log
    
    character(len=60) datafile
    
    write(datafile,'(''prun''i4''.log'')')ist
    Ns_left = int((tend-tnum)/dt)
        
    vx_max = vx_max*v_unit 
    vy_max = vy_max *v_unit
    vz_max = vz_max *v_unit
    cs_max = cs_max *v_unit
    vel_sigma = vel_sigma *v_unit
    
    write(*,'('' step no='',i4,1x&
& ''time = '',f10.8)') nt, tnum
    write(*,'('' cfl ='',f10.8,1x)') dt_cfl
    write(*,'(''left time step No. = '', I8)')Ns_left
    
    write(*,'(''avarage density = '', f15.11,1x,''negative cells # = '', i8)') den_av,ncell_zero
    write(*,'(''minimum/maximum density = '',f13.8,1x,f13.6)')den_min,den_max
    write(*,'(''rms density ='',f12.8)') den_sigma
    write(*,'(''maximum sound velocity = '', f13.6)') cs_max
    write(*,'(''maximum 3d gas velocity = '', 3f13.6)') vx_max,vy_max,vz_max
    write(*,'(''velocity dispersion gas '', f12.8)') vel_sigma
    write(*,'('' '')')
    write(*,'(''---------------------------------------------------------'')')
    write(*,'('' '')')

    open (68,file=datadir(1:lblnk(datadir))//datafile(1:lblnk(datafile)),&
&       status='old',position='append')
    write(68,'('' step no='',i4,1x&
& ''time = '',f10.6)') nt, tnum
    write(68,'('' cfl ='',f10.8,1x,''/ div-vel ='',f8.5,1x)') dt_cfl
    write(68,'(''left time step No. = '', I8)')Ns_left
    
    write(68,'(''avarage density = '', f15.12,1x,''negative cells # = '', i8)') den_av,ncell_zero
    write(68,'(''minimum/maximum density = '',f13.8,1x,f13.6)')den_min,den_max
    write(68,'(''rms density ='',f12.8)') den_sigma
    write(68,'(''maximum sound velocity = '', f13.6)') cs_max
    write(68,'(''maximum 3d gas velocity = '', 3f13.6)') vx_max,vy_max,vz_max
    write(68,'(''velocity dispersion gas '', f12.8)') vel_sigma
    write(68,'(''--------------------------------------------------------'')')
    close(68)
  end subroutine write_log

  integer function lblnk(char)
    character char*(*)
    lblnk=index(char,' ')-1
    return
  end function lblnk

  integer FUNCTION integerread(text)
    Character text*(*)
    integer ix
    write (*,'(A,$)')text
    read (*,*) ix
    integerread=ix
    Return
  End FUNCTION integerread
  
  subroutine CPU_screentime(t_init,t_step,text)
    include 'mpif.h'
    character :: text*(*)
    real :: t_init,t_run,t_step
    
    t_run=mpi_wtime()
    t_step=t_run-t_init
    write(*,'("on",i4,i4,A,''='',f15.6)')myid,Lx,text,t_step
    t_init=t_run
    
    return
  end subroutine CPU_screentime
  
  function feks(rr)
    implicit none
    real:: rr,feks
    if(rr==0) then
       feks=1.0
    else
       feks=log(1+(rr/radc)**2)/2.0+atan(rr/radc)/(rr/radc)
    endif
    
    return
  end function feks


  SUBROUTINE polint(xa,ya,n,x,y,dy)
! polynomial interpolation
    integer n,nmax
    real dy,x,y,xa(n),ya(n)
    parameter (nmax=10)
    integer i,m,ns
    real den,dif,dift,ho,hp,w,c(nmax),d(nmax)
    ns=1
    dif=abs(x-xa(1))
    do i=1,n
       dift=abs(x-xa(i))
       if (dift.lt.dif) then
          ns=i
          dif=dift
       endif
       c(i)=ya(i)
       d(i)=ya(i)
    enddo
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if (den.eq.0.) PAUSE 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       enddo
       if (2*ns.lt.n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       endif
       y=y+dy
    enddo
    RETURN
  end SUBROUTINE polint
  
end module readin_out

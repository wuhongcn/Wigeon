program main

  use datablock
  use weno_kernel
  use readin_out
  use hydro_stat
  use IGM_heating_cooling

  logical output,stop_sign 
  integer ierr,type0
  integer surface,i
  character ucout*20
  real time1,time2
  real init_time,cpu_dt
  character arg_name*80,arg_value*80

!     -------------------------------------------------------------------
!     Name:      weno-ustctw.f
!     Function:  drive routine
!      System to solve: u_t + (f(u)_x + g(u)_y +h(u)_z) / a(t) = forc(u,t)
!     or
!      u_t = RHS = (-f(u)_x - g(u)_y - h(u)_z) / a(t) + forc(u,t)
!	 ---------------------------------------------------------------

! ****** readin parameters ******
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
  time1=mpi_wtime()
  if(myid.eq.0) then
     print*,"There are",nprocs,"processors running this job."
  end if

! read command line parameters
  if (iargc().lt.4) then
     if (myid.eq.0) write(*,*)'Usage: sup_flat -Lx 256 -steps 15'
     call MPI_FINALIZE(ierr)
     stop
  endif
  do i=1,iargc(),2
    call getarg(i,arg_name)
    call getarg(i+1,arg_value)
    if (arg_name.eq.'-lx' .or. arg_name.eq.'-Lx' .or. arg_name.eq.'-LX') then
      read(arg_value,*) Lx
    else if (arg_name.eq.'-ly' .or. arg_name.eq.'-Ly' .or. arg_name.eq.'-LY') then
      read(arg_value,*) Ly
    else if (arg_name.eq.'-lz' .or. arg_name.eq.'-Lz' .or. arg_name.eq.'-LZ') then
      read(arg_value,*) Lz
    else if (arg_name.eq.'-steps') then
      read(arg_value,*) nsteps
    endif
  enddo
  Ly=Lx
  Lz=Lx
  Ls=max(Lx,Ly,Lz)

!  call OMP_set_num_threads(n_proc)
  npx = nprocs**(1./3.)
  do while (nprocs .ne. nprocs/npx * npx)
    npx=npx-1
  end do
  surface = nprocs / npx
  npy = surface**(1./2.)
  do while (surface .ne. surface/npy * npy)
    npy=npy-1
  end do
  npz = surface / npy

! myid = nzid*npx*npy + nyid*npx + nxid
! myid is (nxid,nyid,nzid) in whole procs
  nzid = myid/(npx*npy)
  nyid = (myid-nzid*npx*npy)/npx
  nxid = myid-nzid*npx*npy-nyid*npx

  if(nxid.gt.0) then
     x_left = nxid-1 + nyid*npx + nzid*npx*npy
  else
     x_left=MPI_PROC_NULL
  endif

  if(nxid.lt.npx-1) then
     x_right = nxid+1 + nyid*npx + nzid*npx*npy
  else
     x_right = MPI_PROC_NULL
  endif

  if(nyid.gt.0) then
     y_left = (nyid-1)*npx + nxid + nzid*npx*npy
  else
     y_left=MPI_PROC_NULL
  endif

  if(nyid.lt.npy-1) then
     y_right = (nyid+1)*npx + nxid + nzid*npx*npy
  else
     y_right = MPI_PROC_NULL
  endif

  if(nzid.gt.0) then
     z_left = (nzid-1)*npx*npy + nyid*npx + nxid
  else
     z_left=MPI_PROC_NULL
  endif

  if(nzid.lt.npz-1) then
     z_right = (nzid+1)*npx*npy + nyid*npx + nxid
  else
     z_right = MPI_PROC_NULL
  endif

  nLx=Lx/npx
  nLy=Ly/npy
  nLz=Lz/npz

  nlocal = nLx*nLy*nLz
  nxst = nLx*nxid
  nyst = nLy*nyid
  nzst = nLz*nzid

  allocate(uc(mnp,-4:nLx+5,-4:nLy+5,-4:nLz+5),frc(3,nLx,nLy,nLz))

! 1. commit uc(1:mnp,(-4:0) or (nLx+1:nLx+5),1:nLy,1:nLz) datatype to do boundary transfer
  call mpi_type_vector(nLy,mnp*5,(nLx+10)*mnp,mpi_double_precision,type0,ierr)
  call mpi_type_commit(type0,ierr)
  call mpi_type_extent(mpi_double_precision,lengthDB,ierr)
  call mpi_type_hvector(nLz,1,mnp*(nLx+10)*(nLy+10)*lengthDB,type0,xb_type,ierr)
  call mpi_type_free(type0,ierr)
  call mpi_type_commit(xb_type,ierr)

! 2. commit uc(1:mnp,1:nLx,(-4:0) or (nLy+1:nLy+5),1:nLz) datatype to do boundary transfer
  call mpi_type_vector(5,mnp*nLx,(nLx+10)*mnp,mpi_double_precision,type0,ierr)
  call mpi_type_commit(type0,ierr)
  call mpi_type_hvector(nLz,1,mnp*(nLx+10)*(nLy+10)*lengthDB,type0,yb_type,ierr)
  call mpi_type_free(type0,ierr)
  call mpi_type_commit(yb_type,ierr)

! 3. commit uc(1:mnp,1:nLx,1:nLy,(-4:0) or (nLz+1:nLz+5)) datatype to do boundary transfer
  call mpi_type_vector(nLy,mnp*nLx,(nLx+10)*mnp,mpi_double_precision,type0,ierr)
  call mpi_type_commit(type0,ierr)
  call mpi_type_hvector(5,1,mnp*(nLx+10)*(nLy+10)*lengthDB,type0,zb_type,ierr)
  call mpi_type_free(type0,ierr)
!  call mpi_type_free(lengthDB,ierr)
  call mpi_type_commit(zb_type,ierr)

  call setup
! ****** initialization ******

  call initialize

  allocate(rhs(mnp,nLx,nLy,nLz))
  rhs=0.

  time_cpu=mpi_wtime()

!     'io' is the main controller which rotates between 0 and "mt-1"
!     corresponding to the different stages in Runge-Kutta schemes

  io = 0
  stop_sign = .false.
  output = .false.
! ------------------------------------------------------------------

  do

!     the Runge-Kutta inner stage time --- this works only for the 3rd order
!     Runge-Kutta and must be changed for fourth order
!     note: when io>=1 tnum has already been updated
     if( io==0 .and. stop_sign) exit
     if(nt.gt.nsteps) exit
! ****** compute -f(u)_x, -g(u)_y, -h(u)_z ******

     time_init=mpi_wtime()
     call weno_hydrodynamics(io)
     if(myid.eq.0) call CPU_screentime(time_init,cpu_weno,'CPU time in WENO')

! ****** compute time step size "dt" ******

     if (io==0) then
        time_init=mpi_wtime()

        call cflc(aam)

        dt_cfl= cfl/aam
        if(myid.eq.0) then
           print*,'aam,dt_cfl',aam,dt_cfl
        endif
        dt=dt_cfl

        if(.not.supin) then
           if ((tnum+dt)-t_out(i_out)>=0.) then
              dt = t_out(i_out)-tnum
              if(i_out>=n_out) stop_sign=.true.
              i_out = i_out + 1
              output = .true.
           endif
        else  if(((tnum+dt)-t_out(i_out)>=0.).or.((tnum+dt)-tsupend>=0.)) then
           if(tsupend<t_out(i_out)) then
              dt=tsupend-tnum
           else if(t_out(i_out)<tsupend)  then
              dt=t_out(i_out)-tnum
              if(i_out>=n_out) stop_sign=.true.
              i_out=i_out+1
              output=.true.
           else
              dt=t_out(i_out)-tnum
              if(i_out>=n_out) stop_sign=.true.
              i_out=i_out+1
              output=.true.
           endif
        endif

        tnum = tnum + dt
        nt=nt+1
        if(myid.eq.0) call CPU_screentime(time_init,cpu_dt,'CPU time in dt')
     endif

! ****** Runge-Kutta scheme for advancing in time ******

! ****** compute the gravity and evolve the collisionless particles ******

     time_init=mpi_wtime()
     call RK_forward(io)
     if(myid.eq.0) call CPU_screentime(time_init,cpu_grav,'CPU time in Grav')

     if (io==2) then

        if(myid.eq.0) call CPU_screentime(time_cpu,cpu_step,'CPU time per step')

        if(myid.eq.0) call write_log

        if(tnum==tsupend)  then
           supin=.false.
           write(*,*) 'the end of injection'
        endif
        if (output) then
           !           call data_output(output) 
           output=.false.
        else if (mod(nt,ipdump)==0) then
           !           call data_output(output)
        endif

     endif

     io = mod( io+1, mt )
  enddo

!     --------------------- end  time evolution ---------------------------

  call mpi_type_free(xb_type,ierr)
  call mpi_type_free(yb_type,ierr)
  call mpi_type_free(zb_type,ierr)
  deallocate(uc,frc,rhs,STAT=istat)
  if(myid.eq.0) call CPU_screentime(time1,time2,'^3,time=')
  call MPI_FINALIZE(ierr)
  stop

contains

  subroutine RK_forward(io)
    use weno_kernel,only:rhs_out,data_output2

    real,dimension(5)::fq(5)
    integer,save:: iii=16
    integer:: idbg=0

    qk = ark(io+1)
    rk3 = brk(io+1)
    dtrk = dt*rk3
    
    do m=2,4
       ic=m-1
!$omp parallel do default(shared) private(i,j,k)
       do k=1,nLz
          do j=1,nLy
             do i=1,nLx
                
                rhs(m,i,j,k)=rhs(m,i,j,k)-uc(1,i,j,k)*frc(ic,i,j,k)
                rhs(5,i,j,k)=rhs(5,i,j,k)-uc(m,i,j,k)*frc(ic,i,j,k)
                
             enddo
          enddo
       enddo
!$omp end parallel do
    enddo
!$omp parallel do default(shared) private(i,j,k,den,fq,Q)      
    do k = 1,nLz
       do j = 1,nLy
          do i = 1,nLx
             
             den=uc(1,i,j,k)
             
             if(den>0.) then
                fq(1:5)=uc(1:5,i,j,k)
                call  extQ2(fq,Q)
                rhs(5,i,j,k)=rhs(5,i,j,k)-Q
             endif
          enddo
       enddo
    enddo
!$omp end parallel do
!$omp parallel do default(shared) private(iz)
    do iz=1,nLz
       uc(:,1:nLx,1:nLy,iz) = uc(:,1:nLx,1:nLy,iz) + dtrk * rhs(:,1:nLx,1:nLy,iz)
    enddo
!$omp end parallel do

  end subroutine RK_forward
  
end program main

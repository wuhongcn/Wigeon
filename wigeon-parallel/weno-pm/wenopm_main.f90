!    WENOPM code designed for simulating the coupled system of both 
!    baryonic matter and dark matter. The hydrodynamics is solved using 5th 
!    WENO scheme in 3-dim space incorporated with 3rd Runge-Kutta TVD 
!    time discretization. The gravity solver is using standard particle-mesh 
!    scheme. The code was developed by Feng LongLong (fengll@pmo.ac.cn)and 
!    Shu ChiWang (shu@cfm.brown.edu).  
!    Chi-Wang Shu, 12/21/00, changed at USTC for the astrophysics problem
!    this code is capable to handle 3D cartesian non-uniform but smooth meshes
!    computing the 3D compressible Euler equations. 
!    Long-Long Feng 04/15/01, added the PM scheme for self-gravity system 
!    with periodic boundary condition. 
!    Updated with (1) adding ionization & radiative heating-cooling under optical 
!    thin approximation (2) modified entropy equation for solving high-Mach problem.   
!    The implementation of the WENO-PM code is detailed in astro-ph/0403233
!    "A Hybrid Cosmological Hydrodynamic/N-body Code Based on a Weighted
!     Essentially Non-Oscillatory Scheme" Feng, Shu & Zhang 2004 submitted to ApJS.
!    version : Mar.30, 2004. The code is not publically available for time-being
!    and is limitted for distributing among our research group. 

! -----------------------------------------------------------------

program main

use datablock
use cosmology
use weno_kernel
use gravity
use IGM_heating_cooling
use readin_out
use recipes
use hydro_stat

logical output,stop_sign,negative_sign 

!     -------------------------------------------------------------------
!     Name:      weno-ustctw.f
!     Function:  drive routine
!      System to solve: u_t + (f(u)_x + g(u)_y +h(u)_z) / a(t) = forc(u,t)
!     or
!      u_t = RHS = (-f(u)_x - g(u)_y - h(u)_z) / a(t) + forc(u,t)
!	 ---------------------------------------------------------------

! ****** readin parameters ******

call OMP_set_num_threads(n_proc)

call setup

! ****** initialization ******

time_init=omp_get_wtime()
call initialize

do ic=1,mnp
      uct(ic,1:Lx,1:Ly,1:Lz)=uc(ic,1:Lx,1:Ly,1:Lz)
enddo
!$omp parallel do default(shared) private(i)
do i=1,nobj
   rtemp(1,i) =  r(1,i) 
   rtemp(2,i) =  r(2,i) 
   rtemp(3,i) =  r(3,i) 
   vtemp(1,i) = v(1,i) 
   vtemp(2,i) = v(2,i) 
   vtemp(3,i) = v(3,i) 
enddo
!$omp end parallel do
call CPU_screentime(time_init,cpu_init,'CPU time in initialization')

time_init=omp_get_wtime()
time_cpu=time_init

!     'io' is the main controller which rotates between 0 and "mt-1"
!     corresponding to the different stages in Runge-Kutta schemes

io = 0
stop_sign = .false.
output = .false.
negative_sign= .false.
dt_neg=1.0
if (ntot==0) stop
istoptimes=0
! ------------------------------------------------------------------

do

!     the Runge-Kutta inner stage time --- this works only for the 3rd order
!     Runge-Kutta and must be changed for fourth order
!     note: when io>=1 tnum has already been updated

! stop when several times
istoptimes=istoptimes+1

trk=tnum

if (io==1) trk = tnum+(brk(1)-1.)*dt
if (io==2) trk = tnum+(brk(1)+brk(2)*(ark(2)+1.)-1.)*dt

! ****** computing cosmic age  ******

call redshift_comovingtime(trk,zf)

af = 1./(1.+zf) 
cf = 1./af


!generate cooling heating table
if(zf<=z_reion) xJ=xJ21

if(io==0) call cooling_tab(xJ,zf)
!if(io==0.and.zf==z_reion) call cooling_tab(xJ,zf)

! zf_exact=(sqrt(1.+zinit)-0.5*trk)**2-1. ! flat model 
! print*,'cosmic time',zf_exact,zf

if( io==0 .and. stop_sign) exit
if( istoptimes > 3) exit
! ****** compute -f(u)_x, -g(u)_y, -h(u)_z ******
time_init=omp_get_wtime()
call weno_hydrodynamics(io)
call CPU_screentime(time_init,cpu_weno,'CPU time in WENO')

! ****** compute time step size "dt" ******

! ****** Runge-Kutta scheme for advancing in time ******

! ****** compute the gravity and evolve the collisionless particles ******

time_init=omp_get_wtime()
call RK_forward(io)
call CPU_screentime(time_init,cpu_grav,'CPU time in Grav')

if (io==2) then

   if(.not.negative_sign) then
       do ic=1,mnp
         uc(ic,1:Lx,1:Ly,1:Lz)=uct(ic,1:Lx,1:Ly,1:Lz)
       enddo
       call regulation
!$omp parallel do default(shared) private(i)
        do i=1,nobj
          r(1,i) = rtemp(1,i) 
          r(2,i) = rtemp(2,i) 
          r(3,i) = rtemp(3,i) 
          v(1,i) = vtemp(1,i) 
          v(2,i) = vtemp(2,i) 
          v(3,i) = vtemp(3,i) 
         enddo
!$omp end parallel do
   endif

   call CPU_screentime(time_cpu,cpu_step,'CPU time per step')

   call write_log
   if (output) then
      call data_output(output) ! output at redshift z=z_out(:)
      output=.false.
   else if (mod(nt,ipdump)==0) then
      call data_output(output)
   endif
  if(.not.negative_sign) then   
      dt_neg=1.0
  else
!$omp parallel do default(shared) private(i)
    do i=1,nobj
       rtemp(1,i) = r(1,i) 
       rtemp(2,i) = r(2,i) 
       rtemp(3,i) = r(3,i) 
       vtemp(1,i) = v(1,i) 
       vtemp(2,i) = v(2,i) 
       vtemp(3,i) = v(3,i) 
    enddo
!$omp end parallel do
      tnum=tnum-dt
      dt_neg=0.5*dt_neg
  endif 
    do ic=1,mnp
      uct(ic,1:Lx,1:Ly,1:Lz)=uc(ic,1:Lx,1:Ly,1:Lz)
    enddo
  write(*,*) 'dt_neg',dt_neg
endif

io = mod( io+1, mt )

enddo

!     --------------------- end  time evolution ---------------------------

close(68)

deallocate (uc,uct,STAT=istat)
deallocate (r,v,STAT=istat)
deallocate (rtemp,vtemp,STAT=istat)
deallocate (rs,vs,STAT=istat)

stop

contains

subroutine RK_forward(io)
real,dimension(5)::fq(5)
real,allocatable,dimension(:)::g

allocate (g(nobj),STAT=istat)

rLx = float(Lx)
rLy = float(Ly)
rLz = float(Lz)

qk = ark(io+1)
rk3 = brk(io+1)
!dtrk = dt*rk3

!  mass assignment using CIC or TSC scheme      
call mass_assignment

!$omp parallel do default(shared) private(ix,iy,iz)
do iz=1,Lz
  do iy=1,Ly
    do ix=1,Lx
    !d3(ix,iy,iz) = af * ( DM_baryon * d3(ix,iy,iz) + uc(1,ix,iy,iz))
    d3(ix,iy,iz) = af * ( DM_baryon * d3(ix,iy,iz) + uct(1,ix,iy,iz))
    enddo
  enddo
enddo
!$omp end parallel do

! convolution with Green function and computing gravitational potential  
!time_init=omp_get_wtime()
call cnvl3d(d3,gk)

!call CPU_screentime(time_init,cpu_grav,'CPU time in cnvl')


!time_init=omp_get_wtime()
do m = 2, 4

   ic = m-1
   call gridf10(ic,d3,frc)  ! 4pt difference
   call CIC_interp(frc,g)

!$omp parallel do default(shared) private(i,j,k)
   do k = 1, Lz
      do j = 1, Ly
         do i = 1, Lx
            
            !rhs(m,i,j,k) = rhs(m,i,j,k) + uc(1,i,j,k) * frc(i,j,k)
            !rhs(5,i,j,k) = rhs(5,i,j,k) + uc(m,i,j,k) * frc(i,j,k) 
            rhs(m,i,j,k) = rhs(m,i,j,k) + uct(1,i,j,k) * frc(i,j,k)
            rhs(5,i,j,k) = rhs(5,i,j,k) + uct(m,i,j,k) * frc(i,j,k) 

         enddo 
      enddo
   enddo
!$omp end parallel do


!$omp parallel do default(shared) private(i)
   do i=1,nobj
      !rs(ic,i) = qk * rs(ic,i) + v(ic,i)
      rs(ic,i) = qk * rs(ic,i) + vtemp(ic,i)
      vs(ic,i) = qk * vs(ic,i) + g(i) 
   enddo
!$omp end parallel do

enddo
!call CPU_screentime(time_init,cpu_grav,'CPU time in gridf')

!   the energy equation and entropy equation:

!xJ = 0.
!if (zf <= z_reion) xJ = xJ21 

dt_cool=1.0


time_cool=omp_get_wtime()
!$omp parallel do default(shared) private(i,j,k,den,fq,Q,dt_col) & 
!$ reduction(min:dt_cool)    
do k = 1,Lz
   do j = 1,Ly
      do i = 1,Lx
	     !den=uc(1,i,j,k)
	     den=uct(1,i,j,k)

         if (id_coolheat==1.and.den/=0.) then
            !fq(1:5)=uc(1:5,i,j,k)
            fq(1:5)=uct(1:5,i,j,k)
            call extQ(zf,xJ,fq,Q) 
            rhs(5,i,j,k)=rhs(5,i,j,k)-Q
             rhs(6,i,j,k)=rhs(6,i,j,k)-gm1*Q/entropy_unit/abs(den)**gm1
         endif
         
      enddo
   enddo
enddo
!$omp end parallel do
call CPU_screentime(time_cool,cpu_cool,'CPU time in cool')
if (io==0) then
         
   call cflc(aam)
   dt_cfl=1.0*cfl/aam
   if(zf<=z_reion)   dt_cfl=1.0*cfl/aam/(sqrt(1.+zf))
   dt_exp =1.00*c_exp * cf * cf /  E(zf)
   dt=min(dt_cfl,dt_exp)
   dt=dt_neg*dt
   if ((tnum+dt)-t_out(i_out)>=0.) then
      dt = t_out(i_out)-tnum
      if (i_out>=n_out) stop_sign = .true.
         i_out = i_out + 1
         output = .true.
   endif

      tnum = tnum + dt
      nt = nt + 1
      if (nt>=ntot.or.tnum>=tend) stop_sign = .true.
      if(ncell_zero>7) stop_sign=.true. 
endif

dtrk = dt*rk3
negative_sign=.false.

if(io==2) then
!$omp parallel do default(shared) private(ic,ix,iy,iz)
do iz=1,Lz
  do iy=1,Ly
    do ix=1,Lx
       do ic=1,mnp
!      uc(:,:,:,iz) = uc(:,:,:,iz) + dtrk * rhs(:,:,:,iz)
          uct(ic,ix,iy,iz)=uct(ic,ix,iy,iz)+dtrk*rhs(ic,ix,iy,iz)
       if(uct(1,ix,iy,iz)<0.0) then
           negative_sign=.true.
       endif
     enddo
    enddo
  enddo
enddo
!$omp end parallel do
else
!$omp parallel do default(shared) private(ic,ix,iy,iz)
do iz=1,Lz
  do iy=1,Ly
    do ix=1,Lx
       do ic=1,mnp
!      uc(:,:,:,iz) = uc(:,:,:,iz) + dtrk * rhs(:,:,:,iz)
          uct(ic,ix,iy,iz)=uct(ic,ix,iy,iz)+dtrk*rhs(ic,ix,iy,iz)
     enddo
    enddo
  enddo
enddo
!$omp end parallel do

endif


!$omp parallel do default(shared) private(i)
do i=1,nobj
   rtemp(1,i) = mod(rtemp(1,i) + dtrk * rs(1,i) + rLx-1., rLx) + 1.
   rtemp(2,i) = mod(rtemp(2,i) + dtrk * rs(2,i) + rLy-1., rLy) + 1.
   rtemp(3,i) = mod(rtemp(3,i) + dtrk * rs(3,i) + rLz-1., rLz) + 1. 
   vtemp(1,i) = vtemp(1,i) + dtrk * vs(1,i)
   vtemp(2,i) = vtemp(2,i) + dtrk * vs(2,i)
   vtemp(3,i) = vtemp(3,i) + dtrk * vs(3,i)
enddo
!$omp end parallel do


deallocate (g,STAT=istat)
       
end subroutine RK_forward

end program main


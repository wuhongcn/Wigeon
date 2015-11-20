!     rhs: represents all the stuff on the right hand side of the ODE
!     in time, i.e.  u_t = rhs
module datablock

integer, parameter:: mt=3,mn=5,mnp=6,ipdump=10,n_out=9,n_proc=16
integer, parameter:: Lx=256,Ly=Lx,Lz=Lx,Ls=max(Lx,Ly,Lz)
real, parameter:: pi=3.1415926535, tpi=2.*pi

real(kind=8),parameter:: x_hydrogen=0.76, y_helium=0.24
real(kind=8),parameter:: hydro_mass=1.6735344e-24 ! gram
real(kind=8),parameter:: k_boltzmann=1.380658e-16 ! erg K^-1

!  constants:
real, save:: omega0,xlambda0,omegab,gamma0,h100
real, save:: DM_baryon,omg,rLbox,soft
real, save:: zinit,dt_ave,time_cpu
real, save:: xJ21,z_reion,gamma,gm1  
real, save:: l_unit,t_unit,v_unit,rho_unit,n_unit,e_rate_unit,tem_unit,cooling_coff,entropy_unit
real, save,dimension(n_out)::t_out

integer :: nobj,ntot,nt,ist,i_out,nstop,ncell_zero
integer :: id_coolheat
real :: xJ,den_av,den_min,den_max,vx_max,vy_max,vz_max,cs_max,den_sigma,vel_sigma
real :: vxmax_cdm,vymax_cdm,vzmax_cdm,vsigma_cdm,div_max
real :: cpu_weno,cpu_grav,cpu_step,cpu_cool,cpu_cnvl

character(len=60),parameter::datadir='./'

!     grid related variables:
real :: cfl,c_exp,c_div,c_cdm
real :: tend,dt,tnum,trk,af,cf,zf,dt_cfl,dt_exp,dt_div,dt_cdm,dt_neg

!  necessary vars for solutions:  
real,dimension(3)::ark,brk  
!real,dimension(mnp,Lx,Ly,Lz)::uc
!real,dimension(mnp,Lx,Ly,Lz)::uct
real,allocatable,dimension(:,:,:,:)::uc
real,allocatable,dimension(:,:,:,:)::uct
real,dimension(mnp,Lx,Ly,Lz)::rhs

!  necessary vars for CDM particles
real,allocatable,dimension(:,:)::r,v
real,allocatable,dimension(:,:)::rtemp,vtemp
real,allocatable,dimension(:,:)::rs,vs
integer,allocatable,dimension(:)::ll
real,dimension(Lx,Ly,Lz)::d3,frc
real,save,dimension(Lx/2+1,Ly/2+1,Lz/2+1)::gk ! Green Function: Gravity-Solver

! cooling table 
integer,parameter::Ns_table=1500
real,dimension(Ns_table)::tem_table,den_table
real,dimension(Ns_table,Ns_table)::ch_table

integer :: kd1,kd2
save tem_table,den_table,ch_table,kd1,kd2
data kd1,kd2/0,0/
contains

subroutine TVD_RK_parameter
real (kind=8) ::c2,z1,z2,z3,z4,z5,z6

c2=0.924574D0
z1=dsqrt(36.*c2**4+36.*c2**3-135.*c2**2+84.*c2-12.)
z2=2.*c2**2+c2-2.
z3=12.*c2**4-18.*c2**3+18.*c2**2-11.*c2+2.
z4=36.*c2**4-36.*c2**3+13.*c2**2-8.*c2+4.
z5=69.*c2**3-62.*c2**2+28.*c2-8.
z6=34.*c2**4-46.*c2**3+34.*c2**2-13.*c2+2.
brk(1)=c2
brk(2)=(12.*c2*(c2-1.)*(3.*z2-z1)-(3.*z2-z1)**2) / (144.*c2*(3.*c2-2.)*(c2-1.)**2)
brk(3)=-24.*(3.*c2-2.)*(c2-1.)**2 / ((3.*z2-z1)**2-12.*c2*(c2-1.)*(3.*z2-z1))
ark(1)=0.
ark(2)=(-z1*(6.*c2**2-4.*c2+1.)+3.*z3) / ((2.*c2+1.)*z1-3.*(c2+2)*(2.*c2-1.)**2)
ark(3)=(-z4*z1+108.*(2.*c2-1.)*c2**5-3.*(2.*c2-1)*z5) / &
& (24.*z1*c2*(c2-1.)**4+72.*c2*z6+72.*c2**6*(2.*c2-13.))

end subroutine TVD_RK_parameter

subroutine Phys_Unit

!DM_baryon=omega0/omegab
DM_baryon=(omega0-omegab)/omegab  !/02/04/2009/

gm1=gamma-1.0

rho_crit=1.8788e-29*h100*h100 ! g/cm^3

! length unit: in unit of h^{-1}Mpc    /(?2007/08/01 what does h^{-1}  mean,it seems
! should be deleted       
l_unit=rLbox/float(Ls)
 
! time unit: km^-1 s Mpc
t_unit=0.01/h100
t_unit_second=3.0857e+17/h100

! velocity unit: in unit of km/s 
v_unit=l_unit/t_unit

! cofficinets in Poisson Green function: = G*rho_unit*t_unit^2
omg=3.*omegab/4./tpi

! density unit: in unit of g/cm^3  
rho_unit = omegab*rho_crit  

! H number density unit: 1/cm^3
n_unit = x_hydrogen * rho_unit / hydro_mass

! energy rate unit dE/dt = rho_unit*v_unit^2/t_unit_second: in unit of 10^{-37} erg/cm^3/s
e_rate_unit = 6.08873 * omegab * h100**3 * v_unit * v_unit 

! temperature unit: tem_phys=\mu*(P/rho)*T_unit 
tem_unit= 1.e+10 * hydro_mass/k_boltzmann * v_unit * v_unit           

! cooling-heating cofficients
cooling_coff = 1.e+16 * n_unit * n_unit / e_rate_unit

! initial temperature and entropy unit 
zt=137.
T_cmb=2.73
T0_igm=T_cmb*(1.+zinit)*(1.+zinit)/(1.+zt)
!entropy_unit=T0_igm/tem_unit
entropy_unit=T0_igm/tem_unit*(x_hydrogen+y_helium/4.0) !wszhu 01/24/2009
end subroutine Phys_Unit

end module datablock

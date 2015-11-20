!     rhs: represents all the stuff on the right hand side of the ODE
!     in time, i.e.  u_t = rhs
module datablock
include'mpif.h'
integer, parameter:: mt=3,mn=5,mnp=6,ipdump=50,n_out=2,n_cool=91,n_proc=4
!integer, parameter:: Lx=2560,Ly=2560,Lz=2560,Ls=max(Lx,Ly,Lz)
integer :: nsteps
integer :: Lx,Ly,Lz,Ls
integer :: nprocs,myid,nLz,nLx,nLy,nxid,nyid,nzid,nxst,nyst,nzst,nlocal
integer :: x_left,x_right,y_left,y_right,z_left,z_right
integer :: npx,npy,npz,xb_type,yb_type,zb_type
integer :: istatus(mpi_status_size,4),ireq(4),STATUS(MPI_STATUS_SIZE)
real, parameter:: pi=3.1415926536, tpi=2.*pi

real, parameter::x_hydrogen=0.76,y_helium=0.24
real, parameter::hydro_mass=1.6735344e-24  !gram
real, parameter::k_boltzmann=1.380658e-16  !erg K^-1
real, parameter::g_gravity=6.67e-8   !cm^3 s^-2 g^-1
!  constants:
real, save:: gamma,gm1,rLbox,rhoc,rho0,radc,radd,g0,dedp,dmdp
real, save,dimension(n_out)::t_out
real, save,dimension(n_cool)::coolrate
real, save::l_unit,v_unit,t_unit,rho_unit,n_unit,tem_unit,ene_unit,e_rate_unit,cooling_coff,entropy_unit,p_igm,cs_init
integer :: ist,i_out,ncell_zero,nt
real :: den_av,den_min,den_max,vx_max,vy_max,vz_max,cs_max,den_sigma,vel_sigma
real :: cpu_weno,cpu_step,cpu_weno_nocom
logical::supin
character(len=60),parameter::datadir='./'

!     grid related variables:
real,save :: cfl
real :: tend,tsupend,dt,tnum,dt_cfl

!  necessary vars for solutions:  
real,dimension(3)::ark,brk  
real,allocatable,dimension(:,:,:,:)::uc,rhs,frc
contains

subroutine TVD_RK_parameter
real (kind=8) ::c2,z1,z2,z3,z4,z5,z6

c2=0.924574
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

! length unit: in unit of kpc
l_unit=rLbox/float(Ls-1)

! velocity unit: in unit of km/s
v_unit=1.0

!time unit : in unit of Myr
t_unit=l_unit*3.09e+16/(v_unit*1.e+6*365*86400)

!density  unit : in unit of g/cm^3
rho_unit=1.0e-10*rho0*p_igm/cs_init/cs_init

! H number density unit: 1/cm^3
n_unit=x_hydrogen*rho_unit/hydro_mass

!temperature unit :  in unit of K
tem_unit=1.e+10*hydro_mass/k_boltzmann*v_unit*v_unit

!energy unit:  in unit of erg cm^-3
ene_unit=1.0e+10*rho_unit*v_unit*v_unit

!energy rate unit dE/dt =rho_unit*v_unit^2/t_unit_second: in unit of erg/cm^3/s
e_rate_unit=ene_unit/t_unit/(1.e+6*365*86400)

!cooling-heating cofficients
cooling_coff=n_unit*n_unit/e_rate_unit


! entropy unit
!entropy_unit=

end subroutine Phys_Unit
end module datablock

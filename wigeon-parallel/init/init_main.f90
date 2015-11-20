! The code is to initialize the particle distribution and density field 
! for N-body/hydrodynamic simulation using Zeldovish approximation 
! writen by Feng LongLong in August,2001 (1st version) 
! revised in March, 2003  (fengll@ustc.edu.cn) 
 
module set_parameter

real,parameter::pi=3.14159265358979,tpi=2.*pi

! Cosmological Parameters 
! omega0   : density parameter
! xlambda0 : cosmological constant
! h100: Hubble constant
! omegab: baryonic density parameter 
! ind : spectral index of primodial power spectrum (ind=1 -> scale free spectrum)
! ispec : kinds of the power spectra shaped by transfer function 
!          ispec = 1 -> power law
!          ispec = 2 -> CDM power spectrum from BBKS
!          ispec = 3 -> power spectrum from Hu Wayne 
!          ispec = 4 -> HDM power spectrum from BBKS
! ksp : normalizing the power spectrum
!          ksp = 1 -> structure normalization (cluster abundance)
!          ksp = 2 -> COBE normalization 

! Simulation Parameters
! boxsize: simulation box size in unit of h^-1Mpc
! n_object: particle number  
! Lmesh: mesh size in 1D (any combination of 2^p3^p5^r)
! soft: soft parameter for PM calculation 
! zstart: initial redshift 

! units
! r : grid size [0,Lmesh]
! v : km/s

real ::omega0,xlambda0,omegab,omegah,h100,sigma80 
!integer :: ist,ind,ispec,ksp,nobj,L
integer :: ist,ispec,ksp,nobj,L
!real,save :: rLbox,zinit,alpha,soft,cfl,c_exp,c_div,c_cdm
real,save :: rLbox,zinit,alpha,soft,cfl,c_exp,c_div,c_cdm,ind
real,save :: omega_i,xlambda_i,delta,omegahh,gamma,tstart,h0t0
character(len=60),parameter::datadir='../weno-pm/'

namelist/cosmic/omega0,xlambda0,omegab,omegah,h100,ind,ispec,ksp,sigma80
namelist/simulation/nobj,L,rLbox,zinit,soft,alpha,cfl,c_exp,c_div,c_cdm

end module set_parameter

!   -----------------------------------------------------------------

Program Normalization
use set_parameter
real,allocatable,dimension(:,:)::r,v
real,allocatable,dimension(:,:,:):: d3,p3,den
character fileno*4
      
!ist=integerread(' input simulation no. =')
ist=1000
write(fileno,'(i4)') ist

call reading_parameters(ist)

allocate (r(3,nobj),v(3,nobj),STAT=istat)
allocate (d3(L,L,L),p3(L,L,L),STAT=istat)
allocate (den(L,L,L),STAT=istat)

ai=1./(1.+zinit)      
CALL timeah(tstart,omega0,xlambda0,ai,htime)
CALL timeah(t0,omega0,xlambda0,1.,h0t0)

vd=100.*rLbox/float(L)/(1.+zinit)
vd=vd*E(zinit)

call zeldovich(zinit,nobj,L,rLbox,r,v,d3,p3)

v=vd*alpha*v
p3=vd*alpha*p3
d3=1.+d3

!     ---------------- Saving the Data -------------------------

open (2,file=datadir(1:lnblnk(datadir))//'r'//fileno//'.0000',&
& status='unknown',form='unformatted')
write(2) ((r(ic,i),ic=1,3),i=1,nobj)
close(2)

open (2,file=datadir(1:lnblnk(datadir))//'v'//fileno//'.0000',&
& status='unknown',form='unformatted')
write(2) ((v(ic,i),ic=1,3),i=1,nobj)
close(2)

open (3,file=datadir(1:lnblnk(datadir))//'igm'//fileno//'.0000',&
& status='unknown',form='unformatted')
write(3) d3
close (3)

open (3,file=datadir(1:lnblnk(datadir))//'pot'//fileno//'.0000',&
& status='unknown',form='unformatted')
write(3) p3
close (3)
      
deallocate(r,v,STAT=istat)
deallocate(d3,p3,den,STAT=istat)
     
stop
end

subroutine reading_parameters(iq)
use set_parameter
character fileno*4

write(fileno,'(i4)')iq

open (2,file=datadir(1:lnblnk(datadir))//'weno'//fileno//'.ini',status='old')
read (2,nml=cosmic)
write(*,nml=cosmic)
read (2,nml=simulation)
write(*,nml=simulation)
close (2)

return
end subroutine reading_parameters

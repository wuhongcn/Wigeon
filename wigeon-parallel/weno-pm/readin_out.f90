module readin_out
use datablock
use cosmology
use recipes
use gravity
use hydro_stat
use IGM_heating_cooling
contains

!     ---------------------------------------------------------------------------
!     Name:      data_io.for
!     Function:  read in some parameters
!     ---------------------------------------------------------------------------

subroutine setup
real,dimension(n_out)::z_out
character(len=60):: parafile

!namelist/cosmic/omega0,xlambda0,omegab,omegah,h100,ind,ispec,ksp,sigma80
namelist/cosmic/omega0,xlambda0,omegab,h100
namelist/simulation/nobj,L,rLbox,zinit,soft,alpha,cfl,c_exp,c_div,c_cdm
namelist/hydro/gamma,id_coolheat,z_reion,xJ21,ntot,nstop,z_out
     
!  readin serial no. of the input simulation

!ist=IntegerRead('input fileno =')
ist=1000
write(parafile,'(''cosm''i4''.ini'')')ist

open(55,file=datadir(1:lblnk(datadir))//parafile(1:lblnk(parafile)),status='old')
read(55,nml=cosmic)
write(*,nml=cosmic)
read(55,nml=simulation)
write(*,nml=simulation)
read(55,nml=hydro)
write(*,nml=hydro)
close(55)
write(*,*) 'L=',L,'Ls=',Ls
if (L/=Ls) stop 'wrong simulation parameter'

! cosmic age

tint=0.
tend=age(0.)
dt_ave=(tend-tint)/float(ntot)

! table the output time sequence 

do i=1,n_out
   zi=z_out(i)
   t_out(i)=age(zi)
enddo 

call Phys_Unit

! table the Green function for Poisson equation  
! call greenc(soft,Lx,Ly,Lz,gk)

call greenz
gk=omg*gk
write(*,'(''generating Green function'')')

allocate(uc(mnp,Lx,Ly,Lz),uct(mnp,Lx,Ly,Lz),STAT=istat)
allocate(r(3,nobj),v(3,nobj),STAT=istat)
allocate(rtemp(3,nobj),vtemp(3,nobj),STAT=istat)
allocate(ll(nobj),STAT=istat)

call TVD_RK_parameter

! -----------------------------------------------------------------------------

end subroutine setup

!     -----------------------------------------------------------------------
!     Name:      initialize.f
!     Function:  set up the grid and initial condition ( u(x,y,z) at t= 0 )
!     -----------------------------------------------------------------------

subroutine initialize
character (len=60):: datafile

!write(*,'(''enter the status of the initial condition'')')
!i_run=IntegerRead('initial at t=0. (0) / restart from a restart file (1) = ')
i_run=0

if (i_run==0) then
         
   nt=0  
   tnum = 0.
   i_out=1
   a_init=1./(1+zinit)
   zf=zinit

   write(datafile,'(i4''.0000'')')ist 

!  readin the initial distribution of CDM 

   open(10,FILE=datadir(1:lblnk(datadir))//'r'//datafile(1:lblnk(datafile)), &
&      STATUS='unknown',FORM='UNFORMATTED')
   read(10) ((r(ic,i),ic=1,3),i=1,nobj)
   close(10)

   open(20,FILE=datadir(1:lblnk(datadir))//'v'//datafile(1:lblnk(datafile)), &
&      STATUS='unknown',FORM='UNFORMATTED')
   read(20) ((v(ic,i),ic=1,3),i=1,nobj)
   close(20)

   v(:,:)=a_init*v(:,:)/v_unit

!  initialize the gas component

   open(30,FILE=datadir(1:lblnk(datadir))//'igm'//datafile(1:lblnk(datafile)), &
&      STATUS='unknown',FORM='UNFORMATTED')
   read(30) d3
   close(30)

   uc(1,:,:,:)=d3(:,:,:)

   open(40,FILE=datadir(1:lblnk(datadir))//'pot'//datafile(1:lblnk(datafile)), &
&      STATUS='unknown',FORM='UNFORMATTED')
   read(40) d3
   close(40)

   d3(:,:,:)=d3(:,:,:)/v_unit

   do ic=1,3
     call gridf10(ic,d3,frc)
     do iz=1,Lz
       do iy=1,Ly
         do ix=1,Lx
           uc(ic+1,ix,iy,iz) = a_init * uc(1,ix,iy,iz) * frc(ix,iy,iz)
         enddo
       enddo
     enddo
   enddo

   do iz=1,Lz
      do iy=1,Ly
         do ix=1,Lx
            den=uc(1,ix,iy,iz)
            cden=1./den
            vx=uc(2,ix,iy,iz)*cden
            vy=uc(3,ix,iy,iz)*cden
            vz=uc(4,ix,iy,iz)*cden
            vk=0.5*(vx*vx+vy*vy+vz*vz)
            uk=a_init*a_init*entropy_unit/gm1  
            uc(5,ix,iy,iz)=den*(uk+vk)
            uc(6,ix,iy,iz)=1./den**(gm1-1.)
         enddo
      enddo
   enddo
   
   call data_output(.true.)

else if (i_run==1) then

   open (2,file=datadir(1:lblnk(datadir))//'ucdump.dat',status='old',form='unformatted')
   read (2) tnum,nt,ntot,i_out
   read (2) l_unit,t_unit,v_unit,rho_unit,n_unit,e_rate_unit,tem_unit,cooling_coff,entropy_unit
   do ic=1,mnp
      read (2) d3
      uc(ic,1:Lx,1:Ly,1:Lz)=d3(1:Lx,1:Ly,1:Lz)
   enddo
   read (2) r
   read (2) v
   close(2)

   call redshift_comovingtime(tnum,zf)

endif
ntot=9500
kd1=0
kd2=0

if(zf>z_reion) then
    xJ=0.0 
    call cooling_tab(xJ,zf)
else
    xJ=xJ21
    call cooling_tab(xJ,zf)

endif


allocate(rs(3,nobj),vs(3,nobj),STAT=istat)

rs = 0.
vs = 0.
rhs = 0.

call head_logfile(i_run)

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

end subroutine head_logfile

! -----------------------------------------------------------------------------

subroutine data_output(write_disk)
real,allocatable,dimension(:,:,:)::ut
character ucout*11,sts*8
logical write_disk

allocate(ut(Lx,Ly,Lz),STAT=istat)

!     -------------------------------------------------------------------------
!     Name:      data_output.f
!     Function:  save the solution
!     -------------------------------------------------------------------------


if (write_disk) then 
   write(ucout,'(''uc'',i4.4,''.'',i4.4)')ist,nt
   sts='unknown'
   write(*,'(''save the numerical solution'',''step no='',i4)')nt
else     
   ucout='ucdump.dat'
   sts='unknown'
endif 
     
open(10,file=datadir(1:lblnk(datadir))//ucout,status=sts,form='unformatted')
write(10) tnum,nt,ntot,i_out
write(10) l_unit,t_unit,v_unit,rho_unit,n_unit,e_rate_unit,tem_unit,cooling_coff,entropy_unit    
do ic=1,mnp
   ut(1:Lx,1:Ly,1:Lz)=uc(ic,1:Lx,1:Ly,1:Lz)
   write (10) ut
enddo
write(10) r
write(10) v
close(10)

deallocate(ut,STAT=istat)

return
end subroutine data_output

! ----------------------------- write log file and output on screen -------------------------------

subroutine write_log

call redshift_comovingtime(tnum,znt)
af0 = 1./(1+znt) 
cf0 = 1./af0 
Ns_left = int((tend-tnum)/dt)

call CDM_velocity_stat

v_phys = cf0 * v_unit
vx_max = vx_max * v_phys
vy_max = vy_max * v_phys
vz_max = vz_max * v_phys
cs_max = cs_max * v_phys
vel_sigma = vel_sigma * v_phys
vxmax_cdm = vxmax_cdm * v_phys
vymax_cdm = vymax_cdm * v_phys
vzmax_cdm = vzmax_cdm * v_phys
vsigma_cdm = vsigma_cdm * v_phys

write(*,'('' step no='',i4,1x,''redshift='',f8.4,1x,&
& '' expansion factor = '', f8.4,1x,''time = '',f8.4)') nt, znt, af0, tnum
write(*,'('' cfl ='',f8.5,1x,''/ div-vel ='',f8.5,1x,&
& ''/ exp ='',f8.5,1x ''/ cdm-vel ='',f8.5)') dt_cfl,dt_div,dt_exp,dt_cdm
write(*,'(''left time step No. = '', I8)')Ns_left

write(*,'(''avarage density = '', f15.12,1x,''negative cells # = '', i8)') den_av,ncell_zero
write(*,'(''minimum/maximum density = '',f12.6,1x,f12.6)')den_min,den_max
write(*,'(''rms density ='',f12.8)') den_sigma
write(*,'(''maximum sound velocity = '', f12.6)') cs_max
write(*,'(''maximum 3d gas velocity = '', 3f12.6)') vx_max,vy_max,vz_max
write(*,'(''maximum 3d cdm velocity = '', 3f12.6)') vxmax_cdm,vymax_cdm,vzmax_cdm
write(*,'(''velocity dispersion gas/cdm = '', 2(1x,f12.8))') vel_sigma, vsigma_cdm
write(*,'('' '')')
write(*,'(''---------------------------------------------------------'')')
write(*,'('' '')')

write(68,'('' step no='',i4,2x,''redshift='',f8.4,2x, &
&  '' expansion factor = '', f8.4,1x,''time = '',f8.4)') nt, znt, af0, tnum
write(68,'('' cfl ='',f8.5,1x,''/ div-vel ='',f8.5,1x,&
& ''/ exp ='',f8.5,1x ''/ cdm-vel ='',f8.5)') dt_cfl,dt_div,dt_exp,dt_cdm
write(68,'(''CPU time in weno/grav/step'',3(1x,f12.6))')cpu_weno,cpu_grav,cpu_step
write(68,'(''avarage density = '', f15.12)') den_av
write(68,'(''negative cell # = '', i10)') ncell_zero
write(68,'(''minimum/maximum density = '',f12.6,1x,f12.6)')den_min,den_max
write(68,'(''rms density ='',f12.8)') den_sigma
write(68,'(''maximum sound velocity = '', f12.6)') cs_max
write(68,'(''maximum 3d gas velocity = '', 3f12.6)') vx_max,vy_max,vz_max
write(68,'(''maximum 3d cdm velocity = '', 3f12.6)') vxmax_cdm,vymax_cdm,vzmax_cdm
write(68,'(''velocity dispersion gas/cdm = '', 2(1x,f12.8))') vel_sigma, vsigma_cdm
write(68,'(''--------------------------------------------------------'')')

end subroutine write_log

end module readin_out

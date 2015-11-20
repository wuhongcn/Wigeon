module IGM_heating_cooling
use datablock
use recipes

contains

subroutine extQ(zf,xJ,fq,Q)
real,dimension(5)::fq(5)
real ::zf,n_eff,den,cden,vx,vy,vz,eng,vk,pres,den_H,Q,xJ,tem_H,CH_rate

if (zf>z_reion) then
   n_eff = x_hydrogen+y_helium/4.
else
   n_eff = 2.*x_hydrogen +3. * y_helium /4.
endif

den=fq(1)
cden=1./den
vx=fq(2)*cden
vy=fq(3)*cden
vz=fq(4)*cden
eng=fq(5)
vk=0.5*(vx*vx+vy*vy+vz*vz)
pres=gm1*(eng*cden-vk)*cf*cf
!tem_H = pres * T_unit * n_eff
tem_H = pres * tem_unit / n_eff  ! 01/24/2009 AZ
den_H = n_unit * den * cf * cf * cf 

call GasModeling(zf,xJ,tem_H,den_H,CH_rate)  

Q = cooling_coff * af * den * den * CH_rate
return
end subroutine extQ

! -------------------------------------------------------------
subroutine cooling_tab(xJ,z_redshift)
real, parameter::uv_index=1. ! pow index of the UV spectrum 
real :: CR,z_redshift,xJ,temlog_min,tem_i,den_j
!data tem1,tem2,temlog_max,denlog_min,denlog_max/3.95,0.,8.,-8.,2.0/
data tem1,tem2,temlog_max,denlog_min,denlog_max/3.95,1.,7.70,-5.05,2.55/

!real,dimension(Ns)::tem_table,den_table
!real,dimension(Ns,Ns)::ch_table

!! above vectors are defined in weno_module.f90  wszhu 2009 Tucson   !!

!integer :: kd1,kd2
!save tem_table,den_table,ch_table,kd1,kd2
!data kd1,kd2/0,0/

if (xJ==0.) temlog_min=tem1
if (xJ>0.) temlog_min=tem2
!if ((kd1==0.and.xJ==0.).or.(kd2==0.and.xJ>0.)) then  
   do j=1,Ns_table
         den_table(j)=denlog_min+(denlog_max-denlog_min)*float(j)/float(Ns_table)
         den_j=10.**den_table(j)
      do i=1,Ns_table
         tem_table(i)=temlog_min+(temlog_max-temlog_min)*float(i)/float(Ns_table)
         tem_i=10.**tem_table(i)

         call IGMprocess(z_redshift,tem_i,den_j,xJ,uv_index,CR)         
         ch_table(i,j)=CR

      enddo
   enddo
   if (kd1==0.and.xJ==0.) then
      write(*,'(''generate cooling-heating table without UV '')')
      kd1=1
   endif
   if (kd2==0.and.xJ>0.) then
      write(*,'(''generate cooling-heating table with UV'')')    
      kd2=1
   endif
!endif

end subroutine cooling_tab
!	 --------------------------------------------------------------------------


subroutine GasModeling(z_redshift,xJ,tem_igm,den_igm,CHrate)
integer,parameter::Nx=4
real :: tem_igm,den_igm,CR,CHrate,z_redshit,xJ,temlog_min,tem10,den10
!data tem1,tem2,temlog_max,denlog_min,denlog_max/3.95,0.,8.,-8.,2.0/
data tem1,tem2,temlog_max,denlog_min,denlog_max/3.95,1.,7.70,-5.05,2.55/
integer::i_tem,j_den,io

real,dimension(Nx+1)::temx,denx
real,dimension(Nx+1,Nx+1)::chx


if (xJ==0.) temlog_min=tem1
if (xJ>0.) temlog_min=tem2

!if (den_igm<=0..or.tem_igm<=temlog_min) then
if (den_igm<=1.0e-5.or.den_igm>10.**2.55.or.tem_igm<=10.**temlog_min.or.tem_igm>=10.**temlog_max) then
   CHrate=0.
   return
endif

tem10=log10(tem_igm)
den10=log10(den_igm)


i_tem=int(real(Ns_table)*(tem10-temlog_min)/(temlog_max-temlog_min))
j_den=int(real(Ns_table)*((den10-denlog_min)/(denlog_max-denlog_min)))


Nx2=Nx/2
im=i_tem-Nx2
ip=i_tem+Nx2
jm=j_den-Nx2
jp=j_den+Nx2

if (im>=1.and.ip<=Ns_table.and.jm>=1.and.jp<=Ns_table) then
   temx(1:Nx+1)=tem_table(im:ip)
   denx(1:Nx+1)=den_table(jm:jp)
   chx(1:Nx+1,1:Nx+1)=ch_table(im:ip,jm:jp)
   call polin2d(temx,denx,chx,Nx+1,Nx+1,tem10,den10,CR,df)
   CHrate=CR
else
   CHrate=0.
endif

return
end subroutine GasModeling


! ---------------------------- IGMprocess -------------------------------------

subroutine IGMprocess(z_redshift,T,den_H,xJ,uv_index,CR)
real :: alpha_H1,alpha_He1,alpha_He2,alpha_d
real :: gamma_eH0,gamma_eHe0,gamma_eHe1,gamma_ph0,gamma_pHe0,gamma_pHe1
real :: xH0,xH1,xHe0,xHe1,xHe2,xe
real :: z_redshift,T,den_H,xJ,uv_index,CR


! computing the ionizing cofficients
call ionizing_coff

! computing the fraction of each component in den_H
call fraction_den

! computing the net cooling rate 
call CoolingRate

return

contains

! ---------------------------- sub-subroutine ---------------------------------

! computing the fraction of each components 
subroutine fraction_den
!real, parameter::eps=1.E-10
real, parameter::eps=1.E-6

y = 0.25 * y_helium / x_hydrogen

!write(*,*) 'y',y,'y_helium',y_helium,'x_hydrogen',x_hydrogen
!pause
if (xJ == 0.) then

   xH0 = alpha_H1 / ( alpha_H1 + gamma_eH0 )
   xH1 = 1. - xH0
  if(gamma_eHe0>0.) then
     xHe1 = y /(1. + ( alpha_He1 + alpha_d ) /gamma_eHe0 + gamma_eHe1 /alpha_He2 )
     xHe0 = xHe1 * ( alpha_He1 + alpha_d) /gamma_eHe0
     xHe2 = xHe1 * gamma_eHe1 /alpha_He2
  else
     xHe0=y
     xHe1=0.
     xHe2=0.
  endif
   xe = xH1 + xHe1 + 2. * xHe2

else if (xJ > 0.) then 
   xe0 = 1.D0
10 xH0 = alpha_H1 /( alpha_H1 + gamma_eH0 + gamma_pH0 /xe0 /den_H )
   xH1 = 1. - xH0
   xHe1 = y /( 1.+ ( alpha_He2 + alpha_d ) / ( gamma_eHe0 + gamma_pHe0 /xe0 /den_H ) &
&      + ( gamma_eHe1 + gamma_pHe1 /xe0 /den_H ) /alpha_He2 )
   xHe0 = xHe1 * ( alpha_He1 + alpha_d ) /( gamma_eHe0 + gamma_pHe0 /xe0 /den_H )
   xHe2 = xHe1 * ( gamma_eHe1 + gamma_pHe1 /xe0 /den_H ) /alpha_He2
   xe = xH1 + xHe1 + 2. * xHe2

   if (abs(xe-xe0)/xe0<eps) then
       return
   else
       xe0 = xe
       goto 10
   endif

endif


return
end subroutine fraction_den

! ionizing rate in unit of 10^-10 cm^3 s^-1  
! coeffiecient changed according to Theuns et al. 1998  wszhu /02/03/2009/ UA 
subroutine ionizing_coff
real:: T3,T5,T6,Ts,T5s

T3=T/1000.
T5=T3/100.
T6=T5/10.
Ts=sqrt(T)
T5s=sqrt(T5)

alpha_H1 = 0.63 /Ts * T3**(-0.2) /(1.+T6**0.7)
alpha_He1 = 1.5 * T**(-0.6353)
alpha_He2 = 3.36 /Ts * T3**(-0.2) /(1.+ T6**0.7)
alpha_d = (1.9E+7) * T**(-1.5) * exp(-4.7/T5) * (1.+0.3 * exp(-0.94/T5))

gamma_eH0 = 1.17 * Ts * exp(-157809.1/T)/(1.+T5s)
gamma_eHe0 = 0.476 * Ts * exp(-285335.4/T)/(1.+T5s)
gamma_eHe1 = 0.114 * Ts * exp(-631515.0/T)/(1.+T5s)

if (xJ>0.) then
   xn = uv_index
   gamma_pH0 = 0.126* xJ /(3.+xn)
   gamma_pHe0 = 0.148 * xJ * 0.553**xn * (1.66/(xn+2.05)-0.66/(xn+3.05))
   gamma_pHe1 = 0.0334 * xJ * 0.249**xn /(3.+xn)
endif


return
end subroutine ionizing_coff

! cooling rate in unit of 10^-21 ergs cm^-3 s^-1
subroutine CoolingRate
real :: T3,T5,T6,Ts,T5s,TT
real :: CIC_H0,CIC_He0,CIC_He1,CIC_He2S,DRC_He1,FFC_ALL
real :: RCC_H0,RCC_He1,RCC_He2,CEC_H0,CEC_He1,CMB

T3=T/1000.
T5=T3/100.
T6=T5/10.
Ts=sqrt(T)
T5s=sqrt(T5)
TT=Ts/(1.+T5s)

! collisional ionization
CIC_H0 = 2.54 * TT * exp(-157809.1/T) * xe * xH0
CIC_He0 = 1.88 * TT * exp(-285335.4/T) * xe * xHe0
CIC_He1 = 0.99 * TT * exp(-631515.0/T) * xe * xHe1
CIC_He2S = (5.01D-6) * T**(-0.1687) / (1.+T5s) * exp(-55338.0/T) * xe * xe * xHe1 * den_H
CIC = CIC_H0 + CIC_He0 + CIC_He1 + CIC_He2S


! recombination
RCC_H0 = (8.70D-6) * Ts * T3**(-0.2) / (1.+T6**0.7) * xe * xH1
RCC_He1 = (1.55D-5) * T**0.3647 * xe * xHe1
RCC_He2 = (3.48D-5) * Ts * T3**(-0.2) / (1.+T6**0.7) * xe * xHe2
RCC = RCC_H0 + RCC_He1 + RCC_He2


! collisional excitation
CEC_H0 = 750. * exp(-118348.0/T)/(1.+T5s)*xe*xH0
CEC_He1 = (5.54D+4) * T**(-0.397) * exp(-473638.0/T) / (1.+T5s) * xe * xHe1
CEC = CEC_H0 + CEC_He1


! dielectric recombination
DRC_He1 = (1.24D+8) * T**(-1.5) * exp(-4.7/T5) * (1.+0.3 * exp(-0.94/T5)) * xe * xHe1 
DRC = DRC_He1


! free-free emission 
gff = 1.1 + 0.34 * exp(-(5.5-log10(T))**2/3.0)
FFC_ALL=(1.42E-6) * gff * Ts * ( xH1 + xHe1 + 4. * xHe2 ) * xe
FFC = FFC_ALL


! CMB + electron Compton scattering
CMB = (5.406D-10) * (T5-2.7E-5*(1. + z_redshift))*(1.+ z_redshift)**4 * xe / den_H

! total cooling rate
CR = CIC + RCC + CEC + DRC + FFC + CMB


! UV vackground heating
if (xJ > 0.) then
   xn = uv_index
   Heat_H0 = 0.291 * xJ /(2.+xn) /(3.+xn)
   Heat_He0 = 0.584 * xJ * 0.553**xn * (1.66/(xn+1.05)-2.32/(xn+2.05)+0.66/(xn+3.05))
   Heat_He1 = 0.292 * xJ * 0.249**xn /(2.+xn) /(3.+xn)
   HR = (xH0 * Heat_H0 + xHe1 * Heat_He1 + xHe0 * Heat_He0) /den_H
endif

! the net cooling rate
CR = CR - HR


return
end subroutine CoolingRate 

end subroutine IGMprocess

end module IGM_heating_cooling

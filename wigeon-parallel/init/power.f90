! normalized linear power spectrum at redshift z

function CDMpower(xk,z)
use set_parameter
real,external:: rint
common/tophat/rth
data id/0/
save id

rth=8.0
boxsize=rLbox

if (id==0) then

   Dz=growth(z )
   D0=growth(0.)
write(*,*) 'Dz',Dz,'D0',D0
!pause
   if (ksp==1) then

!	  sigma80=0.55*omega0**(-0.6)
!      sigma80=0.9
      sigma80=0.812  !wmap 5
      sigma8=sigma80*Dz/D0
      delta=sigma8/(pnorm(rth)*boxsize**1.5)
   else if (ksp==2) then
      xn1=ind-1.
      if (xlambda0.eq.0) then
         deltaH=(1.95e-5)*omega0**(-0.35-0.19*alog(omega0)-0.17*xn1)*exp(-xn1-0.14*xn1*xn1)
      else if(xlambda0+omega0.eq.1.) then
         deltaH=(1.94e-5)*omega0**(-0.785-0.05*alog(omega0))*exp(-0.95*xn1-0.169*xn1*xn1)
      endif

      CALL qromb(rint,0.,100./rth,sig2)
      sigma80=deltaH*(3000.)**(1.5+ind/2.)*sqrt(sig2)

      deltabox=(sqrt(2.*pi*pi))*((3000./boxsize)**1.5)*sqrt(3000**ind)
      delta=deltabox*deltaH
      delta=delta*(Dz/D0)
   else if (ksp==3) then 
      fb=0.25
      delta=sqrt(4.*pi*fb*boxsize)/boxsize**1.5
      delta=delta*Dz/D0
   endif

!   write (*,'(''amplitude of power spectrum ='',f15.8)') delta
!   write (*,'(''sigma8 = ''f15.8)') sigma80
!   write (*,'(''bias factor = '',f15.8)') 1./sigma80

   id=1

endif

CDMpower=delta**2*pow(xk)

end function CDMpower
      
function pow(rk)
use set_parameter

omegah=omega0*h100*exp(-omegab*(1.+sqrt(2.*h100)/omega0))
omegahh=omegah

if (rk.le.0.0) then
   pow=0.0
   RETURN
end if

if (ispec==1) then  !  Power Law
   q=rk
   pow=q**ind
else if (ispec==2) then ! Transfer function for CDM (from BBKS)
   q=rk/omegahh
   a1=2.34*q
   a2=3.89*q
   a3=16.1*q
   a4=5.46*q
   a5=6.71*q
   t=1.0+a2+a3*a3+a4*a4*a4+a5*a5*a5*a5
   t=log(1.0+a1)/a1/sqrt(sqrt(t))
   pow=q**ind*t**2
else if (ispec==3) then ! Transfer function for CDM (from Hu)
   fb=omegab/omega0
   call Tifit(rk,omega0,fb,h100,2.73,tf_full,tf_baryon,tf_cdm)      
   t=tf_full
   pow=(rk**ind)*(t**2)
else if (ispec==4) then ! Transfer function for HDM (from BBKS).
   q=rk/omegahh
   a1=2.6*q
   a2=1.6*q
   a3=4.0*q
   a4=0.92*q
   t=exp(-0.16*a1-0.5*a1*a1)/(1.0+a2+a3*sqrt(a3)+a4*a4)
   pow=q*t**2
end if

return
end function pow

subroutine Tifit(xk,omega_total,f_baryon,hubble,Tcmb,tf_full,tf_baryon,tf_cdm)      
real  omega_total,f_baryon,hubble,Tcmb
real  tf_full,tf_baryon,tf_cdm
  
 
!  cosmological parameters

omhh = omega_total*hubble*hubble

! call routine to set fitting parameters

call TFset_parameters(omhh, f_baryon,hubble, Tcmb)
call TFtransfer_function(xk*hubble,omhh,f_baryon,tf_full,tf_baryon,tf_cdm)

return
end subroutine Tifit

! PART I:------------------- FITTING FORMULAE ROUTINES ----------------- 

!   There are two routines and a set of functions.  
!   TFset_parameters() sets all the scalar parameters, while 
!   TFtransfer_function() calculates various transfer functions 

! Global variables -- We've left many of the intermediate results as
! global variables in case you wish to access them, e.g. by declaring
! them as a common block in your main program. 

! Note that all internal scales are in Mpc, without any Hubble constants! 

subroutine TFset_parameters(omhh,f_baryon,hubble,Tcmb)

real y,omhh,obhh,Tcmb
real theta_cmb,z_equality,k_equality,z_drag,R_drag,R_equality,&
& sound_horizon,k_silk,alpha_c,beta_c,alpha_b,beta_b,f_baryon,beta_node
common/GLOBALVARIABLES/theta_cmb,z_equality,k_equality,z_drag,&
& R_drag,R_equality,sound_horizon,k_silk,alpha_c,beta_c,alpha_b,beta_b,beta_node

! Set all the scalars quantities for Eisenstein & Hu 1997 fitting formula */
! Input omhh -- The density of CDM and baryons, in units of critical dens,
!                multiplied by the square of the Hubble constant, in units
!                of 100 km/s/Mpc */
!       f_baryon -- The fraction of baryons to CDM */  ! should be baryons to
!       Tcmb -- The temperature of the CMB in Kelvin, 2.728(4) is COBE and is
!		the default reached by inputing Tcmb=0 -- reset on output. */
! Output nothing, but set many global variables in common block 
!       GLOBALVARIABLES. You can access them yourself, if you want:

!	theta_cmb,	/* Tcmb in units of 2.7 K */ 
!	z_equality,	/* Redshift of matter-radiation equality, really 1+z */
!	k_equality,	/* Scale of equality, in Mpc^-1 */
!	z_drag,		/* Redshift of drag epoch */
!	R_drag,		/* Photon-baryon ratio at drag epoch */
!	R_equality,	/* Photon-baryon ratio at equality epoch */
!	sound_horizon,	/* Sound horizon at drag epoch, in Mpc */
!	k_silk,		/* Silk damping scale, in Mpc^-1 */
!	alpha_c,	/* CDM suppression */
!	beta_c,		/* CDM log shift */
!	alpha_b,	/* Baryon suppression */
!	beta_b,		/* Baryon envelope shift */

! Are inputs reasonable?

if (f_baryon.le.0) f_baryon=1.e-5
if (Tcmb.le.0) Tcmb=2.728
if (omhh.le.0.0) then
   write(6,*) 'TFset_parameters(): Illegal input'  
   pause
end if

if (hubble.gt.10.0) then
   write(6,*) 'TFset_parameters(): WARNING, Hubble constant in 100km/s/Mpc desired'
end if

! Auxiliary variables
obhh = omhh*f_baryon
theta_cmb = Tcmb/2.7

! Main variables
z_equality = 2.50e4*omhh*theta_cmb**(-4.) - 1.D0
k_equality = 0.0746*omhh*theta_cmb**(-2.) 

z_drag = 0.313*omhh**(-0.419)*(1.+0.607*omhh**(0.674))
z_drag = 1e0 + z_drag*obhh**(0.238*omhh**(0.223))
z_drag = 1291e0 * omhh**(0.251)/(1e0 + 0.659*omhh**(0.828)) * z_drag
 
R_drag = 31.5*obhh*theta_cmb**(-4.)*1000e0 / (1e0 + z_drag) 
R_equality = 31.5*obhh*theta_cmb**(-4.) * 1000e0/(1e0 + z_equality) 

sound_horizon = 2./3./k_equality*sqrt(6./R_equality)* &
& log(( sqrt(1.+R_drag)+sqrt(R_drag+R_equality) )/(1.+sqrt(R_equality)))

k_silk = 1.6*obhh**(0.52)*omhh**(0.73)*(1e0 + (10.4*omhh)**(-0.95))

alpha_c = ((46.9*omhh)**(0.670)*(1e0+(32.1*omhh)**(-0.532)))
alpha_c = alpha_c**(-f_baryon) 
alpha_c = alpha_c*((12.0*omhh)**(0.424)*(1e0 + (45.0*omhh)**(-0.582)))**(-f_baryon**3.)

beta_c = 0.944/(1+(458.*omhh)**(-0.708))
beta_c = 1.+beta_c*((1.-f_baryon)**((0.395*omhh)**(-0.0266)) - 1e0)
beta_c = 1./beta_c

y = (1e0+z_equality)/(1e0+z_drag)
alpha_b = y*(-6.*sqrt(1.+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)))
alpha_b = 2.07*k_equality*sound_horizon*(1.+R_drag)**(-0.75)*alpha_b

beta_b = 0.5+f_baryon+(3.-2.*f_baryon)*sqrt((17.2*omhh)**2.+1e0)

beta_node = 8.41*omhh**(0.435)

return

end subroutine TFset_parameters

subroutine TFtransfer_function(k,omhh,f_baryon,tf_full,tf_baryon,tf_cdm)

!  Calculate transfer function from the fitting parameters stored in GLOBALVARIABLES.

!  Input: 
!	 k -- wavenumber in Mpc^{-1}  
!        omhh -- The density of CDM and baryons, in units of critical dens,
!                multiplied by the square of the Hubble constant, in units
!                of 100 km/s/Mpc */
!        f_baryon -- The fraction of baryons to CDM */

!  Output:
!	 tf_full -- The full fitting formula, eq. (16), for the matter
!	            transfer function. 
!	 tf_baryon -- The baryonic piece of the full fitting formula, eq. 21.
!	 tf_cdm -- The CDM piece of the full fitting formula, eq. 17.

real k,tf_full,tf_baryon,tf_cdm,q,ks
real theta_cmb,z_equality,k_equality,z_drag,R_drag,R_equality,&
& sound_horizon,k_silk,alpha_c,beta_c,alpha_b,beta_b,f_baryon,beta_node
common/GLOBALVARIABLES/theta_cmb,z_equality,k_equality,z_drag,&
& R_drag,R_equality,sound_horizon,k_silk,alpha_c,beta_c,alpha_b,beta_b,beta_node

!  Reasonable k?

if (k.le.0) then
   write(6,*) 'TFtransfer_function(): Illegal k'
   pause
end if 

!  Auxiliary Variables

q = k/13.41/k_equality
ks = k*sound_horizon

! Main Variables

tf_cdm = 1./(1.+(ks/5.4)**4.)
tf_cdm = tf_cdm*TF_pressureless(q,1.,beta_c) + &
& (1.-tf_cdm)*TF_pressureless(q,alpha_c,beta_c)

s_tilde = sound_horizon/(1.+(beta_node/ks)**3.)**(1./3.) 
tf_baryon = TF_pressureless(q,1.,1.)/(1.+(ks/5.2)**2.)
tf_baryon = tf_baryon + alpha_b/(1.+(beta_b/ks)**3)*exp(-(k/k_silk)**(1.4))
tf_baryon = tf_baryon *(sin(k*s_tilde)/(k*s_tilde))
tf_full = f_baryon*tf_baryon + (1-f_baryon)*tf_cdm

return
end subroutine TFtransfer_function

! auxiliary function: Pressureless TF

real function TF_pressureless(q,a,b)
real q,a,b

TF_pressureless = Log(exp(1.)+1.8*b*q)
TF_pressureless = TF_pressureless/(TF_pressureless + &
& (14.2/a + 386/(1.+69.9*q**1.08))*q**2)

return
end function TF_pressureless


! PART II:------------------- Scaling Functions ROUTINES ----------------- 
!
!       omhh -- The density of CDM and baryons, in units of critical dens,
!                multiplied by the square of the Hubble constant, in units
!                of 100 km/s/Mpc */
!       f_baryon -- The fraction of baryons to CDM */
!
!
!	TF_zerobaryon:     
!	  Input:  q = k/omhh * (Tcmb/2.7)**2    (k in Mpc^{-1})
!	  Output: zero baryon TF Eq(29)
!	TF_nowiggles:      
!	  Input:  k = wavenumber in Mpc^{-1}, omhh, f_baryon, Tcmb
!	  Output: shape approximation TF  Eq(30-31)
!	  Calls: TF_zerobaryon,sound_horizon_fit,alpha_gamma
! 	sound_horizon_fit: 
!         Input:  omhh,f_baryon	
!	  Output: approximate sound horizon in Mpc	
!	kpeak:		   
!	  Input:  omhh,f_baryon
!         Output: first peak location in Mpc
!	  Calls:  sound_horizon_fit
!	alpha_gamma:	   
!	  Input: omhh,f_baryon
!	  Output: effective small scale suppression

real function TF_zerobaryon(q)
real q

TF_zerobaryon = log(2.0*exp(1.)+1.8*q)
TF_zerobaryon = TF_zerobaryon/(TF_zerobaryon+(14.2 + 731.0/(1+62.5*q))*q**2)

return
end function TF_zerobaryon

real function TF_nowiggles(k,omhh,f_baryon,Tcmb)
real k,omhh,f_baryon,q_eff,a

if (Tcmb.le.0) Tcmb=2.728
a = alpha_gamma(omhh,f_baryon)
q_eff = k/omhh*(Tcmb/2.7)**2
q_eff = q_eff/(a+(1.-a)/(1.+(0.43*k*sound_horizon_fit(omhh,f_baryon))**4))

TF_nowiggles = TF_zerobaryon(q_eff)

return
end function TF_nowiggles

real function sound_horizon_fit(omhh,f_baryon)
real omhh,obhh,f_baryon

obhh = f_baryon*omhh
sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1.+10.0*obhh**(0.75))

return
end function sound_horizon_fit

real function k_peak(omhh,f_baryon)
real omhh,obhh,f_baryon

obhh = f_baryon*omhh
k_peak = 5.*3.14159/2.*(1.+0.217*omhh)/sound_horizon_fit(omhh,f_baryon)

return
end function k_peak

real function alpha_gamma(omhh,f_baryon)
real omhh,f_baryon

alpha_gamma = 1.-0.328*log(431.0*omhh)*f_baryon + 0.38*log(22.3*omhh)*(f_baryon)**2
    
return
end function alpha_gamma

! -------------------------------------------------------------------------------

function pnorm(r)
real,parameter::pi=3.14159265358979,tpi=2.*pi
real,external::rint
common/tophat/rth

rth=r
CALL qromb(rint,0.,100./rth,sig2)
sig2=4.*pi/tpi**3*sig2
pnorm=sqrt(sig2)

return
end function pnorm

function rint(x)
common/tophat/rth

rint=w(rth*x)*pow(x)*x**2

return
end

function w(x)

if (x.le.1e-2) then
   w=1.-0.2*x**2
else
   w=9.*(sin(x)-x*cos(x))**2/x**6
end if

return
end function w

subroutine variance(z,n,Rg,var)
integer,parameter::ns=256
real,dimension(ns)::xs,ws
external CDMpower

x1=0.
x2=sqrt(15.*log(10.))
call gauleg(x1,x2,xs,ws,ns)

inn=-2*(n+1)
pi=acos(-1.)
tpi=2.*pi

Si=0.
do i=1,ns
   xi=xs(i)
   wi=ws(i)
   Si=Si+wi*CDMpower(xi,z)*xi**inn*exp(-xi*xi*Rg*Rg)
enddo
var=Si/tpi/pi

end subroutine variance



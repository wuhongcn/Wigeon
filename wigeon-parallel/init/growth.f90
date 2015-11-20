Function YH(x)
use set_parameter
save id
data id/0/

if (id==0) then
   a0=1.+zinit
   Gn=(omega0*a0+1.-xlambda0-omega0)*a0**2+xlambda0
   omega_i=omega0*a0**3/Gn
   xlambda_i=xlambda0/Gn
   id=1
endif

YH=(x/(x+omega_i*(1.-x)+xlambda_i*x*(x**2-1)))**1.5

return
end

Function Growth(z)
use set_parameter
external YH

a=(1.+zinit)/(1.+z)
call qromb(YH,0.,a,ss)

Ha=sqrt(1.+omega_i*(1./a-1.)+xlambda_i*(a*a-1.))/a
Growth=Ha*ss

return
end

subroutine cofficient_velocity(hst)
use set_parameter
external YH

call qromb(YH,0.,1.,ss)

!write(*,*) alpha,beta,ss

beta=xlambda_i-omega_i/2.-1.+1./ss

write(*,*) 'alpha',alpha,'beta',beta,'ss',ss
hst=beta/alpha

return
end

function E(z)
use set_parameter, only:omega0,xlambda0

zp=1.+z
omegak0=1.-omega0-xlambda0
E=sqrt((omega0*zp+omegak0)*zp*zp+xlambda0)

return
end

SUBROUTINE timeah(t,omega0,xlambda0,a,h)

!  Evaluates the time t (in units of t0) and Hubble parameter (units 1/t0)
!  given the expansion factor a (=1 at t=1)
!  for cosmological parameters omega0 and lambda0 (evaluated at t=1).

h0t0=aintegral(0.,1.,omega0,xlambda0)
h0t=aintegral(0.,a,omega0,xlambda0)
t=h0t/h0t0
h=h0t0*sqrt(1.+omega0*(1./a-1.)+xlambda0*(a**2-1.))/a
aold=a

RETURN
end

! -------------------------------------------------------------------
real FUNCTION aintegral(a1,a2,omega,xlambda)

!  Does open-ended romberg integration of the a integral from a1 to a2
!  Based on a Numerical Recipes routine

parameter (eps=1.e-6, jmax=14, jmaxp=jmax+1, k=5, km=k-1)
real h(jmaxp),s(jmaxp)

arange=a2-a1
h(1)=1.
do j=1,jmax
!  midpoint evaluation
   if (j.eq.1) then
      x=0.5*(a1+a2)
      s(j)=arange*func(x,omega,xlambda)
      it=1
   else
      tnm=it
      del=arange/(3.*tnm)
      ddel=2.*del
      x=a1+0.5*del
      sum=0.
      do jt=1,it
         sum=sum+func(x,omega,xlambda)
         x=x+ddel
         sum=sum+func(x,omega,xlambda)
         x=x+del
      end do
      s(j)=(s(j-1)+arange*sum/tnm)/3.
      it=3*it
   end if
!  extrapolate to zero step size
   if (j.ge.k) then
      CALL polint(h(j-km),s(j-km),k,0.,ss,dss)
      if (abs(dss).le.eps*abs(ss)) then
         aintegral=ss
         RETURN
      end if
   endif
   s(j+1)=s(j)
   h(j+1)=h(j)/9.
end do

STOP 'aintegral: too many steps'

contains

function func(a,om,xla)

func = 1/sqrt(1.+om*(1./a-1.)+xla*(a**2-1.))

end function func

end FUNCTION aintegral

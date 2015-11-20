module cosmology 
use recipes

use datablock, only: omega0, xlambda0, zinit

contains 

! cosmic time t -> redshift z 
subroutine redshift_comovingtime(t,z)
real, parameter:: xacc=1.e-10
! real, external:: Ft
common/cosmic_time/ts

if (t==0.) then
   z=zinit
   return
endif

ts=t
z_off = 0.5
z_min = -z_off
z_max = zinit+z_off

z=zbrent(Ft,z_min,z_max,xacc)

return
end subroutine redshift_comovingtime

!redshift -> cosmic time t
real function age(z)
! real,external:: Ta

CALL qromb(Ta,z,zinit,Ez)
age=Ez

return
end function age

real function Ta(z)

Ta=(1.+z)/E(z)

end function Ta

function Ft(x)
common/cosmic_time/ts
Ft=age(x)-ts
return
end function Ft

function E(z)

zp=1.+z
omegak0=1.-omega0-xlambda0
E=sqrt((omega0*zp+omegak0)*zp*zp+xlambda0)

return
end function

end module cosmology 

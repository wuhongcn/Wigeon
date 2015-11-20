module hydro_stat
use datablock
use gravity
use IGM_heating_cooling

contains

! -----------------------------------------------------------------------------

subroutine cflc(aam)
real,dimension(5)::fq

! computed the largest eigenvalue for cfl

den_av = 0.
den_min = 1.0
den_max = 0.
den_sigma = 0.

vx_max = 0.
vy_max = 0.
vz_max = 0.
cs_max = 0.
vel_sigma = 0.

aam = 0.
ncell_zero = 0
fq=0.0
!dt_cool=1.0
!$omp parallel do default(shared) private(i,j,k,den) & 
!$ reduction(+:den_av,den_sigma,ncell_zero) reduction(min:den_min) reduction(max:den_max)
do k=1,Lz
  do j=1,Ly
    do i=1,Lx
       den=uc(1,i,j,k)
       if (den<1.0e-5) then
          ncell_zero=ncell_zero+1
          write(*,'(3i5,1x,''den='',f14.8,1x,''eng='',f14.5,1x,''ent='',f14.5)')i,j,k,den,uc(5,i,j,k),uc(6,i,j,k)
       endif
       den_min = min(den_min,den)
       den_max = max(den_max,den)
       den_sigma = den_sigma+(den  - 1.)**2
       den_av=den_av+den
    enddo
  enddo
enddo

!$omp end parallel do



!$omp parallel do default(shared)  private(i,j,k,den,eng,ent,cden) &
!$ private(sij,tij,rij,q2,eij,uij,pij,cij,c2,w2) &
!$ reduction(+:vel_sigma)  reduction(max:vx_max,vy_max,vz_max,cs_max,aam)
do k=1,Lz
   do j=1,Ly
      do i=1,Lx
         den = uc(1,i,j,k)
         eng = uc(5,i,j,k)
         ent = uc(6,i,j,k)
         cden = 1./den

         rij = uc(2,i,j,k) * cden
         sij = uc(3,i,j,k) * cden
         tij = uc(4,i,j,k) * cden
         q2 = rij**2 + sij**2 + tij**2
         eij = 0.5 * q2 * den
!         eij=0.5*q2*abs(den)
         pij = gm1 * ( eng - eij )
         c2 = gamma * pij * cden
!         c2=gamma*pij*abs(cden)
         cij = sqrt(abs(c2))
         w2 = abs(rij) + cij + abs(sij) + cij + abs(tij) + cij
!         w2=max(abs(rij) + cij,abs(sij) + cij,abs(tij) + cij)
! /2007/08/21/
         aam = max(aam, w2)
         vel_sigma = vel_sigma + q2

!         den_av = den_av + den
!         den_min = min(den_min,den)
!         den_max = max(den_max,den)
         vx_max = max(vx_max,abs(rij))
         vy_max = max(vy_max,abs(sij))
         vz_max = max(vz_max,abs(tij))
         cs_max = max(cs_max,cij)

!         if(den/=0.)  then
!           fq(1:5)=uc(1:5,i,j,k)
!           call extQ(zf,xJ,fq,Q)
!           dtcol=abs(uc(5,i,j,k)/Q)
!           dt_cool=min(dt_cool,dtcol)
!         endif  

      enddo
   enddo
enddo
!$omp end parallel do

vL = float(Lx*Ly*Lz)
den_av = den_av/vL
den_sigma = sqrt(den_sigma/vL)
vel_sigma = sqrt(vel_sigma/vL)


return
end subroutine cflc

! -----------------------------------------------------------------------------

subroutine regulation

Tmin=1.0
tem_min=Tmin/tem_unit/cf**2
Tmax=10.**7.85
tem_max=Tmax/tem_unit/cf**2
crit=0.001
!crit=0.02

!Need test num_threads() here. /01/24/2009 wszhu AZ

!$omp parallel do default(shared)  private(ix,iy,iz,ixp1,iyp1,izp1,&
!$ ixm1,iym1,izm1,uci,den,vxm,vym,vzm,eng,ent,cden,vx,vy,vz,vk,&
!$ tem,tem_eng,tem_ent,pre,pre_eng,pre_ent) 
do iz=1,Lz
   do iy=1,Ly
      do ix=1,Lx 

!         den=abs(uc(1,ix,iy,iz))
         den=uc(1,ix,iy,iz)
         vxm=uc(2,ix,iy,iz)
         vym=uc(3,ix,iy,iz)
         vzm=uc(4,ix,iy,iz)
         eng=uc(5,ix,iy,iz)
         ent=uc(6,ix,iy,iz)
         cden=1./den
         vx=vxm*cden
         vy=vym*cden
         vz=vzm*cden
         vk=0.5*(vx*vx+vy*vy+vz*vz)

         tem_eng=gm1*(eng*cden-vk)
         
         if (tem_eng<tem_min) then
            tem_eng=tem_min
            pre = tem_eng * den
!            pre=tem_eng*abs(den)
            eng = pre/gm1 + den * vk
            uc(5,ix,iy,iz)=eng
!         else if (tem_eng>tem_max) then
!            tem_eng=tem_max
!            pre = tem_eng * den
!            eng = pre/gm1 + den * vk
!            uc(5,ix,iy,iz)=eng
         endif 

!         tem_ent=abs(den)/den*ent*entropy_unit*abs(den)**(gm1-1.)
         if (tem_ent<tem_min) then
            tem_ent=tem_min
!            ent=abs(den)/den*tem_ent/abs(den)**(gm1-1.)/entropy_unit
            ent=tem_ent/abs(den)**(gm1-1.)/entropy_unit
            uc(6,ix,iy,iz)=ent
!         else if (tem_ent>tem_max) then
!            tem_ent=tem_max
!            ent=tem_ent/abs(den)**(gm1-1.)/entropy_unit
!            uc(6,ix,iy,iz)=ent
         endif
         pre_eng=den*tem_eng/eng

         tem=tem_eng
         if (pre_eng<crit) then
             tem=tem_ent
        endif

        if(uc(5,ix,iy,iz)<0.) then 

            eng= den*(tem/gm1 + vk) 
            uc(5,ix,iy,iz)=eng
         endif

         ent=tem/abs(den)**(gm1-1.)/entropy_unit
         uc(6,ix,iy,iz)=ent
        enddo
    enddo
enddo
!$omp end parallel do

return
end subroutine regulation

subroutine velocity_div1d(div_max)
integer,dimension(3)::nn
!integer,dimension(Lx,Ly,Lz)::ind_div
real,allocatable,dimension(:,:,:)::v3,vd,div
real::den,cden,div_max
nn(1) = Lx
nn(2) = Ly
nn(3) = Lz

!!$omp parallel do default(shared) private(iz)
!do iz = 1, Lz
!   ind_div(:,:,iz) = 1
!enddo   
!!$omp end parallel do 

allocate (v3(Lx,Ly,Lz),STAT=istat)
allocate (vd(Lx,Ly,Lz),STAT=istat)
allocate (div(Lx,Ly,Lz),STAT=istat)

div = 0.
do ic = 1, 3
   if ( nn(ic) > 1 ) then 
!$omp parallel do default(shared) private(ix,iy,iz,den,cden)
   do iz = 1,Lz
      do iy = 1,Ly
         do ix = 1,Lx
            den = uc(1,ix,iy,iz)
            cden = 1./den
            v3(ix,iy,iz) = uc(ic+1,ix,iy,iz) * cden
         enddo
      enddo
   enddo
!$omp end parallel do   
   call gridf10(ic,v3,vd)
   endif
   div = div + vd
enddo


div_max = abs(div(1,1,1))
!$omp parallel do  default(shared) private(ix,iy,iz) reduction(max:div_max)
do iz = 1, Lz
   do iy = 1, Ly
      do ix = 1, Lx
!         if (div(ix,iy,iz) < 0.) ind_div(ix,iy,iz) = - 1
         div_max = max( div_max, abs(div(ix,iy,iz)) )
      enddo
   enddo
enddo
!$omp end parallel do
div_max = div_max / 3.
deallocate(v3,STAT=istat)
deallocate(vd,STAT=istat)
deallocate(div,STAT=istat)

return
end subroutine velocity_div1d   

subroutine CDM_velocity_stat


vxmax_cdm = 0.
vymax_cdm = 0.
vzmax_cdm = 0.
vsigma_cdm = 0.

!$omp parallel do default(shared) private(i,vx,vy,vz) &
!$ reduction(max:vxmax_cdm,vymax_cdm,vzmax_cdm) reduction(+:vsigma_cdm)
do i = 1, nobj

   vx = abs(v(1,i))
   vy = abs(v(2,i))
   vz = abs(v(3,i))
   vsigma_cdm = vsigma_cdm + vx*vx + vy*vy + vz*vz
   vxmax_cdm = max(vxmax_cdm,vx)
   vymax_cdm = max(vymax_cdm,vy)
   vzmax_cdm = max(vzmax_cdm,vz)

enddo
!$omp end parallel do
vsigma_cdm = sqrt(vsigma_cdm/float(nobj))

return
end subroutine CDM_velocity_stat

end module hydro_stat

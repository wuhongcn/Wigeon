! --------------------------------------------------------------------------------------------------

module gravity
use datablock
integer,dimension(Ls,Ls)::hoc
integer,dimension(2,Ls,Ls,n_proc)::htoc

contains

subroutine mass_assignment
integer :: npt, kpt
! integer, parameter :: npt=nobj/n_proc+min(1,mod(nobj,n_proc))
! integer, parameter :: kpt=Ls/n_proc+min(1,mod(Ls,n_proc))
integer ix,iy,iz
integer it,ip,toe

npt=nobj/n_proc + min(1,mod(nobj,n_proc))
kpt=Ls/n_proc + min(1,mod(Ls,n_proc))


!$omp parallel do default(shared) private(iz)
do iz=1,Lz
   d3(:,:,iz)=0.
enddo   
!$omp end parallel do  


!! Construct chaining lists in parallel
!! Particles should be ordered z,x,y to coincide with gas (!!CANCEL)
!$omp parallel do default(shared) private(it,ip,iy,iz)
do it=1,n_proc
   htoc(:,:,:,it)=0
   do ip=1+(it-1)*npt,min(nobj,it*npt)
 
      iy=int(rtemp(2,ip))
      iz=int(rtemp(3,ip))
!       write(*,*) 'iy=',iy,'iz=',iz,'ip=',ip
      iy=mod(iy-1+Ly,Ly)+1
      iz=mod(iz-1+Lz,Lz)+1
!      write(*,*) 'iy=',iy,'iz=',iz
!      pause
      ll(ip)=0
      if (htoc(1,iy,iz,it) == 0) then
          htoc(:,iy,iz,it)=ip
      else
         ll(htoc(2,iy,iz,it))=ip
         htoc(2,iy,iz,it)=ip
      endif
   enddo
enddo
!$omp end parallel do



!! Merge chaining lists
!$omp parallel do default(shared) private(it,iy,iz,toe)
do iz=1,Ls
   do iy=1,Ls
      hoc(iy,iz)=0
      do it=1,n_proc
         if (hoc(iy,iz) == 0) then
             hoc(iy,iz)=htoc(1,iy,iz,it)
             toe=htoc(2,iy,iz,it)
         else
            if (htoc(1,iy,iz,it) /= 0) then
                ll(toe)=htoc(1,iy,iz,it)
                toe=htoc(2,iy,iz,it)
            endif
         endif
      enddo
   enddo
enddo
!$omp end parallel do


!! Add particle density to density field
!$omp parallel do default(shared) private(it,iy,iz,ip)
do it=1,n_proc
   do iz=1+(it-1)*kpt,min(Ls,it*kpt),1
      do iy=1,Ls
         ip=hoc(iy,iz)
         call assign_cic(ip)
!          call assign_tsc(ip)
      enddo
   enddo
enddo
!!$omp end parallel do



contains

subroutine assign_cic(iq)
integer:: iq,iax,iay,iaz

do
   if (iq == 0) exit

   rrx=rtemp(1,iq)
   rry=rtemp(2,iq)
   rrz=rtemp(3,iq)
   
   iax=int(rrx)
   iay=int(rry)
   iaz=int(rrz)
   
   hx=rrx-iax
   hy=rry-iay
   hz=rrz-iaz
   hx0=1.-hx
   hy0=1.-hy
   hz0=1.-hz
   hxp1=hx
   hyp1=hy
   hzp1=hz

   iax =mod(iax-1+Lx,Lx)+1
   iay =mod(iay-1+Ly,Ly)+1
   iaz =mod(iaz-1+Lz,Lz)+1
   ixp1=mod(iax+Lx,Lx)+1
   iyp1=mod(iay+Ly,Ly)+1
   izp1=mod(iaz+Lz,Lz)+1

   d3(iax ,iay ,iaz ) = d3(iax ,iay ,iaz ) + hx0 *hy0 *hz0
   d3(ixp1,iay ,iaz ) = d3(ixp1,iay ,iaz ) + hxp1*hy0 *hz0
   d3(iax ,iyp1,iaz ) = d3(iax ,iyp1,iaz ) + hx0 *hyp1*hz0
   d3(ixp1,iyp1,iaz ) = d3(ixp1,iyp1,iaz ) + hxp1*hyp1*hz0
   d3(iax ,iay ,izp1) = d3(iax ,iay ,izp1) + hx0 *hy0 *hzp1
   d3(ixp1,iay ,izp1) = d3(ixp1,iay ,izp1) + hxp1*hy0 *hzp1
   d3(iax ,iyp1,izp1) = d3(iax ,iyp1,izp1) + hx0 *hyp1*hzp1
   d3(ixp1,iyp1,izp1) = d3(ixp1,iyp1,izp1) + hxp1*hyp1*hzp1

   iq=ll(iq)

enddo

return
end subroutine assign_cic

! TSC mass assignment
subroutine assign_tsc(iq)
integer::iq
integer::iax,iay,iaz
do 
   if(iq==0) exit

   rrx=rtemp(1,iq)
   rry=rtemp(2,iq)
   rrz=rtemp(3,iq)
   iax= int(rrx+0.5)
   iay= int(rry+0.5)
   iaz= int(rrz+0.5)
   hx = rrx -iax
   hy = rry -iay
   hz = rrz -iaz
   iax=mod(iax-1+Lx,Lx)+1
   iay=mod(iay-1+Ly,Ly)+1
   iaz=mod(iaz-1+Lz,Lz)+1

   hx0=0.75 - hx*hx
   hxp1=0.5* (0.5 + hx)**2
   hxm1=0.5* (0.5 - hx)**2
   hy0=0.75 - hy*hy
   hyp1=0.5* (0.5 + hy)**2
   hym1=0.5* (0.5 - hy)**2
   hz0= 0.75 - hz*hz
   hzp1=0.5* (0.5 + hz)**2
   hzm1=0.5* (0.5 - hz)**2

   ixp1=mod(iax+Lx,Lx)+1
   iyp1=mod(iay+Ly,Ly)+1
   izp1=mod(iaz+Lz,Lz)+1
   ixm1=mod(iax-2+Lx,Lx)+1
   iym1=mod(iay-2+Ly,Ly)+1
   izm1=mod(iaz-2+Lz,Lz)+1
   
   d3(ixm1,iym1,izm1)   = d3(ixm1,iym1,izm1)+ hxm1*hym1 *hzm1
   d3(iax ,iym1,izm1)   = d3(iax ,iym1,izm1)+ hx0 *hym1 *hzm1
   d3(ixp1,iym1,izm1)   = d3(ixp1,iym1,izm1)+ hxp1*hym1 *hzm1
   d3(ixm1,iay ,izm1)   = d3(ixm1,iay ,izm1)+ hxm1*hy0  *hzm1
   d3(iax ,iay ,izm1)   = d3(iax ,iay ,izm1)+ hx0 *hy0  *hzm1
   d3(ixp1,iay ,izm1)   = d3(ixp1,iay ,izm1)+ hxp1*hy0  *hzm1
   d3(ixm1,iyp1,izm1)   = d3(ixm1,iyp1,izm1)+ hxm1*hyp1 *hzm1
   d3(iax ,iyp1,izm1)   = d3(iax ,iyp1,izm1)+ hx0 *hyp1 *hzm1
   d3(ixp1,iyp1,izm1)   = d3(ixp1,iyp1,izm1)+ hxp1*hyp1 *hzm1
   d3(ixm1,iym1,iaz )   = d3(ixm1,iym1,iaz )+ hxm1*hym1 *hz0
   d3(iax ,iym1,iaz )   = d3(iax ,iym1,iaz )+ hx0 *hym1 *hz0
   d3(ixp1,iym1,iaz )   = d3(ixp1,iym1,iaz )+ hxp1*hym1 *hz0
   d3(ixm1,iay ,iaz )   = d3(ixm1,iay ,iaz )+ hxm1*hy0  *hz0
   d3(iax ,iay ,iaz )   = d3(iax ,iay ,iaz )+ hx0 *hy0  *hz0
   d3(ixp1,iay ,iaz )   = d3(ixp1,iay ,iaz )+ hxp1*hy0  *hz0
   d3(ixm1,iyp1,iaz )   = d3(ixm1,iyp1,iaz )+ hxm1*hyp1 *hz0
   d3(iax ,iyp1,iaz )   = d3(iax ,iyp1,iaz )+ hx0 *hyp1 *hz0
   d3(ixp1,iyp1,iaz )   = d3(ixp1,iyp1,iaz )+ hxp1*hyp1 *hz0
   d3(ixm1,iym1,izp1)   = d3(ixm1,iym1,izp1)+ hxm1*hym1 *hzp1
   d3(iax ,iym1,izp1)   = d3(iax ,iym1,izp1)+ hx0 *hym1 *hzp1
   d3(ixp1,iym1,izp1)   = d3(ixp1,iym1,izp1)+ hxp1*hym1 *hzp1
   d3(ixm1,iay ,izp1)   = d3(ixm1,iay ,izp1)+ hxm1*hy0  *hzp1
   d3(iax ,iay ,izp1)   = d3(iax ,iay ,izp1)+ hx0 *hy0  *hzp1
   d3(ixp1,iay ,izp1)   = d3(ixp1,iay ,izp1)+ hxp1*hy0  *hzp1
   d3(ixm1,iyp1,izp1)   = d3(ixm1,iyp1,izp1)+ hxm1*hyp1 *hzp1
   d3(iax ,iyp1,izp1)   = d3(iax ,iyp1,izp1)+ hx0 *hyp1 *hzp1
   d3(ixp1,iyp1,izp1)   = d3(ixp1,iyp1,izp1)+ hxp1*hyp1 *hzp1

   iq=ll(iq)
enddo

!write(*,*) 'iz after tsc',iz
return
end subroutine assign_tsc

end subroutine mass_assignment

!  Green function using the 7pt finite difference approximation to Laplacian 
subroutine greenz
real sinx(Lx/2+1),siny(Ly/2+1),sinz(Lz/2+1)

Lcx=Lx/2+1
Lcy=Ly/2+1
Lcz=Lz/2+1
gf=1./float(Lx*Ly*Lz)

do i=1,Lcx
   xi=float(i-1)
   sinx(i)=sin(pi*xi/float(Lx))
   sinx(i)=sinx(i)*sinx(i)
enddo

do i=1,Lcy
   yi=float(i-1)
   siny(i)=sin(pi*yi/float(Ly))
   siny(i)=siny(i)*siny(i)
enddo

do i=1,Lcz
   zi=float(i-1)
   sinz(i)=sin(pi*zi/float(Lz))
   sinz(i)=sinz(i)*sinz(i)
enddo

do igz=1,Lcz
   do igy=1,Lcy
      do igx=1,Lcx
         if (igx==1.and.igy==1.and.igz==1) then
            gk(igx,igy,igz)=0.
         else
            gk(igx,igy,igz)=pi/(sinx(igx)+siny(igy)+sinz(igz))
         endif
      enddo
   enddo
enddo

gk=gf*gk

return
end subroutine greenz

! optimized Green function using 7pt finite difference approximation to Laplacian
! soft is the soften length 
subroutine greenc(soft)
real,parameter:: tpi=6.28318530717959
integer,parameter:: ipn=2
real:: Sprf,soft

gf=12.*tpi/float(Lx*Ly*Lz)

tLx=tpi/float(Lx)
tLy=tpi/float(Ly)
tLz=tpi/float(Lz)
Lx2=Lx/2
Ly2=Ly/2
Lz2=Lz/2
Lcx=Lx2+1
Lcy=Ly2+1
Lcz=Lz2+1
Lex=Lx/4+1
Ley=Ly/4+1
Lez=Lz/4+1

af=2./15.
do igz=1,Lcz
   rkz=tLz*float(igz-1)
   s1z=sin(rkz)
   c1z=cos(rkz)
   s2z=sin(rkz/2.)
   s2z2=s2z**2
   psz=2.*s2z
   sz=sin(rkz/4.)**2
   cz=1.-sz
   do igy=1,Lcy
      rky=tLy*float(igy-1)
      s1y=sin(rky)
      c1y=cos(rky)
      s2y=sin(rky/2.)
      s2y2=s2y**2
      psy=2.*s2y
      sy=sin(rky/4.)**2
      cy=1.-sy
      do igx=1,Lcx
         rkx=tLx*float(igx-1)
         s1x=sin(rkx)
         c1x=cos(rkx)
         s2x=sin(rkx/2.)
         s2x2=s2x**2
         psx=2.*s2x
         sx=sin(rkx/4.)**2
         cx=1.-sx

         gp3t=0.
         itest=mod(igx-1,Lx2)+mod(igy-1,Ly2)+mod(igz-1,Lz2)
         if (itest.ne.0) then
            rk2=rkx*rkx+rky*rky+rkz*rkz

            SUM=0.
            SRx=0.
            SRy=0.
            SRz=0.
            do nz=-ipn,ipn
               rkzn=rkz-tpi*nz
               qz=1.
               if(rkzn.ne.0.)qz=psz/rkzn
               do ny=-ipn,ipn
                  rkyn=rky-tpi*ny
                  qy=1.
                  if(rkyn.ne.0.)qy=psy/rkyn
                  do nx=-ipn,ipn
                     rkxn=rkx-tpi*nx
                     qx=1.
                     if (rkxn.ne.0.)qx=psx/rkxn

                     qt=(qx*qy*qz)**3
                     rk2n=rkxn*rkxn+rkyn*rkyn+rkzn*rkzn
                     Sprf=shp2(soft*sqrt(rk2n)/2.)
                     SprfU=Sprf*qt
                     st=SprfU*SprfU/rk2n
                     SUM=SUM+st
                     SRx=SRx+rkxn*st
                     SRy=SRy+rkyn*st
                     SRz=SRz+rkzn*st
                  enddo
               enddo
            enddo

            SAa = (1.+s2x2*(af*s2x2-1.))*(1.+s2y2*(af*s2y2-1.))*(1.+s2z2*(af*s2z2-1.))
            Sd=SAa**2

            Dkx1=s1x*(1.+c1y+c1z)
            Dky1=s1y*(1.+c1z+c1x)
            Dkz1=s1z*(1.+c1x+c1y)

            Dkx=2.*Dkx1
            Dky=2.*Dky1
            Dkz=2.*Dkz1

            DWR=(SRx*Dkx+SRy*Dky+SRz*Dkz)/(Dkx*Dkx+Dky*Dky+Dkz*Dkz)

            gp3t=gf*DWR/Sd

         endif         

         gk(igx,igy,igz)=gp3t

      enddo
   enddo
enddo

return

contains

function shp2(x)
shp2=12.*(2.-2.*cos(x)-x*sin(x))/x**4
end function shp2

end subroutine greenc

! CIC interpolation 


subroutine CIC_interp(g3,fg)
real,dimension(nobj)::fg
real,dimension(Lx,Ly,Lz)::g3
integer ::iax,iay,iaz
fg=0.

!$omp parallel do default(shared) private(i,rrx,rry,rrz,iax,iay,iaz,hx,hy,hz,hx0,hy0,hz0) &
!$ private(hxp1,hyp1,hzp1,ixp1,iyp1,izp1)

do i=1,nobj

   !rrx=r(1,i)
   !rry=r(2,i)
   !rrz=r(3,i)
   rrx=rtemp(1,i)
   rry=rtemp(2,i)
   rrz=rtemp(3,i)
   
   iax=int(rrx)
   iay=int(rry)
   iaz=int(rrz)
   
   hx=rrx-iax
   hy=rry-iay
   hz=rrz-iaz

   hx0=1.-hx
   hy0=1.-hy
   hz0=1.-hz
   hxp1=hx
   hyp1=hy
   hzp1=hz

   iax =mod(iax-1+Lx,Lx)+1
   iay =mod(iay-1+Ly,Ly)+1
   iaz =mod(iaz-1+Lz,Lz)+1
   ixp1=mod(iax +Lx,Lx)+1
   iyp1=mod(iay +Ly,Ly)+1
   izp1=mod(iaz +Lz,Lz)+1

   ac = g3(iax ,iay ,iaz )* hx0 *hy0  *hz0  &
      + g3(ixp1,iay ,iaz )* hxp1*hy0  *hz0  &
      + g3(iax ,iyp1,iaz )* hx0 *hyp1 *hz0  &
      + g3(ixp1,iyp1,iaz )* hxp1*hyp1 *hz0  &
      + g3(iax ,iay ,izp1)* hx0 *hy0  *hzp1 &
      + g3(ixp1,iay ,izp1)* hxp1*hy0  *hzp1 &
      + g3(iax ,iyp1,izp1)* hx0 *hyp1 *hzp1 &
      + g3(ixp1,iyp1,izp1)* hxp1*hyp1 *hzp1 

   fg(i)=ac

enddo 
!$omp end parallel do

return
end subroutine CIC_interp

! 4pt finite difference 
subroutine gridf04(ic,g3,frc)
real,dimension(Lx,Ly,Lz)::g3,frc
integer::ix,iy,iz
Lx1=10*Lx-1
Ly1=10*Ly-1
Lz1=10*Lz-1

select case(ic)
case(1)
!$omp parallel do default(shared) private(ix,iy,iz,ixp1,ixm1,ixp2,ixm2)
   do iz=1,Lz
      do iy=1,Ly
         do ix=1,Lx
            ixp1=mod(ix+1+Lx1,Lx)+1
            ixm1=mod(ix-1+Lx1,Lx)+1
            ixp2=mod(ix+2+Lx1,Lx)+1
            ixm2=mod(ix-2+Lx1,Lx)+1
            frc(ix,iy,iz)=2.*(g3(ixp1,iy  ,iz  )-g3(ixm1,iy  ,iz  ))/3. &
                            -(g3(ixp2,iy  ,iz  )-g3(ixm2,iy  ,iz  ))/12. 
         enddo
      enddo
   enddo
!$omp end parallel do
case(2)
!$omp parallel do default(shared) private(ix,iy,iz,iyp1,iym1,iyp2,iym2)
   do iz=1,Lz
      do iy=1,Ly
         do ix=1,Lx
            iyp1=mod(iy+1+Ly1,Ly)+1
            iym1=mod(iy-1+Ly1,Ly)+1
            iyp2=mod(iy+2+Ly1,Ly)+1
            iym2=mod(iy-2+Ly1,Ly)+1
            frc(ix,iy,iz)=2.*(g3(ix  ,iyp1,iz  )-g3(ix  ,iym1,iz  ))/3. &
                            -(g3(ix  ,iyp2,iz  )-g3(ix  ,iym2,iz  ))/12. 

         enddo
      enddo
   enddo
!$omp end parallel do
case(3)
!$omp parallel do default(shared) private(ix,iy,iz,izp1,izm1,izp2,izm2)
   do iz=1,Lz
      do iy=1,Ly
         do ix=1,Lx
            izp1=mod(iz+1+Lz1,Lz)+1
            izm1=mod(iz-1+Lz1,Lz)+1
            izp2=mod(iz+2+Lz1,Lz)+1
            izm2=mod(iz-2+Lz1,Lz)+1
            frc(ix,iy,iz)=2.*(g3(ix  ,iy  ,izp1)-g3(ix  ,iy  ,izm1))/3. &
                            -(g3(ix  ,iy  ,izp2)-g3(ix  ,iy  ,izm2))/12. 
         enddo
      enddo
   enddo
!$omp end parallel do
end select

return
end subroutine gridf04

! 10pt finite difference
subroutine gridf10(ic,g3,frc)
real,dimension(Lx,Ly,Lz)::g3,frc
integer::ix,iy,iz
Lx1=Lx-1
Ly1=Ly-1
Lz1=Lz-1

if (ic.eq.1) then
!$omp parallel do default(shared) private(ix,iy,iz,ixp1,ixm1,iyp1,iym1,izp1,izm1)        
   do iz=1,Lz
      do iy=1,Ly
         do ix=1,Lx
            ixp1=mod(ix+1+Lx1,Lx)+1
            ixm1=mod(ix-1+Lx1,Lx)+1
            iyp1=mod(iy+1+Ly1,Ly)+1
            iym1=mod(iy-1+Ly1,Ly)+1
            izp1=mod(iz+1+Lz1,Lz)+1
            izm1=mod(iz-1+Lz1,Lz)+1
            frc(ix,iy,iz)=2.*(g3(ixp1,iy  ,iz  )-g3(ixm1,iy  ,iz  )) &
                         +g3(ixp1,iyp1,iz  )-g3(ixm1,iyp1,iz  ) &
                         +g3(ixp1,iym1,iz  )-g3(ixm1,iym1,iz  ) &
                         +g3(ixp1,iy  ,izp1)-g3(ixm1,iy  ,izp1) &
                         +g3(ixp1,iy  ,izm1)-g3(ixm1,iy  ,izm1)
         enddo
      enddo
   enddo
!$omp end parallel do   
else if(ic.eq.2)then
!$omp parallel do default(shared) private(ix,iy,iz,ixp1,ixm1,iyp1,iym1,izp1,izm1)        
     do iz=1,Lz
        do iy=1,Ly
           do ix=1,Lx
              ixp1=mod(ix+1+Lx1,Lx)+1
              ixm1=mod(ix-1+Lx1,Lx)+1
              iyp1=mod(iy+1+Ly1,Ly)+1
              iym1=mod(iy-1+Ly1,Ly)+1
              izp1=mod(iz+1+Lz1,Lz)+1
              izm1=mod(iz-1+Lz1,Lz)+1
              frc(ix,iy,iz)=2.*(g3(ix  ,iyp1,iz  )-g3(ix  ,iym1,iz  )) &
                           +g3(ixp1,iyp1,iz  )-g3(ixp1,iym1,iz  ) &
                           +g3(ixm1,iyp1,iz  )-g3(ixm1,iym1,iz  ) &
                           +g3(ix  ,iyp1,izp1)-g3(ix  ,iym1,izp1) &
                           +g3(ix  ,iyp1,izm1)-g3(ix  ,iym1,izm1)
         enddo
      enddo
   enddo
!$omp end parallel do   
else if(ic.eq.3)then
!$omp parallel do default(shared) private(ix,iy,iz,ixp1,ixm1,iyp1,iym1,izp1,izm1)        
     do iz=1,Lz
        do iy=1,Ly
           do ix=1,Lx
              ixp1=mod(ix+1+Lx1,Lx)+1
              ixm1=mod(ix-1+Lx1,Lx)+1
              iyp1=mod(iy+1+Ly1,Ly)+1
              iym1=mod(iy-1+Ly1,Ly)+1
              izp1=mod(iz+1+Lz1,Lz)+1
              izm1=mod(iz-1+Lz1,Lz)+1
              frc(ix,iy,iz)=2.*(g3(ix  ,iy  ,izp1)-g3(ix  ,iy  ,izm1)) &
                           +g3(ix  ,iyp1,izp1)-g3(ix  ,iyp1,izm1) &
                           +g3(ix  ,iym1,izp1)-g3(ix  ,iym1,izm1) &
                           +g3(ixm1,iy  ,izp1)-g3(ixm1,iy  ,izm1) &
                           +g3(ixp1,iy  ,izp1)-g3(ixp1,iy  ,izm1) 
         enddo
      enddo
   enddo
!$omp end parallel do   
end if

frc=frc/12. 

return
end subroutine gridf10

! convolution of p3 with g3: p3*g3 -> p3
subroutine cnvl3d(p3,g3)
integer,dimension(3)::ndim
real,dimension(Lx,Ly,Lz)::p3
real,dimension(Lx/2+1,Ly/2+1,Lz/2+1)::g3
real,allocatable,dimension(:,:,:)::q3
      
ndim(1)=Lx
ndim(2)=Ly
ndim(3)=Lz
Lx2=Lx/2+1
Ly2=Ly/2+1
Lz2=Lz/2+1
allocate(q3(Lx,Ly,Lz),STAT=istat)

q3=0.
call fft3d(p3,q3,ndim,1)

!$omp parallel do default(shared) private(ix,iy,iz,icx,icy,icz)
do iz=1,Lz
   do iy=1,Ly
      do ix=1,Lx
         if (ix<=Lx2) then
            icx=ix
         else 
            icx=mod(Lx+1-ix,Lx)+1
         endif
         if (iy<=Ly2) then
            icy=iy
         else 
            icy=mod(Ly+1-iy,Ly)+1
         endif
         if (iz<=Lz2) then
            icz=iz
         else 
            icz=mod(Lz+1-iz,Lz)+1
         endif
         p3(ix,iy,iz) = p3(ix,iy,iz) * g3(icx,icy,icz)
         q3(ix,iy,iz) = q3(ix,iy,iz) * g3(icx,icy,icz)
      enddo
   enddo
enddo
!$omp end parallel do

call fft3d(p3,q3,ndim,-1)

deallocate(q3,STAT=istat)

return
end subroutine cnvl3d

subroutine fft3d(a,b,nn,isn)
use MKL_DFTI
integer,dimension(3)::nn
dimension a(*),b(*)
complex,allocatable,dimension(:)::fft3
type (DFTI_DESCRIPTOR),POINTER:: my_fft
integer ::status

n1=nn(1)
n2=nn(2)
n3=nn(3)
ns=n1*n2*n3

allocate (fft3(n1*n2*n3),STAT=istat)

fft3(1:ns)=cmplx(a(1:ns),b(1:ns))

status=DftiCreateDescriptor(my_fft,DFTI_double,DFTI_complex,3,nn)
status=DftiCommitDescriptor(my_fft)
if(isn==1) then
    status=DftiComputeForward(my_fft,fft3)
else if(isn==-1) then 
    status=DftiComputeBackward(my_fft,fft3)
endif
status=DftiFreeDescriptor(my_fft)

a(1:ns)=real(fft3(1:ns))
b(1:ns)=imag(fft3(1:ns))

deallocate (fft3,STAT=istat)

return
end subroutine fft3d

end module gravity




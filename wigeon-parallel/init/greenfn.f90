subroutine greenOpt(gp3,L)
parameter(tpi=6.28318530717959,ipn=2,isgn=1-2*mod(ipn,2))
real,dimension(L/2+1,L/2+1,L/2+1)::gp3

tpiL=tpi/float(L)
L2=L/2
Lc=L/2+1

do igz=1,Lc
   rkz=tpiL*float(igz-1)
   s1z=sin(rkz)
   c1z=cos(rkz)
   s2z=sin(rkz/2.)
   psz=2.*s2z
   pz=1.
   if(igz.ne.1)pz=psz/rkz
   do igy=1,igz
      rky=tpiL*float(igy-1)
      s1y=sin(rky)
      c1y=cos(rky)
      s2y=sin(rky/2.)
      psy=2.*s2y
      py=1.
      if(igy.ne.1)py=psy/rky
      do igx=1,igy
         rkx=tpiL*float(igx-1)
         s1x=sin(rkx)
         c1x=cos(rkx)
         s2x=sin(rkx/2.)
         psx=2.*s2x
         px=1.
         if (igx.ne.1) px=psx/rkx
         gp3t=0.
         itest=mod(igx-1,L2)+mod(igy-1,L2)+mod(igz-1,L2)
         if(itest.ne.0)then
         rk2=rkx*rkx+rky*rky+rkz*rkz
         Smp=0.

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
                  if(rkxn.ne.0.)qx=psx/rkxn

                  qt=(qx*qy*qz)**3

                  Smp=qt-Smp
               enddo
            enddo
         enddo
         Smp=isgn*Smp

         Dkx=s1x*(1.+c1y+c1z)
         Dky=s1y*(1.+c1z+c1x)
         Dkz=s1z*(1.+c1x+c1y)

         gp3t=3.0/((rkx*Dkx+rky*Dky+rkz*Dkz)*Smp)

         end if
         gp3(igx,igy,igz)=gp3t
         gp3(igx,igz,igy)=gp3t
         gp3(igy,igz,igx)=gp3t
         gp3(igy,igx,igz)=gp3t
         gp3(igz,igx,igy)=gp3t
         gp3(igz,igy,igx)=gp3t
      enddo
   enddo
enddo

return
end subroutine greenOpt

! Green function in k space for -\nabla^2 = 1/k^2 
subroutine greenz(gk,L)
real gk(L/2+1,L/2+1,L/2+1), sins(L/2+1)

pi=acos(-1.)
Lc=L/2+1

do i=1,Lc
   xi=float(i-1)
   sins(i)=sin(pi*xi/float(L))
   sins(i)=sins(i)*sins(i)
enddo

do igz=1,Lc
   do igy=1,igz
      do igx=1,igy
         if (igx==1.and.igy==1.and.igz==1) then
            gt=0.
         else
            gt=0.25/(sins(igx)+sins(igy)+sins(igz))
         endif
         gk(igx,igy,igz)=gt
         gk(igx,igz,igy)=gt
         gk(igy,igz,igx)=gt
         gk(igy,igx,igz)=gt
         gk(igz,igx,igy)=gt
         gk(igz,igy,igx)=gt
      enddo
   enddo
enddo

return
end subroutine greenz

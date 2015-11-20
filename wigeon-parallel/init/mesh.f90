subroutine assign(N,r,L,d3)
real,dimension(3,N)::r
real,dimension(L,L,L)::d3

d3=0.

do i=1,N
   rrx=r(1,i)
   rry=r(2,i)
   rrz=r(3,i)
   ix= int(rrx+0.5)
   iy= int(rry+0.5)
   iz= int(rrz+0.5)
   hx = rrx -ix
   hy = rry -iy
   hz = rrz -iz
   ix=mod(ix-1+L,L)+1
   iy=mod(iy-1+L,L)+1
   iz=mod(iz-1+L,L)+1

   hx0=0.75 - hx*hx
   hxp1=0.5* (0.5 + hx)**2
   hxm1=0.5* (0.5 - hx)**2
   hy0=0.75 - hy*hy
   hyp1=0.5* (0.5 + hy)**2
   hym1=0.5* (0.5 - hy)**2
   hz0= 0.75 - hz*hz
   hzp1=0.5* (0.5 + hz)**2
   hzm1=0.5* (0.5 - hz)**2

   ixp1=mod(ix+L,L)+1
   iyp1=mod(iy+L,L)+1
   izp1=mod(iz+L,L)+1
   ixm1=mod(ix-2+L,L)+1
   iym1=mod(iy-2+L,L)+1
   izm1=mod(iz-2+L,L)+1
   
   d3(ixm1,iym1,izm1)   = d3(ixm1,iym1,izm1)+ hxm1*hym1 *hzm1
   d3(ix  ,iym1,izm1)   = d3(ix  ,iym1,izm1)+ hx0 *hym1 *hzm1
   d3(ixp1,iym1,izm1)   = d3(ixp1,iym1,izm1)+ hxp1*hym1 *hzm1
   d3(ixm1,iy  ,izm1)   = d3(ixm1,iy  ,izm1)+ hxm1*hy0  *hzm1
   d3(ix  ,iy  ,izm1)   = d3(ix  ,iy  ,izm1)+ hx0 *hy0  *hzm1
   d3(ixp1,iy  ,izm1)   = d3(ixp1,iy  ,izm1)+ hxp1*hy0  *hzm1
   d3(ixm1,iyp1,izm1)   = d3(ixm1,iyp1,izm1)+ hxm1*hyp1 *hzm1
   d3(ix  ,iyp1,izm1)   = d3(ix  ,iyp1,izm1)+ hx0 *hyp1 *hzm1
   d3(ixp1,iyp1,izm1)   = d3(ixp1,iyp1,izm1)+ hxp1*hyp1 *hzm1
   d3(ixm1,iym1,iz  )   = d3(ixm1,iym1,iz  )+ hxm1*hym1 *hz0
   d3(ix  ,iym1,iz  )   = d3(ix  ,iym1,iz  )+ hx0 *hym1 *hz0
   d3(ixp1,iym1,iz  )   = d3(ixp1,iym1,iz  )+ hxp1*hym1 *hz0
   d3(ixm1,iy  ,iz  )   = d3(ixm1,iy  ,iz  )+ hxm1*hy0  *hz0
   d3(ix  ,iy  ,iz  )   = d3(ix  ,iy  ,iz  )+ hx0 *hy0  *hz0
   d3(ixp1,iy  ,iz  )   = d3(ixp1,iy  ,iz  )+ hxp1*hy0  *hz0
   d3(ixm1,iyp1,iz  )   = d3(ixm1,iyp1,iz  )+ hxm1*hyp1 *hz0
   d3(ix  ,iyp1,iz  )   = d3(ix  ,iyp1,iz  )+ hx0 *hyp1 *hz0
   d3(ixp1,iyp1,iz  )   = d3(ixp1,iyp1,iz  )+ hxp1*hyp1 *hz0
   d3(ixm1,iym1,izp1)   = d3(ixm1,iym1,izp1)+ hxm1*hym1 *hzp1
   d3(ix  ,iym1,izp1)   = d3(ix  ,iym1,izp1)+ hx0 *hym1 *hzp1
   d3(ixp1,iym1,izp1)   = d3(ixp1,iym1,izp1)+ hxp1*hym1 *hzp1
   d3(ixm1,iy  ,izp1)   = d3(ixm1,iy  ,izp1)+ hxm1*hy0  *hzp1
   d3(ix  ,iy  ,izp1)   = d3(ix  ,iy  ,izp1)+ hx0 *hy0  *hzp1
   d3(ixp1,iy  ,izp1)   = d3(ixp1,iy  ,izp1)+ hxp1*hy0  *hzp1
   d3(ixm1,iyp1,izp1)   = d3(ixm1,iyp1,izp1)+ hxm1*hyp1 *hzp1
   d3(ix  ,iyp1,izp1)   = d3(ix  ,iyp1,izp1)+ hx0 *hyp1 *hzp1
   d3(ixp1,iyp1,izp1)   = d3(ixp1,iyp1,izp1)+ hxp1*hyp1 *hzp1

enddo

xL=float(L)
d3_av=float(N)/xL**3
d3=d3/d3_av

return
end subroutine assign      

subroutine interp(ic,N,r,v,L,frc)
real r(3,N),v(3,N)
real frc(L,L,L)

do i=1,N
   rrx=r(1,i)
   rry=r(2,i)
   rrz=r(3,i)
   ix= int(rrx+0.5)
   iy= int(rry+0.5)
   iz= int(rrz+0.5)
   hx = rrx -ix
   hy = rry -iy
   hz = rrz -iz
   ix=mod(ix-1,L)+1
   iy=mod(iy-1,L)+1
   iz=mod(iz-1,L)+1

   hx0=0.75 - hx*hx
   hxp1=0.5* (0.5 + hx)**2
   hxm1=0.5* (0.5 - hx)**2
   hy0=0.75 - hy*hy
   hyp1=0.5* (0.5 + hy)**2
   hym1=0.5* (0.5 - hy)**2
   hz0=0.75 - hz*hz
   hzp1=0.5* (0.5 + hz)**2
   hzm1=0.5* (0.5 - hz)**2
   ixp1=mod(ix,L)+1
   iyp1=mod(iy,L)+1
   izp1=mod(iz,L)+1
   ixm1=mod(ix-2+L,L)+1
   iym1=mod(iy-2+L,L)+1
   izm1=mod(iz-2+L,L)+1
   ac =   frc(ixm1,iym1,izm1)* hxm1*hym1 *hzm1  &
&       + frc(ix  ,iym1,izm1)* hx0 *hym1 *hzm1  &
&       + frc(ixp1,iym1,izm1)* hxp1*hym1 *hzm1  &
&       + frc(ixm1,iy  ,izm1)* hxm1*hy0  *hzm1  &
&       + frc(ix  ,iy  ,izm1)* hx0 *hy0  *hzm1  &
&       + frc(ixp1,iy  ,izm1)* hxp1*hy0  *hzm1  &
&       + frc(ixm1,iyp1,izm1)* hxm1*hyp1 *hzm1  &
&       + frc(ix  ,iyp1,izm1)* hx0 *hyp1 *hzm1  &
&       + frc(ixp1,iyp1,izm1)* hxp1*hyp1 *hzm1  &
&       + frc(ixm1,iym1,iz  )* hxm1*hym1 *hz0   &
&       + frc(ix  ,iym1,iz  )* hx0 *hym1 *hz0   &
&       + frc(ixp1,iym1,iz  )* hxp1*hym1 *hz0   &
&       + frc(ixm1,iy  ,iz  )* hxm1*hy0  *hz0   &
&       + frc(ix  ,iy  ,iz  )* hx0 *hy0  *hz0
   ac =   ac   &
&       + frc(ixp1,iy  ,iz  )* hxp1*hy0  *hz0   &
&       + frc(ixm1,iyp1,iz  )* hxm1*hyp1 *hz0   &
&       + frc(ix  ,iyp1,iz  )* hx0 *hyp1 *hz0   &
&       + frc(ixp1,iyp1,iz  )* hxp1*hyp1 *hz0   &
&       + frc(ixm1,iym1,izp1)* hxm1*hym1 *hzp1  &
&       + frc(ix  ,iym1,izp1)* hx0 *hym1 *hzp1  &
&       + frc(ixp1,iym1,izp1)* hxp1*hym1 *hzp1  &
&       + frc(ixm1,iy  ,izp1)* hxm1*hy0  *hzp1  &
&       + frc(ix  ,iy  ,izp1)* hx0 *hy0  *hzp1  &
&       + frc(ixp1,iy  ,izp1)* hxp1*hy0  *hzp1  &
&       + frc(ixm1,iyp1,izp1)* hxm1*hyp1 *hzp1  &
&       + frc(ix  ,iyp1,izp1)* hx0 *hyp1 *hzp1  &
&       + frc(ixp1,iyp1,izp1)* hxp1*hyp1 *hzp1  

!   v(i,ic)=v(i,ic)+ac
    v(ic,i)=v(ic,i)+ac
enddo

return
end

subroutine gridf2(ic,d3,frc,L)
real d3(L,L,L),frc(L,L,L)

L1=L-1
if (ic.eq.1)then
   do iz=1,L
      do iy=1,L
         do ix=1,L
            ixp1=mod(ix+1+L1,L)+1
            ixm1=mod(ix-1+L1,L)+1
            frc(ix,iy,iz)=0.5*(d3(ixp1,iy,iz)-d3(ixm1,iy,iz))
         enddo
      enddo
   enddo
else if (ic.eq.2)then
   do iz=1,L
      do iy=1,L
         do ix=1,L
            iyp1=mod(iy+1+L1,L)+1
            iym1=mod(iy-1+L1,L)+1
            frc(ix,iy,iz)=0.5*(d3(ix,iyp1,iz)-d3(ix,iym1,iz))
         enddo
      enddo
   enddo
else if (ic.eq.3)then
   do iz=1,L
      do iy=1,L
         do ix=1,L
            izp1=mod(iz+1+L1,L)+1
            izm1=mod(iz-1+L1,L)+1
            frc(ix,iy,iz)=0.5*(d3(ix,iy,izp1)-d3(ix,iy,izm1))
         enddo
      enddo
   enddo
end if

return
end

! 10pt direction (ic=1,2,3) derivative of a grid field d3  

subroutine gridf10(ic,d3,frc,L)
real,dimension(L,L,L)::d3,frc

L1=L-1
if (ic.eq.1)then
   do iz=1,L
      do iy=1,L
         do ix=1,L
            ixp1=mod(ix+1+L1,L)+1
            ixm1=mod(ix-1+L1,L)+1
            iyp1=mod(iy+1+L1,L)+1
            iym1=mod(iy-1+L1,L)+1
            izp1=mod(iz+1+L1,L)+1
            izm1=mod(iz-1+L1,L)+1
            frc(ix,iy,iz) = 2.*(d3(ixp1,iy  ,iz )-d3(ixm1,iy ,iz  )) &
     &                    + d3(ixp1,iyp1,iz  )-d3(ixm1,iyp1,iz  )    &
     &                    + d3(ixp1,iym1,iz  )-d3(ixm1,iym1,iz  )    &
     &                    + d3(ixp1,iy  ,izp1)-d3(ixm1,iy  ,izp1)    &
     &                    + d3(ixp1,iy  ,izm1)-d3(ixm1,iy  ,izm1)
         enddo
      enddo
   enddo
else if (ic.eq.2)then
   do iz=1,L
      do iy=1,L
         do ix=1,L
            ixp1=mod(ix+1+L1,L)+1
            ixm1=mod(ix-1+L1,L)+1
            iyp1=mod(iy+1+L1,L)+1
            iym1=mod(iy-1+L1,L)+1
            izp1=mod(iz+1+L1,L)+1
            izm1=mod(iz-1+L1,L)+1
            frc(ix,iy,iz) = 2.*(d3(ix  ,iyp1,iz )-d3(ix ,iym1, iz ))  &
     &                    + d3(ixp1,iyp1,iz  )-d3(ixp1,iym1,iz  )     &
     &                    + d3(ixm1,iyp1,iz  )-d3(ixm1,iym1,iz  )     &
     &                    + d3(ix  ,iyp1,izp1)-d3(ix  ,iym1,izp1)     &
     &                    + d3(ix  ,iyp1,izm1)-d3(ix  ,iym1,izm1)     
         enddo
      enddo
   enddo
else if (ic.eq.3)then
   do iz=1,L
      do iy=1,L
         do ix=1,L
            ixp1=mod(ix+1+L1,L)+1
            ixm1=mod(ix-1+L1,L)+1
            iyp1=mod(iy+1+L1,L)+1
            iym1=mod(iy-1+L1,L)+1
            izp1=mod(iz+1+L1,L)+1
            izm1=mod(iz-1+L1,L)+1
            frc(ix,iy,iz) = 2.*(d3(ix  ,iy  ,izp1)-d3(ix ,iy ,izm1))  &
     &                    + d3(ix  ,iyp1,izp1)-d3(ix  ,iyp1,izm1)     &
     &                    + d3(ix  ,iym1,izp1)-d3(ix  ,iym1,izm1)     &
     &                    + d3(ixm1,iy  ,izp1)-d3(ixm1,iy  ,izm1)     &
     &                    + d3(ixp1,iy  ,izp1)-d3(ixp1,iy  ,izm1)
         enddo
      enddo
   enddo
end if

frc=frc/12.

return
end

subroutine gridf10_cell(L,n3,d3,frc)
real,dimension(L,L,L)::d3
integer,dimension(3)::n3
real,dimension(3)::frc

ix=n3(1)
iy=n3(2)
iz=n3(3)

L1=L-1

ixp1=mod(ix+1+L1,L)+1
ixm1=mod(ix-1+L1,L)+1
iyp1=mod(iy+1+L1,L)+1
iym1=mod(iy-1+L1,L)+1
izp1=mod(iz+1+L1,L)+1
izm1=mod(iz-1+L1,L)+1
frc(1) = 2.*(d3(ixp1,iy  ,iz )-d3(ixm1,iy ,iz  )) &
&      + d3(ixp1,iyp1,iz  )-d3(ixm1,iyp1,iz  )    &
&      + d3(ixp1,iym1,iz  )-d3(ixm1,iym1,iz  )    &
&      + d3(ixp1,iy  ,izp1)-d3(ixm1,iy  ,izp1)    &
&      + d3(ixp1,iy  ,izm1)-d3(ixm1,iy  ,izm1)

ixp1=mod(ix+1+L1,L)+1
ixm1=mod(ix-1+L1,L)+1
iyp1=mod(iy+1+L1,L)+1
iym1=mod(iy-1+L1,L)+1
izp1=mod(iz+1+L1,L)+1
izm1=mod(iz-1+L1,L)+1
frc(2) = 2.*(d3(ix  ,iyp1,iz )-d3(ix ,iym1, iz ))  &
&      + d3(ixp1,iyp1,iz  )-d3(ixp1,iym1,iz  )     &
&      + d3(ixm1,iyp1,iz  )-d3(ixm1,iym1,iz  )     &
&      + d3(ix  ,iyp1,izp1)-d3(ix  ,iym1,izp1)     &
&      + d3(ix  ,iyp1,izm1)-d3(ix  ,iym1,izm1)     

ixp1=mod(ix+1+L1,L)+1
ixm1=mod(ix-1+L1,L)+1
iyp1=mod(iy+1+L1,L)+1
iym1=mod(iy-1+L1,L)+1
izp1=mod(iz+1+L1,L)+1
izm1=mod(iz-1+L1,L)+1
frc(3) = 2.*(d3(ix  ,iy  ,izp1)-d3(ix ,iy ,izm1))  &
&      + d3(ix  ,iyp1,izp1)-d3(ix  ,iyp1,izm1)     &
&      + d3(ix  ,iym1,izp1)-d3(ix  ,iym1,izm1)     &
&      + d3(ixm1,iy  ,izp1)-d3(ixm1,iy  ,izm1)     &
&      + d3(ixp1,iy  ,izp1)-d3(ixp1,iy  ,izm1)

frc(1:3)=frc(1:3)/12.

return
end subroutine gridf10_cell

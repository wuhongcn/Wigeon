subroutine smooth3D(L,rad,data_array)
real,dimension(L,L,L)::data_array
real,allocatable,dimension(:,:,:)::gp3

Lc=L/2+1
allocate(gp3(Lc,Lc,Lc),STAT=istat)
      
call wthgreen(rad,L,gp3) 
call cnvl3D(L,L,L,data_array,gp3)

deallocate (gp3,STAT=istat)

return
end subroutine smooth3D

subroutine wthgreen(rad,L,gp3)
!	implicit double precision (a-h,o-z)
real,parameter::tpi=6.28318530717959
real,dimension(L/2+1,L/2+1,L/2+1)::gp3

w(xk)=exp(-0.5*xk*xk)

sf=1./float(L**3)

tpiL=tpi/float(L)
L2=L/2
Lc=L/2+1
      
do 400 igz=1,Lc
   rkz=tpiL*float(igz-1)
   do 400 igy=1,igz
      rky=tpiL*float(igy-1)
      do 400 igx=1,igy
         rkx=tpiL*float(igx-1)
         rk2n=rkx*rkx+rky*rky+rkz*rkz
         gp3t=W(rad*sqrt(rk2n))
         gp3(igx,igy,igz)=gp3t
         gp3(igx,igz,igy)=gp3t
         gp3(igy,igz,igx)=gp3t
         gp3(igy,igx,igz)=gp3t
         gp3(igz,igx,igy)=gp3t
         gp3(igz,igy,igx)=gp3t
400   continue

gp3=gp3*sf

return
end subroutine wthgreen	

subroutine cnvl3d(Lgx,Lgy,Lgz,d3,g3)
integer,dimension(3)::ndim
real,dimension(Lgx,Lgy,Lgz)::d3
real,dimension(Lgx/2+1,Lgy/2+1,Lgz/2+1)::g3
real,allocatable,dimension(:,:,:)::r3
      
ndim(1)=Lgx
ndim(2)=Lgy
ndim(3)=Lgz
Lgx2=Lgx/2+1
Lgy2=Lgy/2+1
Lgz2=Lgz/2+1

allocate(r3(Lgx,Lgy,Lgz),STAT=istat)

r3=0.
call fft3d(d3,r3,ndim,1)

do iz=1,Lgz
   do iy=1,Lgy
      do ix=1,Lgx
         if (ix<=Lgx2) then
            icx=ix
         else 
            icx=mod(Lgx+1-ix,Lgx)+1
         endif
         if (iy<=Lgy2) then
            icy=iy
         else 
            icy=mod(Lgy+1-iy,Lgy)+1
         endif
         if (iz<=Lgz2) then
            icz=iz
         else 
            icz=mod(Lgz+1-iz,Lgz)+1
         endif
         d3(ix,iy,iz) = d3(ix,iy,iz) * g3(icx,icy,icz)
   	     r3(ix,iy,iz) = r3(ix,iy,iz) * g3(icx,icy,icz)
      enddo
   enddo
enddo

call fft3d(d3,r3,ndim,-1)

deallocate(r3,STAT=istat)

return
end subroutine cnvl3d

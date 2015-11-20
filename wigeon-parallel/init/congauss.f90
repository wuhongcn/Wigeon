subroutine contraintGauss(z,rLbox,L,d3)
use set_parameter, only:omega0
real,allocatable,dimension(:,:,:):: den,r3
real d3(L,L,L)
integer,dimension(3):: nn,n3
real,dimension(3)::frc
      
pi=acos(-1.)
tpi=2.*pi

Lc=L/2+1
ndim=3
nn(1)=L
nn(2)=L
nn(3)=L
rL=float(L)
tpiL=tpi/rL
rLb=rL/rLbox

allocate (den(L,L,L),STAT=istat)
allocate (r3(L,L,L),STAT=istat)

den=d3

! ---------------- set constranints --------------------- 

n3(1)=L/2
n3(2)=L/2
n3(3)=L/2

xM15=1.0
Rp=6.115785*(xM15/omega0)**(1./3.) ! h^-1Mpc
Rg=Rp*rLb

! -------------------------------------------------------

ixx=n3(1)
iyy=n3(2)
izz=n3(3)

call variance(z,0,Rg,var)
v0=var
call variance(z,1,Rg,var)
v1=var

call smooth3D(L,rad,den)
den_smooth=den(ixx,iyy,izz)

den=d3

c0=(3.*sqrt(v0)-den_smooth)/v0

call gridf10_cell(L,n3,den,frc)

c1=-3.*frc(1)/v1
c2=-3.*frc(2)/v1
c3=-3.*frc(3)/v1

do iz=1,Lc
   icz=mod(L+1-iz,L)+1
   do iy=1,Lc
      icy=mod(L+1-iy,L)+1
      do ix=1,Lc
         icx=mod(L+1-ix,L)+1
         rk=tpiL*sqrt(float((ix-1)**2+(iy-1)**2+(iz-1)**2))
         rkrl=rk*rLb
         amp=CDMpower(rkrl,z)*exp(-0.5*rkrl*rkrl*Rp*Rp)

         xk=tpiL*float(ix-1)
		 yk=tpiL*float(iy-1)
		 zk=tpiL*float(iz-1)
         
		 ck=(c1*xk+c2*yk+c3*zk)*rLb
		 rkd=xk*float(ixx)+yk*float(iyy)+zk*float(izz)

         d3(ix ,iy ,iz )=amp*(c0*cos(rkd)+ck*sin(rkd))
         r3(ix ,iy ,iz )=amp*(c0*sin(rkd)-ck*cos(rkd))
         d3(icx,icy,icz)= d3(ix ,iy ,iz )
         r3(icx,icy,icz)=-r3(ix ,iy ,iz )
         
		 ck=(c1*xk+c2*yk-c3*zk)*rLb
		 rkd=xk*float(ixx)+yk*float(iyy)-zk*float(izz)

         d3(ix ,iy ,icz)=amp*(c0*cos(rkd)+ck*sin(rkd))
         r3(ix ,iy ,icz)=amp*(c0*sin(rkd)-ck*cos(rkd))
         d3(icx,icy,iz )= d3(ix ,iy ,icz)
         r3(icx,icy,iz )=-r3(ix ,iy ,icz)

		 ck=(c1*xk-c2*yk+c3*zk)*rLb
		 rkd=xk*float(ixx)-yk*float(iyy)+zk*float(izz)

         d3(ix ,icy,iz )=amp*(c0*cos(rkd)+ck*sin(rkd))
         r3(ix ,icy,iz )=amp*(c0*sin(rkd)-ck*cos(rkd))
         d3(icx,iy ,icz)= d3(ix ,icy,iz )
         r3(icx,iy ,icz)=-r3(ix ,icy,iz )

		 ck=(c1*xk-c2*yk-c3*zk)*rLb
		 rkd=xk*float(ixx)-yk*float(iyy)-zk*float(izz)
         d3(ix ,icy,icz)=amp*(c0*cos(rkd)+ck*sin(rkd))
         r3(ix ,icy,icz)=amp*(c0*sin(rkd)-ck*cos(rkd))
         d3(icx,iy ,iz )= d3(ix ,icy,icz)
         r3(icx,iy ,iz )=-r3(ix ,icy,icz)

      enddo
   enddo
enddo

call fft3d(d3,r3,nn,+1)

d3=d3+den

deallocate(r3,STAT=istat)
deallocate(den,STAT=istat)

return
end subroutine contraintGauss

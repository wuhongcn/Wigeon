subroutine zeldovich(z,N,L,rLbox,r,v,d3,p3)
real,dimension(3,N):: r,v
real,dimension(L,L,L)::d3,p3
real,allocatable,dimension(:,:,:)::gk

iseed0=15001
raninit=ran3(iseed0)      
pi=acos(-1.)
tpi=2.*pi
Lc=L/2+1
v3=float(L*L*L)

allocate (gk(Lc,Lc,Lc),STAT=istat)

! generating the initial matter distribution in the linear regime
gk(:,:,:)=1.
!call potinOpt(z,rLbox,L,d3,gk)
call potin3D(z,rLbox,L,d3,gk)

! generating Green function
call greenOpt(gk,L)
!call greenz(gk,L)
write(*,'('' generate green function'')')

! producing uniformly distribution of DM particle
call initr(0,N,L,r)
!call initr(1,N,L,r)
write(*,'('' initialize particle distribution'')')

! tabulating the gravitational potential  
p3=d3
call cnvl3d(L,L,L,p3,gk)
p3=p3/v3
write(*,'('' generate mesh gravitational potential'')')

! Generating the initial DM particle distribution under ZA
call displ(N,L,r,v,p3)

deallocate (gk,STAT=istat)

return
end subroutine zeldovich

! -----------------------------------------------------------------------------

subroutine initr(iru,N,L,r)
real r(3,N)


rL=float(L)

n1d=nint(float(N)**(1./3.))
spc=rL/float(n1d)

if (iru.eq.0) then
   n1d2=n1d*n1d
   do i=1,N
      r(1,i)=mod(i-1,n1d)*spc+0.5
      r(2,i)=mod((i-1)/n1d,n1d)*spc+0.5
      r(3,i)=mod((i-1)/n1d2,n1d)*spc+0.5
   enddo
else if(iru.eq.1) then
   do i=1,N
      r(1,i)=ran3(iseed)*rL+1.
      r(2,i)=ran3(iseed)*rL+1.
      r(3,i)=ran3(iseed)*rL+1.
   enddo
end if

return
end

subroutine displ(N,L,r,v,p3)
real,dimension(L,L,L):: p3
real,allocatable,dimension(:,:,:)::frc
real,dimension(3,N):: r,v
iseed=192701
v=0.

rL=float(L)
rL1=10.*rL-1.
call cofficient_velocity(hst)

write(*,'(''hst = '',f12.5)')hst

allocate (frc(L,L,L),STAT=istat)

do ic=1,3
   call gridf10(ic,p3,frc,L)
   call interp(ic,N,r,v,L,frc)
enddo

deallocate (frc,STAT=istat)

do i=1,N
   do ic=1,3
      r(ic,i)=mod(r(ic,i)+v(ic,i)+rL1,rL)+1.
	  v(ic,i)=hst*v(ic,i)
   enddo
enddo
!v=v*(-1.0)
!p3=hst*p3
!p3=p3*(-1.0)
return
end




subroutine potinOpt(z,rLbox,L,d3,pf)
real,allocatable,dimension(:,:,:):: r3
real d3(L,L,L)
real pf(L/2+1,L/2+1,L/2+1)
integer nn(3)
data iseed/-130007/
      
pi=acos(-1.)
tpi=2.*pi
rke=0.8*pi
Ne=16

Lc=L/2+1
ndim=3
nn(1)=L
nn(2)=L
nn(3)=L
rL=float(L)
tpiL=tpi/rL
rLb=rL/rLbox
iseedt=iseed

allocate (r3(L,L,L),STAT=istat)
do iz=1,L
   if(iz>L/2+1) then
      icz=mod(L+1-iz,L)+1
   else
      icz=iz
   endif
  do iy=1,L
    if(iy>L/2+1) then
      icy=mod(L+1-iy,L)+1
    else 
      icy=iy
    endif
      do ix=1,L
          if(ix>L/2+1) then
             icx=mod(L+1-ix,L)+1
          else 
             icx=ix 
          endif
        
         rk=tpiL*sqrt(float((icx-1)**2+(icy-1)**2+(icz-1)**2))
         rkrl=rk*rLb
         expf=exp(-(rk/rke)**Ne)
         amp=sqrt(CDMpower(rkrl,z))*expf
         amp=amp*pf(icx,icy,icz)
         
         arg=tpi*ran3(iseedt)
         d3(ix,iy,iz)=amp*cos(arg)
         r3(ix,iy,iz)=amp*sin(arg) 
       enddo
    enddo
enddo
call fft3d(d3,r3,nn,+1)

deallocate (r3,STAT=istat)

return
end subroutine potinOpt


subroutine potin3D(z,rLbox,L,d3,pf)
real,allocatable,dimension(:,:,:):: r3
real d3(L,L,L)
real pf(L/2+1,L/2+1,L/2+1)
integer nn(3)
data iseed/-208007/
      
pi=acos(-1.)
tpi=2.*pi
rke=0.8*pi
Ne=16

Lc=L/2+1
ndim=3
nn(1)=L
nn(2)=L
nn(3)=L
rL=float(L)
tpiL=tpi/rL
rLb=rL/rLbox
iseedt=iseed

allocate (r3(L,L,L),STAT=istat)

do iz=1,L
   if(iz>L/2+1) then
      icz=mod(L+1-iz,L)+1
   else
      icz=iz
   endif
  do iy=1,L
    if(iy>L/2+1) then
      icy=mod(L+1-iy,L)+1
    else 
      icy=iy
    endif
      do ix=1,L
          if(ix>L/2+1) then
             icx=mod(L+1-ix,L)+1
          else 
             icx=ix 
          endif
        
         rk=tpiL*sqrt(float((icx-1)**2+(icy-1)**2+(icz-1)**2))
         rkrl=rk*rLb
         expf=exp(-(rk/rke)**Ne)
         amp=sqrt(CDMpower(rkrl,z))*expf
         amp=amp*pf(icx,icy,icz)

         ampran=amp*gasdev(iseedt)
         arg1=tpi*ran3(iseedt)
         d3(ix ,iy ,iz )=ampran*cos(arg1)
         r3(ix ,iy ,iz )=ampran*sin(arg1)
      enddo
   enddo
enddo

call fft3d(d3,r3,nn,+1)

deallocate (r3,STAT=istat)

return
end subroutine potin3D


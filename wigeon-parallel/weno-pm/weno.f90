! ------------------------- WENO Computation --------------------------------

module weno_kernel 
use datablock, only: uc,rhs,mn,mnp,gamma,gm1,Lx,Ly,Lz,ark
integer io

contains

subroutine weno_hydrodynamics(io)


! write(*,'('' start weno'')') 

! time_init=omp_get_wtime()
if (Lx>1) call fx
! write(*,'('' computing fx'')')
! call CPU_screentime(time_init,'CPU time in Fx')
if (Ly>1) call gy
! write(*,'('' computing gy'')')
! call CPU_screentime(time_init,'CPU time in Fy')

if (Lz>1) call hz 
! write(*,'('' computing hz'')')
! call CPU_screentime(time_init,'CPU time in Fz')
 
end subroutine weno_hydrodynamics

subroutine fx
real,dimension(mnp,Lx)::fh,u1d

!     -----------------------------------------------------------------------
!     Name:      fx.f
!     Function:  approximate "-df/dx"
!     -----------------------------------------------------------------------

!     -------------- begin of outer loops in y and z directions -------------

qk=ark(io+1)
!(2007/07/31 the above expression moved from the line below the  OpenMP  directive here )

!$omp parallel do default(shared) private(i,j,k,in1,u1d,fh)

do k=1,Lz
   do j=1,Ly
   
      u1d(:,:)=uc(:,:,j,k)
      call weno(1,Lx,u1d,fh)

      do i = 1, Lx
         in1=mod(i-2+Lx,Lx)+1
         rhs(:,i,j,k) = qk * rhs(:,i,j,k) + fh(:,in1) - fh(:,i)
      enddo

   enddo
enddo
!$omp end parallel do

!     ------------ end  of outer loops in y and z directions -----------

return
end subroutine fx

subroutine gy
real,dimension(mnp,Ly)::fh,u1d

!     --------------------------------------------------------------------
!     Name:      gy.f
!     Function:  approximate "-dg/dy"
!     --------------------------------------------------------------------

!     -------------- begin of outer loops in x and z directions ----------

!$omp parallel do default(shared) private(i,j,k,in1,u1d,fh)

do k=1,Lz
   do j=1,Lx

      u1d((/1,3,2,4,5,6/),:)=uc(:,j,:,k)
      call weno(2,Ly,u1d,fh)

      do i = 1, Ly
         in1=mod(i-2+Ly,Ly)+1
         rhs(:,j,i,k) = rhs(:,j,i,k) + fh(:,in1) - fh(:,i) 
      enddo
   
   enddo
enddo

!$omp end parallel do 

!     --------------  end  of outer loops in x and z directions -----------

return
end subroutine gy

subroutine hz
real,dimension(mnp,Lz)::fh,u1d

!     ----------------------------------------------------------------------
!     Name:      hz.f
!     Function:  approximate "-dh/dz"
!     ----------------------------------------------------------------------

!     --------------- begin of outer loops in x and y directions ----------

!$omp parallel do default(shared) private(i,j,k,in1,u1d,fh)

do j=1,Ly
   do k=1,Lx
    
      u1d((/1,4,3,2,5,6/),:)=uc(:,k,j,:)
      call weno(3,Lz,u1d,fh)

      do i = 1, Lz
         in1=mod(i-2+Lz,Lz)+1
         rhs(:,k,j,i) = rhs(:,k,j,i) + fh(:,in1) - fh(:,i) 
      enddo
   
   enddo
enddo

!$omp end parallel do

!     ---------------  end  of outer loops in x and y directions ------------

return
end subroutine hz

subroutine weno(idim,Lq,u1d,fh)
real,parameter::ama=1.10
real,dimension(mnp,Lq)::u1d,fd,ud,fh
real,dimension(mn,mn,Lq)::evl,evr
real,dimension(Lq)::w,vx,vy,vz,h,cs,cmm
!real,dimension(mn)::am
real,dimension(mn)::amalpha
real,dimension(mn,Lq)::am

!am(1)=1.d-15
!am(2)=1.d-15
!am(5)=1.d-15
amalpha(1)=1.d-15
amalpha(2)=1.d-15
amalpha(5)=1.d-15
am(1,:)=1.d-15
am(2,:)=1.d-15
am(5,:)=1.d-15


do i=1, Lq
   den=u1d(1,i)
   xmt=u1d(2,i)
   ymt=u1d(3,i)
   zmt=u1d(4,i)
   eng=u1d(5,i)
   ent=u1d(6,i)
   t0=1./den
   vex=xmt*t0
   vey=ymt*t0
   vez=zmt*t0
   if(den<1.0d-12) then
       vex=0.0
       vey=0.0
       vez=0.0
  endif
   enk=0.5*(xmt*vex+ymt*vey+zmt*vez)
   pre=gm1*(eng-enk)
!   pre=gm1*(eng-0.5*(xmt*vex+ymt*vey+zmt*vez))
!   pre=gm1*(eng-0.5*abs(xmt*vex+ymt*vey+zmt*vez))
   
   cvel=sqrt(abs(gamma*gm1*(eng-abs(enk))*t0))
   cs(i)=cvel

   fd(1,i)=xmt
   fd(2,i)=xmt*vex+pre
   fd(3,i)=ymt*vex
   fd(4,i)=zmt*vex
   fd(5,i)=vex*(pre+eng)
   fd(6,i)=ent*vex

   ud(1,i)=den
   ud(2,i)=xmt
   ud(3,i)=ymt
   ud(4,i)=zmt
   ud(5,i)=eng
   ud(6,i)=ent

   w(i)=sqrt(abs(den))
   vx(i)=vex
   vy(i)=vey
   vz(i)=vez
!  h(i)=(pre+eng)*t0
   h(i)=cs(i)*cs(i)/gm1+0.5*(vex*vex+vey*vey+vez*vez)

!   ip1=mod(i-3+Lq-1,Lq)+1
!   im1=mod(i+2+Lq-1,Lq)+1
    amlimit=1.0
!   do ia=i-1,i+2
!      ip=mod(ia+Lq-1,Lq)+1
!      am(1,ip)=max(am(1,ip),amlimit*abs(vex-cvel))
!      am(2,ip)=max(am(2,ip),amlimit*abs(vex))
!      am(5,ip)=max(am(5,ip),amlimit*abs(vex+cvel))
!   enddo   
!write(*,*) am(:,ip)
!pause
!   am(1)=max(am(1), abs(vex-cvel))
!   am(2)=max(am(2), abs(vex))
!   am(5)=max(am(5), abs(vex+cvel))
   amalpha(1)=max(amalpha(1), abs(vex-cvel))
   amalpha(2)=max(amalpha(2), abs(vex))
   amalpha(5)=max(amalpha(5), abs(vex+cvel))

enddo

!am(1,:)=am(1,:)*ama
!am(2,:)=am(2,:)*ama
!am(3,:)=am(2,:)
!am(4,:)=am(2,:)
!am(5,:)=am(5,:)*ama
am(1,:)=amalpha(1)*ama
am(2,:)=amalpha(2)*ama
am(3,:)=am(2,:)
am(4,:)=am(2,:)
am(5,:)=amalpha(5)*ama
!am(1)=am(1)*ama
!am(2)=am(2)*ama
!am(3)=am(2)
!am(4)=am(2)
!am(5)=am(5)*ama
do i=1,Lq  
!  Compute e'vectors using Roe's average:
   ip=mod(i+Lq,Lq)+1
   !t0=w(i)/(w(i)+w(ip))
   t0=w(i)/(w(i)+w(ip))+1.0d-15
   t1=1.0-t0
   vxm=t0*vx(i)+t1*vx(ip)
   vym=t0*vy(i)+t1*vy(ip)
   vzm=t0*vz(i)+t1*vz(ip)
   hm=t0*h(i)+t1*h(ip)
   qm=0.5*(vxm*vxm+vym*vym+vzm*vzm)
   vxx=vx(i)-vx(ip)
   vyy=vy(i)-vy(ip)
   vzz=vz(i)-vz(ip)
   vcc=0.5*t0*t1*(vxx*vxx+vyy*vyy+vzz*vzz)
   ccc=t0*cs(i)*cs(i)+t1*cs(ip)*cs(ip)
!  cm=sqrt(abs(gm1*(hm-qm)))
   cm=sqrt(ccc+gm1*vcc)+1.0d-15
   t0=vxm*cm

   evr(1,1,i)=1.0
   evr(1,2,i)=0.0
   evr(1,3,i)=0.0
   evr(1,4,i)=1.0
   evr(1,5,i)=1.0
   evr(2,1,i)=vxm-cm
   evr(2,2,i)=0.0
   evr(2,3,i)=0.0
   evr(2,4,i)=vxm
   evr(2,5,i)=vxm+cm
   evr(3,1,i)=vym
   evr(3,2,i)=1.0
   evr(3,3,i)=0.0
   evr(3,4,i)=vym
   evr(3,5,i)=vym
   evr(4,1,i)=vzm
   evr(4,2,i)=0.0
   evr(4,3,i)=1.0
   evr(4,4,i)=vzm
   evr(4,5,i)=vzm
   evr(5,1,i)=hm-t0
   evr(5,2,i)=vym
   evr(5,3,i)=vzm
   evr(5,4,i)=qm
   evr(5,5,i)=hm+t0
                 
   rcm=1./cm
   cmm(i)=rcm
   b1=gm1
   b2=qm*b1
   t0=vxm*cm
   t1=b1*vxm
   t2=0.5*b1
   t3=b1*vym
   t4=b1*vzm
   cm2=cm*cm

   evl(1,1,i)=0.5*(b2+t0)
   evl(1,2,i)=-0.5*(cm+t1)
   evl(1,3,i)=-0.5*t3
   evl(1,4,i)=-0.5*t4
   evl(1,5,i)=t2
   evl(2,1,i)=-vym*cm2
   evl(2,2,i)=0.0
   evl(2,3,i)=1.0*cm2
   evl(2,4,i)=0.0
   evl(2,5,i)=0.0
   evl(3,1,i)=-vzm*cm2
   evl(3,2,i)=0.0
   evl(3,3,i)=0.0
   evl(3,4,i)=1.0*cm2
   evl(3,5,i)=0.0
   evl(4,1,i)=cm2-b2
   evl(4,2,i)=t1
   evl(4,3,i)=t3
   evl(4,4,i)=t4
   evl(4,5,i)=-b1
   evl(5,1,i)=0.5*(b2-t0)
   evl(5,2,i)=0.5*(cm-t1)
   evl(5,3,i)=-0.5*t3
   evl(5,4,i)=-0.5*t4
   evl(5,5,i)=t2

enddo

select case(idim)
case(1)
    call weno5(Lx,fd,ud,evr,evl,am,cmm,fh)
case(2) 
    call weno5(Ly,fd((/1,3,2,4,5,6/),:),ud((/1,3,2,4,5,6/),:), &
   &       evr((/1,3,2,4,5/),:,:),evl(:,(/1,3,2,4,5/),:),am,cmm,fh ) 
case(3)
    call weno5(Lz,fd((/1,4,3,2,5,6/),:),ud((/1,4,3,2,5,6/),:), &
   &     evr((/1,4,3,2,5/),:,:),evl(:,(/1,4,3,2,5/),:),am,cmm,fh )
end select

return

end subroutine weno


subroutine weno5(LL,f,u,vecr,vecl,am,cmm,fh)
real,parameter::epweno=1.d-6

!     ------------------------------------------------------------------------
!     Name:      wenolf.f
!     Function:  Use WENO-LF-4 or WENO-LF-5 to approximate fluxes
!     ------------------------------------------------------------------------
real,dimension(mnp,LL)::f,u,fh
real,dimension(mn,mn,LL)::vecr,vecl
real,dimension(mnp,LL):: df,du,ff
!real,dimension(mn)::am
real,dimension(mn,LL)::am
real,dimension(LL)::cmm
real :: gg(2,mn,LL), hh(2,4,LL),h4(4)
real :: tt1,tt2,tt3,s1,s2,s3,t12,t13,t31,t32
      
do i = 1,LL 
   ip1=mod(i+LL,LL)+1 
   df(1:mn,i) = f(1:mn,ip1) - f(1:mn,i)
   du(1:mn,i) = u(1:mn,ip1) - u(1:mn,i)
enddo

!     ----------------- loop in "m" starts here  -------------------

do m = 1, mn

!     use Lax-Friedrichs building block to split the fluxes
         
   do i =1, LL
      do m1 = 1, mn
!         gg(1,m1,i) = 0.5 * ( df(m1,i) + am(m1) * du(m1,i) )
         gg(1,m1,i) = 0.5 * ( df(m1,i) + am(m1,i) * du(m1,i) )
         gg(2,m1,i) = gg(1,m1,i) - df(m1,i)
      enddo
   enddo
 

!     Project the positive and negative part of the fluxes to the
!     'm'th characteristic field

   do i = 1, LL
      do m1 = 1, 4
         k0 = m1 - 3
         k1 =  3 - m1
         ipk0=mod(i+k0-1+LL,LL)+1
         ipk1=mod(i+k1-1+LL,LL)+1 
         hh(1,m1,i) = vecl(m,1,i)*gg(1,1,ipk0) &
     &              + vecl(m,2,i)*gg(1,2,ipk0) &
     &              + vecl(m,3,i)*gg(1,3,ipk0) &
     &              + vecl(m,4,i)*gg(1,4,ipk0) &
     &              + vecl(m,5,i)*gg(1,5,ipk0)
         hh(2,m1,i) = vecl(m,1,i)*gg(2,1,ipk1) &
     &              + vecl(m,2,i)*gg(2,2,ipk1) &
     &              + vecl(m,3,i)*gg(2,3,ipk1) &
     &              + vecl(m,4,i)*gg(2,4,ipk1) &
     &              + vecl(m,5,i)*gg(2,5,ipk1)
      enddo
   enddo
!     compute the weights and approximate the fluxes

   do i = 1, LL
      ff(m,i) = 0.0
   enddo

   do i = 1, LL
      do m1 = 1, 2

         t1 = hh(m1,1,i) - hh(m1,2,i)
         t2 = hh(m1,2,i) - hh(m1,3,i)
         t3 = hh(m1,3,i) - hh(m1,4,i)

         tt1 = 13. * t1**2 + 3. * (    hh(m1,1,i) - 3.*hh(m1,2,i) )**2
         tt2 = 13. * t2**2 + 3. * (    hh(m1,2,i) +    hh(m1,3,i) )**2
         tt3 = 13. * t3**2 + 3. * ( 3.*hh(m1,3,i) -    hh(m1,4,i) )**2

         tt1 =  epweno + tt1  
         tt2 =  epweno + tt2  
         tt3 =  epweno + tt3  
 
         t12=tt1/tt2
         t13=tt1/tt3
         t31=1./t13
         t32=tt3/tt2

         t12=t12*t12
         t13=t13*t13
         t31=t31*t31
         t32=t32*t32

         s1=1./(1. + 6.* t12 + 3. * t13 )  
         s3=1./(1. + t31 / 3.  + 2. * t32 )
         ff(m,i)=ff(m,i)+(s1*(t2-t1)+(0.5*s3-0.25)*(t3-t2))/3.
      enddo
   enddo
         
enddo

!     ----------------- loop in "m"  ends  here  -------------------

!     Project the fluxes to the physical space:

do i = 1, LL

   in1=mod(i-2+LL,LL)+1
   ip1=mod(i  +LL,LL)+1
   ip2=mod(i+1+LL,LL)+1 

   do m = 1, mn

      fh(m,i) = ( vecr(m,1,i) * ff(1,i) + vecr(m,2,i) * ff(2,i) &
     &        +   vecr(m,3,i) * ff(3,i) + vecr(m,4,i) * ff(4,i) &
     &        +   vecr(m,5,i) * ff(5,i)) * cmm(i) * cmm(i)      &
     &        + ( -f(m,in1) + 7.*( f(m,i)+f(m,ip1) )-f(m,ip2) )/12.

   enddo

enddo

! **************************** Entropy Equation ******************************      

if (mnp==mn) return

do ic=mn+1,mnp

!---------------------- Computing flux_x,y,z ---------------------------------

! first compute the fp part

do i = 1,LL
   ip1=mod(i+LL,LL)+1
!   df(ic,i) = 0.5*(f(ic,ip1) - f(ic,i) + am(2)*(u(ic,ip1) - u(ic,i)))

   df(ic,i) = 0.5*(f(ic,ip1) - f(ic,i) + am(2,i)*(u(ic,ip1) - u(ic,i)))
enddo

do i = 1,LL

   ip1=mod(i  +LL,LL)+1
   ip2=mod(i+1+LL,LL)+1
   in1=mod(i-2+LL,LL)+1
   in2=mod(i-3+LL,LL)+1

   h4(1) =  df(ic,in2)
   h4(2) =  df(ic,in1)
   h4(3) =  df(ic,  i)
   h4(4) =  df(ic,ip1)

   fh(ic,i)=(-f(ic,in1)+7.*(f(ic,i)+f(ic,ip1))-f(ic,ip2))/12.

   t1 = h4(1) - h4(2)
   t2 = h4(2) - h4(3)
   t3 = h4(3) - h4(4)

   tt1 = 13 * t1 * t1 + 3 * (   h4(1) - 3*h4(2) )**2
   tt2 = 13 * t2 * t2 + 3 * (   h4(2) +   h4(3) )**2
   tt3 = 13 * t3 * t3 + 3 * ( 3*h4(3) -   h4(4) )**2

   tt1 =  epweno + tt1
   tt2 =  epweno + tt2
   tt3 =  epweno + tt3  
   
   t12=tt1/tt2
   t13=tt1/tt3
   t31=1./t13
   t32=tt3/tt2

   t12=t12*t12
   t13=t13*t13
   t31=t31*t31
   t32=t32*t32
   
   s1=1./(1. + 6.* t12 + 3. * t13 )
   s3=1./(1. + t31 / 3.  + 2. * t32 )

   fh(ic,i) = fh(ic,i)+( s1*(t2-t1)+(0.5*s3-0.25)*(t3-t2))/3.

enddo

! then compute the fn part

do i = 1,LL
   ip1=mod(i+LL,LL)+1
!   df(ic,i) = 0.5*(f(ic,ip1) - f(ic,i) - am(2)*(u(ic,ip1) - u(ic,i)))
   df(ic,i) = 0.5*(f(ic,ip1) - f(ic,i) - am(2,i)*(u(ic,ip1) - u(ic,i)))
enddo

do i = 1, LL

   ip1=mod(i  +LL,LL)+1
   ip2=mod(i+1+LL,LL)+1
   in1=mod(i-2+LL,LL)+1
   in2=mod(i-3+LL,LL)+1

   h4(1) = -df(ic,ip2)
   h4(2) = -df(ic,ip1)
   h4(3) = -df(ic,  i)
   h4(4) = -df(ic,in1)

   t1 = h4(1) - h4(2)
   t2 = h4(2) - h4(3)
   t3 = h4(3) - h4(4)

   tt1 = 13 * t1 * t1 + 3 * (   h4(1) - 3*h4(2) )**2
   tt2 = 13 * t2 * t2 + 3 * (   h4(2) +   h4(3) )**2
   tt3 = 13 * t3 * t3 + 3 * ( 3*h4(3) -   h4(4) )**2

   tt1 =  epweno + tt1
   tt2 =  epweno + tt2
   tt3 =  epweno + tt3

   t12=tt1/tt2
   t13=tt1/tt3
   t31=1./t13
   t32=tt3/tt2

   t12=t12*t12
   t13=t13*t13
   t31=t31*t31
   t32=t32*t32
   
   s1=1./(1. + 6.* t12 + 3. * t13 )
   s3=1./(1. + t31 / 3.  + 2. * t32 )
    
   fh(ic,i)=fh(ic,i)+(s1*(t2-t1)+(0.5*s3-0.25)*(t3-t2))/3.

enddo

enddo

return
end subroutine weno5


end module weno_kernel



      

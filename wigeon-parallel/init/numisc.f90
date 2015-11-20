! ----------------------------------------------------------------------------

SUBROUTINE qromb(func,a,b,ss)
!     romberg integrator
!     uses polint,trapzd
integer j,jmax,jmaxp,k,km
real a,b,func,ss,eps
parameter (eps=1.d-6, jmax=20, jmaxp=jmax+1, k=5, km=k-1)
real dss,h(jmaxp),s(jmaxp)
external func

h(1)=1.
do j=1,jmax
   call trapzd(func,a,b,s(j),j)
   if (j.ge.k) then
      call polint(h(j-km),s(j-km),k,0.,ss,dss)
      if (abs(dss).le.eps*abs(ss)) return
   endif
   s(j+1)=s(j)
   h(j+1)=0.25*h(j)
enddo
write(*,*) 'too many steps in qromb'
write(*,*) 'ss = ',ss,', dss = ',dss
PAUSE 'qromb error'
END

SUBROUTINE trapzd(func,a,b,s,n)
! trapezium rule integrator
integer n
real a,b,s,func
external func
integer it,j
real del,sum,tnm,x
if (n.eq.1) then
   s=0.5*(b-a)*(func(a)+func(b))
else
   it=2**(n-2)
   tnm=it
   del=(b-a)/tnm
   x=a+0.5*del
   sum=0.
   do j=1,it
      sum=sum+func(x)
      x=x+del
   enddo
   s=0.5*(s+(b-a)*sum/tnm)
endif
RETURN
END

SUBROUTINE polint(xa,ya,n,x,y,dy)
! polynomial interpolation
integer n,nmax
real dy,x,y,xa(n),ya(n)
parameter (nmax=10)
integer i,m,ns
real den,dif,dift,ho,hp,w,c(nmax),d(nmax)

ns=1
dif=abs(x-xa(1))
do i=1,n
   dift=abs(x-xa(i))
   if (dift.lt.dif) then
      ns=i
      dif=dift
   endif
   c(i)=ya(i)
   d(i)=ya(i)
enddo
y=ya(ns)
ns=ns-1
do m=1,n-1
   do i=1,n-m
      ho=xa(i)-x
      hp=xa(i+m)-x
      w=c(i+1)-d(i)
      den=ho-hp
      if (den.eq.0.) PAUSE 'failure in polint'
      den=w/den
      d(i)=hp*den
      c(i)=ho*den
   enddo
   if (2*ns.lt.n-m) then
      dy=c(ns+1)
   else
      dy=d(ns)
      ns=ns-1
   endif
   y=y+dy
enddo

RETURN
end

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
!status=DftiComputeForward(my_fft,fft3)  modified ,Nov,24th,2007
if(isn==1) then
    status=DftiComputeForward(my_fft,fft3)
else  if(isn==-1) then
    status=DftiComputeBackward(my_fft,fft3)
endif    
status=DftiFreeDescriptor(my_fft)

a(1:ns)=real(fft3(1:ns))
b(1:ns)=imag(fft3(1:ns))

deallocate (fft3,STAT=istat)

return
end subroutine fft3d

FUNCTION RAN3(IDUM)
SAVE
PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
DIMENSION MA(55)
DATA IFF /0/
IF (IDUM.LT.0.OR.IFF.EQ.0)THEN
   IFF=1
   MJ=MSEED-IABS(IDUM)
   MJ=MOD(MJ,MBIG)
   MA(55)=MJ
   MK=1
   DO 11 I=1,54
      II=MOD(21*I,55)
      MA(II)=MK
      MK=MJ-MK
      IF (MK.LT.MZ) MK=MK+MBIG
      MJ=MA(II)
11 CONTINUE
   DO 13 K=1,4
      DO 12 I=1,55
         MA(I)=MA(I)-MA(1+MOD(I+30,55))
         IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
12    CONTINUE
13 CONTINUE
   INEXT=0
   INEXTP=31
   IDUM=1
ENDIF
INEXT=INEXT+1
IF(INEXT.EQ.56)INEXT=1
INEXTP=INEXTP+1
IF(INEXTP.EQ.56)INEXTP=1
MJ=MA(INEXT)-MA(INEXTP)
IF(MJ.LT.MZ)MJ=MJ+MBIG
MA(INEXT)=MJ
RAN3=MJ*FAC
RETURN
END

integer FUNCTION IntegerRead(text)
Character text*(*)
integer ix
write (*,'(A,$)')text
read (*,*) ix
IntegerRead=ix
Return
End

integer function lnblnk(char)
character char*(*)
lnblnk=index(char,' ')-1
return
end


      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
!        do 11 j=NTAB+8,1,-1
        do j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
        enddo
!11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END

FUNCTION gasdev(idum)
INTEGER idum
REAL gasdev
!    USES ran1
INTEGER iset
REAL fac,gset,rsq,v1,v2,ran1
SAVE iset,gset
DATA iset/0/
   if (iset.eq.0) then
1     v1=2.*ran1(idum)-1.
      v2=2.*ran1(idum)-1.
      rsq=v1**2+v2**2
      if(rsq.ge.1..or.rsq.eq.0.)goto 1
      fac=sqrt(-2.*log(rsq)/rsq)
      gset=v1*fac
      gasdev=v2*fac
      iset=1
   else
      gasdev=gset
      iset=0
endif
return
END

SUBROUTINE gauleg(x1,x2,x,w,n)
INTEGER n
REAL x1,x2,x(n),w(n)
DOUBLE PRECISION EPS
PARAMETER (EPS=3.d-14)
INTEGER i,j,m
DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
m=(n+1)/2
xm=0.5d0*(x2+x1)
xl=0.5d0*(x2-x1)
do 12 i=1,m
   z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1  continue
   p1=1.d0
   p2=0.d0
   do 11 j=1,n
      p3=p2
      p2=p1
      p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11    continue
      pp=n*(z*p1-p2)/(z*z-1.d0)
      z1=z
      z=z1-p1/pp
      if(abs(z-z1).gt.EPS)goto 1
      x(i)=xm-xl*z
      x(n+1-i)=xm+xl*z
      w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
      w(n+1-i)=w(i)
12 continue
return
END


module recipes

contains

! ------------------ taken from Numerical Recipes ---------------------

SUBROUTINE qromb(func,a,b,ss)
!     romberg integrator
!     uses polint,trapzd
integer j,jmax,jmaxp,k,km
real a,b,func,ss,eps
parameter (eps=1.e-8, jmax=20, jmaxp=jmax+1, k=5, km=k-1)
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
END SUBROUTINE qromb

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
END SUBROUTINE trapzd

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
end SUBROUTINE polint

subroutine polin2d(X1A,X2A,YA,M,N,X1,X2,Y,DY)
integer,parameter::NMAX=20,MMAX=20
real:: X1A(M),X2A(N),YA(M,N),YNTMP(NMAX),YMTMP(MMAX)
real:: X1,X2,Y,DY

do j=1,M
   do k=1,N
      YNTMP(k)=YA(j,k)
   enddo
   call polint(X2A,YNTMP,N,X2,YMTMP(j),DY)
enddo 
call polint(X1A,YMTMP,M,X1,Y,DY)

RETURN
END subroutine polin2d

FUNCTION zbrent(func,x1,x2,tol)
INTEGER ITMAX
REAL zbrent,tol,x1,x2,func,EPS
EXTERNAL func
PARAMETER (ITMAX=100,EPS=3.e-8)
INTEGER iter
REAL a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
a=x1
b=x2
fa=func(a)
fb=func(b)
if ((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then
   write(*,'(''root must be bracketed for zbrent'')')
   pause
endif
   c=b
   fc=fb
   do iter=1,ITMAX
      if ((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
         c=a
         fc=fa

         d=b-a
         e=d
      endif
      if (abs(fc).lt.abs(fb)) then
         a=b
         b=c
         c=a
         fa=fb
         fb=fc
         fc=fa
      endif
      tol1=2.*EPS*abs(b)+0.5*tol
      xm=.5*(c-b)
      if (abs(xm).le.tol1 .or. fb.eq.0.)then
         zbrent=b
         return
      endif
      if (abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
         s=fb/fa
         if (a.eq.c) then
            p=2.*xm*s
            q=1.-s
         else

            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if (p.gt.0.) q=-q
          p=abs(p)
          if (2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
             e=d
             d=p/q
          else
             d=xm
             e=d
          endif
      else
          d=xm
          e=d
      endif
      a=b
      fa=fb
      if (abs(d) .gt. tol1) then
         b=b+d
      else
         b=b+sign(tol1,xm)
      endif
      fb=func(b)
enddo
stop 'zbrent exceeding maximum iterations'
zbrent=b
return
END FUNCTION zbrent

integer function lblnk(char)
character char*(*)
lblnk=index(char,' ')-1
return
end function lblnk

integer FUNCTION IntegerRead(text)
Character text*(*)
integer ix
write (*,'(A,$)')text
read (*,*) ix
IntegerRead=ix
Return
End FUNCTION IntegerRead

subroutine CPU_screentime(t_init,t_step,text)
character :: text*(*)
real :: t_init,t_run,t_step

t_run=omp_get_wtime()
t_step=t_run-t_init
write(*,'(A,''='',f15.6)') text,t_step
t_init=t_run

return
end subroutine CPU_screentime

end module recipes

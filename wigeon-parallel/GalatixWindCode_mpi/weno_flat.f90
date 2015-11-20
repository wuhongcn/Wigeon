! ------------------------- WENO Computation --------------------------------

module weno_kernel 
  use datablock
  use readin_out
  integer io
  
contains
  subroutine rhs_out(iii)
    character fname*6,ucout*20,sts*8
    integer iii

   write(ucout,'(''rhs'',i2.2,''.'',i2.2)')myid,iii

   open(10,file=ucout,form='unformatted')
   write(10)rhs
   close(10)

return
end subroutine rhs_out

  subroutine data_output2(iii)
    character fname*6,ucout*20,sts*8
    integer iii

if(1)then
   write(ucout,'(''uc'',i2.2,''.'',i2.2)')myid,iii
   if(myid.eq.0) write(*,'(''save the numerical solution'',''step no='',i4)')nt

   open(10,file=datadir(1:lblnk(datadir))//ucout,form='unformatted')

   do ic=1,mnp
      write(10)uc(ic,1:Lx,1:Ly,1:nLz)
   enddo

   close(10)
    if(iii.ge.8) then
       call mpi_barrier(mpi_comm_world,ierr)
       call mpi_abort(mpi_comm_world,ierr)
       call mpi_finalize(ierr)
       stop
    endif
endif
    return
  end subroutine data_output2

  subroutine data_piece(ii)
    include 'mpif.h'
    character fname*4,fname2*4

    write(fname,'(i4.4)')myid
    write(fname2,'(i4.4)')ii

    open(10,file='ppu_x.'//fname//'.'//fname2//'.-4:0')
    write(10,*)uc(1:6,-4:0,1:nLy,1:nLz)
    close(10)

    open(10,file='ppu_x.'//fname//'.'//fname2//'.1:5')
    write(10,*)uc(1:6,1:5,1:nLy,1:nLz)
    close(10)

    open(10,file='ppu_x.'//fname//'.'//fname2//'.nLx-4:nLx')
    write(10,*)uc(1:6,nLx-4:nLx,1:nLy,1:nLz)
    close(10)

    open(10,file='ppu_x.'//fname//'.'//fname2//'.nLx+1:nLx+5')
    write(10,*)uc(1:6,nLx+1:nLx+5,1:nLy,1:nLz)
    close(10)

    open(10,file='ppu_y.'//fname//'.'//fname2//'.-4:0')
    write(10,*)uc(1:6,1:nLx,-4:0,1:nLz)
    close(10)

    open(10,file='ppu_y.'//fname//'.'//fname2//'.1:5')
    write(10,*)uc(1:6,1:nLx,1:5,1:nLz)
    close(10)

    open(10,file='ppu_y.'//fname//'.'//fname2//'.nLy-4:nLy')
    write(10,*)uc(1:6,1:nLx,nLy-4:nLy,1:nLz)
    close(10)

    open(10,file='ppu_y.'//fname//'.'//fname2//'.nLy+1:nLy+5')
    write(10,*)uc(1:6,1:nLx,nLy+1:nLy+5,1:nLz)
    close(10)

    open(10,file='ppu_z.'//fname//'.'//fname2//'.-4:0')
    write(10,*)uc(1:6,1:nLx,1:nLy,-4:0)
    close(10)

    open(10,file='ppu_z.'//fname//'.'//fname2//'.1:5')
    write(10,*)uc(1:6,1:nLx,1:nLy,1:5)
    close(10)

    open(10,file='ppu_z.'//fname//'.'//fname2//'.nLz-4:nLz')
    write(10,*)uc(1:6,1:nLx,1:nLy,nLz-4:nLz)
    close(10)

    open(10,file='ppu_z.'//fname//'.'//fname2//'.nLz+1:nLz+5')
    write(10,*)uc(1:6,1:nLx,1:nLy,nLz+1:nLz+5)
    close(10)

    if(ii.ge.2) then
       call mpi_barrier(mpi_comm_world,ierr)
       call mpi_abort(mpi_comm_world,ierr)
       call mpi_finalize(ierr)
       stop
    endif
    
    return
  end subroutine data_piece

  subroutine weno_hydrodynamics(io)
    include'mpif.h'
    integer,save:: ii=1
    integer::fgh

ts=mpi_wtime()
! 1. Communicate ...

!         uc(:,1:Lx,1:Ly,-4:0)=uc(:,1:Lx,1:Ly,1)
!         uc(:,1:Lx,1:Ly,nLz+1:nLz+5)=uc(:,1:Lx,1:Ly,nLz)

    call mpi_sendrecv(uc(1,nLx-4,1,1),1,xb_type,x_right,101,&
&		uc(1,-4,1,1),1,xb_type,x_left,101,mpi_comm_world,status,ierr)
    call mpi_sendrecv(uc(1,1,1,1),1,xb_type,x_left,100,&
&		uc(1,nLx+1,1,1),1,xb_type,x_right,100,mpi_comm_world,status,ierr)

    call mpi_sendrecv(uc(1,1,nLy-4,1),1,yb_type,y_right,201, &
&		uc(1,1,-4,1),1,yb_type,y_left,201,mpi_comm_world,status,ierr)
    call mpi_sendrecv(uc(1,1,1,1),1,yb_type,y_left,200,&
&		uc(1,1,nLy+1,1),1,yb_type,y_right,200,mpi_comm_world,status,ierr)

    call mpi_sendrecv(uc(1,1,1,nLz-4),1,zb_type,z_right,301,&
&		uc(1,1,1,-4),1,zb_type,z_left,301,mpi_comm_world,status,ierr)
    call mpi_sendrecv(uc(1,1,1,1),1,zb_type,z_left,300,&
&		uc(1,1,1,nLz+1),1,zb_type,z_right,300,mpi_comm_world,status,ierr)

if (x_left.eq.MPI_PROC_NULL) then 
   do i=-4,0
      uc(:,i,1:nLy,1:nLz)=uc(:,1,1:nLy,1:nLz)
   enddo
endif

if (x_right.eq.MPI_PROC_NULL) then
   do i=nLx+1,nLx+5
      uc(:,i,1:nLy,1:nLz)=uc(:,nLx,1:nLy,1:nLz)
   enddo
endif

if (y_left.eq.MPI_PROC_NULL) then
   do i=-4,0
      uc(:,1:nLx,i,1:nLz)=uc(:,1:nLx,1,1:nLz)
   enddo
endif

if (y_right.eq.MPI_PROC_NULL) then
   do i=nLy+1,nLy+5
      uc(:,1:nLx,i,1:nLz)=uc(:,1:nLx,nLy,1:nLz)
   enddo
endif

if (z_left.eq.MPI_PROC_NULL) then
   do i=-4,0
      uc(:,1:nLx,1:nLy,i)=uc(:,1:nLx,1:nLy,1)
   enddo
endif

if (z_right.eq.MPI_PROC_NULL) then
   do i=nLz+1,nLz+5
      uc(:,1:nLx,1:nLy,i)=uc(:,1:nLx,1:nLy,nLz)
   enddo
endif

te=mpi_wtime()
if(myid.eq.0) write (*,*) "comm time:", te-ts

if(0) then
    if(myid.eq.0) then
if(0)then
       ! uc(nLz-4:nLz) -> uc_z_right(-4:0)
       call mpi_send(uc(1,1,1,nLz-4),1,zb_type,z_right,301,mpi_comm_world,ierr)
       ! uc(nLz+1:nLz+5) <- uc_z_right(1:5)
       call mpi_recv(uc(1,1,1,nLz+1),1,zb_type,z_right,300,mpi_comm_world,status,ierr)
       call mpi_send(uc(1,nLx-4,1,1),1,xb_type,x_right,101,mpi_comm_world,ierr)
       call mpi_recv(uc(1,nLx+1,1,1),1,xb_type,x_right,100,mpi_comm_world,status,ierr)
       call mpi_send(uc(1,1,nLy-4,1),1,yb_type,y_right,201,mpi_comm_world,ierr)
       call mpi_recv(uc(1,1,nLy+1,1),1,yb_type,y_right,200,mpi_comm_world,status,ierr)
endif

      do i=-4,0
          uc(:,i,1:nLy,1:nLz)=uc(:,1,1:nLy,1:nLz)
          uc(:,1:nLx,i,1:nLz)=uc(:,1:nLx,1,1:nLz)
          uc(:,1:nLx,1:nLy,i)=uc(:,1:nLx,1:nLy,1)
       enddo

    else if(myid.eq.nprocs-1) then
if(0) then
       ! uc(-4:0) <- uc_z_left(nLz-4:nLz)
       call mpi_recv(uc(1,1,1,-4),1,zb_type,z_left,301,mpi_comm_world,status,ierr)
       ! uc(1:5) -> uc_z_left(nLz+1:nLz+5)
       call mpi_send(uc(1,1,1,1),1,zb_type,z_left,300,mpi_comm_world,ierr)
       
       call mpi_recv(uc(1,-4,1,1),1,xb_type,x_left,101,mpi_comm_world,status,ierr)
       call mpi_send(uc(1,1,1,1),1,xb_type,x_left,100,mpi_comm_world,ierr)

       call mpi_recv(uc(1,1,-4,1),1,yb_type,y_left,201,mpi_comm_world,status,ierr)
       call mpi_send(uc(1,1,1,1),1,yb_type,y_left,200,mpi_comm_world,ierr)
endif

       do i=Ls+1,Ls+5
          uc(:,i,1:nLy,1:nLz)=uc(:,nLx,1:nLy,1:nLz)
          uc(:,1:nLx,i,1:nLz)=uc(:,1:nLx,nLy,1:nLz)
          uc(:,1:nLx,1:nLy,i)=uc(:,1:nLx,1:nLy,nLz)
       enddo
       
    else
if(0) then
       ! uc(-4:0) <- uc_z_left(nLz-4:nLz)
       call mpi_recv(uc(1,1,1,-4),1,zb_type,z_left,301,mpi_comm_world,status,ierr)
       ! uc(nLz-4:nLz) -> uc_z_right(-4:0)
       call mpi_send(uc(1,1,1,nLz-4),1,zb_type,z_right,301,mpi_comm_world,ierr)
       ! uc(1:5) -> uc_z_left(nLz+1:nLz+5)
       call mpi_send(uc(1,1,1,1),1,zb_type,z_left,300,mpi_comm_world,ierr)
       ! uc(nLz+1:nLz+5) <- uc_z_right(1:5)
       call mpi_recv(uc(1,1,1,nLz+1),1,zb_type,z_right,300,mpi_comm_world,status,ierr)

       call mpi_recv(uc(1,-4,1,1),1,xb_type,x_left,101,mpi_comm_world,status,ierr)
       call mpi_send(uc(1,nLx-4,1,1),1,xb_type,x_right,101,mpi_comm_world,ierr)
       call mpi_send(uc(1,1,1,1),1,xb_type,x_left,100,mpi_comm_world,ierr)
       call mpi_recv(uc(1,nLx+1,1,1),1,xb_type,x_right,100,mpi_comm_world,status,ierr)

       call mpi_recv(uc(1,1,-4,1),1,yb_type,y_left,201,mpi_comm_world,status,ierr)
       call mpi_send(uc(1,1,nLy-4,1),1,yb_type,y_right,201,mpi_comm_world,ierr)
       call mpi_send(uc(1,1,1,1),1,yb_type,y_left,200,mpi_comm_world,ierr)
       call mpi_recv(uc(1,1,nLy+1,1),1,yb_type,y_right,200,mpi_comm_world,status,ierr)
endif
    endif
endif

if(0) then
    call data_piece(ii)
    ii=ii+1
endif
 !   if(ii==2) then
    if(ii==-1) then
        open(10,file='./mydata')
        write(10,*)mnp,nLx,nLy,nLz,io,gamma,ark,uc,rhs
        close(10)
    endif

    time1=mpi_wtime() !!mengchen
    if (nLx>1) call fx
    time2=mpi_wtime() !!mengchen
    ptime=time2-time1 !!mengchen
    if(myid.eq.0)then
         write(*,'("-df/dx time =",f15.6)')ptime !!mengchen
    endif

    time1=mpi_wtime() !!mengchen
    if (nLy>1) call gy
    time2=mpi_wtime() !!mengchen
    ptime=time2-time1 !!mengchen
    if(myid.eq.0)then
        write(*,'("-dg/dy time =",f15.6)')ptime !!mengchen
    endif

    time1=mpi_wtime() !!mengchen
    if (nLz>1) call hz
    time2=mpi_wtime() !!mengchen
    ptime=time2-time1 !!mengchen
    if(myid.eq.0)then
        write(*,'("-dh/dz time =",f15.6)')ptime !!mengchen
    endif
    ii=ii+1
  end subroutine weno_hydrodynamics

subroutine fx
real,dimension(0:nLx,mnp)::fh !! zycao_0808
real,dimension(-4:nLx+5,mnp)::ud !! zycao_0808
!     -----------------------------------------------------------------------
!     Name:      fx.f
!     Function:  approximate "-df/dx"
!     -----------------------------------------------------------------------

!     -------------- begin of outer loops in y and z directions -------------

qk=ark(io+1)

!$omp parallel do default(shared) private(i,j,k,in1,u1d,fh)
do k=1,nLz
    do j=1, nLy
        do i=-4,nLx+5
            ud(i,:)=uc(:,i,j,k)
        enddo

        call weno_new(nLx,ud,fh)
        do m=1,mnp
            do i = 1, nLx
                rhs(m,i,j,k) = qk * rhs(m,i,j,k) + fh(i-1,m) - fh(i,m)
            enddo
        enddo

    enddo !! zycao_0808
!! zycao_0808  enddo
enddo
!$omp end parallel do

!     ------------ end  of outer loops in y and z directions -----------

return
end subroutine fx


subroutine gy
real,dimension(0:nLy,mnp)::fh !! zycao_0808
real,dimension(-4:nLy+5,mnp)::ud !! zycao_0808

!     --------------------------------------------------------------------
!     Name:      gy.f
!     Function:  approximate "-dg/dy"
!     --------------------------------------------------------------------

!     -------------- begin of outer loops in x and z directions ----------

!$omp parallel do default(shared) private(i,j,k,in1,u1d,fh)
do k=1,nLz
    do i=1,nLx
        do j=-4,nLy+5
            ud(j,(/1,3,2,4,5,6/))=uc(:,i,j,k)
        enddo
        call weno_new(nLy,ud,fh)
        do j = 1, nLy
                rhs(:,i,j,k) = rhs(:,i,j,k) + fh(j-1,(/1,3,2,4,5,6/)) - fh(j,(/1,3,2,4,5,6/))
        enddo
    enddo
enddo
!$omp end parallel do 

!     --------------  end  of outer loops in x and z directions -----------

return
end subroutine gy

subroutine hz
real,dimension(0:nLz,mnp)::fh !! zycao_0808
real,dimension(-4:nLz+5,mnp)::ud !! zycao_0808
!     ----------------------------------------------------------------------
!     Name:      hz.f
!     Function:  approximate "-dh/dz"
!     ----------------------------------------------------------------------

!     --------------- begin of outer loops in x and y directions ----------
!$omp parallel do default(shared) private(i,j,k,in1,u1d,fh)

do j=1,nLy
    do i=1,nLx
        do k=-4,nLz+5
            ud(k,(/1,4,3,2,5,6/))=uc(:,i,j,k)
        enddo
        call weno_new(nLz,ud,fh)
        do k = 1, nLz
            rhs(:,i,j,k) = rhs(:,i,j,k) + fh(k-1,(/1,4,3,2,5,6/)) - fh(k,(/1,4,3,2,5,6/))
        enddo
    enddo
enddo

!$omp end parallel do

!     ---------------  end  of outer loops in x and y directions ------------

return
end subroutine hz

subroutine weno_new(Lq,ud,fh)
real,parameter::ama=1.1,epweno=1.0d-6
real,dimension(0:Lq,mnp)::fh !! zycao_0808
real,dimension(-4:Lq+5,mnp)::ud !! zycao_0808
real,dimension(-4:Lq+5,mnp)::fd
real,dimension(-2:Lq+2,mnp)::gg1,gg2
real,dimension(-4:Lq+5)::cs,vx,vy,vz,w,h
real::am1,am2,am5
integer::i,j,k,ip3,ip2,ip1,i0,in1,in2

do i=-4,0
    den=ud(i,1)
    xmt=ud(i,2)
    ymt=ud(i,3)
    zmt=ud(i,4)
    eng=ud(i,5)
    tracer=ud(i,6)
    
    tv0=1./den
    vxm=xmt*tv0
    vym=ymt*tv0
    vzm=zmt*tv0

    tv2=0.5*(vxm*vxm+vym*vym+vzm*vzm)

    tv1=gm1*(eng-tv2/tv0)

    tva=xmt
    tvb=xmt*vxm+tv1
    tvc=ymt*vxm
    tvd=zmt*vxm
    tve=vxm*(tv1+eng)
    tvf=tracer*vxm

    fd(i,1)=tva
    fd(i,2)=tvb
    fd(i,3)=tvc
    fd(i,4)=tvd
    fd(i,5)=tve
    fd(i,6)=tvf

    cs(i)=sqrt(abs(gamma*tv1*tv0))
    vx(i)=vxm
    vy(i)=vym
    vz(i)=vzm
    w(i)=sqrt(abs(den))
    h(i)=cs(i)**2/gm1+tv2

enddo

do i=1,4
    ip2=i-5
    ip1=i-4
    i0 =i-3
    in1=i-2
    in2=i-1

    den=ud(i,1)
    xmt=ud(i,2)
    ymt=ud(i,3)
    zmt=ud(i,4)
    eng=ud(i,5)
    tracer=ud(i,6)

    tv0=1./den
    vxm=xmt*tv0
    vym=ymt*tv0
    vzm=zmt*tv0

    tv2=0.5*(vxm*vxm+vym*vym+vzm*vzm)

    tv1=gm1*(eng-tv2/tv0)

    tva=xmt
    tvb=xmt*vxm+tv1
    tvc=ymt*vxm
    tvd=zmt*vxm
    tve=vxm*(tv1+eng)
    tvf=tracer*vxm

    fd(i,1)=tva
    fd(i,2)=tvb
    fd(i,3)=tvc
    fd(i,4)=tvd
    fd(i,5)=tve
    fd(i,6)=tvf

    cs(i)=sqrt(abs(gamma*tv1*tv0))
    vx(i)=vxm
    vy(i)=vym
    vz(i)=vzm
    w(i)=sqrt(abs(den))
    h(i)=cs(i)**2/gm1+tv2

    ta1=max(abs(vx(ip2)-cs(ip2)),abs(vx(ip1)-cs(ip1)))
    ta2=max(abs(vx(ip2))        ,abs(vx(ip1)))
    ta3=max(abs(vx(ip2)+cs(ip2)),abs(vx(ip1)+cs(ip1)))
    ta4=max(abs(vx(i0)-cs(i0)),abs(vx(in1)-cs(in1)))
    ta5=max(abs(vx(i0))        ,abs(vx(in1)))
    ta6=max(abs(vx(i0)+cs(i0)),abs(vx(in1)+cs(in1)))
    ta7=max(abs(vx(in2)-cs(in2)),ta1)
    ta8=max(abs(vx(in2))      ,ta2)
    ta9=max(abs(vx(in2)+cs(in2)),ta3)
    taa=max(abs(vx(i)-cs(i)),ta4)
    tab=max(abs(vx(i))      ,ta5)
    tac=max(abs(vx(i)+cs(i)),ta6)
    am1=max(ta7,taa)*ama
    am2=max(ta8,tab)*ama
    am5=max(ta9,tac)*ama

    tg1 = fd(in1,1) - fd(i0,1)
    tg2 = fd(in1,2) - fd(i0,2)
    tg3 = fd(in1,3) - fd(i0,3)
    tg4 = fd(in1,4) - fd(i0,4)
    tg5 = fd(in1,5) - fd(i0,5)
    tg6 = fd(in1,6) - fd(i0,6)

    sg1 = ud(in1,1) - ud(i0,1)
    sg2 = ud(in1,2) - ud(i0,2)
    sg3 = ud(in1,3) - ud(i0,3)
    sg4 = ud(in1,4) - ud(i0,4)
    sg5 = ud(in1,5) - ud(i0,5)
    sg6 = ud(in1,6) - ud(i0,6)

    ttg1 = 0.5 * ( tg1 + am1 * sg1 )
    ttg2 = 0.5 * ( tg2 + am2 * sg2 )
    ttg3 = 0.5 * ( tg3 + am2 * sg3 )
    ttg4 = 0.5 * ( tg4 + am2 * sg4 )
    ttg5 = 0.5 * ( tg5 + am5 * sg5 )
    ttg6 = 0.5 * ( tg6 + am2 * sg6 )

    gg1(i0,1) = ttg1
    gg1(i0,2) = ttg2
    gg1(i0,3) = ttg3
    gg1(i0,4) = ttg4
    gg1(i0,5) = ttg5
    gg1(i0,6) = ttg6

    gg2(i0,1) = ttg1 - tg1
    gg2(i0,2) = ttg2 - tg2
    gg2(i0,3) = ttg3 - tg3
    gg2(i0,4) = ttg4 - tg4
    gg2(i0,5) = ttg5 - tg5
    gg2(i0,6) = ttg6 - tg6
enddo

do i=5,Lq+5
    ip4=i-7
    ip3=i-6
    ip2=i-5
    ip1=i-4
    i0 =i-3
    in1=i-2
    in2=i-1

    den=ud(i,1)
    xmt=ud(i,2)
    ymt=ud(i,3)
    zmt=ud(i,4)
    eng=ud(i,5)
    tracer=ud(i,6)

    tv0=1./den
    vxm=xmt*tv0
    vym=ymt*tv0
    vzm=zmt*tv0

    tv2=0.5*(vxm*vxm+vym*vym+vzm*vzm)

    tv1=gm1*(eng-tv2/tv0)

    tva=xmt
    tvb=xmt*vxm+tv1
    tvc=ymt*vxm
    tvd=zmt*vxm
    tve=vxm*(tv1+eng)
    tvf=tracer*vxm

    fd(i,1)=tva
    fd(i,2)=tvb
    fd(i,3)=tvc
    fd(i,4)=tvd
    fd(i,5)=tve
    fd(i,6)=tvf

    cs(i)=sqrt(abs(gamma*tv1*tv0))
    vx(i)=vxm
    vy(i)=vym
    vz(i)=vzm
    w(i)=sqrt(abs(den))
    h(i)=cs(i)**2/gm1+tv2

    ta1=max(abs(vx(ip2)-cs(ip2)),abs(vx(ip1)-cs(ip1)))
    ta2=max(abs(vx(ip2))        ,abs(vx(ip1)))
    ta3=max(abs(vx(ip2)+cs(ip2)),abs(vx(ip1)+cs(ip1)))
    ta4=max(abs(vx(i0)-cs(i0)),abs(vx(in1)-cs(in1)))
    ta5=max(abs(vx(i0))        ,abs(vx(in1)))
    ta6=max(abs(vx(i0)+cs(i0)),abs(vx(in1)+cs(in1)))
    ta7=max(abs(vx(in2)-cs(in2)),ta1)
    ta8=max(abs(vx(in2))      ,ta2)
    ta9=max(abs(vx(in2)+cs(in2)),ta3)
    taa=max(abs(vx(i)-cs(i)),ta4)
    tab=max(abs(vx(i))      ,ta5)
    tac=max(abs(vx(i)+cs(i)),ta6)
    am1=max(ta7,taa)*ama
    am2=max(ta8,tab)*ama
    am5=max(ta9,tac)*ama

    tg1 = fd(in1,1) - fd(i0,1)
    tg2 = fd(in1,2) - fd(i0,2)
    tg3 = fd(in1,3) - fd(i0,3)
    tg4 = fd(in1,4) - fd(i0,4)
    tg5 = fd(in1,5) - fd(i0,5)
    tg6 = fd(in1,6) - fd(i0,6)

    sg1 = ud(in1,1) - ud(i0,1)
    sg2 = ud(in1,2) - ud(i0,2)
    sg3 = ud(in1,3) - ud(i0,3)
    sg4 = ud(in1,4) - ud(i0,4)
    sg5 = ud(in1,5) - ud(i0,5)
    sg6 = ud(in1,6) - ud(i0,6)

    ttg1 = 0.5 * ( tg1 + am1 * sg1 )
    ttg2 = 0.5 * ( tg2 + am2 * sg2 )
    ttg3 = 0.5 * ( tg3 + am2 * sg3 )
    ttg4 = 0.5 * ( tg4 + am2 * sg4 )
    ttg5 = 0.5 * ( tg5 + am5 * sg5 )
    ttg6 = 0.5 * ( tg6 + am2 * sg6 )

    gg1(i0,1) = ttg1
    gg1(i0,2) = ttg2
    gg1(i0,3) = ttg3
    gg1(i0,4) = ttg4
    gg1(i0,5) = ttg5
    gg1(i0,6) = ttg6

    gg2(i0,1) = ttg1 - tg1
    gg2(i0,2) = ttg2 - tg2
    gg2(i0,3) = ttg3 - tg3
    gg2(i0,4) = ttg4 - tg4
    gg2(i0,5) = ttg5 - tg5
    gg2(i0,6) = ttg6 - tg6
!if ((myid.eq.12).and.(i.eq.46)) then
!     write(*,"(I3,I4,F14.9,F14.9,F14.9,F14.9,F14.9,F14.9)"), myid, i , am1, am2, am5, gg1(i0,1), gg1(i0,2), gg1(i0,5)
!call mpi_finalize(ierr)
!stop
!endif

    te0=w(ip2)/(w(ip2)+w(ip1))
    te1=1-te0

    te2=vx(ip2)
    te3=vy(ip2)
    te4=vz(ip2)
    te5=vx(ip1)
    te6=vy(ip1)
    te7=vz(ip1)
    
    tea=(te0*te2+te1*te5)
    teb=(te0*te3+te1*te6)
    tec=(te0*te4+te1*te7)

    hm=te0*h(ip2)+te1*h(ip1)
    qm=0.5*(tea*tea+teb*teb+tec*tec)
    te8=gm1*qm

    cm2=te0*(cs(ip2)**2)+te1*(cs(ip1)**2)+ &
     &      0.5*te0*te1*gm1*((te2-te5)**2+(te3-te6)**2+(te4-te7)**2)
    cm=sqrt(cm2)

    cmm=1./(3.*cm2)

    te9=tea*cm
    tea1=gm1*tea 
    teb1=gm1*teb 
    tec1=gm1*tec 

!    evr11=1.0
!    evr21=0.0
!    evr31=0.0
!    evr41=1.0
!    evr51=1.0
!    evr12=tea-cm
!    evr22=0.0
!    evr32=0.0
!    evr42=tea
!    evr52=tea+cm
!    evr13=teb
!    evr23=1.0
!    evr33=0.0
!    evr43=teb
!    evr53=teb
!    evr14=tec
!    evr24=0.0
!    evr34=1.0
!    evr44=tec
!    evr54=tec
    evr15=hm-te9
!    evr25=teb
!    evr35=tec
!    evr45=qm
    evr55=hm+te9

    evl11=te8+te9
    evl21=-(cm+tea1)
!    evl31=-teb1
!    evl41=-tec1
!    evl51=gm1 
    evl12=-teb*cm2
!    evl22=0.0
!    evl32=cm2 
!    evl42=0.0
!    evl52=0.0
    evl13=-tec*cm2
!    evl23=0.0
!    evl33=0.0
!    evl43=cm2
!    evl53=0.0
    evl14=cm2-te8
!    evl24=tea1 
!    evl34=teb1 
!    evl44=tec1 
!    evl54=-gm1
    evl15=te8-te9
    evl25=cm-tea1
!    evl35=-teb1
!    evl45=-tec1
!    evl55=gm1

    h11 = (evl11*gg1(ip4,1)+evl21*gg1(ip4,2)-teb1*gg1(ip4,3)-tec1*gg1(ip4,4)+gm1*gg1(ip4,5))*0.5
    h21 = evl12*gg1(ip4,1)+cm2*gg1(ip4,3)
    h31 = evl13*gg1(ip4,1)+cm2*gg1(ip4,4)
    h41 = evl14*gg1(ip4,1)+tea1*gg1(ip4,2)+teb1*gg1(ip4,3)+tec1*gg1(ip4,4)-gm1*gg1(ip4,5)
    h51 = (evl15*gg1(ip4,1)+evl25*gg1(ip4,2)-teb1*gg1(ip4,3)-tec1*gg1(ip4,4)+gm1*gg1(ip4,5))*0.5
    h61 = gg1(ip4,6)
    h12 = (evl11*gg1(ip3,1)+evl21*gg1(ip3,2)-teb1*gg1(ip3,3)-tec1*gg1(ip3,4)+gm1*gg1(ip3,5))*0.5
    h22 = evl12*gg1(ip3,1)+cm2*gg1(ip3,3) 
    h32 = evl13*gg1(ip3,1)+cm2*gg1(ip3,4)
    h42 = evl14*gg1(ip3,1)+tea1*gg1(ip3,2)+teb1*gg1(ip3,3)+tec1*gg1(ip3,4)-gm1*gg1(ip3,5)
    h52 = (evl15*gg1(ip3,1)+evl25*gg1(ip3,2)-teb1*gg1(ip3,3)-tec1*gg1(ip3,4)+gm1*gg1(ip3,5))*0.5
    h62 = gg1(ip3,6)
    h13 = (evl11*gg1(ip2,1)+evl21*gg1(ip2,2)-teb1*gg1(ip2,3)-tec1*gg1(ip2,4)+gm1*gg1(ip2,5))*0.5
    h23 = evl12*gg1(ip2,1)+cm2*gg1(ip2,3)
    h33 = evl13*gg1(ip2,1)+cm2*gg1(ip2,4)
    h43 = evl14*gg1(ip2,1)+tea1*gg1(ip2,2)+teb1*gg1(ip2,3)+tec1*gg1(ip2,4)-gm1*gg1(ip2,5)
    h53 = (evl15*gg1(ip2,1)+evl25*gg1(ip2,2)-teb1*gg1(ip2,3)-tec1*gg1(ip2,4)+gm1*gg1(ip2,5))*0.5
    h63 = gg1(ip2,6)
    h14 = (evl11*gg1(ip1,1)+evl21*gg1(ip1,2)-teb1*gg1(ip1,3)-tec1*gg1(ip1,4)+gm1*gg1(ip1,5))*0.5
    h24 = evl12*gg1(ip1,1)+cm2*gg1(ip1,3)
    h34 = evl13*gg1(ip1,1)+cm2*gg1(ip1,4)
    h44 = evl14*gg1(ip1,1)+tea1*gg1(ip1,2)+teb1*gg1(ip1,3)+tec1*gg1(ip1,4)-gm1*gg1(ip1,5)
    h54 = (evl15*gg1(ip1,1)+evl25*gg1(ip1,2)-teb1*gg1(ip1,3)-tec1*gg1(ip1,4)+gm1*gg1(ip1,5))*0.5
    h64 = gg1(ip1,6)
    h15 = (evl11*gg2(i0,1)+evl21*gg2(i0,2)-teb1*gg2(i0,3)-tec1*gg2(i0,4)+gm1*gg2(i0,5))*0.5
    h25 = evl12*gg2(i0,1)+cm2*gg2(i0,3)
    h35 = evl13*gg2(i0,1)+cm2*gg2(i0,4)
    h45 = evl14*gg2(i0,1)+tea1*gg2(i0,2)+teb1*gg2(i0,3)+tec1*gg2(i0,4)-gm1*gg2(i0,5)
    h55 = (evl15*gg2(i0,1)+evl25*gg2(i0,2)-teb1*gg2(i0,3)-tec1*gg2(i0,4)+gm1*gg2(i0,5))*0.5
    h65 = gg2(i0,6)
    h16 = (evl11*gg2(ip1,1)+evl21*gg2(ip1,2)-teb1*gg2(ip1,3)-tec1*gg2(ip1,4)+gm1*gg2(ip1,5))*0.5
    h26 = evl12*gg2(ip1,1)+cm2*gg2(ip1,3)
    h36 = evl13*gg2(ip1,1)+cm2*gg2(ip1,4)
    h46 = evl14*gg2(ip1,1)+tea1*gg2(ip1,2)+teb1*gg2(ip1,3)+tec1*gg2(ip1,4)-gm1*gg2(ip1,5)
    h56 = (evl15*gg2(ip1,1)+evl25*gg2(ip1,2)-teb1*gg2(ip1,3)-tec1*gg2(ip1,4)+gm1*gg2(ip1,5))*0.5
    h66 = gg2(ip1,6)
    h17 = (evl11*gg2(ip2,1)+evl21*gg2(ip2,2)-teb1*gg2(ip2,3)-tec1*gg2(ip2,4)+gm1*gg2(ip2,5))*0.5
    h27 = evl12*gg2(ip2,1)+cm2*gg2(ip2,3)
    h37 = evl13*gg2(ip2,1)+cm2*gg2(ip2,4)
    h47 = evl14*gg2(ip2,1)+tea1*gg2(ip2,2)+teb1*gg2(ip2,3)+tec1*gg2(ip2,4)-gm1*gg2(ip2,5)
    h57 = (evl15*gg2(ip2,1)+evl25*gg2(ip2,2)-teb1*gg2(ip2,3)-tec1*gg2(ip2,4)+gm1*gg2(ip2,5))*0.5
    h67 = gg2(ip2,6)
    h18 = (evl11*gg2(ip3,1)+evl21*gg2(ip3,2)-teb1*gg2(ip3,3)-tec1*gg2(ip3,4)+gm1*gg2(ip3,5))*0.5
    h28 = evl12*gg2(ip3,1)+cm2*gg2(ip3,3)
    h38 = evl13*gg2(ip3,1)+cm2*gg2(ip3,4)
    h48 = evl14*gg2(ip3,1)+tea1*gg2(ip3,2)+teb1*gg2(ip3,3)+tec1*gg2(ip3,4)-gm1*gg2(ip3,5)
    h58 = (evl15*gg2(ip3,1)+evl25*gg2(ip3,2)-teb1*gg2(ip3,3)-tec1*gg2(ip3,4)+gm1*gg2(ip3,5))*0.5
    h68 = gg2(ip3,6)

    ts11 = h11 - h12
    ts12 = h12 - h13
    ts13 = h13 - h14

    ts14 = h15 - h16
    ts15 = h16 - h17
    ts16 = h17 - h18

    tt11 = 13.*ts11*ts11 + 3.*(   h11 - 3.*h12)**2+epweno
    tt12 = 13.*ts12*ts12 + 3.*(   h12 +    h13)**2+epweno
    tt13 = 13.*ts13*ts13 + 3.*(3.*h13 -    h14)**2+epweno
    
    tt14 = 13.*ts14*ts14 + 3.*(   h15 - 3.*h16)**2+epweno
    tt15 = 13.*ts15*ts15 + 3.*(   h16 +    h17)**2+epweno
    tt16 = 13.*ts16*ts16 + 3.*(3.*h17 -    h18)**2+epweno

    s11=(ts12-ts11)/(1. + 6.* (tt11/tt12)**2 + 3. * ((tt11/tt13)**2 ))  !!s11
    s12=(ts13-ts12)/(2. + (tt13/tt11)**2 / 1.5  + 4. * (tt13/tt12)**2 )-0.25*(ts13-ts12) !!s12

    s13=(ts15-ts14)/(1. + 6.* (tt14/tt15)**2 + 3. * ((tt14/tt16)**2 ))  !!s11
    s14=(ts16-ts15)/(2. + (tt16/tt14)**2 / 1.5  + 4. * (tt16/tt15)**2 )-0.25*(ts16-ts15) !!s12

    s15=(s11+s12+s13+s14)*cmm

    ts21 = h21 - h22
    ts22 = h22 - h23
    ts23 = h23 - h24

    ts24 = h25 - h26
    ts25 = h26 - h27
    ts26 = h27 - h28

    tt21 = 13.*ts21*ts21 + 3.*(   h21 - 3.*h22)**2+epweno
    tt22 = 13.*ts22*ts22 + 3.*(   h22 +    h23)**2+epweno
    tt23 = 13.*ts23*ts23 + 3.*(3.*h23 -    h24)**2+epweno

    tt24 = 13.*ts24*ts24 + 3.*(   h25 - 3.*h26)**2+epweno
    tt25 = 13.*ts25*ts25 + 3.*(   h26 +    h27)**2+epweno
    tt26 = 13.*ts26*ts26 + 3.*(3.*h27 -    h28)**2+epweno

    s21=(ts22-ts21)/(1. + 6.* (tt21/tt22)**2 + 3. * ((tt21/tt23)**2 ))  !!s21
    s22=(ts23-ts22)/(2. + (tt23/tt21)**2 / 1.5  + 4. * (tt23/tt22)**2 )-0.25*(ts23-ts22) !!s22

    s23=(ts25-ts24)/(1. + 6.* (tt24/tt25)**2 + 3. * ((tt24/tt26)**2 ))  !!s21
    s24=(ts26-ts25)/(2. + (tt26/tt24)**2 / 1.5  + 4. * (tt26/tt25)**2 )-0.25*(ts26-ts25) !!s22

    s25=(s21+s22+s23+s24)*cmm


    ts31 = h31 - h32
    ts32 = h32 - h33
    ts33 = h33 - h34

    ts34 = h35 - h36
    ts35 = h36 - h37
    ts36 = h37 - h38

    tt31 = 13.*ts31*ts31 + 3.*(   h31 - 3.*h32)**2+epweno
    tt32 = 13.*ts32*ts32 + 3.*(   h32 +    h33)**2+epweno
    tt33 = 13.*ts33*ts33 + 3.*(3.*h33 -    h34)**2+epweno
    
    tt34 = 13.*ts34*ts34 + 3.*(   h35 - 3.*h36)**2+epweno
    tt35 = 13.*ts35*ts35 + 3.*(   h36 +    h37)**2+epweno
    tt36 = 13.*ts36*ts36 + 3.*(3.*h37 -    h38)**2+epweno

    s31=(ts32-ts31)/(1. + 6.* (tt31/tt32)**2 + 3. * ((tt31/tt33)**2 ))  !!s31
    s32=(ts33-ts32)/(2. + (tt33/tt31)**2 / 1.5  + 4. * (tt33/tt32)**2 )-0.25*(ts33-ts32) !!s32

    s33=(ts35-ts34)/(1. + 6.* (tt34/tt35)**2 + 3. * ((tt34/tt36)**2 ))  !!s31
    s34=(ts36-ts35)/(2. + (tt36/tt34)**2 / 1.5  + 4. * (tt36/tt35)**2 )-0.25*(ts36-ts35) !!s32

    s35=(s31+s32+s33+s34)*cmm

    ts41 = h41 - h42
    ts42 = h42 - h43
    ts43 = h43 - h44

    ts44 = h45 - h46
    ts45 = h46 - h47
    ts46 = h47 - h48

    tt41 = 13.*ts41*ts41 + 3.*(   h41 - 3.*h42)**2+epweno
    tt42 = 13.*ts42*ts42 + 3.*(   h42 +    h43)**2+epweno
    tt43 = 13.*ts43*ts43 + 3.*(3.*h43 -    h44)**2+epweno
    
    tt44 = 13.*ts44*ts44 + 3.*(   h45 - 3.*h46)**2+epweno
    tt45 = 13.*ts45*ts45 + 3.*(   h46 +    h47)**2+epweno
    tt46 = 13.*ts46*ts46 + 3.*(3.*h47 -    h48)**2+epweno

    s41=(ts42-ts41)/(1. + 6.* (tt41/tt42)**2 + 3. * ((tt41/tt43)**2 ))  !!s41
    s42=(ts43-ts42)/(2. + (tt43/tt41)**2 / 1.5  + 4. * (tt43/tt42)**2 )-0.25*(ts43-ts42) !!s42

    s43=(ts45-ts44)/(1. + 6.* (tt44/tt45)**2 + 3. * ((tt44/tt46)**2 ))  !!s41
    s44=(ts46-ts45)/(2. + (tt46/tt44)**2 / 1.5  + 4. * (tt46/tt45)**2 )-0.25*(ts46-ts45) !!s42

    s45=(s41+s42+s43+s44)*cmm

    ts51 = h51 - h52
    ts52 = h52 - h53
    ts53 = h53 - h54

    ts54 = h55 - h56
    ts55 = h56 - h57
    ts56 = h57 - h58

    tt51 = 13.*ts51*ts51 + 3.*(   h51 - 3.*h52)**2+epweno
    tt52 = 13.*ts52*ts52 + 3.*(   h52 +    h53)**2+epweno
    tt53 = 13.*ts53*ts53 + 3.*(3.*h53 -    h54)**2+epweno
    
    tt54 = 13.*ts54*ts54 + 3.*(   h55 - 3.*h56)**2+epweno
    tt55 = 13.*ts55*ts55 + 3.*(   h56 +    h57)**2+epweno
    tt56 = 13.*ts56*ts56 + 3.*(3.*h57 -    h58)**2+epweno

    s51=(ts52-ts51)/(1. + 6.* (tt51/tt52)**2 + 3. * ((tt51/tt53)**2 ))  !!s51
    s52=(ts53-ts52)/(2. + (tt53/tt51)**2 / 1.5  + 4. * (tt53/tt52)**2 )-0.25*(ts53-ts52) !!s52

    s53=(ts55-ts54)/(1. + 6.* (tt54/tt55)**2 + 3. * ((tt54/tt56)**2 ))  !!s51
    s54=(ts56-ts55)/(2. + (tt56/tt54)**2 / 1.5  + 4. * (tt56/tt55)**2 )-0.25*(ts56-ts55) !!s52

    s55=(s51+s52+s53+s54)*cmm

    ts61 = h61 - h62
    ts62 = h62 - h63
    ts63 = h63 - h64

    ts64 = h65 - h66
    ts65 = h66 - h67
    ts66 = h67 - h68

    tt61 = 13.*ts61*ts61 + 3.*(   h61 - 3.*h62)**2+epweno
    tt62 = 13.*ts62*ts62 + 3.*(   h62 +    h63)**2+epweno
    tt63 = 13.*ts63*ts63 + 3.*(3.*h63 -    h64)**2+epweno
    
    tt64 = 13.*ts64*ts64 + 3.*(   h65 - 3.*h66)**2+epweno
    tt65 = 13.*ts65*ts65 + 3.*(   h66 +    h67)**2+epweno
    tt66 = 13.*ts66*ts66 + 3.*(3.*h67 -    h68)**2+epweno

    s61=(ts62-ts61)/(3. + 18.* (tt61/tt62)**2 + 9. * ((tt61/tt63)**2 ))  !!s61
    s62=(ts63-ts62)/(6. + (tt63/tt61)**2 * 2.  + 12. * (tt63/tt62)**2 )-0.25*(ts63-ts62)/3. !!s62

    s63=(ts65-ts64)/(3. + 18.* (tt64/tt65)**2 + 9. * ((tt64/tt66)**2 ))  !!s61
    s64=(ts66-ts65)/(6. + (tt66/tt64)**2 * 2.  + 12. * (tt66/tt65)**2 )-0.25*(ts66-ts65)/3. !!s62

    s65=s61+s62+s63+s64

    tfh1=s15+s45+s55

    fh(ip2,1) = (-fd(ip3,1) + 7.*(fd(ip2,1)+fd(ip1,1) )-fd(i0,1) )/12.+tfh1
    fh(ip2,2) = (-fd(ip3,2) + 7.*(fd(ip2,2)+fd(ip1,2) )-fd(i0,2) )/12.+tfh1*tea-(s15-s55)*cm
    fh(ip2,3) = (-fd(ip3,3) + 7.*(fd(ip2,3)+fd(ip1,3) )-fd(i0,3) )/12.+tfh1*teb+s25
    fh(ip2,4) = (-fd(ip3,4) + 7.*(fd(ip2,4)+fd(ip1,4) )-fd(i0,4) )/12.+tfh1*tec+s35
    fh(ip2,5) = (-fd(ip3,5) + 7.*(fd(ip2,5)+fd(ip1,5) )-fd(i0,5) )/12.+evr15*s15+teb*s25+tec*s35+qm*s45+evr55*s55
    fh(ip2,6) = (-fd(ip3,6) + 7.*(fd(ip2,6)+fd(ip1,6) )-fd(i0,6) )/12.+s65

enddo
!do i=0,Lq
!     write(*,"(I3,I4,F14.9,F14.9,F14.9,F14.9,F14.9,F14.9)"), myid, i , fh(i,1), fh(i,2), fh(i,3), fh(i,4), fh(i,5), fh(i,6)
!     write(*,"(I3,I4,F14.9,F14.9,F14.9,F14.9,F14.9,F14.9)"), myid, i , gg1(i,1), gg1(i,2), gg1(i,3), gg1(i,4), gg1(i,5), gg1(i,5)
!enddo
!call mpi_finalize(ierr)
!stop

return

end subroutine weno_new

end module weno_kernel

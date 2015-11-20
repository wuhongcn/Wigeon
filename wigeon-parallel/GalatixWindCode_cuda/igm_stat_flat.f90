module hydro_stat
  use datablock
contains
  
! -----------------------------------------------------------------------------
  subroutine cflc(aam)
    double precision value1(6),value2(3),value1b(6),value2b(3),den_min_b,vL
    integer*4 nczb1,nczb2

! computed the largest eigenvalue for cfl
      
    den_av = 0.
    den_min = 1.0
    den_max = 0.
    den_sigma = 0.
    
    vx_max = 0.
    vy_max = 0.
    vz_max = 0.
    cs_max = 0.
    vel_sigma = 0.
    
    aam = 0.
    ncell_zero = 0
    vL = (1.d0/Lx)*(1.d0/Ly)*(1.d0/Lz)

!nnull_zero = 0
!$omp parallel do default(shared) private(i,j,k,den,eng,cden) &
!$ private(sij,tij,rij,q2,eij,uij,pij,cij,c2,w2) &
!$ reduction(+:ncell_zero) reduction(+:den_av) reduction(min:den_min) &
!$ reduction(+:den_sigma) reduction(+:vel_sigma) &
!$ reduction(max:den_max) reduction(max:vx_max) reduction(max:vy_max) &
!$ reduction(max:vz_max) reduction(max:cs_max) reduction(max:aam)
    do k=1,nLz
       do j=1,nLy
          do i=1,nLx
             
             den = uc(1,i,j,k)
             eng = uc(5,i,j,k)
             cden = 1./den
             
             if (den<1.0e-10) then
                ncell_zero=ncell_zero+1
                write(*,'(3i5,1x,''den='',f13.10,1x,''eng='',f13.10)')i,j,k,den,eng
             endif
             
             den_sigma = den_sigma+(den - 30.95)**2
             
             rij = uc(2,i,j,k) * cden
             sij = uc(3,i,j,k) * cden
             tij = uc(4,i,j,k) * cden 
             
             q2 = rij**2 + sij**2 + tij**2
             eij = 0.5 * q2 * den
             pij = gm1 * ( eng - eij )
             c2 = gamma * pij * cden
             cij = sqrt(abs(c2))
             
             w2 = abs(rij) + cij + abs(sij) + cij + abs(tij) + cij
             aam = max(aam, w2)
             vel_sigma = vel_sigma + q2
             
             den_av = den_av + den
             den_min = min(den_min,den)
             den_max = max(den_max,den)
             vx_max = max(vx_max,abs(rij))
             vy_max = max(vy_max,abs(sij))
             vz_max = max(vz_max,abs(tij))
             cs_max = max(cs_max,cij)
          enddo
       enddo
    enddo
!$omp end parallel do
    
    value1(1)=den_max
    value1(2)=vx_max
    value1(3)=vy_max
    value1(4)=vz_max
    value1(5)=cs_max
    value1(6)=aam
    call mpi_allreduce(value1,value1b,6,mpi_double_precision,mpi_max,mpi_comm_world,ierr)
    call mpi_allreduce(den_min,den_min_b,1,mpi_double_precision,mpi_min,mpi_comm_world,ierr)
    den_max=value1b(1)
    vx_max=value1b(2)
    vy_max=value1b(3)
    vz_max=value1b(4)
    cs_max=value1b(5)
    aam=value1b(6)
    den_min=den_min_b

! Convert integer*8 to integer*4, to avoid 64-bit trouble!
    nczb1=ncell_zero
    call mpi_allreduce(nczb1,nczb2,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
    ncell_zero=nczb2

    value2(1)=den_av*vL
    value2(2)=den_sigma*vL
    value2(3)=vel_sigma*vL

!    print*,'id:value2',myid,value2

    call mpi_allreduce(value2,value2b,3,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

!    if(myid.eq.0) print*,'value2b',value2b

    den_av=value2b(1)
    den_sigma=sqrt(value2b(2))
    vel_sigma=sqrt(value2b(3))

    return
  end subroutine cflc
  
! -----------------------------------------------------------------------------
end module hydro_stat


! This is the velocity update. It uses the Boris algorithm as described in f.e.x
! Birdsall and Langdon.


subroutine AdvanceVelocity(real_particles,int_particles,E,B0,Lx_min,Ly_min,Lz_min,dxyz,dt,v0,ve,eta,tB,sB,&
             num_local,max_per_proc,Nspecies,Nx_local,Ny_local,Nz_local,fluid,gyration)

  implicit none


! Parameters  
  logical, intent(in) :: fluid(8), gyration(8)

  integer, intent(in) :: num_local, max_per_proc, Nspecies, Nx_local, Ny_local, Nz_local
  real*8, intent(in) :: Lx_min, Ly_min, Lz_min, dxyz(3), dt, v0(4,8), ve(4), eta(8)

  real*8, intent(in) :: tB(4,8), sB(4,8)

  real*8, intent(inout) :: real_particles(8,max_per_proc) 
  integer, intent(in) :: int_particles(max_per_proc)

  real*8, intent(in) :: B0(4)
  real*8, intent(in) :: E(3,Nx_local+2,Ny_local+2,Nz_local+2)


! Local variables
  real*8 :: Ep(3), Bp(3), vplus(3), vmin(3),&
       vmid(3), etadt2(8)
  integer :: ii, jj, kk, pp, index_local(3)
  real*8 :: rest(3), drest(3), E_eff(3,2,2,2),&
       F1, F2, F3, F4, F5, F6, F7, F8, dV

  dV = dxyz(1)*dxyz(2)*dxyz(3)
  etadt2(:) = 0.5d0*dt*eta(:)


! Loop over all particles
  do pp = num_local, 1, -1

! Start interpolation of fields to particle positions.
! Loop over x, y, and z directions to find the index corresponding to the lowest x, y, z, which the particle should have a force from.

    rest(1) = dmod( real_particles(1,pp)-Lx_min+0.5d0*dxyz(1), dxyz(1) )
    rest(2) = dmod( real_particles(2,pp)-Ly_min+0.5d0*dxyz(2), dxyz(2) )
    rest(3) = dmod( real_particles(3,pp)-Lz_min+0.5d0*dxyz(3), dxyz(3) )
    drest   = dxyz-rest
    
    index_local(1) = 1+idint( (real_particles(1,pp)-Lx_min)/dxyz(1)+0.5d0 )
    index_local(2) = 1+idint( (real_particles(2,pp)-Ly_min)/dxyz(2)+0.5d0 )
    index_local(3) = 1+idint( (real_particles(3,pp)-Lz_min)/dxyz(3)+0.5d0 )


! Find fields Fx, Fy, Fz  at eight closest gridpoints
    do kk = 1, 2
      do jj = 1, 2
        do ii = 1, 2
          E_eff(1,ii,jj,kk) = E(1,index_local(1)+ii-1,index_local(2)+jj-1,index_local(3)+kk-1)
          E_eff(2,ii,jj,kk) = E(2,index_local(1)+ii-1,index_local(2)+jj-1,index_local(3)+kk-1)
          E_eff(3,ii,jj,kk) = E(3,index_local(1)+ii-1,index_local(2)+jj-1,index_local(3)+kk-1)
        end do
      end do
    end do

    F1 = drest(1)*drest(2)*drest(3)
    F2 =  rest(1)*drest(2)*drest(3)
    F3 = drest(1)* rest(2)*drest(3)
    F4 = drest(1)*drest(2)* rest(3)
    F5 =  rest(1)* rest(2)*drest(3)
    F6 = drest(1)* rest(2)* rest(3)
    F7 =  rest(1)*drest(2)* rest(3)
    F8 =  rest(1)* rest(2)* rest(3)

    Ep(:) = E_eff(:,1,1,1)*F1+E_eff(:,2,1,1)*F2+ &
            E_eff(:,1,2,1)*F3+E_eff(:,1,1,2)*F4+ &
            E_eff(:,2,2,1)*F5+E_eff(:,1,2,2)*F6+ &
            E_eff(:,2,1,2)*F7+E_eff(:,2,2,2)*F8
    Ep(:) = Ep(:)/dV


! Now we have the fields and we can start the boris rotation.
! The scheme is identical to the one presented in Birdsall & Langdon.
! I even use the same variable names for ease of reading. 

    if ( fluid(int_particles(pp)) ) then
      if ( gyration(int_particles(pp)) ) then
        Bp(1:3) = B0(1:3)
  
      ! Parallel to B
        vmin(1:3) = ( real_particles(4,pp)*Bp(1) + real_particles(5,pp)*Bp(2) + real_particles(6,pp)*Bp(3) )*Bp(1:3)
      ! Adding E-cross-B drift
        vmid(1)   = vmin(1) + ( Ep(2)*Bp(3)-Ep(3)*Bp(2) )
        vmid(2)   = vmin(2) + ( Ep(3)*Bp(1)-Ep(1)*Bp(3) )
        vmid(3)   = vmin(3) + ( Ep(1)*Bp(2)-Ep(2)*Bp(1) )
        real_particles(4:6,pp) = vmid(:)*sum( Bp(1:3)**2 )**(-1)
      end if
    else
  
    !!!! Calculating 'tB' was done before

    ! Calculate vmin by moving to the frame where the ambient electric field is zero, and then add the first half acceleration
      vmin(:) = real_particles(4:6,pp) - ve(1:3) + etadt2( int_particles(pp) )*Ep(:)

    ! Cross product
      vmid(1) = vmin(1) + vmin(2)*tB(3,int_particles(pp))-vmin(3)*tB(2,int_particles(pp))
      vmid(2) = vmin(2) + vmin(3)*tB(1,int_particles(pp))-vmin(1)*tB(3,int_particles(pp))
      vmid(3) = vmin(3) + vmin(1)*tB(2,int_particles(pp))-vmin(2)*tB(1,int_particles(pp))

    !!!! Calculating 'sB' was done before

    ! Another cross product
      vplus(1) = vmin(1) + vmid(2)*sB(3,int_particles(pp))-vmid(3)*sB(2,int_particles(pp))
      vplus(2) = vmin(2) + vmid(3)*sB(1,int_particles(pp))-vmid(1)*sB(3,int_particles(pp))
      vplus(3) = vmin(3) + vmid(1)*sB(2,int_particles(pp))-vmid(2)*sB(1,int_particles(pp))

    ! Set new velocities (second half acceleration), and add the ambient speed.
      real_particles(4:6,pp) = vplus(:) + ve(1:3) + etadt2( int_particles(pp) )*Ep(:)

    end if    

  end do

  return

end subroutine AdvanceVelocity


subroutine AdvanceVelocityExp(real_particles,int_particles,E,B0,Lx_min,Ly_min,Lz_min,dxyz,deltat,dt,v0,ve,eta,tB,&
             num_local,max_per_proc,Nspecies,Nx_local,Ny_local,Nz_local,fluid,gyration)

! Parameters  
  logical, intent(in) :: fluid(8), gyration(8)

  integer, intent(in) :: num_local, max_per_proc, Nspecies, Nx_local, Ny_local, Nz_local
  real*8, intent(in) :: Lx_min, Ly_min, Lz_min, dxyz(3), deltat, dt, v0(4,8), ve(4), eta(8)

  real*8, intent(in) :: tB(4,8)

  real*8, intent(inout) :: real_particles(8,max_per_proc) 
  integer, intent(in) :: int_particles(max_per_proc)

  real*8, intent(in) :: B0(4)
  real*8, intent(in) :: E(3,Nx_local+2,Ny_local+2,Nz_local+2)


! Local variables
  real*8 :: Ep(3), Bp(3), vplus(3), vmin(3),&
       vmid(3), etadt2(8), tB_temp(4,8)
  integer :: ii, jj, kk, pp, index_local(3)
  real*8 :: rest(3), drest(3), E_eff(3,2,2,2),&
       F1, F2, F3, F4, F5, F6, F7, F8, dV

  dV = dxyz(1)*dxyz(2)*dxyz(3)
  etadt2(:) = 0.5d0*dt*eta(:) * 2.0d0*deltat/dt
  tB_temp = tB * 2.0d0*deltat/dt



! Loop over all particles
  do pp = num_local, 1, -1

! Start interpolation of fields to particle positions.
! Loop over x, y, and z directions to find the index corresponding to the lowest x, y, z, which the particle should have a force from.

    rest(1) = dmod( real_particles(1,pp)-Lx_min+0.5d0*dxyz(1), dxyz(1) )
    rest(2) = dmod( real_particles(2,pp)-Ly_min+0.5d0*dxyz(2), dxyz(2) )
    rest(3) = dmod( real_particles(3,pp)-Lz_min+0.5d0*dxyz(3), dxyz(3) )
    drest   = dxyz-rest
    
    index_local(1) = 1+idint( (real_particles(1,pp)-Lx_min)/dxyz(1)+0.5d0 )
    index_local(2) = 1+idint( (real_particles(2,pp)-Ly_min)/dxyz(2)+0.5d0 )
    index_local(3) = 1+idint( (real_particles(3,pp)-Lz_min)/dxyz(3)+0.5d0 )


! Find fields Fx, Fy, Fz  at eight closest gridpoints
    do kk = 1, 2
      do jj = 1, 2
        do ii = 1, 2
          E_eff(1,ii,jj,kk) = E(1,index_local(1)+ii-1,index_local(2)+jj-1,index_local(3)+kk-1)
          E_eff(2,ii,jj,kk) = E(2,index_local(1)+ii-1,index_local(2)+jj-1,index_local(3)+kk-1)
          E_eff(3,ii,jj,kk) = E(3,index_local(1)+ii-1,index_local(2)+jj-1,index_local(3)+kk-1)
        end do
      end do
    end do

    F1 = drest(1)*drest(2)*drest(3)
    F2 =  rest(1)*drest(2)*drest(3)
    F3 = drest(1)* rest(2)*drest(3)
    F4 = drest(1)*drest(2)* rest(3)
    F5 =  rest(1)* rest(2)*drest(3)
    F6 = drest(1)* rest(2)* rest(3)
    F7 =  rest(1)*drest(2)* rest(3)
    F8 =  rest(1)* rest(2)* rest(3)

    Ep(:) = E_eff(:,1,1,1)*F1+E_eff(:,2,1,1)*F2+ &
            E_eff(:,1,2,1)*F3+E_eff(:,1,1,2)*F4+ &
            E_eff(:,2,2,1)*F5+E_eff(:,1,2,2)*F6+ &
            E_eff(:,2,1,2)*F7+E_eff(:,2,2,2)*F8
    Ep(:) = Ep(:)/dV

    if ( fluid(int_particles(pp)) ) then
      if ( gyration(int_particles(pp)) ) then
        Bp(1:3) = B0(1:3)
  
      ! Parallel to B
        vmin(1:3) = ( real_particles(4,pp)*Bp(1) + real_particles(5,pp)*Bp(2) + real_particles(6,pp)*Bp(3) )*Bp(1:3)
      ! Adding E-cross-B drift
        vmid(1)   = vmin(1) + ( Ep(2)*Bp(3)-Ep(3)*Bp(2) )
        vmid(2)   = vmin(2) + ( Ep(3)*Bp(1)-Ep(1)*Bp(3) )
        vmid(3)   = vmin(3) + ( Ep(1)*Bp(2)-Ep(2)*Bp(1) )
        real_particles(4:6,pp) = vmid(:)*sum( Bp(1:3)**2 )**(-1)
      end if
    else
  
    ! Calculate vmin by moving to the frame where the ambient electric field is zero
      vmin(:) = real_particles(4:6,pp) - ve(1:3)

    ! Add magnetic field (rotation)
      vmid(1) = vmin(1) + vmin(2)*tB_temp(3,int_particles(pp))-vmin(3)*tB_temp(2,int_particles(pp))
      vmid(2) = vmin(2) + vmin(3)*tB_temp(1,int_particles(pp))-vmin(1)*tB_temp(3,int_particles(pp))
      vmid(3) = vmin(3) + vmin(1)*tB_temp(2,int_particles(pp))-vmin(2)*tB_temp(1,int_particles(pp))

    ! normalize
      if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
        vmid = vmid * dsqrt( sum( vmin(:)**2 ) /  sum( vmid(:)**2 ) )
      end if
    ! Add ambient velocity and electric field
      real_particles(4:6,pp) = vmid(:) + ve(1:3) + etadt2( int_particles(pp) )*Ep(:)
    
    end if

  end do

  return

end subroutine AdvanceVelocityExp


subroutine AdvanceVelocityImp(real_particles,int_particles,E,B0,Lx_min,Ly_min,Lz_min,dxyz,deltat,dt,v0,ve,eta,tB,&
             num_local,max_per_proc,Nspecies,Nx_local,Ny_local,Nz_local,fluid,gyration)

! Parameters  
  logical, intent(in) :: fluid(8), gyration(8)

  integer, intent(in) :: num_local, max_per_proc, Nspecies, Nx_local, Ny_local, Nz_local
  real*8, intent(in) :: Lx_min, Ly_min, Lz_min, dxyz(3), deltat, dt, v0(4,8), ve(4), eta(8)

  real*8, intent(in) :: tB(4,8)

  real*8, intent(inout) :: real_particles(8,max_per_proc) 
  integer, intent(in) :: int_particles(max_per_proc)

  real*8, intent(in) :: B0(4)
  real*8, intent(in) :: E(3,Nx_local+2,Ny_local+2,Nz_local+2)


! Local variables
  real*8 :: Ep(3), Bp(3), vplus(3), vmin(3),&
       vmid(3), etadt2(8), tB_temp(4,8)
  integer :: ii, jj, kk, pp, index_local(3)
  real*8 :: rest(3), drest(3), E_eff(3,2,2,2),&
       F1, F2, F3, F4, F5, F6, F7, F8, dV

  dV = dxyz(1)*dxyz(2)*dxyz(3)
  etadt2(:) = 0.5d0*dt*eta(:) * 2.0d0*deltat/dt
  tB_temp = tB * 2.0d0*deltat/dt


! Loop over all particles
  do pp = num_local, 1, -1

! Start interpolation of fields to particle positions.
! Loop over x, y, and z directions to find the index corresponding to the lowest x, y, z, which the particle should have a force from.

    rest(1) = dmod( real_particles(1,pp)-Lx_min+0.5d0*dxyz(1), dxyz(1) )
    rest(2) = dmod( real_particles(2,pp)-Ly_min+0.5d0*dxyz(2), dxyz(2) )
    rest(3) = dmod( real_particles(3,pp)-Lz_min+0.5d0*dxyz(3), dxyz(3) )
    drest   = dxyz-rest
    
    index_local(1) = 1+idint( (real_particles(1,pp)-Lx_min)/dxyz(1)+0.5d0 )
    index_local(2) = 1+idint( (real_particles(2,pp)-Ly_min)/dxyz(2)+0.5d0 )
    index_local(3) = 1+idint( (real_particles(3,pp)-Lz_min)/dxyz(3)+0.5d0 )


! Find fields Fx, Fy, Fz  at eight closest gridpoints
    do kk = 1, 2
      do jj = 1, 2
        do ii = 1, 2
          E_eff(1,ii,jj,kk) = E(1,index_local(1)+ii-1,index_local(2)+jj-1,index_local(3)+kk-1)
          E_eff(2,ii,jj,kk) = E(2,index_local(1)+ii-1,index_local(2)+jj-1,index_local(3)+kk-1)
          E_eff(3,ii,jj,kk) = E(3,index_local(1)+ii-1,index_local(2)+jj-1,index_local(3)+kk-1)
        end do
      end do
    end do

    F1 = drest(1)*drest(2)*drest(3)
    F2 =  rest(1)*drest(2)*drest(3)
    F3 = drest(1)* rest(2)*drest(3)
    F4 = drest(1)*drest(2)* rest(3)
    F5 =  rest(1)* rest(2)*drest(3)
    F6 = drest(1)* rest(2)* rest(3)
    F7 =  rest(1)*drest(2)* rest(3)
    F8 =  rest(1)* rest(2)* rest(3)

    Ep(:) = E_eff(:,1,1,1)*F1+E_eff(:,2,1,1)*F2+ &
            E_eff(:,1,2,1)*F3+E_eff(:,1,1,2)*F4+ &
            E_eff(:,2,2,1)*F5+E_eff(:,1,2,2)*F6+ &
            E_eff(:,2,1,2)*F7+E_eff(:,2,2,2)*F8
    Ep(:) = Ep(:)/dV

    if ( fluid(int_particles(pp)) ) then
      if ( gyration(int_particles(pp)) ) then
        Bp(1:3) = B0(1:3)
  
      ! Parallel to B
        vmin(1:3) = ( real_particles(4,pp)*Bp(1) + real_particles(5,pp)*Bp(2) + real_particles(6,pp)*Bp(3) )*Bp(1:3)
      ! Adding E-cross-B drift
        vmid(1)   = vmin(1) + ( Ep(2)*Bp(3)-Ep(3)*Bp(2) )
        vmid(2)   = vmin(2) + ( Ep(3)*Bp(1)-Ep(1)*Bp(3) )
        vmid(3)   = vmin(3) + ( Ep(1)*Bp(2)-Ep(2)*Bp(1) )
        real_particles(4:6,pp) = vmid(:)*sum( Bp(1:3)**2 )**(-1)
      end if
    else
  
    ! Calculate vmin by moving to the frame where the ambient electric field is zero, and then add the electric field
      vmin(:) = real_particles(4:6,pp) - ve(1:3) + etadt2( int_particles(pp) )*Ep(:)

    ! Add magnetic field (rotation)
      vmid(1) = vmin(1) + vmin(2)*tB_temp(3,int_particles(pp))-vmin(3)*tB_temp(2,int_particles(pp))
      vmid(2) = vmin(2) + vmin(3)*tB_temp(1,int_particles(pp))-vmin(1)*tB_temp(3,int_particles(pp))
      vmid(3) = vmin(3) + vmin(1)*tB_temp(2,int_particles(pp))-vmin(2)*tB_temp(1,int_particles(pp))

    ! normalize
      if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
        vmid = vmid * dsqrt( sum( vmin(:)**2 ) /  sum( vmid(:)**2 ) )
      end if

    ! Add ambient velocity
      real_particles(4:6,pp) = vmid(:) + ve(1:3)
    
    end if

  end do

  return

end subroutine AdvanceVelocityImp




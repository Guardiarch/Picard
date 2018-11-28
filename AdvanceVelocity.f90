! This is the velocity update. It uses the Boris algorithm as described by e.g.
! Birdsall and Langdon.

subroutine AdvanceVelocity(particles,E,B0,Lx_min,Ly_min,Lz_min,dxyz,dt,ve,&
     num_local,max_per_proc,Nspecies,Nx_local,Ny_local,Nz_local,species)

  use SpecificTypes
  implicit none


! Parameters
  integer, intent(in) :: num_local, max_per_proc, Nspecies, &
       Nx_local, Ny_local, Nz_local
  real*8, intent(in) :: Lx_min, Ly_min, Lz_min, dxyz(3), dt, ve(4)

  real*8, intent(in) :: B0(4)
  real*8, intent(in) :: E(3,Nx_local+2,Ny_local+2,Nz_local+2)
  type(particlearrays) particles
  type(particlespecies) species(Nspecies)


! Local variables
  real*8 :: Ep(3), Bp(3), vplus(3), vmin(3),&
       vmid(3), etadt2(Nspecies)
  integer :: ii, jj, kk, pp, index_local(3)
  real*8 :: rest(3), drest(3), E_eff(3,2,2,2),&
       F1, F2, F3, F4, F5, F6, F7, F8, dV

  dV = dxyz(1)*dxyz(2)*dxyz(3)
  etadt2 = 0.5d0*dt*species%eta

! Loop over all particles
  do pp = num_local, 1, -1
     
! Start interpolation of fields to particle positions.
! Loop over x, y, and z directions to find the index corresponding
! to the lowest x, y, z, which the particle should have a force from.

    rest(1) = dmod( particles%coordinates(1,pp)-Lx_min+0.5d0*dxyz(1), dxyz(1) )
    rest(2) = dmod( particles%coordinates(2,pp)-Ly_min+0.5d0*dxyz(2), dxyz(2) )
    rest(3) = dmod( particles%coordinates(3,pp)-Lz_min+0.5d0*dxyz(3), dxyz(3) )
    drest   = dxyz-rest
    
    index_local(1)=1+idint( (particles%coordinates(1,pp)-Lx_min)/dxyz(1)+0.5d0)
    index_local(2)=1+idint( (particles%coordinates(2,pp)-Ly_min)/dxyz(2)+0.5d0)
    index_local(3)=1+idint( (particles%coordinates(3,pp)-Lz_min)/dxyz(3)+0.5d0)


! Find fields Fx, Fy, Fz  at eight closest gridpoints
    do kk = 1, 2
      do jj = 1, 2
        do ii = 1, 2
           E_eff(1,ii,jj,kk) = &
                E(1,index_local(1)+ii-1,index_local(2)+jj-1,index_local(3)+kk-1)
           E_eff(2,ii,jj,kk) = &
                E(2,index_local(1)+ii-1,index_local(2)+jj-1,index_local(3)+kk-1)
           E_eff(3,ii,jj,kk) = &
                E(3,index_local(1)+ii-1,index_local(2)+jj-1,index_local(3)+kk-1)
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

    if ( species(particles%species(pp))%fluid ) then
      if ( species(particles%species(pp))%gyration ) then
        Bp(1:3) = B0(1:3)
  
      ! Parallel to B
        vmin(1:3) = ( particles%coordinates(4,pp)*Bp(1) + &
             particles%coordinates(5,pp)*Bp(2) + &
             particles%coordinates(6,pp)*Bp(3) )*Bp(1:3)
        ! Adding E-cross-B drift
        vmid(1)   = vmin(1) + ( Ep(2)*Bp(3)-Ep(3)*Bp(2) )
        vmid(2)   = vmin(2) + ( Ep(3)*Bp(1)-Ep(1)*Bp(3) )
        vmid(3)   = vmin(3) + ( Ep(1)*Bp(2)-Ep(2)*Bp(1) )
        particles%coordinates(4:6,pp) = vmid(:)*sum( Bp(1:3)**2 )**(-1)
      end if
    else
  
    !!!! Calculating 'tB' was done before

       ! Calculate vmin by moving to the frame where the ambient
       ! electric field is zero, and then add the first half acceleration
       vmin(:) = particles%coordinates(4:6,pp) - ve(1:3) + &
            etadt2( particles%species(pp) )*Ep(:)

       ! Cross product
       vmid(1) = vmin(1) +vmin(2)*species(particles%species(pp))%tB(3) - &
            vmin(3)*species(particles%species(pp))%tB(2)
       vmid(2) = vmin(2) +vmin(3)*species(particles%species(pp))%tB(1) - &
            vmin(1)*species(particles%species(pp))%tB(3)
       vmid(3) = vmin(3) +vmin(1)*species(particles%species(pp))%tB(2) - &
            vmin(2)*species(particles%species(pp))%tB(1)

!!!! Calculating 'sB' was done before

    ! Another cross product
       vplus(1) = vmin(1) +vmid(2)*species(particles%species(pp))%sB(3) - &
            vmid(3)*species(particles%species(pp))%sB(2)
       vplus(2) = vmin(2) +vmid(3)*species(particles%species(pp))%sB(1) - &
            vmid(1)*species(particles%species(pp))%sB(3)
       vplus(3) = vmin(3) +vmid(1)*species(particles%species(pp))%sB(2) - &
            vmid(2)*species(particles%species(pp))%sB(1)

    ! Set new velocities (second half acceleration), and add the ambient speed.
       particles%coordinates(4:6,pp) = vplus(:) + ve(1:3) + &
            etadt2( particles%species(pp) )*Ep(:)

    end if

 end do

 return

end subroutine AdvanceVelocity

!------------------------------------------------------------------------

subroutine AdvanceVelocityExp(particles,E,B0,Lx_min,Ly_min,Lz_min,dxyz, &
     deltat,dt,ve,num_local,max_per_proc,Nspecies, &
     Nx_local,Ny_local,Nz_local,species)

  use SpecificTypes
  implicit none

  ! Parameters
  integer, intent(in) :: num_local, max_per_proc, Nspecies, &
       Nx_local, Ny_local, Nz_local
  real*8, intent(in) :: Lx_min, Ly_min, Lz_min, dxyz(3), deltat, dt, ve(4)

  real*8, intent(in) :: B0(4)
  real*8, intent(in) :: E(3,Nx_local+2,Ny_local+2,Nz_local+2)
  type(particlearrays) particles
  type(particlespecies) species(Nspecies)


! Local variables
  real*8 :: Ep(3), Bp(3), vplus(3), vmin(3),&
       vmid(3), etadt2(8), tB_temp(4,Nspecies)
  integer :: ii, jj, kk, pp, index_local(3)
  real*8 :: rest(3), drest(3), E_eff(3,2,2,2),&
       F1, F2, F3, F4, F5, F6, F7, F8, dV

  dV = dxyz(1)*dxyz(2)*dxyz(3)
  etadt2 = 0.5d0*dt*species%eta * 2.0d0*deltat/dt

  do ii = 1, Nspecies
     tB_temp(:,ii) = species(ii)%tB(:) * 2.0d0*deltat/dt
  end do


! Loop over all particles
  do pp = num_local, 1, -1

! Start interpolation of fields to particle positions.
! Loop over x, y, and z directions to find the index corresponding to the lowest x, y, z, which the particle should have a force from.

    rest(1) = dmod( particles%coordinates(1,pp)-Lx_min+0.5d0*dxyz(1), dxyz(1) )
    rest(2) = dmod( particles%coordinates(2,pp)-Ly_min+0.5d0*dxyz(2), dxyz(2) )
    rest(3) = dmod( particles%coordinates(3,pp)-Lz_min+0.5d0*dxyz(3), dxyz(3) )
    drest   = dxyz-rest
    
    index_local(1)=1+idint( (particles%coordinates(1,pp)-Lx_min)/dxyz(1)+0.5d0)
    index_local(2)=1+idint( (particles%coordinates(2,pp)-Ly_min)/dxyz(2)+0.5d0)
    index_local(3)=1+idint( (particles%coordinates(3,pp)-Lz_min)/dxyz(3)+0.5d0)


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

    if ( species(particles%species(pp))%fluid ) then
      if ( species(particles%species(pp))%gyration ) then
        Bp(1:3) = B0(1:3)
  
      ! Parallel to B
        vmin(1:3) = ( particles%coordinates(4,pp)*Bp(1) + &
             particles%coordinates(5,pp)*Bp(2) + &
             particles%coordinates(6,pp)*Bp(3) )*Bp(1:3)
      ! Adding E-cross-B drift
        vmid(1)   = vmin(1) + ( Ep(2)*Bp(3)-Ep(3)*Bp(2) )
        vmid(2)   = vmin(2) + ( Ep(3)*Bp(1)-Ep(1)*Bp(3) )
        vmid(3)   = vmin(3) + ( Ep(1)*Bp(2)-Ep(2)*Bp(1) )
        particles%coordinates(4:6,pp) = vmid(:)*sum( Bp(1:3)**2 )**(-1)
      end if
    else
  
    ! Calculate vmin by moving to the frame where the ambient electric field is zero
      vmin(:) = particles%coordinates(4:6,pp) - ve(1:3)

    ! Add magnetic field (rotation)
      vmid(1) = vmin(1) + vmin(2)*tB_temp(3,particles%species(pp))- &
           vmin(3)*tB_temp(2,particles%species(pp))
      vmid(2) = vmin(2) + vmin(3)*tB_temp(1,particles%species(pp))- &
           vmin(1)*tB_temp(3,particles%species(pp))
      vmid(3) = vmin(3) + vmin(1)*tB_temp(2,particles%species(pp))- &
           vmin(2)*tB_temp(1,particles%species(pp))

    ! normalize
      if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
        vmid = vmid * dsqrt( sum( vmin(:)**2 ) /  sum( vmid(:)**2 ) )
      end if
    ! Add ambient velocity and electric field
      particles%coordinates(4:6,pp) = vmid(:) + ve(1:3) + &
           etadt2( particles%species(pp) )*Ep(:)
   end if

end do

return

end subroutine AdvanceVelocityExp

!------------------------------------------------------------------------

subroutine AdvanceVelocityImp(particles,E,B0,Lx_min,Ly_min,Lz_min,dxyz, &
     deltat,dt,ve,num_local,max_per_proc,Nspecies, &
     Nx_local,Ny_local,Nz_local,species)

  use SpecificTypes
  implicit none
  
! Parameters  
  integer, intent(in) :: num_local, max_per_proc, Nspecies, &
       Nx_local, Ny_local, Nz_local
  real*8, intent(in) :: Lx_min, Ly_min, Lz_min, dxyz(3), deltat, dt, ve(4)

  real*8, intent(in) :: B0(4)
  real*8, intent(in) :: E(3,Nx_local+2,Ny_local+2,Nz_local+2)
  type(particlearrays) particles
  type(particlespecies) species(Nspecies)


! Local variables
  real*8 :: Ep(3), Bp(3), vplus(3), vmin(3),&
       vmid(3), etadt2(Nspecies), tB_temp(4,Nspecies)
  integer :: ii, jj, kk, pp, index_local(3)
  real*8 :: rest(3), drest(3), E_eff(3,2,2,2),&
       F1, F2, F3, F4, F5, F6, F7, F8, dV

  dV = dxyz(1)*dxyz(2)*dxyz(3)
  etadt2 = 0.5d0*dt*species%eta * 2.0d0*deltat/dt

  do ii = 1, Nspecies
     tB_temp(:,ii) = species(ii)%tB(:) * 2.0d0*deltat/dt
  end do

! Loop over all particles
  do pp = num_local, 1, -1

! Start interpolation of fields to particle positions.
! Loop over x, y, and z directions to find the index corresponding to the lowest x, y, z, which the particle should have a force from.

    rest(1) = dmod( particles%coordinates(1,pp)-Lx_min+0.5d0*dxyz(1), dxyz(1) )
    rest(2) = dmod( particles%coordinates(2,pp)-Ly_min+0.5d0*dxyz(2), dxyz(2) )
    rest(3) = dmod( particles%coordinates(3,pp)-Lz_min+0.5d0*dxyz(3), dxyz(3) )
    drest   = dxyz-rest
    
    index_local(1)=1+idint( (particles%coordinates(1,pp)-Lx_min)/dxyz(1)+0.5d0)
    index_local(2)=1+idint( (particles%coordinates(2,pp)-Ly_min)/dxyz(2)+0.5d0)
    index_local(3)=1+idint( (particles%coordinates(3,pp)-Lz_min)/dxyz(3)+0.5d0)


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

    if ( species(particles%species(pp))%fluid ) then
      if ( species(particles%species(pp))%gyration ) then
        Bp(1:3) = B0(1:3)
  
      ! Parallel to B
        vmin(1:3) = ( particles%coordinates(4,pp)*Bp(1) + &
             particles%coordinates(5,pp)*Bp(2) + &
             particles%coordinates(6,pp)*Bp(3) )*Bp(1:3)
      ! Adding E-cross-B drift
        vmid(1)   = vmin(1) + ( Ep(2)*Bp(3)-Ep(3)*Bp(2) )
        vmid(2)   = vmin(2) + ( Ep(3)*Bp(1)-Ep(1)*Bp(3) )
        vmid(3)   = vmin(3) + ( Ep(1)*Bp(2)-Ep(2)*Bp(1) )
        particles%coordinates(4:6,pp) = vmid(:)*sum( Bp(1:3)**2 )**(-1)
      end if
    else
  
    ! Calculate vmin by moving to the frame where the ambient electric field is zero, and then add the electric field
       vmin(:) = particles%coordinates(4:6,pp) - ve(1:3) + &
            etadt2( particles%species(pp) )*Ep(:)

       ! Add magnetic field (rotation)       
      vmid(1) = vmin(1) + vmin(2)*tB_temp(3,particles%species(pp))-vmin(3)*tB_temp(2,particles%species(pp))
      vmid(2) = vmin(2) + vmin(3)*tB_temp(1,particles%species(pp))-vmin(1)*tB_temp(3,particles%species(pp))
      vmid(3) = vmin(3) + vmin(1)*tB_temp(2,particles%species(pp))-vmin(2)*tB_temp(1,particles%species(pp))

    ! normalize
      if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
        vmid = vmid * dsqrt( sum( vmin(:)**2 ) /  sum( vmid(:)**2 ) )
      end if

    ! Add ambient velocity
      particles%coordinates(4:6,pp) = vmid(:) + ve(1:3)
    
    end if

  end do

  return

end subroutine AdvanceVelocityImp




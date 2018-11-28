subroutine BCinflow(particles,Lx_min,Ly_min,Lz_min,Lx_max,Ly_max,Lz_max, &
     dxyz,ve,max_per_proc,num_local,Nspecies,Nx_local,Ny_local,Nz_local,&
     xmin,xmax,ymin,ymax,zmin,zmax,xflow,yflow,zflow,species)


  use SpecificTypes
  implicit none

! Parameters
  logical, intent(in) :: xflow, yflow, zflow
  integer, intent(in) :: max_per_proc, Nspecies, Nx_local, Ny_local, Nz_local
  real*8, intent(in) :: Lx_min, Ly_min, Lz_min, Lx_max, Ly_max, Lz_max, &
       dxyz(3), ve(4), xmin, xmax, ymin, ymax, zmin, zmax
  integer, intent(inout) :: num_local
  type(particlearrays) particles
  type(particlespecies) species(Nspecies)

! Local variables
  integer ee, pp, ss, ii, jj, kk, yflow_neg, zflow_neg, yflow_pos, zflow_pos
  real*8 rndp(3), x0, y0, z0, vmin(3), vmid(3)


! Do we have the +y boundary?
if ( (yflow) .and. (Ly_max+0.5d0*dxyz(2) .gt. ymax) ) then
  yflow_pos = 1
else
  yflow_pos = 0
end if

! Do we have the -y boundary?
if ( (yflow) .and. (Ly_min-0.5d0*dxyz(2) .lt. ymin) ) then
  yflow_neg = 1
else
  yflow_neg = 0
end if

! Do we have the +z boundary?
if ( (zflow) .and. (Lz_max+0.5d0*dxyz(3) .gt. zmax) ) then
  zflow_pos = 1
else
  zflow_pos = 0
end if

! Do we have the -z boundary?
if ( (zflow) .and. (Lz_min-0.5d0*dxyz(3) .lt. zmin) ) then
  zflow_neg = 1
else
  zflow_neg = 0
end if


! Do we have the +x boundary (main inflow boundary)?
if ( (xflow) .and. (Lx_max+0.5d0*dxyz(1) .gt. xmax) ) then

! Produce undisturbed particles in guard cells at +x boundary
do kk = 2-zflow_neg, Nz_local+1+zflow_pos
  do jj = 2-yflow_neg, Ny_local+1+yflow_pos
      ! Lower edge of cells
      x0 = Lx_max
      y0 = Ly_min + dxyz(2)*(jj-2)
      z0 = Lz_min + dxyz(3)*(kk-2)
      ! Loop through all ions
      do ss = 1, Nspecies

        ! Loop through all particles to be created
        do pp = 1, species(ss)%ipp
          particles%species(num_local+1)        = ss
          call random_number(rndp)
          particles%coordinates(1,num_local+1)     = x0 + dxyz(1)*rndp(1)
          particles%coordinates(2,num_local+1)     = y0 + dxyz(2)*rndp(2)
          particles%coordinates(3,num_local+1)     = z0 + dxyz(3)*rndp(3)
          call Gaussian( particles%coordinates(4:6,num_local+1), &
               species(ss)%v0(1:3),species(ss)%upstreamkelvin,species(ss)%mass )
          ! Advance velocity dt/2
          vmin(1:3) = particles%coordinates(4:6,num_local+1) - ve(1:3)
          vmid(1) = vmin(1) + &
               vmin(2)*species(particles%species(num_local+1))%tB(3) - &
               vmin(3)*species(particles%species(num_local+1))%tB(2)
          vmid(2) = vmin(2) + &
               vmin(3)*species(particles%species(num_local+1))%tB(1) - &
               vmin(1)*species(particles%species(num_local+1))%tB(3)
          vmid(3) = vmin(3) + &
               vmin(1)*species(particles%species(num_local+1))%tB(2) - &
               vmin(2)*species(particles%species(num_local+1))%tB(1)
          if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
            vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
          end if
          particles%coordinates(4:6,num_local+1) = vmid(1:3) + ve(1:3)

          ! For each ion, create as many electrons as its charged state.
          do ee = 1, species(ss)%epp
            particles%species(num_local+1+ee)         = 1
            particles%coordinates(:,num_local+1+ee)      = particles%coordinates(:,num_local+1)
            call Gaussian( particles%coordinates(4:6,num_local+1+ee),species(ss)%v0(1:3),species(1)%upstreamkelvin,species(1)%mass )
            ! Advance velocity dt/2
            vmin(1:3) = particles%coordinates(4:6,num_local+1+ee) - ve(1:3)
            vmid(1) = vmin(1) +vmin(2)*species(1)%tB(3)-vmin(3)*species(1)%tB(2)
            vmid(2) = vmin(2) +vmin(3)*species(1)%tB(1)-vmin(1)*species(1)%tB(3)
            vmid(3) = vmin(3) +vmin(1)*species(1)%tB(2)-vmin(2)*species(1)%tB(1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            particles%coordinates(4:6,num_local+1+ee) = vmid(1:3) + ve(1:3)
          end do
          
          do ee = 1, species(ss)%hpp
            particles%species(num_local+1+species(ss)%epp+ee)    = 2
            particles%coordinates(:,num_local+1+species(ss)%epp+ee) = particles%coordinates(:,num_local+1)
            call Gaussian( particles%coordinates(4:6,num_local+1+species(ss)%epp+ee),species(ss)%v0(1:3),species(2)%upstreamkelvin,species(2)%mass )
            ! Advance velocity dt/2
            vmin(1:3) = particles%coordinates(4:6,num_local+1+species(ss)%epp+ee) - ve(1:3)
            vmid(1) = vmin(1) +vmin(2)*species(2)%tB(3)-vmin(3)*species(2)%tB(2)
            vmid(2) = vmin(2) +vmin(3)*species(2)%tB(1)-vmin(1)*species(2)%tB(3)
            vmid(3) = vmin(3) +vmin(1)*species(2)%tB(2)-vmin(2)*species(2)%tB(1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            particles%coordinates(4:6,num_local+1+species(ss)%epp+ee) = vmid(1:3) + ve(1:3)
          end do

          num_local = num_local+1+species(ss)%epp+species(ss)%hpp

        end do
      end do
  end do
end do
end if


! Do we have the -x boundary (main outflow boundary)?
if ( (xflow) .and. (Lx_min-0.5d0*dxyz(1) .lt. xmin) ) then

! Produce undisturbed particles in guard cells at -x boundary
do kk = 2-zflow_neg, Nz_local+1+zflow_pos
  do jj = 2-yflow_neg, Ny_local+1+yflow_pos
      ! Lower edge of cells
      x0 = Lx_min - dxyz(1)
      y0 = Ly_min + dxyz(2)*(jj-2)
      z0 = Lz_min + dxyz(3)*(kk-2)
      ! Loop through all ions
      do ss = 1, Nspecies

        ! Loop through all particles to be created
        do pp = 1, species(ss)%ipp
          particles%species(num_local+1)        = ss
          call random_number(rndp)
          particles%coordinates(1,num_local+1)     = x0 + dxyz(1)*rndp(1)
          particles%coordinates(2,num_local+1)     = y0 + dxyz(2)*rndp(2)
          particles%coordinates(3,num_local+1)     = z0 + dxyz(3)*rndp(3)
          call Gaussian( particles%coordinates(4:6,num_local+1),species(ss)%v0(1:3),species(ss)%upstreamkelvin,species(ss)%mass )
          ! Advance velocity dt/2
          vmin(1:3) = particles%coordinates(4:6,num_local+1) - ve(1:3)
          vmid(1) = vmin(1) + &
               vmin(2)*species(particles%species(num_local+1))%tB(3) - &
               vmin(3)*species(particles%species(num_local+1))%tB(2)
          vmid(2) = vmin(2) + &
               vmin(3)*species(particles%species(num_local+1))%tB(1) - &
               vmin(1)*species(particles%species(num_local+1))%tB(3)
          vmid(3) = vmin(3) + &
               vmin(1)*species(particles%species(num_local+1))%tB(2) - &
               vmin(2)*species(particles%species(num_local+1))%tB(1)
          if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
            vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
          end if
          particles%coordinates(4:6,num_local+1) = vmid(1:3) + ve(1:3)

          ! For each ion, create as many electrons as its charged state.
          do ee = 1, species(ss)%epp
            particles%species(num_local+1+ee)         = 1
            particles%coordinates(:,num_local+1+ee)      = particles%coordinates(:,num_local+1)
            call Gaussian( particles%coordinates(4:6,num_local+1+ee),species(ss)%v0(1:3),species(1)%upstreamkelvin,species(1)%mass )
            ! Advance velocity dt/2
            vmin(1:3) = particles%coordinates(4:6,num_local+1+ee) - ve(1:3)
            vmid(1) = vmin(1) +vmin(2)*species(1)%tB(3)-vmin(3)*species(1)%tB(2)
            vmid(2) = vmin(2) +vmin(3)*species(1)%tB(1)-vmin(1)*species(1)%tB(3)
            vmid(3) = vmin(3) +vmin(1)*species(1)%tB(2)-vmin(2)*species(1)%tB(1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            particles%coordinates(4:6,num_local+1+ee) = vmid(1:3) + ve(1:3)
          end do
          
          do ee = 1, species(ss)%hpp
            particles%species(num_local+1+species(ss)%epp+ee)    = 2
            particles%coordinates(:,num_local+1+species(ss)%epp+ee) = particles%coordinates(:,num_local+1)
            call Gaussian( particles%coordinates(4:6,num_local+1+species(ss)%epp+ee),species(ss)%v0(1:3),species(2)%upstreamkelvin,species(2)%mass )
            ! Advance velocity dt/2
            vmin(1:3) = particles%coordinates(4:6,num_local+1+species(ss)%epp+ee) - ve(1:3)
            vmid(1) = vmin(1) +vmin(2)*species(2)%tB(3)-vmin(3)*species(2)%tB(2)
            vmid(2) = vmin(2) +vmin(3)*species(2)%tB(1)-vmin(1)*species(2)%tB(3)
            vmid(3) = vmin(3) +vmin(1)*species(2)%tB(2)-vmin(2)*species(2)%tB(1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            particles%coordinates(4:6,num_local+1+species(ss)%epp+ee) = vmid(1:3) + ve(1:3)
          end do

          num_local = num_local+1+species(ss)%epp+species(ss)%hpp

        end do
      end do
  end do
end do
end if


! Produce undisturbed particles in guard cells at +y boundary
if ( yflow_pos .gt. 0 ) then

do kk = 2, Nz_local+1+zflow_pos
  do ii = 2, Nx_local+1
      ! Lower edge of cells
      x0 = Lx_min + dxyz(1)*(ii-2)
      y0 = Ly_max
      z0 = Lz_min + dxyz(3)*(kk-2)
      ! Loop through all ions
      do ss = 2, Nspecies

        ! Loop through all particles to be created
        do pp = 1, species(ss)%ipp
          particles%species(num_local+1)        = ss
          call random_number(rndp)
          particles%coordinates(1,num_local+1)     = x0 + dxyz(1)*rndp(1)
          particles%coordinates(2,num_local+1)     = y0 + dxyz(2)*rndp(2)
          particles%coordinates(3,num_local+1)     = z0 + dxyz(3)*rndp(3)
          call Gaussian( particles%coordinates(4:6,num_local+1),species(ss)%v0(1:3),species(ss)%upstreamkelvin,species(ss)%mass )
          ! Advance velocity dt/2
          vmin(1:3) = particles%coordinates(4:6,num_local+1) - ve(1:3)
          vmid(1) = vmin(1) + &
               vmin(2)*species(particles%species(num_local+1))%tB(3) - &
               vmin(3)*species(particles%species(num_local+1))%tB(2)
          vmid(2) = vmin(2) + &
               vmin(3)*species(particles%species(num_local+1))%tB(1) - &
               vmin(1)*species(particles%species(num_local+1))%tB(3)
          vmid(3) = vmin(3) + &
               vmin(1)*species(particles%species(num_local+1))%tB(2) - &
               vmin(2)*species(particles%species(num_local+1))%tB(1)
          if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
            vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
          end if
          particles%coordinates(4:6,num_local+1) = vmid(1:3) + ve(1:3)

          ! For each ion, create as many electrons as its charged state.
          do ee = 1, species(ss)%epp
            particles%species(num_local+1+ee)         = 1
            particles%coordinates(:,num_local+1+ee)      = particles%coordinates(:,num_local+1)
            call Gaussian( particles%coordinates(4:6,num_local+1+ee),species(ss)%v0(1:3),species(1)%upstreamkelvin,species(1)%mass )
            ! Advance velocity dt/2
            vmin(1:3) = particles%coordinates(4:6,num_local+1+ee) - ve(1:3)
            vmid(1) = vmin(1) +vmin(2)*species(1)%tB(3)-vmin(3)*species(1)%tB(2)
            vmid(2) = vmin(2) +vmin(3)*species(1)%tB(1)-vmin(1)*species(1)%tB(3)
            vmid(3) = vmin(3) +vmin(1)*species(1)%tB(2)-vmin(2)*species(1)%tB(1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            particles%coordinates(4:6,num_local+1+ee) = vmid(1:3) + ve(1:3)
          end do
          
          do ee = 1, species(ss)%hpp
            particles%species(num_local+1+species(ss)%epp+ee)    = 2
            particles%coordinates(:,num_local+1+species(ss)%epp+ee) = particles%coordinates(:,num_local+1)
            call Gaussian( particles%coordinates(4:6,num_local+1+species(ss)%epp+ee),species(ss)%v0(1:3),species(2)%upstreamkelvin,species(2)%mass )
            ! Advance velocity dt/2
            vmin(1:3) = particles%coordinates(4:6,num_local+1+species(ss)%epp+ee) - ve(1:3)
            vmid(1) = vmin(1) +vmin(2)*species(2)%tB(3)-vmin(3)*species(2)%tB(2)
            vmid(2) = vmin(2) +vmin(3)*species(2)%tB(1)-vmin(1)*species(2)%tB(3)
            vmid(3) = vmin(3) +vmin(1)*species(2)%tB(2)-vmin(2)*species(2)%tB(1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            particles%coordinates(4:6,num_local+1+species(ss)%epp+ee) = vmid(1:3) + ve(1:3)
          end do

          num_local = num_local+1+species(ss)%epp+species(ss)%hpp

        end do
      end do
  end do
end do
end if


! Produce undisturbed particles in guard cells at -y boundary
if ( yflow_neg .gt. 0 ) then

do kk = 2-zflow_neg, Nz_local+1
  do ii = 2, Nx_local+1
      ! Lower edge of cells
      x0 = Lx_min + dxyz(1)*(ii-2)
      y0 = Ly_min - dxyz(2)
      z0 = Lz_min + dxyz(3)*(kk-2)
      ! Loop through all ions
      do ss = 2, Nspecies

        ! Loop through all particles to be created
        do pp = 1, species(ss)%ipp
          particles%species(num_local+1)        = ss
          call random_number(rndp)
          particles%coordinates(1,num_local+1)     = x0 + dxyz(1)*rndp(1)
          particles%coordinates(2,num_local+1)     = y0 + dxyz(2)*rndp(2)
          particles%coordinates(3,num_local+1)     = z0 + dxyz(3)*rndp(3)
          call Gaussian( particles%coordinates(4:6,num_local+1),species(ss)%v0(1:3),species(ss)%upstreamkelvin,species(ss)%mass )
          ! Advance velocity dt/2
          vmin(1:3) = particles%coordinates(4:6,num_local+1) - ve(1:3)
          vmid(1) = vmin(1) + &
               vmin(2)*species(particles%species(num_local+1))%tB(3) - &
               vmin(3)*species(particles%species(num_local+1))%tB(2)
          vmid(2) = vmin(2) + &
               vmin(3)*species(particles%species(num_local+1))%tB(1) - &
               vmin(1)*species(particles%species(num_local+1))%tB(3)
          vmid(3) = vmin(3) + &
               vmin(1)*species(particles%species(num_local+1))%tB(2) - &
               vmin(2)*species(particles%species(num_local+1))%tB(1)
          if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
            vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
          end if
          particles%coordinates(4:6,num_local+1) = vmid(1:3) + ve(1:3)

          ! For each ion, create as many electrons as its charged state.
          do ee = 1, species(ss)%epp
            particles%species(num_local+1+ee)         = 1
            particles%coordinates(:,num_local+1+ee)      = particles%coordinates(:,num_local+1)
            call Gaussian( particles%coordinates(4:6,num_local+1+ee),species(ss)%v0(1:3),species(1)%upstreamkelvin,species(1)%mass )
            ! Advance velocity dt/2
            vmin(1:3) = particles%coordinates(4:6,num_local+1+ee) - ve(1:3)
            vmid(1) = vmin(1) +vmin(2)*species(1)%tB(3)-vmin(3)*species(1)%tB(2)
            vmid(2) = vmin(2) +vmin(3)*species(1)%tB(1)-vmin(1)*species(1)%tB(3)
            vmid(3) = vmin(3) +vmin(1)*species(1)%tB(2)-vmin(2)*species(1)%tB(1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            particles%coordinates(4:6,num_local+1+ee) = vmid(1:3) + ve(1:3)
          end do
          
          do ee = 1, species(ss)%hpp
            particles%species(num_local+1+species(ss)%epp+ee)    = 2
            particles%coordinates(:,num_local+1+species(ss)%epp+ee) = particles%coordinates(:,num_local+1)
            call Gaussian( particles%coordinates(4:6,num_local+1+species(ss)%epp+ee),species(ss)%v0(1:3),species(2)%upstreamkelvin,species(2)%mass )
            ! Advance velocity dt/2
            vmin(1:3) = particles%coordinates(4:6,num_local+1+species(ss)%epp+ee) - ve(1:3)
            vmid(1) = vmin(1) +vmin(2)*species(2)%tB(3)-vmin(3)*species(2)%tB(2)
            vmid(2) = vmin(2) +vmin(3)*species(2)%tB(1)-vmin(1)*species(2)%tB(3)
            vmid(3) = vmin(3) +vmin(1)*species(2)%tB(2)-vmin(2)*species(2)%tB(1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            particles%coordinates(4:6,num_local+1+species(ss)%epp+ee) = vmid(1:3) + ve(1:3)
          end do

          num_local = num_local+1+species(ss)%epp+species(ss)%hpp

        end do
      end do
  end do
end do
end if


! Produce undisturbed particles in guard cells at +z boundary
if ( zflow_pos .gt. 0 ) then

do jj = 2-yflow_neg, Ny_local+1
  do ii = 2, Nx_local+1
      ! Lower edge of cells
      x0 = Lx_min + dxyz(1)*(ii-2)
      y0 = Ly_min + dxyz(2)*(jj-2)
      z0 = Lz_max
      ! Loop through all ions
      do ss = 2, Nspecies

        ! Loop through all particles to be created
        do pp = 1, species(ss)%ipp
          particles%species(num_local+1)        = ss
          call random_number(rndp)
          particles%coordinates(1,num_local+1)     = x0 + dxyz(1)*rndp(1)
          particles%coordinates(2,num_local+1)     = y0 + dxyz(2)*rndp(2)
          particles%coordinates(3,num_local+1)     = z0 + dxyz(3)*rndp(3)
          call Gaussian( particles%coordinates(4:6,num_local+1), &
               species(ss)%v0(1:3),species(ss)%upstreamkelvin,species(ss)%mass )
          ! Advance velocity dt/2
          vmin(1:3) = particles%coordinates(4:6,num_local+1) - ve(1:3)
          vmid(1) = vmin(1) + &
               vmin(2)*species(particles%species(num_local+1))%tB(3) - &
               vmin(3)*species(particles%species(num_local+1))%tB(2)
          vmid(2) = vmin(2) + &
               vmin(3)*species(particles%species(num_local+1))%tB(1) - &
               vmin(1)*species(particles%species(num_local+1))%tB(3)
          vmid(3) = vmin(3) + &
               vmin(1)*species(particles%species(num_local+1))%tB(2) - &
               vmin(2)*species(particles%species(num_local+1))%tB(1)
          if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
             vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
          end if
          particles%coordinates(4:6,num_local+1) = vmid(1:3) + ve(1:3)

          ! For each ion, create as many electrons as its charged state.
          do ee = 1, species(ss)%epp
            particles%species(num_local+1+ee)         = 1
            particles%coordinates(:,num_local+1+ee)      = particles%coordinates(:,num_local+1)
            call Gaussian( particles%coordinates(4:6,num_local+1+ee),species(ss)%v0(1:3),species(1)%upstreamkelvin,species(1)%mass )
            ! Advance velocity dt/2
            vmin(1:3) = particles%coordinates(4:6,num_local+1+ee) - ve(1:3)
            vmid(1) = vmin(1) +vmin(2)*species(1)%tB(3)-vmin(3)*species(1)%tB(2)
            vmid(2) = vmin(2) +vmin(3)*species(1)%tB(1)-vmin(1)*species(1)%tB(3)
            vmid(3) = vmin(3) +vmin(1)*species(1)%tB(2)-vmin(2)*species(1)%tB(1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            particles%coordinates(4:6,num_local+1+ee) = vmid(1:3) + ve(1:3)
          end do
          
          do ee = 1, species(ss)%hpp
            particles%species(num_local+1+species(ss)%epp+ee)    = 2
            particles%coordinates(:,num_local+1+species(ss)%epp+ee) = particles%coordinates(:,num_local+1)
            call Gaussian( particles%coordinates(4:6,num_local+1+species(ss)%epp+ee),species(ss)%v0(1:3),species(2)%upstreamkelvin,species(2)%mass )
            ! Advance velocity dt/2
            vmin(1:3) = particles%coordinates(4:6,num_local+1+species(ss)%epp+ee) - ve(1:3)
            vmid(1) = vmin(1) +vmin(2)*species(2)%tB(3)-vmin(3)*species(2)%tB(2)
            vmid(2) = vmin(2) +vmin(3)*species(2)%tB(1)-vmin(1)*species(2)%tB(3)
            vmid(3) = vmin(3) +vmin(1)*species(2)%tB(2)-vmin(2)*species(2)%tB(1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            particles%coordinates(4:6,num_local+1+species(ss)%epp+ee) = vmid(1:3) + ve(1:3)
          end do

          num_local = num_local+1+species(ss)%epp+species(ss)%hpp

        end do
      end do
  end do
end do
end if


! Produce undisturbed particles in guard cells at -z boundary
if ( zflow_neg .gt. 0 ) then

do jj = 2, Ny_local+1+yflow_pos
  do ii = 2, Nx_local+1
      ! Lower edge of cells
      x0 = Lx_min + dxyz(1)*(ii-2)
      y0 = Ly_min + dxyz(2)*(jj-2)
      z0 = Lz_min - dxyz(3)
      ! Loop through all ions
      do ss = 2, Nspecies

        ! Loop through all particles to be created
        do pp = 1, species(ss)%ipp
          particles%species(num_local+1)        = ss
          call random_number(rndp)
          particles%coordinates(1,num_local+1)     = x0 + dxyz(1)*rndp(1)
          particles%coordinates(2,num_local+1)     = y0 + dxyz(2)*rndp(2)
          particles%coordinates(3,num_local+1)     = z0 + dxyz(3)*rndp(3)
          call Gaussian( particles%coordinates(4:6,num_local+1), &
               species(ss)%v0(1:3),species(ss)%upstreamkelvin,species(ss)%mass )
          ! Advance velocity dt/2
          vmin(1:3) = particles%coordinates(4:6,num_local+1) - ve(1:3)
          vmid(1) = vmin(1) + &
               vmin(2)*species(particles%species(num_local+1))%tB(3) - &
               vmin(3)*species(particles%species(num_local+1))%tB(2)
          vmid(2) = vmin(2) + &
               vmin(3)*species(particles%species(num_local+1))%tB(1) - &
               vmin(1)*species(particles%species(num_local+1))%tB(3)
          vmid(3) = vmin(3) + &
               vmin(1)*species(particles%species(num_local+1))%tB(2) - &
               vmin(2)*species(particles%species(num_local+1))%tB(1)
          if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
            vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
          end if
          particles%coordinates(4:6,num_local+1) = vmid(1:3) + ve(1:3)

          ! For each ion, create as many electrons as its charged state.
          do ee = 1, species(ss)%epp
            particles%species(num_local+1+ee)         = 1
            particles%coordinates(:,num_local+1+ee)      = particles%coordinates(:,num_local+1)
            call Gaussian( particles%coordinates(4:6,num_local+1+ee),species(ss)%v0(1:3),species(1)%upstreamkelvin,species(1)%mass )
            ! Advance velocity dt/2
            vmin(1:3) = particles%coordinates(4:6,num_local+1+ee) - ve(1:3)
            vmid(1) = vmin(1) +vmin(2)*species(1)%tB(3)-vmin(3)*species(1)%tB(2)
            vmid(2) = vmin(2) +vmin(3)*species(1)%tB(1)-vmin(1)*species(1)%tB(3)
            vmid(3) = vmin(3) +vmin(1)*species(1)%tB(2)-vmin(2)*species(1)%tB(1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            particles%coordinates(4:6,num_local+1+ee) = vmid(1:3) + ve(1:3)
          end do
          
          do ee = 1, species(ss)%hpp
            particles%species(num_local+1+species(ss)%epp+ee)    = 2
            particles%coordinates(:,num_local+1+species(ss)%epp+ee) = particles%coordinates(:,num_local+1)
            call Gaussian( particles%coordinates(4:6,num_local+1+species(ss)%epp+ee),species(ss)%v0(1:3),species(2)%upstreamkelvin,species(2)%mass )
            ! Advance velocity dt/2
            vmin(1:3) = particles%coordinates(4:6,num_local+1+species(ss)%epp+ee) - ve(1:3)
            vmid(1) = vmin(1) +vmin(2)*species(2)%tB(3)-vmin(3)*species(2)%tB(2)
            vmid(2) = vmin(2) +vmin(3)*species(2)%tB(1)-vmin(1)*species(2)%tB(3)
            vmid(3) = vmin(3) +vmin(1)*species(2)%tB(2)-vmin(2)*species(2)%tB(1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            particles%coordinates(4:6,num_local+1+species(ss)%epp+ee) = vmid(1:3) + ve(1:3)
          end do

          num_local = num_local+1+species(ss)%epp+species(ss)%hpp

        end do
      end do
  end do
end do
end if



return
end subroutine BCinflow


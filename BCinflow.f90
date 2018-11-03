
subroutine BCinflow(real_particles,int_particles,Lx_min,Ly_min,Lz_min,Lx_max,Ly_max,Lz_max,dxyz,kelvin,mass,v0,ve,&
                               max_per_proc,num_local,ppc,epp,hpp,ipp,Nspecies,Nx_local,Ny_local,Nz_local,&
                               xmin,xmax,ymin,ymax,zmin,zmax,tB,xflow,yflow,zflow)


  implicit none

! Parameters
  logical, intent(in) :: xflow, yflow, zflow
  integer, intent(in) :: max_per_proc, ppc(8), epp(8), hpp(8), ipp(8), Nspecies, Nx_local, Ny_local, Nz_local
  real*8, intent(in) :: Lx_min, Ly_min, Lz_min, Lx_max, Ly_max, Lz_max, dxyz(3), kelvin(8), mass(8), v0(4,8), ve(4), &
                        xmin, xmax, ymin, ymax, zmin, zmax, tB(4,8)

  real*8, intent(inout) :: real_particles(8,max_per_proc)
  integer, intent(inout) :: int_particles(max_per_proc), num_local

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
        do pp = 1, ipp(ss)
          int_particles(num_local+1)        = ss
          call random_number(rndp)
          real_particles(1,num_local+1)     = x0 + dxyz(1)*rndp(1)
          real_particles(2,num_local+1)     = y0 + dxyz(2)*rndp(2)
          real_particles(3,num_local+1)     = z0 + dxyz(3)*rndp(3)
          call Gaussian( real_particles(4:6,num_local+1),v0(1:3,ss),kelvin(ss),mass(ss) )
          ! Advance velocity dt/2
          vmin(1:3) = real_particles(4:6,num_local+1) - ve(1:3)
          vmid(1) = vmin(1) + vmin(2)*tB(3,int_particles(num_local+1))-vmin(3)*tB(2,int_particles(num_local+1))
          vmid(2) = vmin(2) + vmin(3)*tB(1,int_particles(num_local+1))-vmin(1)*tB(3,int_particles(num_local+1))
          vmid(3) = vmin(3) + vmin(1)*tB(2,int_particles(num_local+1))-vmin(2)*tB(1,int_particles(num_local+1))
          if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
            vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
          end if
          real_particles(4:6,num_local+1) = vmid(1:3) + ve(1:3)

          ! For each ion, create as many electrons as its charged state.
          do ee = 1, epp(ss)
            int_particles(num_local+1+ee)         = 1
            real_particles(:,num_local+1+ee)      = real_particles(:,num_local+1)
            call Gaussian( real_particles(4:6,num_local+1+ee),v0(1:3,ss),kelvin(1),mass(1) )
            ! Advance velocity dt/2
            vmin(1:3) = real_particles(4:6,num_local+1+ee) - ve(1:3)
            vmid(1) = vmin(1) + vmin(2)*tB(3,1)-vmin(3)*tB(2,1)
            vmid(2) = vmin(2) + vmin(3)*tB(1,1)-vmin(1)*tB(3,1)
            vmid(3) = vmin(3) + vmin(1)*tB(2,1)-vmin(2)*tB(1,1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            real_particles(4:6,num_local+1+ee) = vmid(1:3) + ve(1:3)
          end do
          
          do ee = 1, hpp(ss)
            int_particles(num_local+1+epp(ss)+ee)    = 2
            real_particles(:,num_local+1+epp(ss)+ee) = real_particles(:,num_local+1)
            call Gaussian( real_particles(4:6,num_local+1+epp(ss)+ee),v0(1:3,ss),kelvin(2),mass(2) )
            ! Advance velocity dt/2
            vmin(1:3) = real_particles(4:6,num_local+1+epp(ss)+ee) - ve(1:3)
            vmid(1) = vmin(1) + vmin(2)*tB(3,2)-vmin(3)*tB(2,2)
            vmid(2) = vmin(2) + vmin(3)*tB(1,2)-vmin(1)*tB(3,2)
            vmid(3) = vmin(3) + vmin(1)*tB(2,2)-vmin(2)*tB(1,2)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            real_particles(4:6,num_local+1+epp(ss)+ee) = vmid(1:3) + ve(1:3)
          end do

          num_local = num_local+1+epp(ss)+hpp(ss)

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
        do pp = 1, ipp(ss)
          int_particles(num_local+1)        = ss
          call random_number(rndp)
          real_particles(1,num_local+1)     = x0 + dxyz(1)*rndp(1)
          real_particles(2,num_local+1)     = y0 + dxyz(2)*rndp(2)
          real_particles(3,num_local+1)     = z0 + dxyz(3)*rndp(3)
          call Gaussian( real_particles(4:6,num_local+1),v0(1:3,ss),kelvin(ss),mass(ss) )
          ! Advance velocity dt/2
          vmin(1:3) = real_particles(4:6,num_local+1) - ve(1:3)
          vmid(1) = vmin(1) + vmin(2)*tB(3,int_particles(num_local+1))-vmin(3)*tB(2,int_particles(num_local+1))
          vmid(2) = vmin(2) + vmin(3)*tB(1,int_particles(num_local+1))-vmin(1)*tB(3,int_particles(num_local+1))
          vmid(3) = vmin(3) + vmin(1)*tB(2,int_particles(num_local+1))-vmin(2)*tB(1,int_particles(num_local+1))
          if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
            vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
          end if
          real_particles(4:6,num_local+1) = vmid(1:3) + ve(1:3)

          ! For each ion, create as many electrons as its charged state.
          do ee = 1, epp(ss)
            int_particles(num_local+1+ee)         = 1
            real_particles(:,num_local+1+ee)      = real_particles(:,num_local+1)
            call Gaussian( real_particles(4:6,num_local+1+ee),v0(1:3,ss),kelvin(1),mass(1) )
            ! Advance velocity dt/2
            vmin(1:3) = real_particles(4:6,num_local+1+ee) - ve(1:3)
            vmid(1) = vmin(1) + vmin(2)*tB(3,1)-vmin(3)*tB(2,1)
            vmid(2) = vmin(2) + vmin(3)*tB(1,1)-vmin(1)*tB(3,1)
            vmid(3) = vmin(3) + vmin(1)*tB(2,1)-vmin(2)*tB(1,1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            real_particles(4:6,num_local+1+ee) = vmid(1:3) + ve(1:3)
          end do
          
          do ee = 1, hpp(ss)
            int_particles(num_local+1+epp(ss)+ee)    = 2
            real_particles(:,num_local+1+epp(ss)+ee) = real_particles(:,num_local+1)
            call Gaussian( real_particles(4:6,num_local+1+epp(ss)+ee),v0(1:3,ss),kelvin(2),mass(2) )
            ! Advance velocity dt/2
            vmin(1:3) = real_particles(4:6,num_local+1+epp(ss)+ee) - ve(1:3)
            vmid(1) = vmin(1) + vmin(2)*tB(3,2)-vmin(3)*tB(2,2)
            vmid(2) = vmin(2) + vmin(3)*tB(1,2)-vmin(1)*tB(3,2)
            vmid(3) = vmin(3) + vmin(1)*tB(2,2)-vmin(2)*tB(1,2)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            real_particles(4:6,num_local+1+epp(ss)+ee) = vmid(1:3) + ve(1:3)
          end do

          num_local = num_local+1+epp(ss)+hpp(ss)

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
        do pp = 1, ipp(ss)
          int_particles(num_local+1)        = ss
          call random_number(rndp)
          real_particles(1,num_local+1)     = x0 + dxyz(1)*rndp(1)
          real_particles(2,num_local+1)     = y0 + dxyz(2)*rndp(2)
          real_particles(3,num_local+1)     = z0 + dxyz(3)*rndp(3)
          call Gaussian( real_particles(4:6,num_local+1),v0(1:3,ss),kelvin(ss),mass(ss) )
          ! Advance velocity dt/2
          vmin(1:3) = real_particles(4:6,num_local+1) - ve(1:3)
          vmid(1) = vmin(1) + vmin(2)*tB(3,int_particles(num_local+1))-vmin(3)*tB(2,int_particles(num_local+1))
          vmid(2) = vmin(2) + vmin(3)*tB(1,int_particles(num_local+1))-vmin(1)*tB(3,int_particles(num_local+1))
          vmid(3) = vmin(3) + vmin(1)*tB(2,int_particles(num_local+1))-vmin(2)*tB(1,int_particles(num_local+1))
          if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
            vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
          end if
          real_particles(4:6,num_local+1) = vmid(1:3) + ve(1:3)

          ! For each ion, create as many electrons as its charged state.
          do ee = 1, epp(ss)
            int_particles(num_local+1+ee)         = 1
            real_particles(:,num_local+1+ee)      = real_particles(:,num_local+1)
            call Gaussian( real_particles(4:6,num_local+1+ee),v0(1:3,ss),kelvin(1),mass(1) )
            ! Advance velocity dt/2
            vmin(1:3) = real_particles(4:6,num_local+1+ee) - ve(1:3)
            vmid(1) = vmin(1) + vmin(2)*tB(3,1)-vmin(3)*tB(2,1)
            vmid(2) = vmin(2) + vmin(3)*tB(1,1)-vmin(1)*tB(3,1)
            vmid(3) = vmin(3) + vmin(1)*tB(2,1)-vmin(2)*tB(1,1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            real_particles(4:6,num_local+1+ee) = vmid(1:3) + ve(1:3)
          end do
          
          do ee = 1, hpp(ss)
            int_particles(num_local+1+epp(ss)+ee)    = 2
            real_particles(:,num_local+1+epp(ss)+ee) = real_particles(:,num_local+1)
            call Gaussian( real_particles(4:6,num_local+1+epp(ss)+ee),v0(1:3,ss),kelvin(2),mass(2) )
            ! Advance velocity dt/2
            vmin(1:3) = real_particles(4:6,num_local+1+epp(ss)+ee) - ve(1:3)
            vmid(1) = vmin(1) + vmin(2)*tB(3,2)-vmin(3)*tB(2,2)
            vmid(2) = vmin(2) + vmin(3)*tB(1,2)-vmin(1)*tB(3,2)
            vmid(3) = vmin(3) + vmin(1)*tB(2,2)-vmin(2)*tB(1,2)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            real_particles(4:6,num_local+1+epp(ss)+ee) = vmid(1:3) + ve(1:3)
          end do

          num_local = num_local+1+epp(ss)+hpp(ss)

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
        do pp = 1, ipp(ss)
          int_particles(num_local+1)        = ss
          call random_number(rndp)
          real_particles(1,num_local+1)     = x0 + dxyz(1)*rndp(1)
          real_particles(2,num_local+1)     = y0 + dxyz(2)*rndp(2)
          real_particles(3,num_local+1)     = z0 + dxyz(3)*rndp(3)
          call Gaussian( real_particles(4:6,num_local+1),v0(1:3,ss),kelvin(ss),mass(ss) )
          ! Advance velocity dt/2
          vmin(1:3) = real_particles(4:6,num_local+1) - ve(1:3)
          vmid(1) = vmin(1) + vmin(2)*tB(3,int_particles(num_local+1))-vmin(3)*tB(2,int_particles(num_local+1))
          vmid(2) = vmin(2) + vmin(3)*tB(1,int_particles(num_local+1))-vmin(1)*tB(3,int_particles(num_local+1))
          vmid(3) = vmin(3) + vmin(1)*tB(2,int_particles(num_local+1))-vmin(2)*tB(1,int_particles(num_local+1))
          if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
            vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
          end if
          real_particles(4:6,num_local+1) = vmid(1:3) + ve(1:3)

          ! For each ion, create as many electrons as its charged state.
          do ee = 1, epp(ss)
            int_particles(num_local+1+ee)         = 1
            real_particles(:,num_local+1+ee)      = real_particles(:,num_local+1)
            call Gaussian( real_particles(4:6,num_local+1+ee),v0(1:3,ss),kelvin(1),mass(1) )
            ! Advance velocity dt/2
            vmin(1:3) = real_particles(4:6,num_local+1+ee) - ve(1:3)
            vmid(1) = vmin(1) + vmin(2)*tB(3,1)-vmin(3)*tB(2,1)
            vmid(2) = vmin(2) + vmin(3)*tB(1,1)-vmin(1)*tB(3,1)
            vmid(3) = vmin(3) + vmin(1)*tB(2,1)-vmin(2)*tB(1,1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            real_particles(4:6,num_local+1+ee) = vmid(1:3) + ve(1:3)
          end do
          
          do ee = 1, hpp(ss)
            int_particles(num_local+1+epp(ss)+ee)    = 2
            real_particles(:,num_local+1+epp(ss)+ee) = real_particles(:,num_local+1)
            call Gaussian( real_particles(4:6,num_local+1+epp(ss)+ee),v0(1:3,ss),kelvin(2),mass(2) )
            ! Advance velocity dt/2
            vmin(1:3) = real_particles(4:6,num_local+1+epp(ss)+ee) - ve(1:3)
            vmid(1) = vmin(1) + vmin(2)*tB(3,2)-vmin(3)*tB(2,2)
            vmid(2) = vmin(2) + vmin(3)*tB(1,2)-vmin(1)*tB(3,2)
            vmid(3) = vmin(3) + vmin(1)*tB(2,2)-vmin(2)*tB(1,2)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            real_particles(4:6,num_local+1+epp(ss)+ee) = vmid(1:3) + ve(1:3)
          end do

          num_local = num_local+1+epp(ss)+hpp(ss)

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
        do pp = 1, ipp(ss)
          int_particles(num_local+1)        = ss
          call random_number(rndp)
          real_particles(1,num_local+1)     = x0 + dxyz(1)*rndp(1)
          real_particles(2,num_local+1)     = y0 + dxyz(2)*rndp(2)
          real_particles(3,num_local+1)     = z0 + dxyz(3)*rndp(3)
          call Gaussian( real_particles(4:6,num_local+1),v0(1:3,ss),kelvin(ss),mass(ss) )
          ! Advance velocity dt/2
          vmin(1:3) = real_particles(4:6,num_local+1) - ve(1:3)
          vmid(1) = vmin(1) + vmin(2)*tB(3,int_particles(num_local+1))-vmin(3)*tB(2,int_particles(num_local+1))
          vmid(2) = vmin(2) + vmin(3)*tB(1,int_particles(num_local+1))-vmin(1)*tB(3,int_particles(num_local+1))
          vmid(3) = vmin(3) + vmin(1)*tB(2,int_particles(num_local+1))-vmin(2)*tB(1,int_particles(num_local+1))
          if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
            vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
          end if
          real_particles(4:6,num_local+1) = vmid(1:3) + ve(1:3)

          ! For each ion, create as many electrons as its charged state.
          do ee = 1, epp(ss)
            int_particles(num_local+1+ee)         = 1
            real_particles(:,num_local+1+ee)      = real_particles(:,num_local+1)
            call Gaussian( real_particles(4:6,num_local+1+ee),v0(1:3,ss),kelvin(1),mass(1) )
            ! Advance velocity dt/2
            vmin(1:3) = real_particles(4:6,num_local+1+ee) - ve(1:3)
            vmid(1) = vmin(1) + vmin(2)*tB(3,1)-vmin(3)*tB(2,1)
            vmid(2) = vmin(2) + vmin(3)*tB(1,1)-vmin(1)*tB(3,1)
            vmid(3) = vmin(3) + vmin(1)*tB(2,1)-vmin(2)*tB(1,1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            real_particles(4:6,num_local+1+ee) = vmid(1:3) + ve(1:3)
          end do
          
          do ee = 1, hpp(ss)
            int_particles(num_local+1+epp(ss)+ee)    = 2
            real_particles(:,num_local+1+epp(ss)+ee) = real_particles(:,num_local+1)
            call Gaussian( real_particles(4:6,num_local+1+epp(ss)+ee),v0(1:3,ss),kelvin(2),mass(2) )
            ! Advance velocity dt/2
            vmin(1:3) = real_particles(4:6,num_local+1+epp(ss)+ee) - ve(1:3)
            vmid(1) = vmin(1) + vmin(2)*tB(3,2)-vmin(3)*tB(2,2)
            vmid(2) = vmin(2) + vmin(3)*tB(1,2)-vmin(1)*tB(3,2)
            vmid(3) = vmin(3) + vmin(1)*tB(2,2)-vmin(2)*tB(1,2)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            real_particles(4:6,num_local+1+epp(ss)+ee) = vmid(1:3) + ve(1:3)
          end do

          num_local = num_local+1+epp(ss)+hpp(ss)

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
        do pp = 1, ipp(ss)
          int_particles(num_local+1)        = ss
          call random_number(rndp)
          real_particles(1,num_local+1)     = x0 + dxyz(1)*rndp(1)
          real_particles(2,num_local+1)     = y0 + dxyz(2)*rndp(2)
          real_particles(3,num_local+1)     = z0 + dxyz(3)*rndp(3)
          call Gaussian( real_particles(4:6,num_local+1),v0(1:3,ss),kelvin(ss),mass(ss) )
          ! Advance velocity dt/2
          vmin(1:3) = real_particles(4:6,num_local+1) - ve(1:3)
          vmid(1) = vmin(1) + vmin(2)*tB(3,int_particles(num_local+1))-vmin(3)*tB(2,int_particles(num_local+1))
          vmid(2) = vmin(2) + vmin(3)*tB(1,int_particles(num_local+1))-vmin(1)*tB(3,int_particles(num_local+1))
          vmid(3) = vmin(3) + vmin(1)*tB(2,int_particles(num_local+1))-vmin(2)*tB(1,int_particles(num_local+1))
          if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
            vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
          end if
          real_particles(4:6,num_local+1) = vmid(1:3) + ve(1:3)

          ! For each ion, create as many electrons as its charged state.
          do ee = 1, epp(ss)
            int_particles(num_local+1+ee)         = 1
            real_particles(:,num_local+1+ee)      = real_particles(:,num_local+1)
            call Gaussian( real_particles(4:6,num_local+1+ee),v0(1:3,ss),kelvin(1),mass(1) )
            ! Advance velocity dt/2
            vmin(1:3) = real_particles(4:6,num_local+1+ee) - ve(1:3)
            vmid(1) = vmin(1) + vmin(2)*tB(3,1)-vmin(3)*tB(2,1)
            vmid(2) = vmin(2) + vmin(3)*tB(1,1)-vmin(1)*tB(3,1)
            vmid(3) = vmin(3) + vmin(1)*tB(2,1)-vmin(2)*tB(1,1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            real_particles(4:6,num_local+1+ee) = vmid(1:3) + ve(1:3)
          end do
          
          do ee = 1, hpp(ss)
            int_particles(num_local+1+epp(ss)+ee)    = 2
            real_particles(:,num_local+1+epp(ss)+ee) = real_particles(:,num_local+1)
            call Gaussian( real_particles(4:6,num_local+1+epp(ss)+ee),v0(1:3,ss),kelvin(2),mass(2) )
            ! Advance velocity dt/2
            vmin(1:3) = real_particles(4:6,num_local+1+epp(ss)+ee) - ve(1:3)
            vmid(1) = vmin(1) + vmin(2)*tB(3,2)-vmin(3)*tB(2,2)
            vmid(2) = vmin(2) + vmin(3)*tB(1,2)-vmin(1)*tB(3,2)
            vmid(3) = vmin(3) + vmin(1)*tB(2,2)-vmin(2)*tB(1,2)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            real_particles(4:6,num_local+1+epp(ss)+ee) = vmid(1:3) + ve(1:3)
          end do

          num_local = num_local+1+epp(ss)+hpp(ss)

        end do
      end do
  end do
end do
end if



return
end subroutine BCinflow


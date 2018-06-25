
! Removes particles that are outside the grid

      
subroutine BCoutflow(real_particles,int_particles,Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,&
                       xmin,xmax,ymin,ymax,zmin,zmax,max_per_proc,num_local,xflow,yflow,zflow)

  implicit none


! Parameters

  logical, intent(in) :: xflow,yflow,zflow
  real*8, intent(in) :: Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,xmin,xmax,ymin,ymax,zmin,zmax
  integer, intent(in) :: max_per_proc

  real*8, intent(inout) :: real_particles(8,max_per_proc)
  integer, intent(inout) :: int_particles(max_per_proc), num_local
      
! Local variables
  integer :: pp
  logical :: delete

  do pp = num_local, 1, -1
    
    delete = .false.

    if (xflow) then
     
      ! Remove particles outside x-boundaries
      if ( (real_particles(1,pp) .lt. xmin) .or. (real_particles(1,pp) .gt. xmax) ) then
        delete = .true.
      end if

    end if

    if (yflow) then
     
      ! Remove particles outside y-boundaries
      if ( (real_particles(2,pp) .lt. ymin) .or. (real_particles(2,pp) .gt. ymax) ) then
        delete = .true.
      end if

    end if

    if (zflow) then
     
      ! Remove particles outside z-boundaries
      if ( (real_particles(3,pp) .lt. zmin) .or. (real_particles(3,pp) .gt. zmax) ) then
        delete = .true.
      end if

    end if

    ! Remove particles farther outside the local domain than thy neighbours.
    if ( (real_particles(1,pp) .lt. 2.0d0*Lx_min-Lx_max) .or. (real_particles(1,pp) .gt. 2.0d0*Lx_max-Lx_min) ) then
      delete = .true.
    end if    

    if ( (real_particles(2,pp) .lt. 2.0d0*Ly_min-Ly_max) .or. (real_particles(2,pp) .gt. 2.0d0*Ly_max-Ly_min) ) then
      delete = .true.
    end if    

    if ( (real_particles(3,pp) .lt. 2.0d0*Lz_min-Lz_max) .or. (real_particles(3,pp) .gt. 2.0d0*Lz_max-Lz_min) ) then
      delete = .true.
    end if    

!    if ( any( isnan(real_particles(:,pp)) ) ) then
!      delete = .true.
!    end if    

    ! Actual removal
    if (delete) then
      if (pp .eq. num_local) then
        num_local = num_local - 1
      else
        real_particles(:,pp) = real_particles(:,num_local)
        int_particles(pp) = int_particles(num_local)
        num_local = num_local - 1
      end if
    end if

  end do
    

  return
end subroutine BCoutflow

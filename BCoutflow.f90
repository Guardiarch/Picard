! Removes particles that are outside the grid

subroutine BCoutflow(particles,Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max, &
     xmin,xmax,ymin,ymax,zmin,zmax,max_per_proc,num_local,xflow,yflow,zflow)

  use SpecificTypes
  implicit none


! Parameters

  logical, intent(in) :: xflow,yflow,zflow
  real*8, intent(in) :: Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max, &
       xmin,xmax,ymin,ymax,zmin,zmax
  integer, intent(in) :: max_per_proc
  integer, intent(inout) :: num_local
  type(particlearrays) particles

! Local variables
  integer :: pp
  logical :: delete

  do pp = num_local, 1, -1
    
    delete = .false.

    if (xflow) then
     
      ! Remove particles outside x-boundaries
       if ( (particles%coordinates(1,pp) .lt. xmin) .or. &
            (particles%coordinates(1,pp) .gt. xmax) ) then
          delete = .true.
       end if

    end if

    if (yflow) then
     
      ! Remove particles outside y-boundaries
       if ( (particles%coordinates(2,pp) .lt. ymin) .or. &
            (particles%coordinates(2,pp) .gt. ymax) ) then
          delete = .true.
       end if

    end if

    if (zflow) then
     
      ! Remove particles outside z-boundaries
       if ( (particles%coordinates(3,pp) .lt. zmin) .or. &
            (particles%coordinates(3,pp) .gt. zmax) ) then
          delete = .true.
       end if

    end if

    ! Remove particles farther outside the local domain than thy neighbours.
    if ( (particles%coordinates(1,pp) .lt. 2.0d0*Lx_min-Lx_max) .or. &
         (particles%coordinates(1,pp) .gt. 2.0d0*Lx_max-Lx_min) ) then
       delete = .true.
    end if

    if ( (particles%coordinates(2,pp) .lt. 2.0d0*Ly_min-Ly_max) .or. &
         (particles%coordinates(2,pp) .gt. 2.0d0*Ly_max-Ly_min) ) then
       delete = .true.
    end if

    if ( (particles%coordinates(3,pp) .lt. 2.0d0*Lz_min-Lz_max) .or. &
         (particles%coordinates(3,pp) .gt. 2.0d0*Lz_max-Lz_min) ) then
       delete = .true.
    end if

    ! Actual removal
    if (delete) then
      if (pp .eq. num_local) then
        num_local = num_local - 1
      else
        particles%coordinates(:,pp) = particles%coordinates(:,num_local)
        particles%species(pp) = particles%species(num_local)
        num_local = num_local - 1
      end if
    end if

  end do
    

  return
end subroutine BCoutflow

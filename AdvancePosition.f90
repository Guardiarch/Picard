! Particle pusher.

subroutine AdvancePosition(particles, dt, max_per_proc, num_local)

  use SpecificTypes
  implicit none


! Parameters
  integer, intent(in) :: max_per_proc, num_local
  real*8, intent(in) :: dt
  type(particlearrays) particles

  integer ii
  
  ! Push!
  do ii = 1, num_local
     particles%coordinates(1:3, ii) = particles%coordinates(1:3, ii) + &
          particles%coordinates(4:6, ii)*dt
  end do


  return
end subroutine AdvancePosition

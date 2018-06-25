
! Particle pusher.

subroutine AdvancePosition(real_particles, dt, max_per_proc, num_local)

  implicit none


! Parameters
  integer, intent(in) :: max_per_proc, num_local
  real*8, intent(inout) :: real_particles(8,max_per_proc)
  real*8, intent(in) :: dt

! Push!
  real_particles(1:3, 1:num_local) = real_particles(1:3, 1:num_local) + real_particles(4:6, 1:num_local)*dt


  return
end subroutine AdvancePosition

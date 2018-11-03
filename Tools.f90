
subroutine Gaussian(u, v, t, m)
    ! U is a random velocity from the Maxwellian with bulk velocity V, 
    ! and temperature T for a specie of mass M. 
    implicit none
    real*8, intent(out) :: u(3)         ! [m/s]
    real*8, intent(in)  :: v(3), t, m   ! [m/s], [K], [kg]

    real*8 :: a, pw
    real*8, parameter :: kb   = 1.38064852d-23

    if (t .gt. 0.0d0 .and. m .gt. 0.0d0) then

      a = m/(2.0d0*kb*t)
      pw = dsqrt( a**(-1) )         ! approx. width of distribution

      call normal(u(1))
      call normal(u(2))
      call normal(u(3))
      u = u*pw/dsqrt(2.0d0)

      u = u + v

    else
      u = v
    end if

end subroutine Gaussian

subroutine normal(u)
    ! Gaussian distribution by Box-Muller. 
    ! We are wasting one random nmber, see 001219/boxmuller.c
    ! Also, it can be done without trig.  See NR. 
    implicit none
    real*8, intent(out) :: u
    real*8, parameter   :: tau  = 6.2831853071795864769252867665590d0  ! 2*pi
    real*8              :: r1, r2
    call random_number(r1)
    do
       call random_number(r2)
       if (r2 .gt. 0.0d0) then
         if (-2.0d0*dlog(r2) .ge. 0.0d0) then
           exit
         end if
       end if
    end do
    u = dsqrt(-2.0d0*dlog(r2))*dsin(tau*r1)
end subroutine normal


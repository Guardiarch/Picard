subroutine CalcDistance(D,Nx,Ny,Nz,dxyz,xflow,yflow,zflow)

  implicit none

! Parameters
  logical, intent(in) :: xflow, yflow, zflow
  integer, intent(in) :: Nx, Ny, Nz
  real*8, intent(in) :: dxyz(3)
  real*8, intent(out) :: D(2*Nx+1,2*Ny+1,2*Nz+1)

! Local variables
  integer :: ii, jj, kk, pi, pj, pk, ni, nj, nk
  real*8 :: di, dj, dk, temp

  temp = 0.0d0
  D = 0.0d0

  do kk = 1, Nz+1
    do jj = 1, Ny+1
      do ii = 1, Nx+1
        pi = Nx+1+(ii-1)
        pj = Ny+1+(jj-1)
        pk = Nz+1+(kk-1)
        ni = Nx+1-(ii-1)
        nj = Ny+1-(jj-1)
        nk = Nz+1-(kk-1)
        di = dble(ii-1)
        dj = dble(jj-1)
        dk = dble(kk-1)

! Don't save zero distance
        if ( ii .eq. 1 .and. jj .eq. 1 .and. kk .eq. 1 ) then
          temp = 0.0d0
        else
          temp = dsqrt( (di*dxyz(1))**2 + (dj*dxyz(2))**2 + (dk*dxyz(3))**2 )**(-1)
        end if

        if (xflow) then
          if (.not. yflow .and. .not. zflow) then
            temp = temp + dsqrt( (di*dxyz(1))**2 + ((Ny-dj)*dxyz(2))**2 + (dk*dxyz(3))**2 )**(-1)
            temp = temp + dsqrt( (di*dxyz(1))**2 + (dj*dxyz(2))**2      + ((Nz-dk)*dxyz(3))**2 )**(-1)
            temp = temp + dsqrt( (di*dxyz(1))**2 + ((Ny-dj)*dxyz(2))**2 + ((Nz-dk)*dxyz(3))**2 )**(-1)
          else if (.not. yflow) then
            temp = temp + dsqrt( (di*dxyz(1))**2 + ((Ny-dj)*dxyz(2))**2 + (dk*dxyz(3))**2 )**(-1)
          else if (.not. zflow) then
            temp = temp + dsqrt( (di*dxyz(1))**2 + (dj*dxyz(2))**2      + ((Nz-dk)*dxyz(3))**2 )**(-1)
          end if
        else if (yflow) then
          if (.not. zflow) then
            temp = temp + dsqrt( ((Nx-di)*dxyz(1))**2 + (dj*dxyz(2))**2 + (dk*dxyz(3))**2 )**(-1)
            temp = temp + dsqrt( (di*dxyz(1))**2      + (dj*dxyz(2))**2 + ((Nz-dk)*dxyz(3))**2 )**(-1)
            temp = temp + dsqrt( ((Nx-di)*dxyz(1))**2 + (dj*dxyz(2))**2 + ((Nz-dk)*dxyz(3))**2 )**(-1)
          else
            temp = temp + dsqrt( ((Nx-di)*dxyz(1))**2 + (dj*dxyz(2))**2 + (dk*dxyz(3))**2 )**(-1)
          end if
        else if (zflow) then
            temp = temp + dsqrt( ((Nx-di)*dxyz(1))**2 + (dj*dxyz(2))**2      + (dk*dxyz(3))**2 )**(-1)
            temp = temp + dsqrt( (di*dxyz(1))**2      + ((Ny-dj)*dxyz(2))**2 + (dk*dxyz(3))**2 )**(-1)
            temp = temp + dsqrt( ((Nx-di)*dxyz(1))**2 + ((Ny-dj)*dxyz(2))**2 + (dk*dxyz(3))**2 )**(-1)
        else
            temp = temp + dsqrt( ((Nx-di)*dxyz(1))**2 + (dj*dxyz(2))**2      + (dk*dxyz(3))**2 )**(-1)
            temp = temp + dsqrt( (di*dxyz(1))**2      + ((Ny-dj)*dxyz(2))**2 + (dk*dxyz(3))**2 )**(-1)
            temp = temp + dsqrt( (di*dxyz(1))**2      + (dj*dxyz(2))**2      + ((Nz-dk)*dxyz(3))**2 )**(-1)
            temp = temp + dsqrt( ((Nx-di)*dxyz(1))**2 + ((Ny-dj)*dxyz(2))**2 + (dk*dxyz(3))**2 )**(-1)      
            temp = temp + dsqrt( (di*dxyz(1))**2      + ((Ny-dj)*dxyz(2))**2 + ((Nz-dk)*dxyz(3))**2 )**(-1) 
            temp = temp + dsqrt( ((Nx-di)*dxyz(1))**2 + (dj*dxyz(2))**2      + ((Nz-dk)*dxyz(3))**2 )**(-1)  
            temp = temp + dsqrt( ((Nx-di)*dxyz(1))**2 + ((Ny-dj)*dxyz(2))**2 + ((Nz-dk)*dxyz(3))**2 )**(-1)
        end if

        D(pi,pj,pk) = 1.0d0*temp
        D(ni,pj,pk) = 1.0d0*temp
        D(pi,nj,pk) = 1.0d0*temp
        D(pi,pj,nk) = 1.0d0*temp
        D(ni,nj,pk) = 1.0d0*temp
        D(pi,nj,nk) = 1.0d0*temp
        D(ni,pj,nk) = 1.0d0*temp
        D(ni,nj,nk) = 1.0d0*temp

      end do
    end do
  end do


  return
end subroutine CalcDistance

subroutine ParticlesProduce(real_particles,int_particles, &
     Lx_min,Ly_min,Lz_min,dxyz,kelvin,mass,nu_d,nu_i,Ek,Qn,vn,v0,ve,&
     max_per_proc,num_local,ppc,epp,hpp,ipp,Nspecies, &
     Nx_local,Ny_local,Nz_local,dt,n2p,flatradius) 


  implicit none

  ! Parameters
  integer, intent(in) :: max_per_proc, ppc(8), epp(8), hpp(8), ipp(8), &
       Nspecies, Nx_local, Ny_local, Nz_local
  real*8, intent(in) :: Lx_min, Ly_min, Lz_min, dxyz(3), kelvin(8), &
       mass(8), v0(4,8), ve(4), vn(8), Qn(8), nu_i(8), nu_d(8), &
       dt, n2p(8), Ek(8), flatradius

  real*8, intent(inout) :: real_particles(8,max_per_proc)
  integer, intent(inout) :: int_particles(max_per_proc), num_local

  real*8, parameter :: tau  = 6.2831853071795864769252867665590d0  ! 2*pi

  ! Local variables
  integer ee, pp, ss, ii, jj, kk, num_create
  real*8 rnd, rndp(7), x0, y0, z0, vmin(3), vmid(3), r1, r2, nden, dV, &
       ui(8), ue(8), uh(8)

  dV = dxyz(1)*dxyz(2)*dxyz(3)
  ue = 0.0d0
  uh = 0.0d0
  ui = 0.0d0

  do ss = 1, Nspecies
     ue(ss)=dsqrt(2.0d0*mass(ss)*Ek(ss)*(mass(1)*(mass(ss)+mass(1)))**(-1.0d0))
     uh(ss)=dsqrt(2.0d0*mass(ss)*Ek(ss)*(mass(2)*(mass(ss)+mass(2)))**(-1.0d0))
     if (epp(ss) .gt. 0) then
        ui(ss) = dsqrt( 2.0d0*mass(1)*Ek(ss)* &
             (mass(ss)*(mass(ss)+mass(1)))**(-1.0d0) )
     else if (hpp(ss) .gt. 0) then
        ui(ss) = dsqrt( 2.0d0*mass(2)*Ek(ss) * &
             (mass(ss)*(mass(ss)+mass(2)))**(-1.0d0) )
     else
        ui(ss)  = 0.0d0
     end if
  end do

  ! Loop through each cell which is not a guard cell
  do kk = 2, Nz_local+1
     do jj = 2, Ny_local+1
        do ii = 2, Nx_local+1
           ! Lower edge of cells
           x0 = Lx_min + dxyz(1)*(ii-2)
           y0 = Ly_min + dxyz(2)*(jj-2)
           z0 = Lz_min + dxyz(3)*(kk-2)
           r2 =(x0+0.5d0*dxyz(1))**2 + (y0+0.5d0*dxyz(2))**2.0d0 + &
                (z0+0.5d0*dxyz(3))**2
           r1 = dsqrt(r2)
           if ( (4.0d0*r1)**3.0d0 .lt. dV) then
              r1 = ( dsqrt((x0+0.25d0*dxyz(1))**2.0d0 + &
                   (y0+0.25d0*dxyz(2))**2.0d0 + &
                   (z0+0.25d0*dxyz(3))**2.0d0) &
                   + dsqrt((x0+0.25d0*dxyz(1))**2.0d0 + &
                   (y0+0.25d0*dxyz(2))**2.0d0 + &
                   (z0+0.75d0*dxyz(3))**2) &
                   + dsqrt((x0+0.25d0*dxyz(1))**2.0d0 + &
                   (y0+0.75d0*dxyz(2))**2.0d0 + &
                   (z0+0.25d0*dxyz(3))**2) + &
                   dsqrt((x0+0.75d0*dxyz(1))**2.0d0 + &
                   (y0+0.25d0*dxyz(2))**2.0d0 + &
                   (z0+0.25d0*dxyz(3))**2) + &
                   dsqrt((x0+0.25d0*dxyz(1))**2.0d0 + &
                   (y0+0.75d0*dxyz(2))**2.0d0 + &
                   (z0+0.75d0*dxyz(3))**2) + &
                   dsqrt((x0+0.75d0*dxyz(1))**2.0d0 + &
                   (y0+0.25d0*dxyz(2))**2.0d0 + &
                   (z0+0.75d0*dxyz(3))**2) + &
                   dsqrt((x0+0.75d0*dxyz(1))**2.0d0 + &
                   (y0+0.75d0*dxyz(2))**2.0d0 + &
                   (z0+0.25d0*dxyz(3))**2) + &
                   dsqrt((x0+0.75d0*dxyz(1))**2.0d0 + &
                   (y0+0.75d0*dxyz(2))**2.0d0 + &
                   (z0+0.75d0*dxyz(3))**2) )*0.125d0
              r2 = r1**2.0d0
           end if

           do ss = 1, Nspecies
              if ( vn(ss) .gt. 0.0d0 ) then
                 if ( r1 > flatradius ) then
                    nden = Qn(ss)/( 2.0d0*tau*r2*vn(ss) ) * &
                         dexp( - nu_d(ss)*r1/vn(ss) )
                 else
                    nden = Qn(ss)/( 2.0d0*tau*flatradius**2.0d0*vn(ss) ) * &
                         dexp( - nu_d(ss)*flatradius/vn(ss) )
                 end if
              else
                 nden = 0.0d0
              end if

              if ( n2p(ss) .gt. 0.0d0 ) then
                 call random_number(rnd)
                 num_create = idint( nu_i(ss)*nden*dV/n2p(ss) + rnd )
              else
                 num_create = 0
              end if

              do pp = 1, num_create
                 call random_number(rndp)
                 int_particles(num_local+1)      = ss
                 real_particles(1,num_local+1)   = x0 + dxyz(1)*rndp(1)
                 real_particles(2,num_local+1)   = y0 + dxyz(2)*rndp(2)
                 real_particles(3,num_local+1)   = z0 + dxyz(3)*rndp(3)
                 real_particles(4:6,num_local+1) = vn(ss) * &
                      real_particles(1:3,num_local+1) / &
                      dsqrt(sum(real_particles(1:3,num_local+1)**2.0d0)) + &
                      ui(ss) * rndp(4:6)/dsqrt(sum(rndp(4:6)**2.0d0))
                 real_particles(1:3,num_local+1) = &
                      real_particles(1:3,num_local+1) + &
                      real_particles(4:6,num_local+1)*dt*rndp(7)
                 ! For each ion, create as many electrons as its charged state.
                 do ee = 1, epp(ss)
                    int_particles(num_local+1+ee)        = 1 
                    real_particles(4:6,num_local+1+ee)   = &
                         real_particles(4:6,num_local+1) - &
                         (ue(ss)+ui(ss))*rndp(4:6)/dsqrt(sum(rndp(4:6)**2.0d0))
                    real_particles(1:3,num_local+1+ee)   = &
                         real_particles(1:3,num_local+1) &
                         - real_particles(4:6,num_local+1)*dt*rndp(7) &
                         + real_particles(4:6,num_local+1+ee)*dt*rndp(7)
                 end do
                 do ee = 1, hpp(ss)
                    int_particles(num_local+1+epp(ss)+ee)       = 2
                    real_particles(4:6,num_local+1+epp(ss)+ee)  = &
                         real_particles(4:6,num_local+1) - &
                         (uh(ss)+ui(ss))*rndp(4:6)/dsqrt(sum(rndp(4:6)**2.0d0))
                    real_particles(1:3,num_local+1+epp(ss)+ee)  = &
                         real_particles(1:3,num_local+1) &
                         - real_particles(4:6,num_local+1)*dt*rndp(7) &
                         + real_particles(4:6,num_local+1+epp(ss)+ee)*dt*rndp(7)
                 end do
                 num_local = num_local+1+epp(ss)+hpp(ss)
              end do

           end do
        end do
     end do
  end do



  return
end subroutine ParticlesProduce

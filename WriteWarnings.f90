!
! Calculates plasma parameters and issues warnings
!
subroutine WriteWarnings(mass,charge,eta,kelvin,density,v0,ve,dt,n2p,dxyz,B0,Nx,Ny,Nz,xmin,xmax,ymin,ymax,zmin,zmax,&
                          iprocs,jprocs,kprocs,gyration,xflow,yflow,zflow,Qn,vn,Ek,nu_i,nu_d,ppc)


! No implicit variables
  implicit none

! Parameters
  logical, intent(in) :: gyration(8), xflow, yflow, zflow
  integer, intent(in) :: Nx, Ny, Nz, iprocs, jprocs, kprocs, ppc(8)
  real*8, intent(in) :: mass(8), charge(8), eta(8), density(8), kelvin(8), v0(4,8), ve(4), dt, &
       n2p(8), dxyz(3), B0(4), xmin, xmax, ymin, ymax, zmin, zmax, Qn(8), vn(8), Ek(8), nu_i(8), nu_d(8)


! Local variables
  real*8 :: freq, wp(3), wc(3), va, vc(3), Ld(3), rL(3), rPU(3), Li(3), rM(3), CFL, sd, vth(8)
  real*8 :: dt_opt_part, dt_opt_wave, dt_opt, dt_c
  
  integer :: Nprocs

! courant = dx/dt, if courant < vth(1) we should warn. 
! Ld = sqrt(kT/e0n) , warn if Ld > dx/2.
! Rg = mv/qB, warn if Rg < dx. 
! fp = sqrt(q*q*n/m), warn if fpe > 1/dt
! fg = mv/qB, warn if fg > 1/dt

  real*8, parameter :: tau       = 6.2831853071795864769252867665590d0 ! 2*pi
  real*8, parameter :: ec        = 1.60217646d-19
  real*8, parameter :: mp        = 1.6726216d-27
  real*8, parameter :: me        = 9.1093819d-31
  real*8, parameter :: c0        = 2.99792458d8
  real*8, parameter :: mu0       = 1.25663706143591729538505735331180d-6
  real*8, parameter :: eps0      = 8.854187817d-12
  real*8, parameter :: kb        = 1.380650d-23


  vth       = dsqrt(2.0d0*kb*kelvin/(mass+1.0d-99))
! Calculate local variables
  Nprocs    = iprocs*jprocs*kprocs
  freq      = dt**(-1)
! Plasma frequency
  wp(1)     = dsqrt( density(1)*charge(1)*charge(1)/mass(1)/eps0 )
  wp(2)     = maxval( dsqrt( density(2:)*charge(2:)*charge(2:)/mass(2:)/eps0 ) )
  wp(3)     = dsqrt( density(1)*ec**2/me/eps0 )
! Cyclotron frequency
  wc(1)     = dabs(charge(1))*B0(4)/mass(1)
  wc(2)     = maxval( dabs(charge(2:))*B0(4)/mass(2:) )
  wc(3)     = ec*B0(4)/me
! Alfven speed
  va        = B0(4)/dsqrt( mu0*sum( mass(:)*density(:) ) )
! Sound speedish
  vc(1)     = vth(1)*dsqrt(0.5d0)
  vc(2)     = maxval( vth(2:)*dsqrt(0.5d0) )
  vc(3)     = vth(1)*dsqrt( 0.5d0*mass(1)/me )

  if (.not. gyration(1)) then
    wc(1) = wc(2)
  end if

! Debye length
  Ld(:)     = vc(:)/wp(:)
! Gyroradius
  rL(:)     = vc(:)/wc(:)
! Gyroradius (pick-up)
  rPU(:)    = v0(4,1)/wc(:)
! Inertial length
  Li(:)     = va/wc(:)
  sd        = c0/wp(1)
! Mystery length
  rM(:)     = va/wp(:)
! CFL condition
  CFL       = dxyz(1)/dt

! Optimal dt
  dt_c   = dsqrt(3.0d0)**(-1)*dxyz(1)/c0
  dt_opt_part = dsqrt(3.0d0)**(-1)*dxyz(1)*( max(2.0d0*v0(4,1),2.0d0*v0(4,2)) + dsqrt(12.0d0*tau**(-1))*max(vc(1),vc(2)) )**(-1)
  dt_opt_wave = tau*dsqrt(12.0d0)**(-1)*max( wc(1),wc(2),wp(1),wp(2) )**(-1)
  dt_opt = max( min( dt_opt_part, dt_opt_wave ), dt_c )


! Begin by giving summary of parameters
  write(*,*) ' ' 
  write(*,*) '****************************************************************'
  write(*,*) 'Summary of physical quantities'
  write(*,*) 'Plasma freq., wpe (physical) = ', wp(3)
  write(*,*) 'Plasma freq., wpe (simulation) = ', wp(1)
  write(*,*) 'Plasma freq., wpi = ', wp(2)
  write(*,*) 'Cyclotron freq. (physical) , wce = ', wc(3)
  write(*,*) 'Cyclotron freq. (simulation), wce = ', wc(1)
  write(*,*) 'Cyclotron freq., wci = ', wc(2)
  write(*,*) 'Alfven speed, va = ', va
  write(*,*) 'Cyclotron speed (physical) , vce = ', vc(3)
  write(*,*) 'Cyclotron speed (simulation), vce = ', vc(1)
  write(*,*) 'Cyclotron speed, vci = ', vc(2)
  write(*,*) 'Debye length (physical) , Lde = ', Ld(3)
  write(*,*) 'Debye length (simulation), Lde = ', Ld(1)
  write(*,*) 'Debye length, Ldi = ', Ld(2)
  write(*,*) 'Gyroradius (physical) , rLe = ', rL(3)
  write(*,*) 'Gyroradius (simulation), rLe = ', rL(1)
  write(*,*) 'Gyroradius, rLi = ', rL(2)
  write(*,*) 'Gyroradius pick-up (physical) , rLe = ', rPU(3)
  write(*,*) 'Gyroradius pick-up (simulation), rLe = ', rPU(1)
  write(*,*) 'Gyroradius pick-up, rLi = ', rPU(2)
  write(*,*) 'Inertial length (physical) , Lie = ', Li(3)
  write(*,*) 'Inertial length (simulation), Lie = ', Li(1)
  write(*,*) 'Inertial length, Lii = ', Li(2)
  write(*,*) 'Skin depth, sd = ', sd
  write(*,*) 'Mystery length (physical) , rMe = ', rM(3)
  write(*,*) 'Mystery length (simulation), rMe = ', rM(1)
  write(*,*) 'Mystery length, rMi = ', rM(2)
  write(*,*) 'dx = ', dxyz(1)
  write(*,*) 'dy = ', dxyz(2)
  write(*,*) 'dz = ', dxyz(3)
  write(*,*) 'CFL speed, dx/dt = ', CFL
  write(*,*) 'CFL freq., 1/dt = ', freq
  write(*,*) 'dt = ', dt
  write(*,*) 'Speed of light time step, dt_c = ', dt_c
  write(*,*) 'Speed of particles time step, dt_opt_part = ', dt_opt_part
  write(*,*) 'Speed of waves time step, dt_opt_wave = ', dt_opt_wave
  write(*,*) 'Optimal time step, dt_opt = ', dt_opt
  write(*,*) 'Recommended time step, dt_rec = ', dt_opt*0.75d0
  write(*,*) '****************************************************************'
  write(*,*) ' '
  write(*,*) 'WARNINGS!'
  
! Write the parameters to an m-file
      open(unit=1,file='save/parameters.dat')

!       write (1,fmt='(E10.4)') (Field(p,i,j,k), j=1, Ny)

      write(1,fmt='(a,E12.6,a)') 'n2p_1 = ', n2p(1)
      write(1,fmt='(a,E12.6,a)') 'n2p_2 = ', n2p(2)
      write(1,fmt='(a,E12.6,a)') 'n2p_3 = ', n2p(3)
      write(1,fmt='(a,E12.6,a)') 'n2p_4 = ', n2p(4)
      write(1,fmt='(a,E12.6,a)') 'n2p_5 = ', n2p(5)
      write(1,fmt='(a,E12.6,a)') 'n2p_6 = ', n2p(6)
      write(1,fmt='(a,E12.6,a)') 'n2p_7 = ', n2p(7)
      write(1,fmt='(a,E12.6,a)') 'n2p_8 = ', n2p(8)
      write(1,fmt='(a,E12.6,a)') 'mass_1 = ', mass(1)
      write(1,fmt='(a,E12.6,a)') 'mass_2 = ', mass(2)
      write(1,fmt='(a,E12.6,a)') 'mass_3 = ', mass(3)
      write(1,fmt='(a,E12.6,a)') 'mass_4 = ', mass(4)
      write(1,fmt='(a,E12.6,a)') 'mass_5 = ', mass(5)
      write(1,fmt='(a,E12.6,a)') 'mass_6 = ', mass(6)
      write(1,fmt='(a,E12.6,a)') 'mass_7 = ', mass(7)
      write(1,fmt='(a,E12.6,a)') 'mass_8 = ', mass(8)
      write(1,fmt='(a,E12.6,a)') 'charge_1 = ', charge(1)
      write(1,fmt='(a,E12.6,a)') 'charge_2 = ', charge(2)
      write(1,fmt='(a,E12.6,a)') 'charge_3 = ', charge(3)
      write(1,fmt='(a,E12.6,a)') 'charge_4 = ', charge(4)
      write(1,fmt='(a,E12.6,a)') 'charge_5 = ', charge(5)
      write(1,fmt='(a,E12.6,a)') 'charge_6 = ', charge(6)
      write(1,fmt='(a,E12.6,a)') 'charge_7 = ', charge(7)
      write(1,fmt='(a,E12.6,a)') 'charge_8 = ', charge(8)
      write(1,fmt='(a,E12.6,a)') 'density_1 = ', density(1)
      write(1,fmt='(a,E12.6,a)') 'density_2 = ', density(2)
      write(1,fmt='(a,E12.6,a)') 'density_3 = ', density(3)
      write(1,fmt='(a,E12.6,a)') 'density_4 = ', density(4)
      write(1,fmt='(a,E12.6,a)') 'density_5 = ', density(5)
      write(1,fmt='(a,E12.6,a)') 'density_6 = ', density(6)
      write(1,fmt='(a,E12.6,a)') 'density_7 = ', density(7)
      write(1,fmt='(a,E12.6,a)') 'density_8 = ', density(8)
      write(1,fmt='(a,E12.6,a)') 'kelvin_1 = ', kelvin(1)
      write(1,fmt='(a,E12.6,a)') 'kelvin_2 = ', kelvin(2)
      write(1,fmt='(a,E12.6,a)') 'kelvin_3 = ', kelvin(3)
      write(1,fmt='(a,E12.6,a)') 'kelvin_4 = ', kelvin(4)
      write(1,fmt='(a,E12.6,a)') 'kelvin_5 = ', kelvin(5)
      write(1,fmt='(a,E12.6,a)') 'kelvin_6 = ', kelvin(6)
      write(1,fmt='(a,E12.6,a)') 'kelvin_7 = ', kelvin(7)
      write(1,fmt='(a,E12.6,a)') 'kelvin_8 = ', kelvin(8)
      write(1,fmt='(a,E12.6,a)') 'Qn_1 = ', Qn(1)
      write(1,fmt='(a,E12.6,a)') 'Qn_2 = ', Qn(2)
      write(1,fmt='(a,E12.6,a)') 'Qn_3 = ', Qn(3)
      write(1,fmt='(a,E12.6,a)') 'Qn_4 = ', Qn(4)
      write(1,fmt='(a,E12.6,a)') 'Qn_5 = ', Qn(5)
      write(1,fmt='(a,E12.6,a)') 'Qn_6 = ', Qn(6)
      write(1,fmt='(a,E12.6,a)') 'Qn_7 = ', Qn(7)
      write(1,fmt='(a,E12.6,a)') 'Qn_8 = ', Qn(8)
      write(1,fmt='(a,E12.6,a)') 'vn_1 = ', vn(1)
      write(1,fmt='(a,E12.6,a)') 'vn_2 = ', vn(2)
      write(1,fmt='(a,E12.6,a)') 'vn_3 = ', vn(3)
      write(1,fmt='(a,E12.6,a)') 'vn_4 = ', vn(4)
      write(1,fmt='(a,E12.6,a)') 'vn_5 = ', vn(5)
      write(1,fmt='(a,E12.6,a)') 'vn_6 = ', vn(6)
      write(1,fmt='(a,E12.6,a)') 'vn_7 = ', vn(7)
      write(1,fmt='(a,E12.6,a)') 'vn_8 = ', vn(8)
      write(1,fmt='(a,E12.6,a)') 'Ek_1 = ', Ek(1)
      write(1,fmt='(a,E12.6,a)') 'Ek_2 = ', Ek(2)
      write(1,fmt='(a,E12.6,a)') 'Ek_3 = ', Ek(3)
      write(1,fmt='(a,E12.6,a)') 'Ek_4 = ', Ek(4)
      write(1,fmt='(a,E12.6,a)') 'Ek_5 = ', Ek(5)
      write(1,fmt='(a,E12.6,a)') 'Ek_6 = ', Ek(6)
      write(1,fmt='(a,E12.6,a)') 'Ek_7 = ', Ek(7)
      write(1,fmt='(a,E12.6,a)') 'Ek_8 = ', Ek(8)
      write(1,fmt='(a,E12.6,a)') 'nui_1 = ', nu_i(1)
      write(1,fmt='(a,E12.6,a)') 'nui_2 = ', nu_i(2)
      write(1,fmt='(a,E12.6,a)') 'nui_3 = ', nu_i(3)
      write(1,fmt='(a,E12.6,a)') 'nui_4 = ', nu_i(4)
      write(1,fmt='(a,E12.6,a)') 'nui_5 = ', nu_i(5)
      write(1,fmt='(a,E12.6,a)') 'nui_6 = ', nu_i(6)
      write(1,fmt='(a,E12.6,a)') 'nui_7 = ', nu_i(7)
      write(1,fmt='(a,E12.6,a)') 'nui_8 = ', nu_i(8)
      write(1,fmt='(a,E12.6,a)') 'nud_1 = ', nu_d(1)
      write(1,fmt='(a,E12.6,a)') 'nud_2 = ', nu_d(2)
      write(1,fmt='(a,E12.6,a)') 'nud_3 = ', nu_d(3)
      write(1,fmt='(a,E12.6,a)') 'nud_4 = ', nu_d(4)
      write(1,fmt='(a,E12.6,a)') 'nud_5 = ', nu_d(5)
      write(1,fmt='(a,E12.6,a)') 'nud_6 = ', nu_d(6)
      write(1,fmt='(a,E12.6,a)') 'nud_7 = ', nu_d(7)
      write(1,fmt='(a,E12.6,a)') 'nud_8 = ', nu_d(8)
      write(1,fmt='(a,I4,a)') 'ppc_1 = ', ppc(1)
      write(1,fmt='(a,I4,a)') 'ppc_2 = ', ppc(2)
      write(1,fmt='(a,I4,a)') 'ppc_3 = ', ppc(3)
      write(1,fmt='(a,I4,a)') 'ppc_4 = ', ppc(4)
      write(1,fmt='(a,I4,a)') 'ppc_5 = ', ppc(5)
      write(1,fmt='(a,I4,a)') 'ppc_6 = ', ppc(6)
      write(1,fmt='(a,I4,a)') 'ppc_7 = ', ppc(7)
      write(1,fmt='(a,I4,a)') 'ppc_8 = ', ppc(8)
      write(1,*) 'gyration_1 = ', gyration(1)
      write(1,*) 'gyration_2 = ', gyration(2)
      write(1,*) 'gyration_3 = ', gyration(3)
      write(1,*) 'gyration_4 = ', gyration(4)
      write(1,*) 'gyration_5 = ', gyration(5)
      write(1,*) 'gyration_6 = ', gyration(6)
      write(1,*) 'gyration_7 = ', gyration(7)
      write(1,*) 'gyration_8 = ', gyration(8)
      write(1,fmt='(a,E12.6,a)') 'v0_1 = ', ve(1)
      write(1,fmt='(a,E12.6,a)') 'v0_2 = ', ve(2)
      write(1,fmt='(a,E12.6,a)') 'v0_3 = ', ve(3)
      write(1,fmt='(a,E12.6,a)') 'B0_1 = ', B0(1)
      write(1,fmt='(a,E12.6,a)') 'B0_2 = ', B0(2)
      write(1,fmt='(a,E12.6,a)') 'B0_3 = ', B0(3)
      write(1,fmt='(a,E12.6,a)') 'dxyz_1 = ', dxyz(1)
      write(1,fmt='(a,E12.6,a)') 'dxyz_2 = ', dxyz(2)
      write(1,fmt='(a,E12.6,a)') 'dxyz_3 = ', dxyz(3)
      write(1,fmt='(a,E12.6,a)') 'xmin = ', xmin
      write(1,fmt='(a,E12.6,a)') 'xmax = ', xmax
      write(1,fmt='(a,E12.6,a)') 'ymin = ', ymin
      write(1,fmt='(a,E12.6,a)') 'ymax = ', ymax
      write(1,fmt='(a,E12.6,a)') 'zmin = ', zmin
      write(1,fmt='(a,E12.6,a)') 'zmax = ', zmax
      write(1,fmt='(a,E12.6,a)') 'dt = ', dt
      write(1,fmt='(a,I4,a)') 'Nx = ', Nx
      write(1,fmt='(a,I4,a)') 'Ny = ', Ny
      write(1,fmt='(a,I4,a)') 'Nz = ', Nz
      write(1,fmt='(a,I4,a)') 'iprocs = ', iprocs
      write(1,fmt='(a,I4,a)') 'jprocs = ', jprocs
      write(1,fmt='(a,I4,a)') 'kprocs = ', kprocs
      write(1,*) 'xflow = ', xflow
      write(1,*) 'yflow = ', yflow
      write(1,*) 'zflow = ', zflow

      close(1)

! write all warnings to a file too
     open(unit=1,file='save/warnings.dat')

! Warn for plasma freq.
  if (wp(1) .GT. freq) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'wp_e > 1/dt'
     write(*,*) wp(1), freq
     write(*,*) '****************************************************************'
     write(1,*) ' '   
     write(1,*) '****************************************************************'
     write(1,*) 'wp_e > 1/dt'
     write(1,*) wp(1), freq
     write(1,*) '****************************************************************'
  end if

  if (wp(2) .GT. freq) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'wp_i > 1/dt'
     write(*,*) wp(2), freq
     write(*,*) '****************************************************************'
     write(1,*) ' '   
     write(1,*) '****************************************************************'
     write(1,*) 'wp_i > 1/dt'
     write(1,*) wp(2), freq
     write(1,*) '****************************************************************'
  end if

! Warn for cyclotron freq.
  if (wc(1) .GT. freq) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'wc_e > 1/dt'
     write(*,*) wc(1), freq
     write(*,*) '****************************************************************'
     write(1,*) ' '   
     write(1,*) '*******************************************'
     write(1,*) 'wc_e > 1/dt'
     write(1,*) wc(1), freq
     write(1,*) '*******************************************'
  end if

  if (wc(2) .GT. freq) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'wc_i > 1/dt'
     write(*,*) wc(2), freq
     write(*,*) '****************************************************************'
     write(1,*) ' '   
     write(1,*) '*******************************************'
     write(1,*) 'wc_i > 1/dt'
     write(1,*) wc(2), freq
     write(1,*) '*******************************************'
  end if

! Warn for Alfven speed
  if (va .GT. CFL) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'va > dx/dt'
     write(*,*) va, CFL
     write(*,*) '****************************************************************'
     write(1,*) ' '   
     write(1,*) '*******************************************'
     write(1,*) 'va > dx/dt'
     write(1,*) va, CFL
     write(1,*) '*******************************************'
  end if


! Warn for courant condition
  if (v0(4,1)+vth(1) .GT. CFL ) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'v_e > dx/dt'
     write(*,*) v0(4,1), vth(1), CFL
     write(*,*) '****************************************************************'
     write(1,*) ' '   
     write(1,*) '*******************************************'
     write(1,*) 'v_e > dx/dt'
     write(1,*) v0(4,1), vth(1), CFL
     write(1,*) '*******************************************'
  end if

  if (v0(4,2)+vth(2) .GT. CFL ) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'v_i > dx/dt'
     write(*,*) v0(4,2), vth(2), CFL
     write(*,*) '****************************************************************'
     write(1,*) ' '   
     write(1,*) '*******************************************'
     write(1,*) 'v_i > dx/dt'
     write(1,*) v0(4,2), vth(2), CFL
     write(1,*) '*******************************************'
  end if

! Warn for Ld_e
  if (Ld(1) .LT. dxyz(1)) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'Ld_e < dx'
     write(*,*) Ld(1), dxyz(1)
     write(*,*) '****************************************************************'
     write(1,*) ' '   
     write(1,*) '****************************************************************'
     write(1,*) 'Ld_e < dx'
     write(1,*) Ld(1), dxyz(1)
     write(1,*) '****************************************************************'
  end if


! Warn for Ld_i
  if (Ld(2) .LT. dxyz(1)) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'Ld_i < dx'
     write(*,*) Ld(2), dxyz(1)
     write(*,*) '****************************************************************'
     write(1,*) ' '   
     write(1,*) '****************************************************************'
     write(1,*) 'Ld_i < dx'
     write(1,*) Ld(2), dxyz(1)
     write(1,*) '****************************************************************'
  end if
  

! Warn for gyroradius
  if (rL(1) .LT. dxyz(1)) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'rL_e < dx'
     write(*,*) rL(1), dxyz(1)
     write(*,*) '****************************************************************'
     write(1,*) ' '   
     write(1,*) '****************************************************************'
     write(1,*) 'rL_e < dx'
     write(1,*) rL(1), dxyz(1)
     write(1,*) '****************************************************************'
  end if

  if (rL(2) .LT. dxyz(1)) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'rL_i < dx'
     write(*,*) rL(2), dxyz(1)
     write(*,*) '****************************************************************'
     write(1,*) ' '   
     write(1,*) '****************************************************************'
     write(1,*) 'rL_i < dx'
     write(1,*) rL(2), dxyz(1)
     write(1,*) '****************************************************************'
  end if

  if (rPU(1) .LT. dxyz(1)) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'rPU_e < dx'
     write(*,*) rPU(1), dxyz(1)
     write(*,*) '****************************************************************'
     write(1,*) ' '   
     write(1,*) '****************************************************************'
     write(1,*) 'rPU_e < dx'
     write(1,*) rPU(1), dxyz(1)
     write(1,*) '****************************************************************'
  end if

  if (rPU(2) .LT. dxyz(1)) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'rPU_i < dx'
     write(*,*) rPU(2), dxyz(1)
     write(*,*) '****************************************************************'
     write(1,*) ' '   
     write(1,*) '****************************************************************'
     write(1,*) 'rPU_i < dx'
     write(1,*) rPU(2), dxyz(1)
     write(1,*) '****************************************************************'
  end if

! Warn for inertial length
  if (Li(1) .LT. dxyz(1)) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'Li_e < dx'
     write(*,*) Li(1), dxyz(1)
     write(*,*) '****************************************************************'
     write(1,*) ' '   
     write(1,*) '****************************************************************'
     write(1,*) 'Li_e < dx'
     write(1,*) Li(1), dxyz(1)
     write(1,*) '****************************************************************'
  end if

  if (Li(2) .LT. dxyz(1)) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'Li_i < dx'
     write(*,*) Li(2), dxyz(1)
     write(*,*) '****************************************************************'
     write(1,*) ' '   
     write(1,*) '****************************************************************'
     write(1,*) 'Li_i < dx'
     write(1,*) Li(2), dxyz(1)
     write(1,*) '****************************************************************'
  end if

! Warn if dx NE dy NE dz
  if (dxyz(1) .NE. dxyz(2) .OR. dxyz(2) .NE. dxyz(3) ) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'Not a regular grid!'
     write(*,*) 'dx=', dxyz(1), ' dy=', dxyz(2), ' dz=', dxyz(3)
     write(*,*) '****************************************************************'
     write(1,*) ' '   
     write(1,*) '****************************************************************'
     write(1,*) 'Not a regular grid!'
     write(1,*) 'dx=', dxyz(1), ' dy=', dxyz(2), ' dz=', dxyz(3)
     write(1,*) '****************************************************************'
  end if

! Warn if you create particles that should be ionized from Haser
  if ( sum(Qn*dble(ppc)) .gt. 0.0d0 ) then
     write(*,*) ' '   
     write(*,*) '****************************************************************'
     write(*,*) 'You are ionizing incorrectly!'
     write(*,*) 'Qn=', Qn, ' ppc', ppc
     write(*,*) '****************************************************************' 
     write(1,*) ' '   
     write(1,*) '****************************************************************'
     write(1,*) 'You are ionizing incorrectly!'
     write(1,*) 'Qn=', Qn, ' ppc', ppc
     write(1,*) '****************************************************************'
  end if

  close(1)


  write(*,*) ' '   
  

  return
end subroutine WriteWarnings

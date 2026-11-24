module axisymmetric_plume
  use axisymmetric_kinetic
  use iso_fortran_env, only: output_unit
  use ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: run_axisymmetric

  real(8), parameter :: pi8 = acos(-1.d0), temperature_ref = 300.d0, density_ref = 1.d-3
  integer :: nx = 100, nr = 50, nu = 40, ns = 12, na = 12, log_stride = 100
  real(8) :: length = 10.d0, radius = 5.d0, x_front = 0.8d0, x_aft = 1.8d0, r_body = 0.5d0
  real(8) :: gas_r = 403.5744316182839d0, mu_at_300 = 4.05d-5, omega = 1.237507354328596d0
  real(8) :: u_scale = 1.75d0, u_center = 0.8d0, s_max = 10.d0, cfl = 0.45d0
  real(8) :: end_time = 0.006d0, average_start = 0.003d0, rho_initial = 0.d0
  real(8) :: cold_wall = 300.d0, aft_wall = 600.d0
  logical :: uniform_test = .false.
  character(512) :: inlet_file = 'inlet_panels.dat'
  character(512) :: species_file = ''
  real(8) :: vhs_reference_temperature = 273.15d0
  integer :: spatial_order = 1, slope_limiter = 1, x_ref_split(2) = 1, r_ref_split(2) = 1
  ! An open x=length boundary truncates the domain for near-nozzle studies.
  logical :: open_far_boundary = .false.
  real(8) :: x_ref_lo(2) = 0.d0, x_ref_hi(2) = 0.d0, r_ref_lo(2) = 0.d0, r_ref_hi(2) = 0.d0
  real(8) :: timestep_limit_s = 0.d0, translational_prandtl = 1.d0
  logical :: axial_half_range = .false., transport_matching = .false., enable_vibration = .false.
  integer, parameter :: mode_species_limit = 64, vibration_band_limit = 8
  real(8) :: rotation_tstar(mode_species_limit) = 0.d0, vibration_probability(mode_species_limit) = 0.d0
  real(8) :: vibration_bands(vibration_band_limit,mode_species_limit) = 0.d0
  namelist /axisym/ nx, nr, nu, ns, na, length, radius, x_front, x_aft, r_body, gas_r, mu_at_300, omega, &
                    u_scale, u_center, s_max, cfl, end_time, average_start, rho_initial, cold_wall, aft_wall, &
                    uniform_test, inlet_file, log_stride, species_file, vhs_reference_temperature, &
                    spatial_order, slope_limiter, open_far_boundary, x_ref_split, x_ref_lo, x_ref_hi, r_ref_split, r_ref_lo, r_ref_hi, timestep_limit_s, &
                    axial_half_range, transport_matching, translational_prandtl, enable_vibration, rotation_tstar, &
                    vibration_probability, vibration_bands

  type(velocity_ring) :: q
  logical, allocatable :: fluid(:,:)
  real(8), allocatable :: f(:,:,:,:), fx(:,:,:,:), fr(:,:,:,:), profile(:,:,:), average(:,:,:), macro(:,:,:)
  real(8), allocatable :: cold(:,:), hot(:,:), wall_average(:,:), wall_now(:,:)
  real(8) :: dx, dr, c0, mu, dt

  interface
    module subroutine run_mixture()
    end subroutine run_mixture
  end interface

contains

  subroutine run_axisymmetric(config_file)
    character(*), intent(in) :: config_file
    real(8) :: total_mass, initial_mass, cumulative_flux, mass_change, mass_error, dt_limit
    real(8) :: t, t_end, t_average, sample_time, overlap, fallback_mass, peak_fallback
    integer :: unit, info, i, j, step, i_front, i_aft, j_body
    character(512) :: message

    open(newunit=unit,file=config_file,status='old',action='read',iostat=info,iomsg=message)
    if (info /= 0) then
      write(*,*) trim(message)
      error stop 'cannot open axisymmetric namelist'
    end if
    read(unit,nml=axisym,iostat=info,iomsg=message)
    close(unit)
    if (info /= 0) then
      write(*,*) trim(message)
      error stop 'invalid axisymmetric namelist'
    end if
    if (min(nx,nr,log_stride) < 1 .or. min(nu,ns,na) < 4 .or. mod(na,2) /= 0) error stop 'invalid grid'
    if (min(length,radius,gas_r,mu_at_300,end_time,u_scale,s_max,cold_wall,aft_wall) <= 0.d0) &
      error stop 'invalid axisymmetric physical parameters'
    if (cfl <= 0.d0 .or. cfl > 0.85d0 .or. rho_initial < 0.d0 .or. average_start < 0.d0) &
      error stop 'invalid axisymmetric initialization or CFL'
    c0 = sqrt(2.d0*gas_r*temperature_ref)
    mu = mu_at_300/(density_ref*c0)
    dx = length/nx
    dr = radius/nr
    if (len_trim(species_file) > 0) then
      call run_mixture()
      return
    end if
    allocate(fluid(nx,nr))
    fluid = .true.
    if (.not. uniform_test) then
      i_front = nint(x_front/dx)
      i_aft = nint(x_aft/dx)
      j_body = nint(r_body/dr)
      if (i_front < 1 .or. i_front >= i_aft .or. i_aft >= nx .or. j_body < 1 .or. j_body >= nr) &
        error stop 'lander outside the axisymmetric grid'
      if (max(abs(x_front-i_front*dx),abs(x_aft-i_aft*dx),abs(r_body-j_body*dr)) > 1.d-9) &
        error stop 'choose an axisymmetric mesh with the lander walls on faces'
      fluid(i_front+1:i_aft,1:j_body) = .false.
    end if

    call make_velocity_ring(q, nu, ns, na, u_scale, u_center, s_max)
    allocate(f(nu,q%np,nx,nr), fx(nu,q%np,0:nx,nr), fr(nu,q%np,nx,0:nr))
    allocate(profile(nu,q%np,nr), average(4,nx,nr), macro(4,nx,nr))
    allocate(cold(nu,q%np), hot(nu,q%np), wall_average(2,nr), wall_now(2,nr))
    call prescribed(rho=1.d0, ux=0.d0, ur=0.d0, temperature=cold_wall/temperature_ref, result=cold)
    call prescribed(rho=1.d0, ux=0.d0, ur=0.d0, temperature=aft_wall/temperature_ref, result=hot)
    profile = 0.d0
    if (.not. uniform_test) call read_inlet()

    !$omp parallel do collapse(2)
    do j = 1, nr
      do i = 1, nx
        f(:,:,i,j) = 0.d0
        if (fluid(i,j)) f(:,:,i,j) = (rho_initial/density_ref)*cold
        macro(:,i,j) = ring_moments(q,f(:,:,i,j))
      end do
    end do
    !$omp end parallel do
    average = 0.d0
    wall_average = 0.d0
    t = 0.d0
    t_end = end_time*c0
    t_average = average_start*c0
    sample_time = 0.d0
    step = 0
    initial_mass = domain_mass()
    cumulative_flux = 0.d0
    peak_fallback = 0.d0
    dt_limit = cfl / (maxval(abs(q%u))/dx + 2.d0*maxval(abs(q%v))/dr + 2.d0*s_max/(dr*q%angle_width))
    write(*,'(a,4i6,a,es14.6)') 'Axisymmetric UGKS (x,r,u,s,theta), nx nr nu nperp: ', nx,nr,nu,q%np, &
      '  dt(s): ',dt_limit/c0
    write(*,'(a)') 'step time_s mass_kg mass_balance_kg fallbacks_kg'
    flush(output_unit)

    do while (t < t_end)
      dt = min(dt_limit, t_end-t)
      step = step+1
      call advance(fallback_mass)
      mass_change = dt * boundary_mass_flow()
      cumulative_flux = cumulative_flux + mass_change
      peak_fallback = max(peak_fallback,fallback_mass)
      overlap = max(0.d0, t+dt-max(t,t_average))
      if (overlap > 0.d0) then
        average = average + overlap*macro
        wall_average = wall_average + overlap*wall_now
        sample_time = sample_time + overlap
      end if
      t = t+dt
      if (mod(step,log_stride) == 0 .or. t >= t_end) then
        total_mass = domain_mass()
        mass_error = total_mass-initial_mass-cumulative_flux
        if (.not. all(ieee_is_finite(macro))) error stop 'nonfinite axisymmetric solution'
        write(*,'(i7,4es17.8)') step,t/c0,total_mass,mass_error,fallback_mass
        flush(output_unit)
      end if
    end do
    call write_macro('axisymmetric_final.dat',macro)
    if (sample_time > 0.d0) then
      call write_macro('axisymmetric_average.dat',average/sample_time)
      call write_wall(wall_average/sample_time)
    else
      call write_wall(wall_now)
    end if
    write(*,'(a,es16.7,a,es16.7)') 'Completed axisymmetric run. Averaging seconds: ', sample_time/c0, &
      '   Peak mass using collision fallback, kg: ',peak_fallback
  end subroutine run_axisymmetric

  subroutine prescribed(rho, ux, ur, temperature, result)
    real(8), intent(in) :: rho, ux, ur, temperature
    real(8), intent(out) :: result(q%nu,q%np)
    logical :: fits

    call ring_equilibrium(q,[rho,rho*ux,rho*ur,rho*(ux**2+ur**2+1.5d0*temperature)],result,fits)
    if (.not. fits) error stop 'inlet or wall Maxwellian needs a better resolved velocity grid'
  end subroutine prescribed

  subroutine read_inlet()
    integer :: unit, io, j
    real(8) :: panel(q%nu,q%np), lo, hi, a, b, rho, temperature, ux, ur, weighted
    real(8) :: massflow
    character(512) :: header

    open(newunit=unit,file=inlet_file,status='old',action='read')
    read(unit,'(a)') header
    do
      read(unit,*,iostat=io) lo,hi,rho,temperature,ux,ur
      if (io < 0) exit
      if (io /= 0 .or. lo < 0.d0 .or. hi > r_body .or. hi <= lo .or. min(rho,temperature) <= 0.d0) &
        error stop 'invalid axisymmetric inlet panel'
      call prescribed(rho/density_ref,ux/c0,ur/c0,temperature/temperature_ref,panel)
      do j = 1, nr
        a = max(lo,(j-1)*dr)
        b = min(hi,j*dr)
        if (b <= a) cycle
        weighted = (b*b-a*a)/((2*j-1)*dr*dr)
        profile(:,:,j) = profile(:,:,j) + weighted*panel
      end do
    end do
    close(unit)
    massflow = 0.d0
    do j = 1, nr
      massflow = massflow + pi8*(2*j-1)*dr*dr*sum(matmul(max(q%u,0.d0),profile(:,:,j)))
    end do
    write(*,'(a,es16.7)') 'Prescribed positive nozzle mass flux, kg/s: ', massflow*density_ref*c0
  end subroutine read_inlet

  subroutine advance(fallback_mass)
    real(8), intent(out) :: fallback_mass
    integer :: i,j
    real(8) :: r

    !$omp parallel do collapse(2)
    do j = 1, nr
      do i = 0, nx
        call x_face(i,j)
      end do
    end do
    !$omp end parallel do
    !$omp parallel do collapse(2)
    do j = 1, nr
      do i = 1, nx
        call r_face(i,j)
      end do
    end do
    !$omp end parallel do

    fallback_mass = 0.d0
    !$omp parallel do collapse(2) private(r) reduction(+:fallback_mass)
    do j = 1, nr
      do i = 1, nx
        if (.not. fluid(i,j)) cycle
        r = (j-0.5d0)*dr
        call update_cell(i,j,r,fallback_mass)
      end do
    end do
    !$omp end parallel do
    do j = 1, nr
      wall_now(1,j) = density_ref*c0**2*sum(matmul(q%u,fx(:,:,nx,j)))
      wall_now(2,j) = 0.5d0*density_ref*c0**3*(sum(matmul(q%u2,fx(:,:,nx,j)))+ &
                       sum(matmul(fx(:,:,nx,j),q%s2)))
    end do
  end subroutine advance

  subroutine x_face(i,j)
    integer, intent(in) :: i,j
    integer :: k
    logical :: left, right

    if (i == 0) then
      if (uniform_test) then
        call axial_wall(f(:,:,1,j),cold,1,fx(:,:,i,j))
      else
        do k = 1, q%np
          fx(:,k,i,j) = min(0.d0,q%u)*f(:,k,1,j)
        end do
      end if
    else if (i == nx) then
      if (open_far_boundary .and. .not. uniform_test) then
        do k = 1, q%np
          fx(:,k,i,j) = max(0.d0,q%u)*f(:,k,nx,j)
        end do
      else
        call axial_wall(f(:,:,nx,j),cold,-1,fx(:,:,i,j))
      end if
    else
      left = fluid(i,j)
      right = fluid(i+1,j)
      if (left .and. right) then
        call internal_flux(f(:,:,i,j),f(:,:,i+1,j),.true.,fx(:,:,i,j))
      else if (left) then
        call axial_wall(f(:,:,i,j),cold,-1,fx(:,:,i,j))
      else if (right) then
        call axial_wall(f(:,:,i+1,j),hot,1,fx(:,:,i,j))
        do k = 1, q%np
          fx(:,k,i,j) = fx(:,k,i,j) + max(0.d0,q%u)*profile(:,k,j)
        end do
      else
        fx(:,:,i,j) = 0.d0
      end if
    end if
  end subroutine x_face

  subroutine r_face(i,j)
    integer, intent(in) :: i,j
    integer :: k
    logical :: below, above

    if (j == nr) then
      if (uniform_test) then
        call radial_wall(f(:,:,i,j),fr(:,:,i,j))
      else
        do k = 1, q%np
          fr(:,k,i,j) = max(0.d0,q%v(k))*f(:,k,i,j)
        end do
      end if
    else
      below = fluid(i,j)
      above = fluid(i,j+1)
      if (below .and. above) then
        call internal_flux(f(:,:,i,j),f(:,:,i,j+1),.false.,fr(:,:,i,j))
      else if (above) then
        call radial_wall(f(:,:,i,j+1),fr(:,:,i,j),1)
      else
        fr(:,:,i,j) = 0.d0
      end if
    end if
  end subroutine r_face

  subroutine internal_flux(left,right,axial,flux)
    real(8), intent(in) :: left(q%nu,q%np), right(q%nu,q%np)
    logical, intent(in) :: axial
    real(8), intent(out) :: flux(q%nu,q%np)
    real(8) :: upwind(q%nu,q%np), eq(q%nu,q%np), w(4), collisional, cap
    logical :: fits
    integer :: k

    if (axial) then
      do k = 1, q%np
        where (q%u > 0.d0)
          upwind(:,k) = left(:,k)
        elsewhere
          upwind(:,k) = right(:,k)
        end where
      end do
    else
      do k = 1, q%np
        if (q%v(k) > 0.d0) then
          upwind(:,k) = left(:,k)
        else
          upwind(:,k) = right(:,k)
        end if
      end do
    end if
    w = ring_moments(q,upwind)
    collisional = 1.d0-flight_fraction(relaxation_rate(w))
    if (collisional > 1.d-6) then
      call ring_equilibrium(q,w,eq,fits)
      if (fits) then
        ! The conservative interface distribution is the time integral of the
        ! BGK characteristic with frozen, piecewise-constant states. Limit the
        ! flux only where its correction would exhaust an outgoing bin.
        cap = 0.95d0/cfl
        upwind = min(upwind+collisional*(eq-upwind),cap*upwind)
      end if
    end if
    if (axial) then
      do k = 1, q%np
        flux(:,k) = q%u*upwind(:,k)
      end do
    else
      do k = 1, q%np
        flux(:,k) = q%v(k)*upwind(:,k)
      end do
    end if
  end subroutine internal_flux

  subroutine axial_wall(cell,equilibrium,normal,flux)
    real(8), intent(in) :: cell(q%nu,q%np),equilibrium(q%nu,q%np)
    integer, intent(in) :: normal
    real(8), intent(out) :: flux(q%nu,q%np)
    real(8) :: coef
    integer :: k

    coef = -sum(matmul(min(normal*q%u,0.d0),cell)) / sum(matmul(max(normal*q%u,0.d0),equilibrium))
    do k = 1, q%np
      where (normal*q%u < 0.d0)
        flux(:,k) = q%u*cell(:,k)
      elsewhere
        flux(:,k) = coef*q%u*equilibrium(:,k)
      end where
    end do
  end subroutine axial_wall

  subroutine radial_wall(cell,flux,normal)
    real(8), intent(in) :: cell(q%nu,q%np)
    real(8), intent(out) :: flux(q%nu,q%np)
    integer, optional, intent(in) :: normal
    integer :: sign, k
    real(8) :: coef

    sign = -1
    if (present(normal)) sign = normal
    coef = -sum(matmul(cell,min(sign*q%v,0.d0))) / sum(matmul(cold,max(sign*q%v,0.d0)))
    do k = 1, q%np
      if (sign*q%v(k) < 0.d0) then
        flux(:,k) = q%v(k)*cell(:,k)
      else
        flux(:,k) = coef*q%v(k)*cold(:,k)
      end if
    end do
  end subroutine radial_wall

  subroutine update_cell(i,j,r,fallback)
    integer, intent(in) :: i,j
    real(8), intent(in) :: r
    real(8), intent(inout) :: fallback
    real(8) :: transported(q%nu,q%np), eq(q%nu,q%np), angular(q%nu), lower(q%nu), upper(q%nu)
    real(8) :: w(4), rate, radius_lo, radius_hi, factor
    integer :: k
    logical :: fits

    radius_lo = (j-1)*dr
    radius_hi = j*dr
    do k = 1, q%np
      factor = -sqrt(q%s2(k))/(r*q%angle_width)
      lower = factor*q%sin_lower(k)*f(:,k,i,j)
      upper = 0.d0
      if (mod(k,q%na) /= 0) upper = factor*q%sin_upper(k)*f(:,k+1,i,j)
      angular = lower-upper
      transported(:,k) = f(:,k,i,j) + dt*((fx(:,k,i-1,j)-fx(:,k,i,j))/dx + angular)
      transported(:,k) = transported(:,k) - (dt*radius_hi/(r*dr))*fr(:,k,i,j)
      if (j > 1) transported(:,k) = transported(:,k) + (dt*radius_lo/(r*dr))*fr(:,k,i,j-1)
    end do
    if (minval(transported) < -1.d-13*maxval(abs(transported))) error stop 'negative axisymmetric distribution'
    transported = max(transported,0.d0)
    w = ring_moments(q,transported)
    rate = relaxation_rate(w)
    if (rate > 1.d-7) then
      call ring_equilibrium(q,w,eq,fits)
      if (fits) then
        transported = (transported+rate*eq)/(1.d0+rate)
      else
        fallback = fallback + w(1)*(rate/(1.d0+rate))*density_ref*2.d0*pi8*r*dx*dr
      end if
    end if
    f(:,:,i,j) = transported
    macro(:,i,j) = w
  end subroutine update_cell

  real(8) function relaxation_rate(w) result(rate)
    real(8), intent(in) :: w(4)
    real(8) :: thermal

    rate = 0.d0
    if (w(1) <= 1.d-20) return
    thermal = max((2.d0/3.d0)*(w(4)/w(1)-sum((w(2:3)/w(1))**2)),1.d-4)
    rate = (dt/(2.d0*mu))*w(1)*thermal**(1.d0-omega)
  end function relaxation_rate

  real(8) function domain_mass() result(mass)
    integer :: j

    mass = 0.d0
    do j = 1, nr
      mass = mass + (j-0.5d0)*dr*sum(macro(1,:,j))
    end do
    mass = mass * density_ref * 2.d0*pi8*dx*dr
  end function domain_mass

  real(8) function boundary_mass_flow() result(flow)
    integer :: i, j
    real(8) :: rd

    flow = -radius*dx*sum(fr(:,:,1:nx,nr))
    do j = 1, nr
      rd = (j-0.5d0)*dr*dr
      flow = flow + rd*(sum(fx(:,:,0,j))-sum(fx(:,:,nx,j)))
      do i = 1, nx-1
        if (fluid(i,j) .neqv. fluid(i+1,j)) then
          if (fluid(i+1,j)) then
            flow = flow + rd*sum(fx(:,:,i,j))
          else
            flow = flow - rd*sum(fx(:,:,i,j))
          end if
        end if
      end do
    end do
    flow = flow * density_ref * 2.d0*pi8
  end function boundary_mass_flow

  subroutine write_macro(filename,w)
    character(*), intent(in) :: filename
    real(8), intent(in) :: w(4,nx,nr)
    real(8) :: mass, u, v, temp
    integer :: unit, i,j

    open(newunit=unit,file=filename,status='replace',action='write')
    write(unit,'(a)') '# x_m r_m rho_kg_m3 ux_m_s ur_m_s T_K p_Pa fluid'
    do j = 1, nr
      do i = 1, nx
        mass = w(1,i,j)
        u = 0.d0
        v = 0.d0
        temp = 0.d0
        if (mass > 1.d-30) then
          u = w(2,i,j)/mass
          v = w(3,i,j)/mass
          temp = (2.d0*temperature_ref/3.d0)*max(0.d0, w(4,i,j)/mass-u*u-v*v)
        end if
        mass = density_ref*mass
        write(unit,'(7(es23.15,1x),i1)') (i-0.5d0)*dx, (j-0.5d0)*dr, mass, u*c0, v*c0, temp, &
          gas_r*mass*temp, merge(1,0,fluid(i,j))
      end do
    end do
    close(unit)
  end subroutine write_macro

  subroutine write_wall(wall)
    real(8), intent(in) :: wall(2,nr)
    integer :: unit, j

    open(newunit=unit,file='axisymmetric_wall.dat',status='replace',action='write')
    write(unit,'(a)') '# r_m lunar_wall_normal_pressure_Pa lunar_wall_translational_heat_flux_W_m2'
    do j = 1, nr
      write(unit,'(3(es23.15,1x))') (j-0.5d0)*dr,wall(:,j)
    end do
    close(unit)
  end subroutine write_wall
end module axisymmetric_plume

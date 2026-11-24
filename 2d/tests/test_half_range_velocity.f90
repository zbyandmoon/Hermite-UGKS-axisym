program test_half_range_velocity
  use axisymmetric_kinetic, only: velocity_ring, make_velocity_ring_half_range
  implicit none

  integer, parameter :: orders(*) = [2, 4, 14, 28, 40, 64, 128]
  real(8), parameter :: pi = acos(-1.d0), u_scale = 1.6d0
  type(velocity_ring) :: grid
  real(8) :: error
  integer :: i

  do i = 1, size(orders)
    call check_half_moments(orders(i), error)
    write(*,'(a,i3,a,es12.4)') 'Axial nodes ', orders(i), ': maximum relative moment error ', error
  end do
  call check_wall_fluxes(1.6d0, 1.d-7)
  call check_wall_fluxes(2.2d0, 1.d-5)
  write(*,*) 'PASS: mirrored half-range Gauss-Hermite velocity quadrature'

contains

  subroutine check_half_moments(nu, max_error)
    integer, intent(in) :: nu
    real(8), intent(out) :: max_error
    real(8) :: x(nu/2), log_weight(nu/2), log_expected, relative_moment
    integer :: nh, k

    call make_velocity_ring_half_range(grid, nu, 8, 8, u_scale, 8.d0)
    nh = nu/2
    call require(grid%nu == nu .and. grid%np == 64 .and. grid%na == 8, 'velocity ring layout')
    x = grid%u(nh+1:)/u_scale
    call require(all(x > 0.d0) .and. all(x(2:) > x(:nh-1)), 'positive roots are ordered')
    call require(maxval(abs(grid%u(:nh)+grid%u(nu:nh+1:-1))) < 1.d-14, 'nodes are exact reflections')
    call require(maxval(abs(grid%log_wu(:nh)-grid%log_wu(nu:nh+1:-1))) < 1.d-14, 'weights are reflected')
    log_weight = grid%log_wu(nh+1:)-log(u_scale)-x**2
    max_error = 0.d0
    do k = 0, 2*nh-1
      ! The kth half-Gaussian moment is 0.5*Gamma((k+1)/2).
      log_expected = log(0.5d0)+log_gamma(0.5d0*(k+1))
      relative_moment = sum(exp(log_weight+k*log(x)-log_expected))
      call require(abs(relative_moment-1.d0) < 1.5d-12, 'half-Gaussian polynomial moment')
      max_error = max(max_error, abs(relative_moment-1.d0))
    end do
  end subroutine check_half_moments

  subroutine check_wall_fluxes(axial_scale, tolerance)
    real(8), intent(in) :: axial_scale, tolerance
    integer, parameter :: nu = 28
    real(8), parameter :: gas_r = 403.5744316182838d0
    real(8), parameter :: molar_mass(*) = [28.013d0, 18.016d0, 44.011d0, 28.011d0, 2.016d0]
    integer, parameter :: rot_dof(*) = [2, 3, 2, 2, 2]
    real(8) :: u(nu/2), f(nu/2), scale, theta, internal_energy, flux(3), exact(3), max_error(3)
    integer :: s, wall

    ! Integrate transverse and rotational factors analytically to isolate axial error.
    call make_velocity_ring_half_range(grid, nu, 8, 8, axial_scale, 8.d0)
    max_error = 0.d0
    do s = 1, size(molar_mass)
      scale = sqrt(1.380649d-23*6.02214076d23*1000.d0/(gas_r*molar_mass(s)))
      u = scale*grid%u(nu/2+1:)
      do wall = 1, 2
        theta = wall*scale**2
        internal_energy = 0.25d0*rot_dof(s)*theta
        f = exp(grid%log_wu(nu/2+1:)+log(scale)-u**2/theta)/sqrt(pi*theta)
        flux(1) = sum(u*f)
        flux(2) = sum(u**2*f)
        flux(3) = sum(u*(0.5d0*(u**2+theta)+internal_energy)*f)
        exact(1) = 0.5d0*sqrt(theta/pi)
        exact(2) = 0.25d0*theta
        exact(3) = (theta+internal_energy)*exact(1)
        call require(maxval(abs(flux/exact-1.d0)) < tolerance, 'stationary species wall half-flux')
        max_error = max(max_error, abs(flux/exact-1.d0))
      end do
    end do
    write(*,'(a,f3.1,a,3es12.4)') 'Axial scale ', axial_scale, ': relative wall mass, pressure and energy errors ', max_error
  end subroutine check_wall_fluxes

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (condition) return
    write(*,*) 'FAIL: ', message
    error stop 1
  end subroutine require
end program test_half_range_velocity

module axisymmetric_kinetic
  use config, only: nv, Tref
  use alldata, only: uspace, weight
  use velocity, only: init_velocity_gauss_hermite_2d
  implicit none
  private
  public :: velocity_ring, make_velocity_ring, make_velocity_ring_half_range
  public :: ring_moments, ring_equilibrium, flight_fraction

  real(8), parameter :: pi8 = acos(-1.d0)

  type velocity_ring
    integer :: nu, na, np
    real(8) :: angle_width
    real(8), allocatable :: u(:), u2(:), log_wu(:)
    real(8), allocatable :: v(:), s2(:), sin_lower(:), sin_upper(:), log_wp(:)
  end type velocity_ring

contains

  subroutine make_velocity_ring(q, nu, ns, na, u_scale, u_center, s_max)
    type(velocity_ring), intent(out) :: q
    integer, intent(in) :: nu, ns, na
    real(8), intent(in) :: u_scale, u_center, s_max
    integer :: i

    ! Reuse the verified axial Gauss-Hermite nodes and unweighted integral.
    ! In the transverse velocity plane use annuli and angle sectors.
    nv = nu
    Tref = 1.d0
    call init_velocity_gauss_hermite_2d()
    q%nu = nu
    allocate(q%u(nu), q%u2(nu), q%log_wu(nu))
    q%u = u_center + u_scale * uspace(:,1)
    q%u2 = q%u**2
    do i = 1, nu
      q%log_wu(i) = log(u_scale * sqrt(weight(i,i)))
    end do
    call make_transverse_ring(q, ns, na, s_max)
  end subroutine make_velocity_ring

  subroutine make_velocity_ring_half_range(q, nu, ns, na, u_scale, s_max)
    type(velocity_ring), intent(out) :: q
    integer, intent(in) :: nu, ns, na
    real(8), intent(in) :: u_scale, s_max
    integer :: nh

    if (nu < 2 .or. mod(nu,2) /= 0) error stop 'half-range Gauss-Hermite requires an even axial node count'
    if (min(ns,na) < 1 .or. min(u_scale,s_max) <= 0.d0) error stop 'invalid half-range velocity grid'
    nh = nu/2
    q%nu = nu
    allocate(q%u(nu), q%u2(nu), q%log_wu(nu))
    call half_gauss_hermite(nh, q%u(nh+1:), q%log_wu(nh+1:))

    ! Remove the Gaussian measure and apply the velocity Jacobian. Zero is
    ! the sign boundary for the axial wall flux on both mirrored intervals.
    q%log_wu(nh+1:) = q%log_wu(nh+1:) + q%u(nh+1:)**2 + log(u_scale)
    q%u(nh+1:) = u_scale*q%u(nh+1:)
    q%u(:nh) = -q%u(nu:nh+1:-1)
    q%log_wu(:nh) = q%log_wu(nu:nh+1:-1)
    q%u2 = q%u**2
    call make_transverse_ring(q, ns, na, s_max)
  end subroutine make_velocity_ring_half_range

  subroutine make_transverse_ring(q, ns, na, s_max)
    type(velocity_ring), intent(inout) :: q
    integer, intent(in) :: ns, na
    real(8), intent(in) :: s_max
    integer :: i, j, k
    real(8) :: s, s_lo, s_hi, lo, hi

    q%na = na
    q%np = ns * na
    q%angle_width = pi8 / na
    allocate(q%v(q%np), q%s2(q%np), q%sin_lower(q%np), q%sin_upper(q%np), q%log_wp(q%np))
    do j = 1, ns
      s_lo = s_max * (dble(j-1) / ns)**1.65d0
      s_hi = s_max * (dble(j) / ns)**1.65d0
      s = sqrt(0.5d0 * (s_lo**2 + s_hi**2))
      do i = 1, na
        k = (j-1)*na+i
        lo = (i-1)*q%angle_width
        hi = i*q%angle_width
        q%s2(k) = s*s
        q%sin_lower(k) = sin(lo)
        q%sin_upper(k) = sin(hi)
        if (i == na) q%sin_upper(k) = 0.d0
        ! Average cos(theta) within a sector to balance r and theta fluxes
        ! exactly. The other azimuthal half-plane follows by reflection.
        q%v(k) = s * (q%sin_upper(k)-q%sin_lower(k)) / q%angle_width
        q%log_wp(k) = log(q%angle_width * (s_hi**2-s_lo**2))
      end do
    end do
  end subroutine make_transverse_ring

  subroutine half_gauss_hermite(n, nodes, log_weights)
    integer, intent(in) :: n
    real(8), intent(out) :: nodes(n), log_weights(n)
    real(8), dimension(4*n+48) :: x, previous, current, next_values
    real(8) :: diagonal(n), coupling(n-1), cutoff, off, left, right, mid, pivot, pivot_min
    real(8) :: p0, p1, p2, factor, log_norm
    integer :: i, j, below

    ! Stieltjes recurrence for exp(-x*x) on the positive half-line. The wide
    ! integration interval makes the Gaussian tail negligible for this order.
    ! Multiplying each sample by sqrt(weight) avoids products of tiny and huge powers.
    cutoff = 2.d0*sqrt(dble(n))+8.d0
    call half_gaussian_samples(cutoff, x, current)
    current = current / sqrt(sum(current**2))
    previous = 0.d0
    off = 0.d0
    do i = 1, n
      diagonal(i) = sum(x*current**2)
      if (i == n) exit
      next_values = (x-diagonal(i))*current-off*previous
      off = sqrt(sum(next_values**2))
      coupling(i) = off
      previous = current
      current = next_values / off
    end do

    ! The Gaussian nodes are eigenvalues of the Jacobi matrix. A Sturm count
    ! from its pivots brackets each node and preserves the smallest roots.
    pivot_min = tiny(1.d0)*cutoff**2
    do i = 1, n
      left = 0.d0
      right = cutoff
      do
        mid = left+0.5d0*(right-left)
        if (mid <= left .or. mid >= right) exit
        pivot = diagonal(1)-mid
        if (abs(pivot) < pivot_min) pivot = -pivot_min
        below = merge(1,0,pivot < 0.d0)
        do j = 2, n
          pivot = diagonal(j)-mid-coupling(j-1)**2/pivot
          if (abs(pivot) < pivot_min) pivot = -pivot_min
          if (pivot < 0.d0) below = below+1
        end do
        if (below < i) then
          left = mid
        else
          right = mid
        end if
      end do
      nodes(i) = mid

      ! The weight is the zeroth moment times the squared first component
      ! of the unit eigenvector. Renormalize as the components are generated.
      p0 = 0.d0
      p1 = 1.d0
      log_norm = 0.d0
      off = 0.d0
      do j = 1, n-1
        p2 = ((mid-diagonal(j))*p1-off*p0)/coupling(j)
        factor = hypot(1.d0,p2)
        p0 = p1/factor
        p1 = p2/factor
        log_norm = log_norm+log(factor)
        off = coupling(j)
      end do
      log_weights(i) = log(0.5d0*sqrt(pi8))-2.d0*log_norm
    end do
  end subroutine half_gauss_hermite

  subroutine half_gaussian_samples(cutoff, x, sqrt_weight)
    real(8), intent(in) :: cutoff
    real(8), intent(out) :: x(:), sqrt_weight(size(x))
    real(8) :: z, p0, p1, p2, derivative, change, sqrt_legendre_weight
    integer :: i, j, k, iteration, nq

    ! Use a high-order Legendre rule on [0,cutoff] for the Stieltjes inner products.
    nq = size(x)
    do i = 1, nq/2
      z = cos(pi8*(i-0.25d0)/(nq+0.5d0))
      do iteration = 1, 30
        p0 = 1.d0
        p1 = z
        do j = 2, nq
          p2 = ((2.d0*j-1.d0)*z*p1-(j-1)*p0)/j
          p0 = p1
          p1 = p2
        end do
        derivative = nq*(p0-z*p1)/((1.d0-z)*(1.d0+z))
        change = p1/derivative
        if (abs(change) <= epsilon(z)) exit
        z = z-change
      end do
      if (iteration > 30) error stop 'half-range Gauss-Hermite integration did not converge'
      k = nq+1-i
      x(i) = 0.5d0*cutoff*(1.d0-z)
      x(k) = 0.5d0*cutoff*(1.d0+z)
      sqrt_legendre_weight = sqrt(cutoff/((1.d0-z)*(1.d0+z)))/abs(derivative)
      sqrt_weight(i) = sqrt_legendre_weight*exp(-0.5d0*x(i)**2)
      sqrt_weight(k) = sqrt_legendre_weight*exp(-0.5d0*x(k)**2)
    end do
  end subroutine half_gaussian_samples

  pure function ring_moments(q, f) result(w)
    type(velocity_ring), intent(in) :: q
    real(8), intent(in) :: f(q%nu,q%np)
    real(8) :: w(4), shell, m1, m2
    integer :: k

    ! f already includes its velocity integration weight. The fourth
    ! invariant is twice the total energy in the monatomic gas model.
    w = 0.d0
    do k = 1, q%np
      shell = sum(f(:,k))
      m1 = dot_product(q%u, f(:,k))
      m2 = dot_product(q%u2, f(:,k))
      w = w + [shell, m1, q%v(k)*shell, m2+q%s2(k)*shell]
    end do
  end function ring_moments

  pure subroutine ring_equilibrium(q, w, f, converged)
    type(velocity_ring), intent(in) :: q
    real(8), intent(in) :: w(4)
    real(8), intent(out) :: f(q%nu,q%np)
    logical, intent(out) :: converged
    real(8) :: target(3), multiplier(3), change(3), model(3), cov(5), residual(3)
    real(8) :: wu(q%nu), wp(q%np), try_u(q%nu), try_p(q%np)
    real(8) :: dual, trial_dual, temperature, remain, var_u, var_v, step
    integer :: n, ls, j

    f = 0.d0
    converged = .true.
    if (w(1) <= tiny(1.d0)) return
    target = w(2:4) / w(1)
    temperature = max((2.d0/3.d0)*(target(3)-sum(target(1:2)**2)), 0.02d0)
    multiplier = [2.d0*target(1), 2.d0*target(2), -1.d0] / temperature
    converged = .false.

    ! A positive discrete Maxwellian with matching moments makes the collision
    ! integral conservative on both a coarse and a fine velocity grid.
    do n = 1, 35
      call probabilities(q, multiplier, target, wu, wp, dual)
      call expectation(q, wu, wp, model, cov)
      residual = target - model
      if (maxval(abs(residual)/(1.d0+abs(target))) <= 5.d-12) then
        converged = .true.
        do j = 1, q%np
          f(:,j) = (w(1)*wp(j))*wu
        end do
        return
      end if

      var_u = max(cov(1), 1.d-15)
      var_v = max(cov(2), 1.d-15)
      remain = max(cov(3)-cov(4)**2/var_u-cov(5)**2/var_v, 1.d-15)
      change(3) = (residual(3)-cov(4)*residual(1)/var_u-cov(5)*residual(2)/var_v)/remain
      change(1) = (residual(1)-cov(4)*change(3))/var_u
      change(2) = (residual(2)-cov(5)*change(3))/var_v
      step = min(1.d0, 16.d0/max(16.d0, maxval(abs(change))))
      do ls = 1, 18
        call probabilities(q, multiplier+step*change, target, try_u, try_p, trial_dual)
        if (trial_dual <= dual+epsilon(dual)*32.d0*max(1.d0,abs(dual))) exit
        step = 0.5d0*step
      end do
      if (ls > 18) return
      multiplier = multiplier + step*change
    end do
  end subroutine ring_equilibrium

  pure subroutine probabilities(q, m, target, wu, wp, dual)
    type(velocity_ring), intent(in) :: q
    real(8), intent(in) :: m(3), target(3)
    real(8), intent(out) :: wu(q%nu), wp(q%np), dual
    real(8) :: lu, lp, su, sp

    wu = q%log_wu + m(1)*q%u + m(3)*q%u2
    wp = q%log_wp + m(2)*q%v + m(3)*q%s2
    lu = maxval(wu)
    lp = maxval(wp)
    wu = exp(wu-lu)
    wp = exp(wp-lp)
    su = sum(wu)
    sp = sum(wp)
    wu = wu / su
    wp = wp / sp
    dual = lu+lp+log(su)+log(sp)-dot_product(m,target)
  end subroutine probabilities

  pure subroutine expectation(q, wu, wp, m, cov)
    type(velocity_ring), intent(in) :: q
    real(8), intent(in) :: wu(q%nu), wp(q%np)
    real(8), intent(out) :: m(3), cov(5)
    real(8) :: kx, kr, du(q%nu), dv(q%np), eu(q%nu), ep(q%np)

    m(1) = dot_product(wu,q%u)
    m(2) = dot_product(wp,q%v)
    kx = dot_product(wu,q%u2)
    kr = dot_product(wp,q%s2)
    m(3) = kx + kr
    du = q%u-m(1)
    dv = q%v-m(2)
    eu = q%u2-kx
    ep = q%s2-kr
    cov = [sum(wu*du**2), sum(wp*dv**2), sum(wu*eu**2)+sum(wp*ep**2), &
           sum(wu*du*eu), sum(wp*dv*ep)]
  end subroutine expectation

  pure real(8) function flight_fraction(dt_over_tau)
    real(8), intent(in) :: dt_over_tau

    if (dt_over_tau < 1.d-4) then
      flight_fraction = 1.d0-dt_over_tau*(0.5d0-dt_over_tau/6.d0)
    else
      flight_fraction = (1.d0-exp(-dt_over_tau))/dt_over_tau
    end if
  end function flight_fraction
end module axisymmetric_kinetic

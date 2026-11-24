module axisymmetric_heat_flux
  use axisymmetric_kinetic, only: velocity_ring
  implicit none
  private
  public :: apply_heat_relaxation

contains

  pure subroutine apply_heat_relaxation(q,w,population,equilibrium,prandtl)
    type(velocity_ring), intent(in) :: q
    real(8), intent(in) :: w(4),population(q%nu,q%np),prandtl
    real(8), intent(inout) :: equilibrium(q%nu,q%np)
    integer, parameter :: x_power(6) = [0,1,0,0,1,0], r_power(6) = [0,0,1,0,0,1]
    integer, parameter :: e_power(6) = [0,0,0,1,1,1]
    real(8) :: scale, ux, ur, x(q%nu), x2(q%nu), x3(q%nu), r(q%np), z(q%np), &
               damp_x(q%nu), damp_r(q%np), along_x(q%nu), along_r(q%np), delta(q%nu), &
               moment_x(0:6), moment_r(0:2,0:2), coefficient(6), metric(6,6), &
               modification(q%nu,q%np), denominator, minimum
    integer :: i,j,k,exponent_x,exponent_r,exponent_e
    logical :: fits

    if (prandtl >= 1.d0 .or. w(1) <= 1.d-30) return
    ux = w(2)/w(1)
    ur = w(3)/w(1)
    scale = sqrt(max(0.d0,(2.d0/3.d0)*(w(4)/w(1)-ux**2-ur**2)))
    if (scale <= 1.d-7) return
    x = (q%u-ux)/scale
    x2 = x*x
    x3 = x2*x
    r = (q%v-ur)/scale
    z = max(q%s2-2.d0*ur*q%v+ur**2,0.d0)/scale**2

    coefficient = 0.d0
    do k = 1,q%np
      delta = population(:,k)-equilibrium(:,k)
      coefficient(5) = coefficient(5)+dot_product(delta,x3)+z(k)*dot_product(delta,x)
      coefficient(6) = coefficient(6)+r(k)*(dot_product(delta,x2)+z(k)*sum(delta))
    end do
    coefficient = coefficient*((1.d0-prandtl)/w(1))
    if (maxval(abs(coefficient(5:6))) < 1.d-10) return

    ! The fitted ring Maxwellian separates into axial and transverse weights.
    ! A separable damping factor bounds the correction in the far velocity tails.
    damp_x = exp(-0.125d0*x2)
    damp_r = exp(-0.125d0*z)
    along_x = sum(equilibrium,dim=2)*damp_x/w(1)
    along_r = sum(equilibrium,dim=1)*damp_r/w(1)
    do i = 0,6
      moment_x(i) = sum(along_x*x**i)
    end do
    do j = 0,2
      do i = 0,2
        moment_r(i,j) = sum(along_r*r**i*z**j)
      end do
    end do

    ! Enforce mass, both momenta and translational energy exactly on the discrete
    ! grid while imposing two heat-flux moments on the equilibrium correction.
    do j = 1,6
      do i = j,6
        exponent_x = x_power(i)+x_power(j)
        exponent_r = r_power(i)+r_power(j)
        exponent_e = e_power(i)+e_power(j)
        metric(i,j) = moment_x(exponent_x+2*exponent_e)*moment_r(exponent_r,0)
        if (exponent_e >= 1) metric(i,j) = metric(i,j)+moment_x(exponent_x)*moment_r(exponent_r,exponent_e)
        if (exponent_e == 2) metric(i,j) = metric(i,j)+2.d0*moment_x(exponent_x+2)*moment_r(exponent_r,1)
      end do
    end do
    call solve_metric(metric,coefficient,fits)
    if (.not. fits) return
    minimum = 0.d0
    do k = 1,q%np
      modification(:,k) = damp_x*damp_r(k)*(coefficient(1)+coefficient(2)*x+coefficient(3)*r(k) &
                            +(x2+z(k))*(coefficient(4)+coefficient(5)*x+coefficient(6)*r(k)))
      where (equilibrium(:,k) <= tiny(1.d0)) modification(:,k) = 0.d0
      minimum = min(minimum,minval(modification(:,k)))
    end do
    denominator = max(1.d0,-minimum/0.98d0)
    equilibrium = equilibrium*(1.d0+modification/denominator)
  end subroutine apply_heat_relaxation

  pure subroutine solve_metric(matrix,rhs,fits)
    real(8), intent(inout) :: matrix(6,6),rhs(6)
    logical, intent(out) :: fits
    real(8) :: diagonal(6), pivot
    integer :: i,j

    fits = .false.
    do i = 1,6
      if (matrix(i,i) <= tiny(1.d0)) return
      diagonal(i) = sqrt(matrix(i,i))
    end do
    rhs = rhs/diagonal
    do j = 1,6
      do i = j,6
        matrix(i,j) = matrix(i,j)/(diagonal(i)*diagonal(j))
      end do
    end do
    do j = 1,6
      pivot = matrix(j,j)-sum(matrix(j,1:j-1)**2)
      if (pivot <= 1.d-10) return
      matrix(j,j) = sqrt(pivot)
      do i = j+1,6
        matrix(i,j) = (matrix(i,j)-dot_product(matrix(i,1:j-1),matrix(j,1:j-1)))/matrix(j,j)
      end do
      rhs(j) = (rhs(j)-dot_product(matrix(j,1:j-1),rhs(1:j-1)))/matrix(j,j)
    end do
    do j = 6,1,-1
      rhs(j) = (rhs(j)-dot_product(matrix(j+1:6,j),rhs(j+1:6)))/matrix(j,j)
    end do
    rhs = rhs/diagonal
    fits = .true.
  end subroutine solve_metric
end module axisymmetric_heat_flux

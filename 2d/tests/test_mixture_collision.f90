program test_mixture_collision
  use axisymmetric_mixture_collision, only: mixture_species, mixture_model, make_mixture_model, mixture_targets
  implicit none

  integer, parameter :: n = 5
  real(8), parameter :: gas_r = 403.5744316182838d0, tref = 300.d0, rho_ref = 1.d-3
  real(8), parameter :: boltzmann = 1.380649d-23, avo = 6.02214076d23
  type(mixture_species) :: species(n)
  type(mixture_model) :: model
  real(8) :: w(5,n), advanced(5,n), equilibrium(5,n), frequency(n), mass_ratio(n), rot(n)
  integer :: i

  species%mass_kg = [28.013d0,18.016d0,44.011d0,28.011d0,2.016d0]*1.d-3/avo
  species%mass_fraction = [0.35d0,0.25d0,0.20d0,0.12d0,0.08d0]
  species%diameter_m = [4.00d0,3.00d0,4.40d0,3.80d0,2.70d0]*1.d-10
  species%omega = [0.74d0,0.80d0,0.85d0,0.76d0,0.70d0]
  species%rot_dof = [2,3,2,2,2]
  species%rot_probability = 0.2d0
  mass_ratio = species%mass_kg*gas_r/boltzmann
  rot = dble(species%rot_dof)
  call make_mixture_model(model,species,gas_r,tref,rho_ref,273.15d0)

  do i = 1, n
    w(:,i) = state(i,0.4d0,-0.1d0,1.8d0,1.8d0)
  end do
  call mixture_targets(model,w,0.d0,equilibrium,frequency)
  if (minval(frequency) <= 0.d0) error stop 'equilibrium pair rates must be positive'
  if (maxval(abs(equilibrium-w)) > 1.d-11) error stop 'common equilibrium is not stationary at a face'
  call collide(0.2d0)
  if (maxval(abs(advanced-w)) > 1.d-11) error stop 'common equilibrium is not stationary in a cell'

  do i = 1, n
    w(:,i) = state(i,0.8d0-0.25d0*i,0.08d0*i,2.5d0+0.3d0*i,0.35d0)
  end do
  call mixture_targets(model,w,0.d0,equilibrium,frequency)
  call check_invariants(matmul(equilibrium-w,frequency),1.d-8)
  call collide(2.d-4)
  if (sum(advanced(5,:)) <= sum(w(5,:))) error stop 'hot translation should heat the rotational modes'
  if (maxval(abs(advanced(2,:)/advanced(1,:)-w(2,:)/w(1,:))) < 1.d-3) error stop 'missing momentum exchange'
  if (any(advanced(5,:) < 0.d0)) error stop 'rotational energy must remain nonnegative'
  do i = 1, n
    if (advanced(4,i)*advanced(1,i) <= sum(advanced(2:3,i)**2)) error stop 'negative translational temperature'
  end do

  w(:,n) = 0.d0
  call collide(2.d-4)
  if (any(abs(advanced(:,n)) > 0.d0)) error stop 'collision created an absent species'

  do i = 1,n
    w(:,i) = state(i,0.d0,0.d0,0.7d0,6.d0)
  end do
  call collide(0.2d0)
  if (sum(advanced(5,:)) >= sum(w(5,:))) error stop 'hot rotation should heat translation in the stiff limit'

  species%rot_probability = 0.d0
  call make_mixture_model(model,species,gas_r,tref,rho_ref,273.15d0)
  call collide(2.d-4)
  if (maxval(abs(advanced(5,:)-w(5,:))) > 1.d-11) error stop 'frozen rotational energy changed'
  print '(a)', 'PASS mixture equilibrium, momentum exchange, rotational relaxation and invariants'

contains

  function state(i,u,v,temperature,rot_temperature) result(values)
    integer, intent(in) :: i
    real(8), intent(in) :: u,v,temperature,rot_temperature
    real(8) :: values(5)
    values = species(i)%mass_fraction * [1.d0,u,v,u*u+v*v+1.5d0*temperature/mass_ratio(i), &
                                        0.25d0*rot(i)*rot_temperature/mass_ratio(i)]
  end function state

  subroutine collide(dt)
    real(8), intent(in) :: dt
    real(8) :: factor
    integer :: k
    call mixture_targets(model,w,dt,equilibrium,frequency)
    do k = 1, n
      factor = dt*frequency(k)/(1.d0+dt*frequency(k))
      advanced(:,k) = w(:,k)+factor*(equilibrium(:,k)-w(:,k))
    end do
    if (maxval(abs(advanced(1,:)-w(1,:))) > 1.d-13) error stop 'species mass changed'
    call check_invariants(sum(advanced-w,dim=2),1.d-10)
  end subroutine collide

  subroutine check_invariants(delta,tolerance)
    real(8), intent(in) :: delta(5), tolerance
    if (maxval(abs(delta(2:3))) > tolerance) error stop 'mixture momentum changed'
    if (abs(0.5d0*delta(4)+delta(5)) > tolerance) error stop 'total energy changed'
  end subroutine check_invariants
end program test_mixture_collision

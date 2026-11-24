program test_thermal_transport
  use axisymmetric_kinetic, only: velocity_ring,make_velocity_ring_half_range,ring_equilibrium,ring_moments
  use axisymmetric_heat_flux, only: apply_heat_relaxation
  implicit none
  integer, parameter :: nu=32, ns=24, na=16, np=ns*na
  real(8), parameter :: pr=2.d0/3.d0, epsilon_t=1.d-5, lambda=1.d0
  type(velocity_ring) :: q
  real(8) :: eq(nu,np), minus(nu,np), plus(nu,np), source(nu,np), response(nu,np), next(nu,np), speed(nu,np)
  real(8) :: weight(nu,np), f(nu,np), modified(nu,np), perturbation(nu,np), w(4), before(2), after(2)
  real(8) :: viscosity,kappa,expected
  integer :: j,k
  logical :: fits

  call make_velocity_ring_half_range(q,nu,ns,na,1.65d0,8.d0)
  w = [1.d0,0.d0,0.d0,1.5d0]
  call ring_equilibrium(q,w,eq,fits)
  if (.not. fits) error stop 'cannot construct reference Maxwellian'
  do j=1,np
    speed(:,j)=q%u
    weight(:,j)=0.5d0*(q%u2+q%s2(j))
    perturbation(:,j)=q%u*(weight(:,j)-1.25d0)+q%v(j)*(weight(:,j)-1.25d0)
  end do
  f = eq*(1.d0+min(max(1.d-4*perturbation,-0.4d0),0.4d0))
  modified=eq
  before=heat_moment(f-eq)
  call apply_heat_relaxation(q,w,f,modified,pr)
  after=heat_moment(modified-eq)
  if (maxval(abs(ring_moments(q,modified)-w)) > 2.d-11) error stop 'heat correction changed collision invariants'
  if (maxval(abs(after-(1.d0-pr)*before)) > 1.d-12) error stop 'incorrect translational heat moments'
  if (minval(modified) < 0.d0) error stop 'negative heat correction population'

  ! Solve the steady linearized kinetic equation for a transverse Couette shear.
  call ring_equilibrium(q,[1.d0,0.d0,epsilon_t,1.5d0+epsilon_t**2],plus,fits)
  if (.not. fits) error stop 'cannot fit Couette test state'
  call ring_equilibrium(q,[1.d0,0.d0,-epsilon_t,1.5d0+epsilon_t**2],minus,fits)
  if (.not. fits) error stop 'cannot fit Couette test state'
  response=-speed*(plus-minus)/(2.d0*epsilon_t*lambda)
  viscosity=0.d0
  do j=1,np
    viscosity=viscosity-q%v(j)*dot_product(q%u,response(:,j))
  end do
  expected = dot_product(q%u2,sum(eq,dim=2))/lambda
  if (abs(viscosity-expected) > 1.d-9) error stop 'incorrect discrete Couette stress coefficient'
  if (abs(viscosity/0.5d0-1.d0) > 0.01d0) error stop 'grid does not resolve the analytic Couette coefficient'

  ! Isobaric Fourier conduction. Fixed point of the local linearized collision
  ! operator driven by dM/dx, with zero macroscopic drift and unit dT/dx.
  call ring_equilibrium(q,[1.d0/(1.d0+epsilon_t),0.d0,0.d0,1.5d0],plus,fits)
  if (.not. fits) error stop 'cannot fit Fourier test state'
  call ring_equilibrium(q,[1.d0/(1.d0-epsilon_t),0.d0,0.d0,1.5d0],minus,fits)
  if (.not. fits) error stop 'cannot fit Fourier test state'
  source=-speed*(plus-minus)/(2.d0*epsilon_t*lambda)
  response=source
  do k=1,25
    f=eq+epsilon_t*response
    modified=eq
    call apply_heat_relaxation(q,w,f,modified,pr)
    next=source+(modified-eq)/epsilon_t
    response=next
  end do
  after=heat_moment(source)
  expected=-after(1)/pr
  after=heat_moment(response)
  kappa=-after(1)
  if (abs(kappa/expected-1.d0) > 1.d-8) error stop 'failed thermal response calibration'
  if (abs(kappa/(1.25d0*viscosity/pr)-1.d0) > 0.015d0) error stop 'thermal coefficient differs from target'
  print '(a,2es16.7)', 'PASS linearized Couette and Fourier response, viscosity and conductivity ',viscosity,kappa
contains
  function heat_moment(values) result(out)
    real(8), intent(in) :: values(nu,np)
    real(8) :: out(2)
    integer :: k

    out(1)=sum(values*weight*speed)
    out(2)=0.d0
    do k=1,np
      out(2)=out(2)+q%v(k)*dot_product(values(:,k),weight(:,k))
    end do
  end function heat_moment
end program test_thermal_transport

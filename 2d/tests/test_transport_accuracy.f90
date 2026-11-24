program test_transport_accuracy
  use axisymmetric_mixture_collision
  implicit none
  real(8), parameter :: kb=1.380649d-23, avo=6.02214076d23, t=300.d0, tv=273.15d0
  real(8), parameter :: gas_r=403.5744316182838d0, rho_ref=1.d-3, c0=sqrt(2.d0*gas_r*t), dt=2.d-4
  type(mixture_species) :: s(2)
  type(mixture_model) :: base, corrected
  real(8) :: w(5,2), a(5,2), b(5,2), fa(2), fb(2), freq(2), post_a(5,2), post_b(5,2)
  real(8) :: factor, drag_a, drag_b, stress_mu, wilke_mu, temperature, ratio, phi(2,2), pure_mu(2)
  real(8) :: number(2), reduced, diameter, omega, encounter, diffusion, benchmark, resistance
  integer :: i,j

  s%mass_kg = [28.013d0,44.011d0]*1.d-3/avo
  s%mass_fraction = [0.6d0,0.4d0]
  s%diameter_m = [4.170d0,4.606d0]*1.d-10
  s%omega = [0.74d0,1.61d0]
  s%rot_dof = 2
  s%rot_probability = 0.2d0
  call make_mixture_model(base,s,gas_r,t,rho_ref,tv)
  call make_mixture_model(corrected,s,gas_r,t,rho_ref,tv,correct_vhs_drag=.true.)
  w(1,:) = s%mass_fraction
  w(2,:) = 1.d-5*[w(1,1),-w(1,2)]
  w(3,:) = 0.d0
  w(4,:) = 1.5d0*w(1,:)/base%mass+w(2,:)**2/w(1,:)
  w(5,:) = 0.5d0*w(1,:)/base%mass
  call mixture_targets(base,w,0.d0,a,fa)
  call mixture_targets(corrected,w,0.d0,b,fb)
  drag_a = fa(1)*(a(2,1)-w(2,1))
  drag_b = fb(1)*(b(2,1)-w(2,1))
  factor = 2.d0*(2.5d0-sum(s%omega)/2.d0)/3.d0
  if (abs(drag_b/drag_a-factor) > 1.d-10) error stop 'incorrect binary drag or diffusion'
  call check_invariants(matmul(b-w,fb))

  ! Infer the Maxwell-Stefan coefficient from the measured linear drag.
  number = rho_ref*w(1,:)/s%mass_kg
  resistance = -rho_ref*c0*drag_b/(2.d-5)
  diffusion = kb*t*product(number)/(sum(number)*resistance)
  reduced = product(s%mass_kg)/sum(s%mass_kg)
  diameter = sum(s%diameter_m)/2.d0
  omega = sum(s%omega)/2.d0
  encounter = acos(-1.d0)*diameter**2*sqrt(8.d0*kb*tv/(acos(-1.d0)*reduced))*(t/tv)**(1.d0-omega)
  benchmark = kb*t/(reduced*sum(number)*encounter*factor)
  if (abs(diffusion/benchmark-1.d0) > 1.d-8) error stop 'dimensional diffusion check failed'

  w(2,:) = 0.d0
  w(4,:) = 1.5d0*w(1,:)/base%mass
  call mixture_targets(corrected,w,0.d0,b,fb)
  if (maxval(abs(b-w)) > 1.d-12) error stop 'corrected closure perturbed equilibrium'
  call shear_frequencies(corrected,w,freq)
  temperature = t
  pure_mu = 15.d0*sqrt(s%mass_kg*kb*temperature/acos(-1.d0))/ &
     (2.d0*(5.d0-2.d0*s%omega)*(7.d0-2.d0*s%omega)*s%diameter_m**2)*(temperature/tv)**(s%omega-0.5d0)
  do i=1,2
    do j=1,2
      ratio=s(j)%mass_kg/s(i)%mass_kg
      phi(i,j)=(1.d0+sqrt(pure_mu(i)/pure_mu(j))*sqrt(sqrt(ratio)))**2/sqrt(8.d0*(1.d0+1.d0/ratio))
    end do
  end do
  wilke_mu=sum(w(1,:)/base%mass*pure_mu/matmul(phi,w(1,:)/base%mass))
  stress_mu=rho_ref*c0/2.d0*sum(w(1,:)/base%mass/freq)
  if (abs(stress_mu/wilke_mu-1.d0) > 1.d-12) error stop 'incorrect shear viscosity'

  w(2,:) = [0.5d0,-0.22d0]
  w(4,:) = 2.3d0*w(4,:)+w(2,:)**2/w(1,:)
  w(5,:) = [0.002d0,0.4d0]
  call mixture_targets(corrected,w,dt,a,fa)
  call retune_target(corrected,w,dt,a,fa,b,fb)
  do i=1,2
    post_a(:,i) = w(:,i)+dt*fa(i)/(1.d0+dt*fa(i))*(a(:,i)-w(:,i))
    post_b(:,i) = w(:,i)+dt*fb(i)/(1.d0+dt*fb(i))*(b(:,i)-w(:,i))
    if (b(4,i)*b(1,i) <= sum(b(2:3,i)**2) .or. b(5,i) < 0.d0) error stop 'unphysical retuned target'
  end do
  if (maxval(abs(post_a-post_b)) > 1.d-12) error stop 'stress rate changed conservative pair updates'
  call check_invariants(sum(post_b-w,dim=2))
  print '(a,3es16.7)', 'PASS VHS drag factor, Wilke and discrete BGK shear viscosity ',factor,wilke_mu,stress_mu
  print '(a,3es16.7)', 'PASS N2-CO2 diffusion (old, retuned, VHS benchmark), m2/s ',diffusion*factor,diffusion,benchmark
contains
  subroutine check_invariants(delta)
    real(8), intent(in) :: delta(5)
    if (maxval(abs(delta(1:3))) > 1.d-12) error stop 'species collision changes mass or mixture momentum'
    if (abs(0.5d0*delta(4)+delta(5)) > 1.d-12) error stop 'species collision changes total energy'
  end subroutine check_invariants
end program test_transport_accuracy

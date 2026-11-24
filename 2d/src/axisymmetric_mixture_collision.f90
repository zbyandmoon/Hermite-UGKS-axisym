module axisymmetric_mixture_collision
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use axisymmetric_internal_modes, only: parker_probability
  implicit none
  private
  public :: mixture_species, mixture_model, make_mixture_model, mixture_targets
  public :: shear_frequencies, retune_target

  real(8), parameter :: kb = 1.380649d-23, pi8 = acos(-1.d0)

  type mixture_species
    character(12) :: name = ''
    real(8) :: mass_kg, mass_fraction, diameter_m, omega, rot_probability
    integer :: rot_dof
    real(8) :: rotation_tstar = 0.d0
  end type mixture_species

  type mixture_model
    integer :: n
    real(8), allocatable :: mass(:), mass_share(:,:), reduced_mass(:,:)
    real(8), allocatable :: rate_prefactor(:,:), rate_exponent(:,:)
    real(8), allocatable :: energy_map(:,:,:,:), drift_map(:,:,:)
    real(8), allocatable :: drag_factor(:,:), pure_mu(:), mu_exponent(:), phi_weight(:,:), phi_denominator(:,:)
    real(8), allocatable :: rot_probability(:), rot_tstar(:)
    integer, allocatable :: rot_dof(:)
    real(8) :: temperature_ratio, temperature_ref, stress_scale
    logical :: variable_rotation
  end type mixture_model

contains

  subroutine make_mixture_model(model,species,gas_r,tref,rho_ref,vhs_tref,correct_vhs_drag)
    type(mixture_model), intent(out) :: model
    type(mixture_species), intent(in) :: species(:)
    real(8), intent(in) :: gas_r,tref,rho_ref,vhs_tref
    logical, intent(in), optional :: correct_vhs_drag
    integer :: i,j,n
    real(8) :: mref,c0,a,b,d,omega,chi,ratio
    logical :: correct_drag

    n = size(species)
    if (n < 1 .or. min(gas_r,tref,rho_ref,vhs_tref) <= 0.d0) error stop 'invalid mixture scales'
    if (.not. all(ieee_is_finite(species%mass_kg)) .or. minval(species%mass_kg) <= 0.d0 .or. &
        .not. all(ieee_is_finite(species%diameter_m)) .or. minval(species%diameter_m) <= 0.d0 .or. &
        .not. all(ieee_is_finite(species%mass_fraction)) .or. minval(species%mass_fraction) < 0.d0 .or. &
        abs(sum(species%mass_fraction)-1.d0) > 1.d-9) error stop 'invalid mixture species data'
    if (.not. all(ieee_is_finite(species%omega)) .or. minval(species%omega) < 0.5d0 .or. &
        maxval(species%omega) >= 2.5d0) error stop 'invalid VHS exponent'
    if (.not. all(ieee_is_finite(species%rot_probability)) .or. minval(species%rot_probability) < 0.d0 .or. &
        maxval(species%rot_probability) > 1.d0 .or. minval(species%rot_dof) < 0 .or. &
        maxval(species%rot_dof) > 3 .or. .not. all(ieee_is_finite(species%rotation_tstar)) .or. &
        minval(species%rotation_tstar) < 0.d0) error stop 'invalid rotational relaxation data'

    mref = kb/gas_r
    c0 = sqrt(2.d0*gas_r*tref)
    model%n = n
    model%temperature_ratio = tref/vhs_tref
    model%temperature_ref = tref
    model%stress_scale = rho_ref*c0/2.d0
    model%mass = species%mass_kg/mref
    model%rot_dof = species%rot_dof
    model%rot_probability = species%rot_probability
    model%rot_tstar = species%rotation_tstar
    model%variable_rotation = any(model%rot_tstar > 0.d0)
    model%mu_exponent = species%omega
    model%pure_mu = 15.d0*sqrt(species%mass_kg*kb*tref/pi8)*model%temperature_ratio**(species%omega-0.5d0) / &
      (2.d0*(5.d0-2.d0*species%omega)*(7.d0-2.d0*species%omega)*species%diameter_m**2)
    correct_drag = .false.
    if (present(correct_vhs_drag)) correct_drag = correct_vhs_drag
    allocate(model%mass_share(n,n),model%reduced_mass(n,n), &
             model%rate_prefactor(n,n),model%rate_exponent(n,n),model%energy_map(4,4,n,n),model%drift_map(4,n,n))
    allocate(model%drag_factor(n,n),model%phi_weight(n,n),model%phi_denominator(n,n))
    do j = 1, n
      do i = 1, n
        a = model%mass(i)/(model%mass(i)+model%mass(j))
        b = 1.d0-a
        model%mass_share(i,j) = a
        model%reduced_mass(i,j) = b*model%mass(i)
        d = 0.5d0*(species(i)%diameter_m+species(j)%diameter_m)
        omega = 0.5d0*(species(i)%omega+species(j)%omega)
        chi = 2.5d0-omega
        model%rate_exponent(i,j) = 1.d0-omega
        model%drag_factor(i,j) = 1.d0
        if (correct_drag) model%drag_factor(i,j) = (2.d0/3.d0)*chi
        ratio = model%mass(j)/model%mass(i)
        model%phi_weight(i,j) = sqrt(sqrt(ratio))
        model%phi_denominator(i,j) = sqrt(8.d0*(1.d0+1.d0/ratio))
        ! VHS rate averaged over a Maxwellian relative velocity at the reference temperature.
        model%rate_prefactor(i,j) = pi8*d*d*sqrt(8.d0*kb*vhs_tref/(pi8*mref*model%reduced_mass(i,j))) &
                                   *rho_ref/(mref*c0)
        call collision_maps(a,model%reduced_mass(i,j),chi,model%drag_factor(i,j), &
          species(i)%rot_dof,species(j)%rot_dof,species(i)%rot_probability,species(j)%rot_probability, &
          model%energy_map(:,:,i,j),model%drift_map(:,i,j))
      end do
    end do
  end subroutine make_mixture_model

  pure subroutine collision_maps(a,reduced,chi,drag_factor,dof_i,dof_j,prob_i,prob_j,energy_map,drift_map)
    real(8), intent(in) :: a,reduced,chi,drag_factor,prob_i,prob_j
    integer, intent(in) :: dof_i,dof_j
    real(8), intent(out) :: energy_map(4,4),drift_map(4)
    real(8) :: b, redistribute, into_relative(3,4),out_of_relative(4,3),center_energy(4,4)
    real(8) :: exchange_i(3,3),exchange_j(3,3),exchange(3,3),identity(3,3)
    integer :: k

    call rotational_exchange(chi,dof_i,prob_i,2,exchange_i)
    call rotational_exchange(chi,dof_j,prob_j,3,exchange_j)
    exchange = 0.5d0*(matmul(exchange_i,exchange_j)+matmul(exchange_j,exchange_i))
    if (abs(drag_factor-1.d0) > epsilon(1.d0)) then
      identity = 0.d0
      do k = 1,3
        identity(k,k) = 1.d0
      end do
      ! Change the momentum transfer rate while keeping the rotational source
      ! at its original rate whenever the rescaled exchange map is positive.
      redistribute = max(drag_factor,1.d0-exchange(1,1),1.d0-exchange(2,2),1.d0-exchange(3,3))
      exchange = identity+(exchange-identity)/redistribute
    end if

    ! Separate center-of-mass energy from the energy available to the collision.
    b = 1.d0-a
    into_relative = 0.d0
    into_relative(1,1:2) = [b,a]
    into_relative(2,3) = 1.d0
    into_relative(3,4) = 1.d0
    out_of_relative = 0.d0
    out_of_relative(1:2,1) = [b,a]
    out_of_relative(3,2) = 1.d0
    out_of_relative(4,3) = 1.d0
    center_energy = 0.d0
    center_energy(1,1:2) = [a*a,a*b]
    center_energy(2,1:2) = [a*b,b*b]
    energy_map = center_energy+matmul(out_of_relative,matmul(exchange,into_relative))
    drift_map = matmul(out_of_relative,exchange(:,1))*(0.5d0*reduced)
  end subroutine collision_maps

  pure subroutine rotational_exchange(chi,dof,probability,mode,exchange)
    real(8), intent(in) :: chi,probability
    integer, intent(in) :: dof,mode
    real(8), intent(out) :: exchange(3,3)
    real(8) :: rot_to_tr,tr_to_rot

    exchange = 0.d0
    exchange(1,1) = 1.d0
    exchange(2,2) = 1.d0
    exchange(3,3) = 1.d0
    if (dof == 0) return
    ! Mean Larsen-Borgnakke exchange with a VHS-weighted relative energy.
    ! Convert the collision-weighted relative energy back to the three
    ! unweighted translational degrees so a common temperature is stationary.
    rot_to_tr = probability*chi/(chi+0.5d0*dof)
    tr_to_rot = rot_to_tr*dof/3.d0
    exchange(1,mode) = rot_to_tr
    exchange(mode,mode) = 1.d0-rot_to_tr
    exchange(mode,1) = tr_to_rot
    exchange(1,1) = 1.d0-tr_to_rot
  end subroutine rotational_exchange

  pure subroutine mixture_targets(model,w,dt,target,frequency,encounters)
    type(mixture_model), intent(in) :: model
    real(8), intent(in) :: w(5,model%n),dt
    real(8), intent(out) :: target(5,model%n),frequency(model%n)
    real(8), intent(out), optional :: encounters(model%n)
    integer :: i,j,n,indices(4)
    real(8) :: density(model%n), u_old(model%n,2), u(model%n,2), energy(2*model%n), rates(model%n,model%n)
    real(8) :: relative_energy, drift_sq, a, factor, ui(2), energy_rhs(2*model%n,1), local_energy(4), sum_energy(2)
    real(8) :: momentum_system(model%n,model%n), energy_system(2*model%n,2*model%n), sum_velocity(2), spread
    real(8) :: energy_maps(4,4,model%n,model%n), drift_maps(4,model%n,model%n), pair_temp, prob_i, prob_j

    n = model%n
    target = w
    frequency = 0.d0
    if (present(encounters)) encounters = 0.d0
    density = w(1,:)
    where (density < 1.d-30) density = 0.d0
    if (sum(density) <= 0.d0) return
    u_old = 0.d0
    energy = 0.d0
    do i = 1,n
      if (density(i) <= 0.d0) cycle
      u_old(i,:) = w(2:3,i)/density(i)
      energy(i) = 0.5d0*model%mass(i)*max(w(4,i)/density(i)-sum(u_old(i,:)**2),0.d0)
      energy(n+i) = model%mass(i)*max(w(5,i)/density(i),0.d0)
    end do
    energy_maps = model%energy_map
    drift_maps = model%drift_map
    do j = 1,n
      do i = 1,n
        relative_energy = model%reduced_mass(i,j)* &
          (energy(i)/model%mass(i)+energy(j)/model%mass(j)+0.5d0*sum((u_old(i,:)-u_old(j,:))**2))
        rates(i,j) = model%rate_prefactor(i,j)*(density(j)/model%mass(j))* &
          max(1.d-4,(4.d0/3.d0)*relative_energy*model%temperature_ratio)**model%rate_exponent(i,j)
        if (model%variable_rotation .and. density(i) > 0.d0 .and. density(j) > 0.d0) then
          pair_temp = max(1.d-4,4.d0*relative_energy/3.d0)*model%temperature_ref
          prob_i = parker_probability(model%rot_probability(i),model%rot_tstar(i),pair_temp,model%temperature_ref)
          prob_j = parker_probability(model%rot_probability(j),model%rot_tstar(j),pair_temp,model%temperature_ref)
          call collision_maps(model%mass_share(i,j),model%reduced_mass(i,j),1.5d0+model%rate_exponent(i,j), &
            model%drag_factor(i,j),model%rot_dof(i),model%rot_dof(j),prob_i,prob_j, &
            energy_maps(:,:,i,j),drift_maps(:,i,j))
        end if
      end do
    end do
    if (present(encounters)) encounters = sum(rates,dim=2)
    rates = rates*model%drag_factor
    frequency = sum(rates,dim=2)
    if (dt > 0.d0 .and. dt*maxval(frequency) < 1.d-8) then
      frequency = 0.d0
      return
    end if

    u = u_old
    if (dt > 0.d0) then
      ! Solve the frozen-rate pair momentum and energy equations implicitly.
      ! The implicit matrices have nonpositive off-diagonals. Combining
      ! positive thermal variances avoids cancellation of the drift energy.
      momentum_system = 0.d0
      do i = 1,n
        momentum_system(i,i) = 1.d0
        do j = 1,n
          if (i == j) cycle
          factor = dt*rates(i,j)*(1.d0-model%mass_share(i,j))
          momentum_system(i,i) = momentum_system(i,i)+factor
          momentum_system(i,j) = -factor
        end do
      end do
      call positive_linear_solve(momentum_system,u)

      energy_system = 0.d0
      energy_rhs(:,1) = energy
      do i = 1,n
        energy_system(i,i) = 1.d0+dt*frequency(i)
        energy_system(n+i,n+i) = energy_system(i,i)
        energy_rhs(i,1) = energy_rhs(i,1)+0.5d0*model%mass(i)*sum((u_old(i,:)-u(i,:))**2)
        do j = 1,n
          factor = dt*rates(i,j)
          drift_sq = sum((u(i,:)-u(j,:))**2)
          a = model%mass_share(i,j)
          indices = [i,j,n+i,n+j]
          call subtract_pair(energy_system(i,:),indices,factor*energy_maps(1,:,i,j))
          call subtract_pair(energy_system(n+i,:),indices,factor*energy_maps(3,:,i,j))
          energy_rhs(i,1) = energy_rhs(i,1)+factor*drift_sq* &
            (drift_maps(1,i,j)+0.5d0*model%mass(i)*(1.d0-a)**2)
          energy_rhs(n+i,1) = energy_rhs(n+i,1)+factor*drift_sq*drift_maps(3,i,j)
        end do
      end do
      call positive_linear_solve(energy_system,energy_rhs)
      energy = energy_rhs(:,1)
    end if

    do i = 1,n
      if (density(i) <= 0.d0 .or. frequency(i) <= 0.d0) cycle
      sum_velocity = 0.d0
      sum_energy = 0.d0
      spread = 0.d0
      do j = 1,n
        factor = rates(i,j)/frequency(i)
        a = model%mass_share(i,j)
        ui = a*u(i,:)+(1.d0-a)*u(j,:)
        local_energy = matmul(energy_maps(:,:,i,j),[energy(i),energy(j),energy(n+i),energy(n+j)]) &
                       +drift_maps(:,i,j)*sum((u(i,:)-u(j,:))**2)
        sum_velocity = sum_velocity+factor*ui
        sum_energy = sum_energy+factor*local_energy([1,3])
        spread = spread+factor*sum((ui-u(i,:))**2)
      end do
      spread = max(spread-sum((sum_velocity-u(i,:))**2),0.d0)
      target(2:3,i) = density(i)*sum_velocity
      target(4,i) = density(i)*(sum(sum_velocity**2)+spread+2.d0*sum_energy(1)/model%mass(i))
      target(5,i) = density(i)*sum_energy(2)/model%mass(i)
    end do
  end subroutine mixture_targets

  pure subroutine shear_frequencies(model,w,frequency)
    type(mixture_model), intent(in) :: model
    real(8), intent(in) :: w(5,model%n)
    real(8), intent(out) :: frequency(model%n)
    real(8) :: number(model%n), root_mu(model%n), temperature, variance, phi
    integer :: i,j

    number = max(w(1,:),0.d0)/model%mass
    frequency = 0.d0
    if (sum(number) <= 1.d-30) return
    variance = 0.d0
    do i = 1,model%n
      if (w(1,i) > 1.d-30) variance = variance+max(w(4,i)-sum(w(2:3,i)**2)/w(1,i),0.d0)
    end do
    temperature = max((2.d0/3.d0)*variance/sum(number),1.d-4)
    root_mu = sqrt(model%pure_mu)*temperature**(0.5d0*model%mu_exponent)
    do i = 1,model%n
      do j = 1,model%n
        phi = (1.d0+model%phi_weight(i,j)*root_mu(i)/root_mu(j))**2/model%phi_denominator(i,j)
        frequency(i) = frequency(i)+number(j)*phi
      end do
    end do
    frequency = (model%stress_scale*temperature)*frequency/root_mu**2
  end subroutine shear_frequencies

  pure subroutine retune_target(model,w,dt,pair_target,pair_frequency,target,frequency)
    type(mixture_model), intent(in) :: model
    real(8), intent(in) :: w(5,model%n),dt,pair_target(5,model%n),pair_frequency(model%n)
    real(8), intent(out) :: target(5,model%n),frequency(model%n)
    real(8) :: advanced(5), proposal(5), ratio, baseline_variance
    integer :: i,k

    call shear_frequencies(model,w,frequency)
    target = pair_target
    do i = 1,model%n
      if (w(1,i) <= 1.d-30 .or. min(frequency(i),pair_frequency(i)) <= 0.d0) then
        frequency(i) = pair_frequency(i)
        cycle
      end if
      advanced = pair_target(:,i)+(w(:,i)-pair_target(:,i))/(1.d0+dt*pair_frequency(i))
      ratio = pair_frequency(i)/frequency(i)
      baseline_variance = pair_target(4,i)*w(1,i)-sum(pair_target(2:3,i)**2)
      do k = 1,12
        proposal = advanced+ratio*(pair_target(:,i)-advanced)
        if (proposal(5) >= 0.d0 .and. proposal(4)*w(1,i)-sum(proposal(2:3)**2) > 0.02d0*baseline_variance) exit
        ratio = 0.5d0*(1.d0+ratio)
      end do
      if (k <= 12) then
        target(:,i) = proposal
        frequency(i) = pair_frequency(i)/ratio
      else
        frequency(i) = pair_frequency(i)
      end if
    end do
  end subroutine retune_target

  pure subroutine subtract_pair(row,indices,value)
    real(8), intent(inout) :: row(:)
    integer, intent(in) :: indices(4)
    real(8), intent(in) :: value(4)
    integer :: k
    ! At self collision the indices repeat and each contribution is accumulated.
    do k = 1,4
      row(indices(k)) = row(indices(k))-value(k)
    end do
  end subroutine subtract_pair

  pure subroutine positive_linear_solve(matrix,rhs)
    real(8), intent(inout) :: matrix(:,:),rhs(:,:)
    integer :: i,k,n
    real(8) :: factor
    n = size(matrix,1)
    do k = 1,n-1
      do i = k+1,n
        factor = matrix(i,k)/matrix(k,k)
        matrix(i,k+1:n) = matrix(i,k+1:n)-factor*matrix(k,k+1:n)
        rhs(i,:) = rhs(i,:)-factor*rhs(k,:)
      end do
    end do
    do i = n,1,-1
      do k = i+1,n
        rhs(i,:) = rhs(i,:)-matrix(i,k)*rhs(k,:)
      end do
      rhs(i,:) = rhs(i,:)/matrix(i,i)
    end do
  end subroutine positive_linear_solve
end module axisymmetric_mixture_collision

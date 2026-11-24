program test_internal_modes
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use axisymmetric_internal_modes, only: parker_probability, vibration_energy, vibration_target
  implicit none

  real(8), parameter :: probability_ref = 0.2d0, tstar = 91.5d0, tref = 300.d0
  real(8), parameter :: gas_r = 300.d0, theta = 3000.d0, translational_equilibrium = 1350000.d0
  real(8), parameter :: vibrational_equilibrium = 523779.0361823938d0
  real(8), parameter :: polyatomic_theta(4) = [960.d0,960.d0,1900.d0,3380.d0]
  real(8), parameter :: polyatomic_translation = 340020.d0, polyatomic_vibration = 429268.0906044184d0
  real(8) :: probability(3), target_tr, target_vib, energy_total

  probability = parker_probability(probability_ref,tstar,[100.d0,tref,600.d0],tref)
  if (maxval(abs(probability-[0.4140639669598913d0,0.2d0,0.1385301958061796d0])) > 1.d-14) &
    error stop 'Parker probability must follow the normalized collision number'
  if (abs(parker_probability(probability_ref,50.d0,600.d0,tref)-0.1478862316696196d0) > 1.d-14) &
    error stop 'Parker response must depend on the species parameter'
  if (abs(parker_probability(probability_ref,0.d0,600.d0,tref)-probability_ref) > 1.d-15) &
    error stop 'zero Parker temperature must recover the constant probability'
  if (abs(parker_probability(probability_ref,tstar,tiny(1.d0),tref)-1.d0) > 1.d-15) &
    error stop 'Parker probability must saturate at one in the cold limit'
  if (abs(parker_probability(0.d0,tstar,1.d0,tref)) > 1.d-15) &
    error stop 'zero rotational probability must remain elastic'

  if (abs(vibration_energy(theta,600.d0,gas_r)-6105.289415673808d0) > 1.d-8) &
    error stop 'harmonic oscillator energy without the zero-point term is incorrect'
  if (abs(vibration_energy(theta,1.d0,gas_r)) > tiny(1.d0)) &
    error stop 'the harmonic oscillator must freeze out in the cold limit'
  if (abs(vibration_energy(theta,1.d16,gas_r)/(gas_r*1.d16)-1.d0) > 1.d-12) &
    error stop 'the harmonic oscillator must recover the high-temperature limit'
  if (abs(vibration_energy(0.d0,600.d0,gas_r)) > tiny(1.d0)) &
    error stop 'zero characteristic temperature must disable the mode'

  energy_total = translational_equilibrium+vibrational_equilibrium
  call vibration_target(energy_total,0.d0,theta,gas_r,target_tr,target_vib)
  call check_target()

  call vibration_target(0.d0,energy_total,theta,gas_r,target_tr,target_vib)
  call check_target()

  call vibration_target(translational_equilibrium,vibrational_equilibrium,theta,gas_r,target_tr,target_vib)
  call check_target()

  call vibration_target(450.d0,0.d0,theta,gas_r,target_tr,target_vib)
  if (abs(target_tr-450.d0) > 1.d-12 .or. abs(target_vib) > tiny(1.d0)) &
    error stop 'a cold translational mode must stay cold without numerical energy loss'

  call vibration_target(1000.d0,20.d0,0.d0,gas_r,target_tr,target_vib)
  if (abs(target_tr-1000.d0) > 1.d-12 .or. abs(target_vib-20.d0) > 1.d-12) &
    error stop 'an inactive vibrational mode must preserve both input energies'

  energy_total = polyatomic_translation+polyatomic_vibration
  call vibration_target(energy_total,0.d0,[polyatomic_theta,0.d0],188.9d0,target_tr,target_vib)
  if (abs(target_tr-polyatomic_translation) > 1.d-8 .or. abs(target_vib-polyatomic_vibration) > 1.d-8) &
    error stop 'the polyatomic target must include each oscillator at the shared temperature'
  print '(a)', 'PASS Parker probability, vibrational energy limits, and conservative T-V equilibrium'

contains

  subroutine check_target()
    if (.not. all(ieee_is_finite([target_tr,target_vib]))) error stop 'the T-V target must be finite'
    if (min(target_tr,target_vib) <= 0.d0) error stop 'the T-V target must be positive'
    if (abs((target_tr+target_vib)-energy_total) > 4.d0*epsilon(1.d0)*energy_total) &
      error stop 'T-V exchange must conserve energy'
    if (abs(target_tr-translational_equilibrium) > 1.d-8 .or. &
        abs(target_vib-vibrational_equilibrium) > 1.d-8) &
      error stop 'the T-V target energies must correspond to the same temperature'
  end subroutine check_target
end program test_internal_modes

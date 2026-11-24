module axisymmetric_internal_modes
  implicit none
  private
  public :: parker_probability, vibration_energy, vibration_target

  real(8), parameter :: pi = acos(-1.d0)
  real(8), parameter :: parker_sqrt = 0.5d0*pi*sqrt(pi), parker_linear = pi*(1.d0+0.25d0*pi)
  real(8), parameter :: cv_translation = 1.5d0

  interface vibration_target
    module procedure single_vibration_target, band_vibration_target
  end interface vibration_target

contains

  pure elemental function parker_probability(probability_ref,tstar_k,temperature_k,tref_k) result(probability)
    real(8), intent(in) :: probability_ref,tstar_k,temperature_k,tref_k
    real(8) :: probability, scale, tref_scaled, tstar_scaled, root_term, linear_term, tref_weight, temperature_weight
    real(8) :: scaled_probability

    probability = max(0.d0,min(1.d0,probability_ref))
    if (probability <= 0.d0 .or. tstar_k <= 0.d0 .or. tref_k <= 0.d0) return
    if (temperature_k <= 0.d0) then
      probability = 1.d0
      return
    end if

    ! P(T) = P(Tref) * Z(Tref)/Z(T) using Parker's collision number (doi:10.1063/1.1724417).
    ! Scale the terms so the low-temperature saturation never overflows before clamping.
    scale = max(tref_k,tstar_k)
    tref_scaled = tref_k/scale
    tstar_scaled = tstar_k/scale
    root_term = parker_sqrt*sqrt(tref_scaled)*sqrt(tstar_scaled)
    linear_term = parker_linear*tstar_scaled
    scale = max(tref_k,temperature_k)
    tref_weight = tref_k/scale
    temperature_weight = temperature_k/scale
    scaled_probability = probability*((tref_scaled*temperature_weight + &
      root_term*sqrt(tref_weight)*sqrt(temperature_weight) + linear_term*tref_weight) / &
      (tref_scaled+root_term+linear_term))
    if (scaled_probability >= temperature_weight) then
      probability = 1.d0
    else
      probability = scaled_probability/temperature_weight
    end if
  end function parker_probability

  pure elemental function vibration_energy(theta_k,temperature_k,gas_r) result(energy)
    real(8), intent(in) :: theta_k,temperature_k,gas_r
    real(8) :: energy, energy_k, heat_capacity

    call oscillator(theta_k,temperature_k,energy_k,heat_capacity)
    energy = gas_r*energy_k
  end function vibration_energy

  pure elemental subroutine oscillator(theta_k,temperature_k,energy_k,heat_capacity)
    real(8), intent(in) :: theta_k,temperature_k
    real(8), intent(out) :: energy_k,heat_capacity
    real(8) :: ratio, ratio_sq, boltzmann, excitation

    energy_k = 0.d0
    heat_capacity = 0.d0
    if (theta_k <= 0.d0 .or. temperature_k <= 0.d0) return
    if (theta_k/700.d0 > temperature_k) return

    ! Energy and its temperature derivative are divided by R. The small-ratio expansion
    ! avoids loss of precision in 1-exp(-theta/T). Zero-point energy is omitted.
    ratio = theta_k/temperature_k
    if (ratio < 1.d-3) then
      ratio_sq = ratio*ratio
      excitation = 1.d0-0.5d0*ratio+ratio_sq*(1.d0/12.d0-ratio_sq/720.d0)
    else
      boltzmann = exp(-ratio)
      excitation = ratio*boltzmann/(1.d0-boltzmann)
    end if
    energy_k = temperature_k*excitation
    heat_capacity = excitation*(excitation+ratio)
  end subroutine oscillator

  pure subroutine single_vibration_target(translation_j_kg,vibration_j_kg,theta_k,gas_r,target_translation,target_vibration)
    real(8), intent(in) :: translation_j_kg,vibration_j_kg,theta_k,gas_r
    real(8), intent(out) :: target_translation,target_vibration

    call band_vibration_target(translation_j_kg,vibration_j_kg,[theta_k],gas_r,target_translation,target_vibration)
  end subroutine single_vibration_target

  pure subroutine band_vibration_target(translation_j_kg,vibration_j_kg,theta_k,gas_r,target_translation,target_vibration)
    real(8), intent(in) :: translation_j_kg,vibration_j_kg,theta_k(:),gas_r
    real(8), intent(out) :: target_translation,target_vibration
    real(8) :: energy_total, total_k, temperature, residual, tolerance, energy_k(size(theta_k)), capacity(size(theta_k))
    integer :: iteration

    target_translation = translation_j_kg
    target_vibration = vibration_j_kg
    if (gas_r <= 0.d0 .or. .not. any(theta_k > 0.d0)) return
    energy_total = translation_j_kg+vibration_j_kg
    if (energy_total <= 0.d0) return

    ! Solve E_tr(T)+E_v(T)=E_tr+E_v. For positive temperatures the left-hand side is
    ! increasing and convex, so Newton iterations from the all-translational state
    ! approach equilibrium from above. Repeat a band temperature for a degenerate mode.
    total_k = energy_total/gas_r
    temperature = total_k/cv_translation
    tolerance = 4.d0*epsilon(1.d0)*total_k
    do iteration = 1,20
      call oscillator(theta_k,temperature,energy_k,capacity)
      residual = (cv_translation*temperature-total_k)+sum(energy_k)
      if (abs(residual) <= tolerance) exit
      temperature = temperature-residual/(cv_translation+sum(capacity))
    end do

    target_vibration = gas_r*sum(energy_k)
    target_translation = energy_total-target_vibration
  end subroutine band_vibration_target

end module axisymmetric_internal_modes

program test_axisymmetric_mesh
  use axisymmetric_mesh, only: build_axis, limited_slope, flight_slope_fraction
  implicit none

  real(8), parameter :: tolerance = 32.d0*epsilon(1.d0), pi = acos(-1.d0)
  real(8), allocatable :: faces(:), centers(:), widths(:), radial_faces(:), radial_centers(:), radial_widths(:)
  real(8) :: annular_area, body_volume, left_state, right_state, slope
  character(24) :: rejection_case
  integer :: i, j

  ! Run each invalid input in a separate process so ERROR STOP can be checked.
  call get_command_argument(1,rejection_case)
  select case (trim(rejection_case))
  case ('misaligned')
    call build_axis(10,1.d0,[0.15d0],[0.5d0],[2],faces,centers,widths)
  case ('overlap')
    call build_axis(10,1.d0,[0.2d0,0.4d0],[0.5d0,0.6d0],[2,3],faces,centers,widths)
  case ('reversed')
    call build_axis(10,1.d0,[0.6d0],[0.5d0],[2],faces,centers,widths)
  case ('outside')
    call build_axis(10,1.d0,[0.5d0],[1.1d0],[2],faces,centers,widths)
  case ('bad_split')
    call build_axis(10,1.d0,[0.2d0],[0.5d0],[0],faces,centers,widths)
  case ('bad_length')
    call build_axis(10,0.d0,[0.2d0],[0.5d0],[2],faces,centers,widths)
  case ('shape')
    call build_axis(10,1.d0,[0.2d0],[0.5d0,0.6d0],[2],faces,centers,widths)
  case ('')
    call check_uniform_axis()
    call check_local_refinement()
    call check_reconstruction()
    call check_flight_coefficient()
    print '(a)', 'PASS axisymmetric mesh partition, body geometry, reconstruction and free-flight coefficient'
  case default
    error stop 'unrecognized axisymmetric mesh test case'
  end select

contains

  subroutine check_uniform_axis()
    call build_axis(4,1.d0,[-1.d0],[2.d0],[1],faces,centers,widths)
    if (lbound(faces,1) /= 0 .or. ubound(faces,1) /= 4) error stop 'uniform face bounds'
    if (lbound(centers,1) /= 1 .or. size(centers) /= 4) error stop 'uniform cell count'
    if (maxval(abs(faces-[0.d0,0.25d0,0.5d0,0.75d0,1.d0])) > tolerance) error stop 'uniform faces'
    if (maxval(abs(centers-[0.125d0,0.375d0,0.625d0,0.875d0])) > tolerance) error stop 'uniform centers'
    if (maxval(abs(widths-0.25d0)) > tolerance) error stop 'uniform widths'
  end subroutine check_uniform_axis

  subroutine check_local_refinement()
    call build_axis(82,8.2d0,[1.8d0,0.8d0],[3.d0,1.8d0],[3,2],faces,centers,widths)
    if (ubound(faces,1) /= 116 .or. size(widths) /= 116) error stop 'local refinement cell count'
    if (any(widths <= 0.d0) .or. abs(sum(widths)-8.2d0) > tolerance) error stop 'invalid axial partition'
    if (maxval(abs(widths(9:28)-0.05d0)) > tolerance) error stop 'first refinement width'
    if (maxval(abs(widths(29:64)-0.1d0/3.d0)) > tolerance) error stop 'second refinement width'
    if (maxval(abs(widths(65:)-0.1d0)) > tolerance) error stop 'width outside refinement'
    do i = 0, 82
      if (.not. any(abs(faces-0.1d0*i) <= tolerance)) error stop 'missing baseline face'
    end do
    if (abs(faces(8)-0.8d0) > tolerance .or. abs(faces(28)-1.8d0) > tolerance) error stop 'axial body faces'

    call build_axis(50,5.d0,[0.d0],[0.8d0],[4],radial_faces,radial_centers,radial_widths)
    if (ubound(radial_faces,1) /= 74) error stop 'radial refinement cell count'
    if (abs(radial_faces(20)-0.5d0) > tolerance) error stop 'radial body face'
    annular_area = 0.d0
    body_volume = 0.d0
    do j = 1, size(radial_widths)
      annular_area = annular_area + pi*(radial_faces(j)**2-radial_faces(j-1)**2)
      if (j <= 20) body_volume = body_volume + pi*(radial_faces(j)**2-radial_faces(j-1)**2)*sum(widths(9:28))
    end do
    if (abs(annular_area-25.d0*pi) > tolerance*annular_area) error stop 'annuli do not cover full section'
    if (abs(body_volume-0.25d0*pi) > tolerance) error stop 'refinement changes body volume'
  end subroutine check_local_refinement

  subroutine check_reconstruction()
    real(8) :: linear(size(centers)), slopes(size(centers)-2)
    integer :: n

    n = size(centers)
    linear = 3.d0 - 0.25d0*centers
    slopes = limited_slope(linear(1:n-2),linear(2:n-1),linear(3:n), &
                           centers(2:n-1)-centers(1:n-2),centers(3:n)-centers(2:n-1))
    if (maxval(abs(slopes+0.25d0)) > 1.d-13) error stop 'linear gradient across nonuniform cells'

    ! The middle cell spans [-0.5, 0.5] and has very small neighboring cells.
    slope = limited_slope(0.d0,1.d0,100.d0,0.51d0,0.51d0)
    left_state = 1.d0-0.5d0*slope
    right_state = 1.d0+0.5d0*slope
    if (left_state < 0.d0 .or. left_state > 1.d0) error stop 'increasing reconstruction left bound'
    if (right_state < 1.d0 .or. right_state > 100.d0) error stop 'increasing reconstruction right bound'
    slope = limited_slope(100.d0,1.d0,0.d0,0.51d0,0.51d0)
    left_state = 1.d0-0.5d0*slope
    right_state = 1.d0+0.5d0*slope
    if (left_state < 1.d0 .or. left_state > 100.d0) error stop 'decreasing reconstruction left bound'
    if (right_state < 0.d0 .or. right_state > 1.d0) error stop 'decreasing reconstruction right bound'
    if (abs(limited_slope(0.d0,2.d0,1.d0,1.d0,2.d0)) > 0.d0) error stop 'local maximum slope must vanish'

    ! Van Leer keeps linear data exact, exceeds minmod on uneven smooth data
    ! and keeps the faces of a cell beside a much narrower cell bounded.
    slopes = limited_slope(linear(1:n-2),linear(2:n-1),linear(3:n), &
                           centers(2:n-1)-centers(1:n-2),centers(3:n)-centers(2:n-1),2, &
                           centers(2:n-1)-faces(1:n-2),faces(2:n-1)-centers(2:n-1))
    if (maxval(abs(slopes+0.25d0)) > 1.d-13) error stop 'van Leer linear gradient across nonuniform cells'
    slope = limited_slope(0.d0,1.d0,4.d0,1.d0,1.d0,2,0.5d0,0.5d0)
    if (abs(slope-1.5d0) > 1.d-14) error stop 'van Leer harmonic slope'
    if (abs(limited_slope(0.d0,2.d0,1.d0,1.d0,2.d0,2,0.5d0,0.5d0)) > 0.d0) error stop 'van Leer extremum'
    slope = limited_slope(0.d0,1.d0,100.d0,0.51d0,0.51d0,2,0.5d0,0.5d0)
    if (1.d0+0.5d0*slope > 100.d0 .or. 1.d0-0.5d0*slope < 0.d0) error stop 'van Leer face bound'
  end subroutine check_reconstruction

  subroutine check_flight_coefficient()
    if (abs(flight_slope_fraction(0.d0)-0.5d0) > tolerance) error stop 'collisionless free-flight limit'
    if (abs(flight_slope_fraction(1.d-12)-(0.5d0-1.d-12/3.d0)) > epsilon(1.d0)) &
      error stop 'small relaxation argument loses precision'
    if (abs(flight_slope_fraction(0.1d0)-0.46788401604444695d0) > tolerance) error stop 'intermediate series value'
    if (abs(flight_slope_fraction(1.d0)-0.26424111765711536d0) > tolerance) error stop 'finite free-flight fraction'
    if (abs(flight_slope_fraction(100.d0)-1.d-4) > 1.d-18) error stop 'stiff relaxation limit'
    if (abs(flight_slope_fraction(1.d200)) > tiny(1.d0)) error stop 'very stiff relaxation overflows'
  end subroutine check_flight_coefficient
end program test_axisymmetric_mesh

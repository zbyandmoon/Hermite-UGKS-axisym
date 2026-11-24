module axisymmetric_mesh
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  private
  public :: build_axis, limited_slope, flight_slope_fraction

contains

  subroutine build_axis(n_base,length,ref_lo,ref_hi,ref_split,faces,centers,widths)
    integer, intent(in) :: n_base, ref_split(:)
    real(8), intent(in) :: length, ref_lo(:), ref_hi(:)
    real(8), allocatable, intent(out) :: faces(:), centers(:), widths(:)
    integer, allocatable :: split(:)
    integer(int64) :: cell_count
    integer :: patch, first_face, last_face, base_cell, part, face, n
    real(8) :: spacing, tolerance, lower, upper

    if (n_base < 1 .or. .not. ieee_is_finite(length)) error stop 'invalid axis size'
    if (length <= 0.d0) error stop 'invalid axis size'
    if (size(ref_lo) /= size(ref_split) .or. size(ref_hi) /= size(ref_split)) &
      error stop 'axis refinement arrays must have the same length'

    spacing = length/n_base
    if (spacing <= 0.d0) error stop 'axis spacing must be positive'
    tolerance = 16.d0*epsilon(length)*length
    allocate(split(n_base))
    split = 1
    do patch = 1, size(ref_split)
      if (ref_split(patch) < 1) error stop 'axis refinement subdivision must be positive'
      if (ref_split(patch) == 1) cycle
      if (.not. ieee_is_finite(ref_lo(patch)) .or. .not. ieee_is_finite(ref_hi(patch))) &
        error stop 'axis refinement bounds must be finite'
      if (ref_lo(patch) < -tolerance .or. ref_lo(patch) >= ref_hi(patch)) &
        error stop 'invalid axis refinement bounds'
      if (ref_hi(patch) > length) then
        if (ref_hi(patch)-length > tolerance) error stop 'axis refinement outside domain'
      end if

      first_face = nint(max(0.d0,ref_lo(patch))/spacing)
      last_face = nint(min(length,ref_hi(patch))/spacing)
      if (first_face >= last_face) error stop 'axis refinement must span a baseline cell'
      if (abs(ref_lo(patch)-first_face*spacing) > tolerance .or. &
          abs(ref_hi(patch)-last_face*spacing) > tolerance) error stop 'axis refinement bounds must be baseline faces'
      if (any(split(first_face+1:last_face) /= 1)) error stop 'overlapping axis refinement patches'
      split(first_face+1:last_face) = ref_split(patch)
    end do

    cell_count = sum(int(split,int64))
    if (cell_count >= huge(n)) error stop 'refined axis has too many cells'
    n = int(cell_count)
    allocate(faces(0:n),centers(n),widths(n))
    faces(0) = 0.d0
    face = 0
    do base_cell = 1, n_base
      lower = faces(face)
      upper = base_cell*spacing
      if (base_cell == n_base) upper = length
      do part = 1, split(base_cell)-1
        faces(face+part) = lower+(upper-lower)*(dble(part)/dble(split(base_cell)))
      end do
      face = face+split(base_cell)
      faces(face) = upper
    end do

    widths = faces(1:n)-faces(0:n-1)
    if (any(widths <= 0.d0)) error stop 'axis refinement produced a nonpositive cell width'
    centers = faces(0:n-1)+0.5d0*widths
  end subroutine build_axis

  pure elemental real(8) function limited_slope(left,center,right,dist_left,dist_right, &
                                                limiter,face_left,face_right) result(slope)
    real(8), intent(in) :: left, center, right, dist_left, dist_right
    integer, intent(in), optional :: limiter
    real(8), intent(in), optional :: face_left, face_right
    real(8) :: slope_left, slope_right, bound

    slope = 0.d0
    if (min(dist_left,dist_right) <= 0.d0) return
    slope_left = (center-left)/dist_left
    slope_right = (right-center)/dist_right
    if (slope_left*slope_right <= 0.d0) return

    ! Each cell face lies inside the corresponding center-to-center interval.
    ! Minmod keeps both reconstructed faces between their adjacent values.
    slope = sign(min(abs(slope_left),abs(slope_right)),slope_left)
    if (.not. present(limiter)) return
    if (limiter /= 2 .or. .not. (present(face_left) .and. present(face_right))) return

    ! Van Leer: the harmonic mean is second order at smooth extrema-free
    ! points and reaches twice the minmod slope. Cap it so neither face
    ! value leaves the interval spanned by the adjacent cell averages.
    slope = 2.d0*slope_left*slope_right/(slope_left+slope_right)
    bound = huge(1.d0)
    if (face_left > 0.d0) bound = min(bound,abs(center-left)/face_left)
    if (face_right > 0.d0) bound = min(bound,abs(right-center)/face_right)
    slope = sign(min(abs(slope),bound),slope_left)
  end function limited_slope

  pure real(8) function flight_slope_fraction(z) result(fraction)
    real(8), intent(in) :: z
    ! Integral of s*exp(-z*s) for 0 <= s <= 1. Near zero its kth Taylor
    ! coefficient is (-1)**k / (k! * (k+2)).
    real(8), parameter :: series(0:11) = [0.5d0,-1.d0/3.d0,1.d0/8.d0,-1.d0/30.d0,1.d0/144.d0, &
                                        -1.d0/840.d0,1.d0/5760.d0,-1.d0/45360.d0,1.d0/403200.d0, &
                                        -1.d0/3991680.d0,1.d0/43545600.d0,-1.d0/518918400.d0]
    integer :: k

    if (abs(z) <= 0.25d0) then
      fraction = series(11)
      do k = 10, 0, -1
        fraction = series(k)+z*fraction
      end do
    else if (z >= 50.d0) then
      fraction = (1.d0/z)/z
    else
      fraction = ((1.d0-(1.d0+z)*exp(-z))/z)/z
    end if
  end function flight_slope_fraction

end module axisymmetric_mesh

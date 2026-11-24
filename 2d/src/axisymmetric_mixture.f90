submodule (axisymmetric_plume) axisymmetric_mixture
  use axisymmetric_mixture_collision, only: mixture_species, mixture_model, make_mixture_model, mixture_targets, retune_target
  use axisymmetric_mesh, only: build_axis, limited_slope, flight_slope_fraction
  use axisymmetric_internal_modes, only: vibration_energy, vibration_target
  use axisymmetric_heat_flux, only: apply_heat_relaxation
  implicit none

  real(8), parameter :: boltzmann = 1.380649d-23, avogadro = 6.02214076d23
  integer :: ng, np, active_modes
  type(mixture_species), allocatable :: species(:)
  type(mixture_model) :: mixture
  type(velocity_ring), allocatable :: grid(:)
  real(8), allocatable :: population(:,:,:,:,:,:), x_flux(:,:,:,:,:,:), r_flux(:,:,:,:,:,:)
  real(8), allocatable :: stream_x(:,:,:), stream_r(:,:,:), angle_lo(:,:), angle_hi(:,:)
  real(8), allocatable :: inlet_population(:,:,:,:,:), cold_shape(:,:,:,:), hot_shape(:,:,:,:)
  real(8), allocatable :: cell_moments(:,:,:,:), mean_moments(:,:,:,:), annulus(:), radial_center(:)
  real(8), allocatable :: axial_faces(:), radial_faces(:), axial_center(:), axial_width(:), radial_width(:), radial_centroid(:)
  real(8), allocatable :: species_vib_probability(:), species_vib_bands(:,:)

contains

  module procedure run_mixture
    real(8) :: scale,t,t_end,t_average,interval,sampled,dt_limit,loss,fallback,peak_fallback,shortest_axial
    real(8), allocatable :: initial(:), supplied(:), current(:), mass_inlet(:), wall(:,:)
    logical :: fits
    integer :: s,i,j,step,unit

    call load_species()
    call make_mixture_model(mixture,species,gas_r,temperature_ref,density_ref,vhs_reference_temperature, &
                            correct_vhs_drag=transport_matching)
    call initialize_mixture_mesh()
    if (axial_half_range) then
      if (abs(u_center) > epsilon(1.d0)) error stop 'half-range axial quadrature requires u_center=0'
      call make_velocity_ring_half_range(q,nu,ns,na,u_scale,s_max)
    else
      call make_velocity_ring(q,nu,ns,na,u_scale,u_center,s_max)
    end if
    ng = mixture%n
    np = q%np
    active_modes = merge(3,2,enable_vibration)
    allocate(grid(ng),stream_x(nu,np,ng),stream_r(nu,np,ng),angle_lo(np,ng),angle_hi(np,ng))
    allocate(cold_shape(nu,np,active_modes,ng),hot_shape(nu,np,active_modes,ng), &
             inlet_population(nu,np,active_modes,ng,nr))
    do s = 1,ng
      scale = 1.d0/sqrt(mixture%mass(s))
      grid(s) = q
      grid(s)%u = u_center+scale*(q%u-u_center)
      grid(s)%u2 = grid(s)%u**2
      grid(s)%v = scale*q%v
      grid(s)%s2 = scale*scale*q%s2
      grid(s)%log_wu = q%log_wu+log(scale)
      grid(s)%log_wp = q%log_wp+2.d0*log(scale)
      stream_x(:,:,s) = spread(grid(s)%u,2,np)
      stream_r(:,:,s) = spread(grid(s)%v,1,nu)
      angle_lo(:,s) = sqrt(grid(s)%s2)*q%sin_lower/q%angle_width
      angle_hi(:,s) = sqrt(grid(s)%s2)*q%sin_upper/q%angle_width
      call fit_species(s,thermal_state(s,1.d0,0.d0,0.d0,cold_wall/temperature_ref),cold_shape(:,:,:,s),fits)
      if (.not. fits) error stop 'mixture cold wall requires a better resolved velocity grid'
      call fit_species(s,thermal_state(s,1.d0,0.d0,0.d0,aft_wall/temperature_ref),hot_shape(:,:,:,s),fits)
      if (.not. fits) error stop 'mixture aft wall requires a better resolved velocity grid'
    end do
    allocate(mass_inlet(ng),initial(ng),supplied(ng),current(ng),wall(2,nr))
    inlet_population = 0.d0
    if (.not. uniform_test) call load_mixture_inlet()
    mass_inlet = 0.d0
    do s = 1,ng
      do j = 1,nr
        mass_inlet(s) = mass_inlet(s)+density_ref*annulus(j)* &
                     sum(max(stream_x(:,:,s),0.d0)*inlet_population(:,:,1,s,j))
      end do
      write(*,'(a,a,es17.8)') 'Nozzle mass flux kg/s ',trim(species(s)%name)//' ',c0*mass_inlet(s)
    end do
    write(*,'(a,es17.8)') 'Prescribed positive nozzle mass flux, kg/s: ',c0*sum(mass_inlet)

    allocate(population(nu,np,active_modes,ng,nx,nr),x_flux(nu,np,active_modes,ng,0:nx,nr), &
             r_flux(nu,np,active_modes,ng,nx,0:nr))
    allocate(cell_moments(6,ng,nx,nr),mean_moments(6,ng,nx,nr))
    !$omp parallel do collapse(2) private(s)
    do j = 1,nr
      do i = 1,nx
        population(:,:,:,:,i,j) = 0.d0
        if (fluid(i,j) .and. rho_initial > 0.d0) then
          do s = 1,ng
            population(:,:,:,s,i,j) = (rho_initial/density_ref)*species(s)%mass_fraction*cold_shape(:,:,:,s)
          end do
        end if
        call population_moments(population(:,:,:,:,i,j),cell_moments(:,:,i,j))
      end do
    end do
    !$omp end parallel do
    mean_moments = 0.d0
    wall = 0.d0
    initial = mixture_mass()
    supplied = 0.d0
    t = 0.d0
    sampled = 0.d0
    step = 0
    peak_fallback = 0.d0
    t_end = c0*end_time
    t_average = c0*average_start
    ! Bound the sum of outgoing physical and angular transport rates for every bin.
    shortest_axial = minval(axial_width)
    loss = 0.d0
    do s = 1,ng
      do j = 1,nr
        loss = max(loss,maxval(abs(stream_x(:,:,s))/shortest_axial+spread(angle_lo(:,s),1,nu)/radial_center(j) &
          +(max(stream_r(:,:,s),0.d0)*radial_faces(j)-min(stream_r(:,:,s),0.d0)*radial_faces(j-1)) / &
            (radial_center(j)*radial_width(j))))
      end do
    end do
    dt_limit = cfl/loss
    if (timestep_limit_s > 0.d0) dt_limit = min(dt_limit,c0*timestep_limit_s)
    write(*,'(a,5i6,a,es14.6)') 'Axisymmetric mixture UGKS, nx nr nu nperp species: ',nx,nr,nu,np,ng, &
                               '  dt(s): ',dt_limit/c0
    write(*,'(a,i2,4(a,l1))') 'Physical reconstruction order ',spatial_order,' half-range axial ',axial_half_range, &
      ' transport matching ',transport_matching,' variable rotation ',mixture%variable_rotation,' vibration ',enable_vibration
    write(*,'(a,i2,a,l1)') 'Slope limiter (1 minmod, 2 van Leer) ',slope_limiter,'  open far boundary ',open_far_boundary
    write(*,'(a)') 'step time_s mass_kg max_species_balance_kg collision_fallback_kg'
    flush(output_unit)

    do while (t < t_end)
      dt = min(dt_limit,t_end-t)
      step = step+1
      call advance_mixture(fallback)
      peak_fallback = max(peak_fallback,fallback)
      supplied = supplied + dt*mixture_boundary_flow(mass_inlet)
      interval = max(0.d0,t+dt-max(t,t_average))
      if (interval > 0.d0) then
        call accumulate_mixture(interval,wall)
        sampled = sampled+interval
      end if
      t = t+dt
      if (mod(step,log_stride) == 0 .or. t >= t_end) then
        current = mixture_mass()
        if (.not. all(ieee_is_finite(cell_moments))) error stop 'nonfinite axisymmetric mixture solution'
        write(*,'(i8,4es18.9)') step,t/c0,sum(current),maxval(abs(current-initial-supplied)),fallback
        flush(output_unit)
      end if
    end do
    call write_mixture('axisymmetric_mixture_final.dat',cell_moments)
    if (sampled > 0.d0) then
      call write_mixture('axisymmetric_mixture_average.dat',mean_moments/sampled)
      open(newunit=unit,file='axisymmetric_mixture_wall.dat',status='replace',action='write')
      write(unit,'(a)') '# r_m normal_pressure_Pa total_heat_flux_W_m2'
      do j = 1,nr
        write(unit,'(3(es24.15e3,1x))') radial_center(j),wall(:,j)/sampled
      end do
      close(unit)
    end if
    current = mixture_mass()
    open(newunit=unit,file='axisymmetric_mixture_mass_balance.dat',status='replace',action='write')
    write(unit,'(a)') '# species initial_kg net_boundary_kg final_kg residual_kg'
    do s = 1,ng
      write(unit,'(a,4(1x,es24.15e3))') trim(species(s)%name),initial(s),supplied(s),current(s), &
                                      current(s)-initial(s)-supplied(s)
    end do
    close(unit)
    write(*,'(a,es16.7,a,es16.7)') 'Completed axisymmetric mixture run. Averaging seconds: ',sampled/c0, &
      '   Peak mass using collision fallback, kg: ',peak_fallback
  end procedure run_mixture

  subroutine initialize_mixture_mesh()
    integer :: i_front,i_aft,j_body,unit

    if (spatial_order < 1 .or. spatial_order > 2) error stop 'spatial_order must be 1 or 2'
    if (slope_limiter < 1 .or. slope_limiter > 2) error stop 'slope_limiter must be 1 (minmod) or 2 (van Leer)'
    if (.not. ieee_is_finite(timestep_limit_s) .or. timestep_limit_s < 0.d0) error stop 'invalid timestep limit'
    if (.not. ieee_is_finite(translational_prandtl) .or. translational_prandtl <= 0.d0 .or. &
        translational_prandtl > 1.d0) error stop 'translational_prandtl must be in (0,1]'
    call build_axis(nx,length,x_ref_lo,x_ref_hi,x_ref_split,axial_faces,axial_center,axial_width)
    call build_axis(nr,radius,r_ref_lo,r_ref_hi,r_ref_split,radial_faces,radial_center,radial_width)
    nx = size(axial_width)
    nr = size(radial_width)
    annulus = pi8*(radial_faces(1:nr)+radial_faces(0:nr-1))*radial_width
    radial_centroid = (2.d0/3.d0)*(radial_faces(1:nr)**2+radial_faces(1:nr)*radial_faces(0:nr-1) &
                     +radial_faces(0:nr-1)**2)/(radial_faces(1:nr)+radial_faces(0:nr-1))
    allocate(fluid(nx,nr))
    fluid = .true.
    if (.not. uniform_test) then
      i_front = minloc(abs(axial_faces-x_front),dim=1)-1
      i_aft = minloc(abs(axial_faces-x_aft),dim=1)-1
      j_body = minloc(abs(radial_faces-r_body),dim=1)-1
      if (i_front < 1 .or. i_front >= i_aft .or. i_aft >= nx .or. j_body < 1 .or. j_body >= nr) &
        error stop 'lander outside the axisymmetric mixture grid'
      if (max(abs(x_front-axial_faces(i_front)),abs(x_aft-axial_faces(i_aft)),abs(r_body-radial_faces(j_body))) &
          > 1.d-9) error stop 'choose an axisymmetric mixture mesh with the lander walls on faces'
      fluid(i_front+1:i_aft,1:j_body) = .false.
    end if
    open(newunit=unit,file='axisymmetric_x_faces.dat',status='replace',action='write')
    write(unit,'(es24.15e3)') axial_faces
    close(unit)
    open(newunit=unit,file='axisymmetric_r_faces.dat',status='replace',action='write')
    write(unit,'(es24.15e3)') radial_faces
    close(unit)
  end subroutine initialize_mixture_mesh

  subroutine load_species()
    integer :: unit,io,count,k,configured
    character(512) :: line

    open(newunit=unit,file=species_file,status='old',action='read')
    count = 0
    do
      read(unit,'(a)',iostat=io) line
      if (io < 0) exit
      if (io > 0) error stop 'cannot read mixture species file'
      line = adjustl(line)
      if (len_trim(line) > 0 .and. line(1:1) /= '#') count = count+1
    end do
    if (count < 1) error stop 'empty mixture species file'
    if (.not. all(ieee_is_finite(vibration_probability)) .or. .not. all(ieee_is_finite(vibration_bands)) .or. &
        any(vibration_probability < 0.d0) .or. any(vibration_probability > 1.d0) .or. any(vibration_bands < 0.d0)) &
      error stop 'invalid vibration parameters'
    configured = min(count,mode_species_limit)
    if (enable_vibration .and. .not. any(vibration_bands(:,1:configured) > 0.d0)) error stop 'no active vibration bands'
    allocate(species(count))
    if (enable_vibration) then
      allocate(species_vib_probability(count),species_vib_bands(vibration_band_limit,count))
      species_vib_probability = 0.d0
      species_vib_bands = 0.d0
      species_vib_probability(1:configured) = vibration_probability(1:configured)
      species_vib_bands(:,1:configured) = vibration_bands(:,1:configured)
    end if
    rewind(unit)
    k = 0
    do while (k < count)
      read(unit,'(a)') line
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      k = k+1
      read(line,*,iostat=io) species(k)%name,species(k)%mass_kg,species(k)%mass_fraction, &
           species(k)%diameter_m,species(k)%omega,species(k)%rot_dof,species(k)%rot_probability
      if (io /= 0) error stop 'invalid mixture species row'
      species(k)%mass_kg = species(k)%mass_kg*(1.d-3/avogadro)
      if (k <= configured) species(k)%rotation_tstar = rotation_tstar(k)
      if (enable_vibration) then
        if (species_vib_probability(k) > 0.d0 .and. .not. any(species_vib_bands(:,k) > 0.d0)) &
          error stop 'missing vibration bands for a relaxing species'
      end if
      if (any(species(1:k-1)%name == species(k)%name)) error stop 'duplicate mixture species name'
    end do
    close(unit)
  end subroutine load_species

  pure function thermal_state(s,rho,u,v,temperature) result(w)
    integer, intent(in) :: s
    real(8), intent(in) :: rho,u,v,temperature
    real(8) :: w(6), thermal_scale
    thermal_scale = temperature/mixture%mass(s)
    w(1:5) = rho*[1.d0,u,v,u*u+v*v+1.5d0*thermal_scale,0.25d0*species(s)%rot_dof*thermal_scale]
    w(6) = 0.d0
    if (enable_vibration) w(6) = rho*sum(vibration_energy(species_vib_bands(:,s),temperature*temperature_ref, &
                                                       gas_r/mixture%mass(s)))/c0**2
  end function thermal_state

  subroutine load_mixture_inlet()
    integer :: unit,io,s,j
    real(8) :: lo,hi,a,b,rho,temperature,ux,ur,area_weight,panel(nu,np,active_modes)
    logical :: fits
    character(512) :: header

    open(newunit=unit,file=inlet_file,status='old',action='read')
    read(unit,'(a)') header
    do
      read(unit,*,iostat=io) lo,hi,rho,temperature,ux,ur
      if (io < 0) exit
      if (io /= 0 .or. lo < 0.d0 .or. hi > r_body .or. hi <= lo .or. min(rho,temperature) <= 0.d0) &
        error stop 'invalid axisymmetric mixture inlet panel'
      do s = 1,ng
        call fit_species(s,thermal_state(s,rho*species(s)%mass_fraction/density_ref,ux/c0,ur/c0, &
                         temperature/temperature_ref),panel,fits)
        if (.not. fits) error stop 'mixture nozzle requires a better resolved velocity grid'
        do j = 1,nr
          a = max(lo,radial_faces(j-1))
          b = min(hi,radial_faces(j))
          if (b <= a) cycle
          area_weight = pi8*(b*b-a*a)/annulus(j)
          inlet_population(:,:,:,s,j) = inlet_population(:,:,:,s,j)+area_weight*panel
        end do
      end do
    end do
    close(unit)
  end subroutine load_mixture_inlet

  pure subroutine fit_species(s,w,result,fits)
    integer, intent(in) :: s
    real(8), intent(in) :: w(6)
    real(8), intent(out) :: result(nu,np,active_modes)
    logical, intent(out) :: fits
    if (w(1) <= 0.d0) then
      result = 0.d0
      fits = .true.
      return
    end if
    call ring_equilibrium(grid(s),w(1:4),result(:,:,1),fits)
    result(:,:,2) = (w(5)/w(1))*result(:,:,1)
    if (enable_vibration) result(:,:,3) = (w(6)/w(1))*result(:,:,1)
  end subroutine fit_species

  pure subroutine population_moments(values,w)
    real(8), intent(in) :: values(nu,np,active_modes,ng)
    real(8), intent(out) :: w(6,ng)
    integer :: s
    do s = 1,ng
      w(1:4,s) = ring_moments(grid(s),values(:,:,1,s))
      w(5,s) = sum(values(:,:,2,s))
      w(6,s) = 0.d0
      if (enable_vibration) w(6,s) = sum(values(:,:,3,s))
    end do
  end subroutine population_moments

  subroutine advance_mixture(fallback)
    real(8), intent(out) :: fallback
    integer :: i,j
    fallback = 0.d0
    !$omp parallel default(shared) private(i,j) reduction(+:fallback)
    !$omp do collapse(2) schedule(dynamic,16)
    do j = 1,nr
      do i = 0,nx
        call mixture_x_face(i,j)
      end do
    end do
    !$omp end do nowait
    !$omp do collapse(2) schedule(dynamic,16)
    do j = 1,nr
      do i = 1,nx
        call mixture_r_face(i,j)
      end do
    end do
    !$omp end do
    !$omp do collapse(2) schedule(dynamic,16)
    do j = 1,nr
      do i = 1,nx
        if (fluid(i,j)) call mixture_cell(i,j,fallback)
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine advance_mixture

  subroutine mixture_x_face(i,j)
    integer, intent(in) :: i,j
    integer :: s,mode
    if (i == 0) then
      if (uniform_test) then
        call mixture_reflect(population(:,:,:,:,1,j),cold_shape,stream_x,1,x_flux(:,:,:,:,i,j))
      else
        do s = 1,ng
          do mode = 1,active_modes
            x_flux(:,:,mode,s,i,j) = min(0.d0,stream_x(:,:,s))*population(:,:,mode,s,1,j)
          end do
        end do
      end if
    else if (i == nx) then
      if (open_far_boundary .and. .not. uniform_test) then
        do s = 1,ng
          do mode = 1,active_modes
            x_flux(:,:,mode,s,i,j) = max(0.d0,stream_x(:,:,s))*population(:,:,mode,s,nx,j)
          end do
        end do
      else
        call mixture_reflect(population(:,:,:,:,nx,j),cold_shape,stream_x,-1,x_flux(:,:,:,:,i,j))
      end if
    else if (fluid(i,j) .and. fluid(i+1,j)) then
      if (spatial_order == 2) then
        call reconstructed_flux(i,j,1,0,axial_center,axial_faces,stream_x,x_flux(:,:,:,:,i,j))
      else
        call mixture_flux(population(:,:,:,:,i,j),population(:,:,:,:,i+1,j),stream_x,x_flux(:,:,:,:,i,j))
      end if
    else if (fluid(i,j)) then
      call mixture_reflect(population(:,:,:,:,i,j),cold_shape,stream_x,-1,x_flux(:,:,:,:,i,j))
    else if (fluid(i+1,j)) then
      call mixture_reflect(population(:,:,:,:,i+1,j),hot_shape,stream_x,1,x_flux(:,:,:,:,i,j))
      do s = 1,ng
        do mode = 1,active_modes
          x_flux(:,:,mode,s,i,j) = x_flux(:,:,mode,s,i,j)+max(0.d0,stream_x(:,:,s))*inlet_population(:,:,mode,s,j)
        end do
      end do
    else
      x_flux(:,:,:,:,i,j) = 0.d0
    end if
  end subroutine mixture_x_face

  subroutine mixture_r_face(i,j)
    integer, intent(in) :: i,j
    integer :: s,mode
    if (j == nr) then
      if (uniform_test) then
        call mixture_reflect(population(:,:,:,:,i,j),cold_shape,stream_r,-1,r_flux(:,:,:,:,i,j))
      else
        do s = 1,ng
          do mode = 1,active_modes
            r_flux(:,:,mode,s,i,j) = max(0.d0,stream_r(:,:,s))*population(:,:,mode,s,i,j)
          end do
        end do
      end if
    else if (fluid(i,j) .and. fluid(i,j+1)) then
      if (spatial_order == 2) then
        call reconstructed_flux(i,j,0,1,radial_centroid,radial_faces,stream_r,r_flux(:,:,:,:,i,j))
      else
        call mixture_flux(population(:,:,:,:,i,j),population(:,:,:,:,i,j+1),stream_r,r_flux(:,:,:,:,i,j))
      end if
    else if (fluid(i,j+1)) then
      call mixture_reflect(population(:,:,:,:,i,j+1),cold_shape,stream_r,1,r_flux(:,:,:,:,i,j))
    else
      r_flux(:,:,:,:,i,j) = 0.d0
    end if
  end subroutine mixture_r_face

  pure subroutine mixture_reflect(cell,shape,speed,normal,flux)
    real(8), intent(in) :: cell(nu,np,active_modes,ng),shape(nu,np,active_modes,ng),speed(nu,np,ng)
    integer, intent(in) :: normal
    real(8), intent(out) :: flux(nu,np,active_modes,ng)
    real(8) :: factor
    integer :: s,mode
    do s = 1,ng
      factor = -sum(min(normal*speed(:,:,s),0.d0)*cell(:,:,1,s))/ &
                sum(max(normal*speed(:,:,s),0.d0)*shape(:,:,1,s))
      do mode = 1,active_modes
        where (normal*speed(:,:,s) < 0.d0)
          flux(:,:,mode,s) = speed(:,:,s)*cell(:,:,mode,s)
        elsewhere
          flux(:,:,mode,s) = factor*speed(:,:,s)*shape(:,:,mode,s)
        end where
      end do
    end do
  end subroutine mixture_reflect

  subroutine reconstructed_flux(i,j,di,dj,centers,faces,speed,flux)
    integer, intent(in) :: i,j,di,dj
    real(8), intent(in) :: centers(:),faces(0:),speed(nu,np,ng)
    real(8), intent(out) :: flux(nu,np,active_modes,ng)
    real(8) :: face_values(nu,np,active_modes,ng), gradient(nu,np,active_modes,ng), &
               left_slope(nu,np,active_modes,ng), right_slope(nu,np,active_modes,ng)
    integer :: face,s,mode

    face = i*di+j*dj
    call local_gradient(i,j,face,left_slope)
    call local_gradient(i+di,j+dj,face+1,right_slope)
    do s = 1,ng
      do mode = 1,active_modes
        where (speed(:,:,s) > 0.d0)
          gradient(:,:,mode,s) = left_slope(:,:,mode,s)
          face_values(:,:,mode,s) = population(:,:,mode,s,i,j)+(faces(face)-centers(face))*left_slope(:,:,mode,s)
        elsewhere
          gradient(:,:,mode,s) = right_slope(:,:,mode,s)
          face_values(:,:,mode,s) = population(:,:,mode,s,i+di,j+dj) &
                                     +(faces(face)-centers(face+1))*right_slope(:,:,mode,s)
        end where
      end do
    end do
    call mixture_flux(population(:,:,:,:,i,j),population(:,:,:,:,i+di,j+dj),speed,flux,face_values,gradient)

  contains
    subroutine local_gradient(ic,jc,k,slope)
      integer, intent(in) :: ic,jc,k
      real(8), intent(out) :: slope(nu,np,active_modes,ng)
      slope = 0.d0
      if (k <= 1 .or. k >= size(centers)) return
      if (.not. (fluid(ic-di,jc-dj) .and. fluid(ic+di,jc+dj))) return
      slope = limited_slope(population(:,:,:,:,ic-di,jc-dj),population(:,:,:,:,ic,jc), &
                  population(:,:,:,:,ic+di,jc+dj),centers(k)-centers(k-1),centers(k+1)-centers(k), &
                  slope_limiter,centers(k)-faces(k-1),faces(k)-centers(k))
    end subroutine local_gradient
  end subroutine reconstructed_flux

  subroutine mixture_flux(left,right,speed,flux,face_values,gradient)
    real(8), intent(in) :: left(nu,np,active_modes,ng),right(nu,np,active_modes,ng),speed(nu,np,ng)
    real(8), intent(out) :: flux(nu,np,active_modes,ng)
    real(8), intent(in), optional :: face_values(nu,np,active_modes,ng),gradient(nu,np,active_modes,ng)
    real(8) :: upwind(nu,np,active_modes,ng),equilibrium(nu,np,active_modes,ng),w(6,ng),target(6,ng), &
               applied(ng),frequency(ng),encounters(ng),slope_weight,cap
    integer :: s,mode
    logical :: all_fit

    do s = 1,ng
      do mode = 1,active_modes
        where (speed(:,:,s) > 0.d0)
          upwind(:,:,mode,s) = left(:,:,mode,s)
        elsewhere
          upwind(:,:,mode,s) = right(:,:,mode,s)
        end where
      end do
    end do
    if (present(face_values)) then
      call population_moments(face_values,w)
    else
      call population_moments(upwind,w)
    end if
    call fitted_collision(w,0.d0,target,frequency,encounters,applied,equilibrium,all_fit)

    ! Integrate the upwind and frozen equilibrium contributions along the
    ! characteristics. Bound outgoing bins using the donor cell averages.
    cap = 0.95d0/cfl
    do s = 1,ng
      if (.not. all_fit) then
        applied(s) = 0.d0
        frequency(s) = 0.d0
      end if
      slope_weight = 0.d0
      if (present(gradient)) slope_weight = dt*flight_slope_fraction(dt*frequency(s))
      if (applied(s) > 0.d0 .and. translational_prandtl < 1.d0) then
        if (present(face_values)) then
          call apply_heat_relaxation(grid(s),target(1:4,s),face_values(:,:,1,s), &
                                     equilibrium(:,:,1,s),translational_prandtl)
        else
          call apply_heat_relaxation(grid(s),target(1:4,s),upwind(:,:,1,s),equilibrium(:,:,1,s),translational_prandtl)
        end if
      end if
      do mode = 1,active_modes
        if (present(face_values)) then
          flux(:,:,mode,s) = face_values(:,:,mode,s)
        else
          flux(:,:,mode,s) = upwind(:,:,mode,s)
        end if
        if (applied(s) > 0.d0) flux(:,:,mode,s) = flux(:,:,mode,s) &
                                    +applied(s)*(equilibrium(:,:,mode,s)-flux(:,:,mode,s))
        if (present(gradient)) flux(:,:,mode,s) = flux(:,:,mode,s)-slope_weight*speed(:,:,s)*gradient(:,:,mode,s)
        flux(:,:,mode,s) = speed(:,:,s)*min(max(flux(:,:,mode,s),0.d0),cap*upwind(:,:,mode,s))
      end do
    end do
  end subroutine mixture_flux

  subroutine fitted_collision(w,cell_dt,target,frequency,encounters,applied,equilibrium,all_fit)
    real(8), intent(in) :: w(6,ng),cell_dt
    real(8), intent(out) :: target(6,ng),frequency(ng),encounters(ng),applied(ng)
    real(8), intent(out) :: equilibrium(nu,np,active_modes,ng)
    logical, intent(out) :: all_fit
    real(8) :: pair_target(5,ng),pair_frequency(ng)
    integer :: attempt,s
    logical :: fits

    call mixture_targets(mixture,w(1:5,:),cell_dt,pair_target,pair_frequency,encounters)
    do attempt = 1,merge(2,1,transport_matching)
      target(1:5,:) = pair_target
      target(6,:) = w(6,:)
      frequency = pair_frequency
      if (transport_matching .and. attempt == 1) &
        call retune_target(mixture,w(1:5,:),cell_dt,pair_target,pair_frequency,target(1:5,:),frequency)
      all_fit = .true.
      do s = 1,ng
        if (cell_dt > 0.d0) then
          applied(s) = cell_dt*frequency(s)/(1.d0+cell_dt*frequency(s))
        else
          applied(s) = 1.d0-flight_fraction(dt*frequency(s))
          if (applied(s) <= 1.d-7) applied(s) = 0.d0
        end if
        if (w(1,s) <= 1.d-30) applied(s) = 0.d0
        if (applied(s) <= 0.d0) cycle
        call fit_species(s,target(:,s),equilibrium(:,:,:,s),fits)
        all_fit = all_fit .and. fits
      end do
      if (all_fit) exit
    end do
  end subroutine fitted_collision

  subroutine mixture_cell(i,j,fallback)
    integer, intent(in) :: i,j
    real(8), intent(inout) :: fallback
    real(8) :: transported(nu,np,active_modes,ng),equilibrium(nu,np,active_modes,ng),w(6,ng),target(6,ng), &
               frequency(ng),weight(ng),encounters(ng)
    real(8) :: center,hi,lo,width,volume_factor
    integer :: s,k
    logical :: all_fit

    center = radial_center(j)
    hi = radial_faces(j)
    lo = radial_faces(j-1)
    width = radial_width(j)
    volume_factor = density_ref*annulus(j)*axial_width(i)
    transported = population(:,:,:,:,i,j)+dt*((x_flux(:,:,:,:,i-1,j)-x_flux(:,:,:,:,i,j))/axial_width(i) &
                    -(hi/(center*width))*r_flux(:,:,:,:,i,j))
    if (j > 1) transported = transported+(dt*lo/(center*width))*r_flux(:,:,:,:,i,j-1)
    do s = 1,ng
      do k = 1,np
        transported(:,k,:,s) = transported(:,k,:,s)-(dt*angle_lo(k,s)/center)*population(:,k,:,s,i,j)
        if (mod(k,na) /= 0) transported(:,k,:,s) = transported(:,k,:,s) &
                                              +(dt*angle_hi(k,s)/center)*population(:,k+1,:,s,i,j)
      end do
    end do
    if (minval(transported) < -1.d-13*maxval(abs(transported))) error stop 'negative axisymmetric mixture population'
    transported = max(transported,0.d0)
    call population_moments(transported,w)
    call fitted_collision(w,dt,target,frequency,encounters,weight,equilibrium,all_fit)
    if (all_fit) then
      do s = 1,ng
        if (weight(s) <= 0.d0) cycle
        if (translational_prandtl < 1.d0) &
          call apply_heat_relaxation(grid(s),target(1:4,s),transported(:,:,1,s),equilibrium(:,:,1,s),translational_prandtl)
        transported(:,:,:,s) = (1.d0-weight(s))*transported(:,:,:,s)+weight(s)*equilibrium(:,:,:,s)
        w(:,s) = w(:,s)+weight(s)*(target(:,s)-w(:,s))
      end do
    else
      ! Skip the entire coupled collision if one velocity grid cannot
      ! represent the target and retain the conservative transport step.
      fallback = fallback+volume_factor*dot_product(weight,w(1,:))
    end if
    if (enable_vibration) then
      do s = 1,ng
        call relax_vibration(s,encounters(s),transported(:,:,:,s),w(:,s),volume_factor,fallback)
      end do
    end if
    population(:,:,:,:,i,j) = transported
    cell_moments(:,:,i,j) = w
  end subroutine mixture_cell

  subroutine relax_vibration(s,encounter,values,w,volume_factor,fallback)
    integer, intent(in) :: s
    real(8), intent(in) :: encounter,volume_factor
    real(8), intent(inout) :: values(nu,np,active_modes),w(6),fallback
    real(8) :: target(6), equilibrium(nu,np), rate, weight, velocity_sq, translation, vibration, new_tr, new_vib
    logical :: fits

    if (w(1) <= 1.d-30) return
    rate = dt*encounter*species_vib_probability(s)
    if (rate <= 1.d-8) return
    weight = rate/(1.d0+rate)
    velocity_sq = sum((w(2:3)/w(1))**2)
    translation = 0.5d0*c0**2*max(w(4)/w(1)-velocity_sq,0.d0)
    vibration = w(6)*c0**2/w(1)
    call vibration_target(translation,vibration,species_vib_bands(:,s),gas_r/mixture%mass(s),new_tr,new_vib)
    if (abs(new_vib-vibration) <= 1.d-10*max(1.d0,translation+vibration)) return
    target = w
    target(4) = w(1)*(velocity_sq+2.d0*new_tr/c0**2)
    target(6) = w(1)*new_vib/c0**2
    call ring_equilibrium(grid(s),target(1:4),equilibrium,fits)
    if (.not. fits) then
      fallback = fallback+volume_factor*weight*w(1)
      return
    end if
    values(:,:,1) = (1.d0-weight)*values(:,:,1)+weight*equilibrium
    values(:,:,3) = (1.d0-weight)*values(:,:,3)+weight*(target(6)/w(1))*equilibrium
    w = w+weight*(target-w)
  end subroutine relax_vibration

  function mixture_mass() result(mass)
    real(8) :: mass(ng)
    integer :: j
    mass = 0.d0
    do j = 1,nr
      mass = mass+density_ref*annulus(j)*matmul(cell_moments(1,:,:,j),axial_width)
    end do
  end function mixture_mass

  function mixture_boundary_flow(inlet_rate) result(flow)
    real(8), intent(in) :: inlet_rate(ng)
    real(8) :: flow(ng)
    integer :: j,i,s
    flow = inlet_rate
    do s = 1,ng
      do j = 1,nr
        flow(s) = flow(s)+density_ref*annulus(j)*(sum(x_flux(:,:,1,s,0,j))-sum(x_flux(:,:,1,s,nx,j)))
      end do
      do i = 1,nx
        flow(s) = flow(s)-density_ref*(2.d0*pi8*radius*axial_width(i))*sum(r_flux(:,:,1,s,i,nr))
      end do
    end do
  end function mixture_boundary_flow

  subroutine accumulate_mixture(interval,wall)
    real(8), intent(in) :: interval
    real(8), intent(inout) :: wall(2,nr)
    integer :: j,s,k
    real(8) :: force,heat
    mean_moments = mean_moments+interval*cell_moments
    do j = 1,nr
      force = 0.d0
      heat = 0.d0
      do s = 1,ng
        do k = 1,np
          force = force+dot_product(grid(s)%u,x_flux(:,k,1,s,nx,j))
          heat = heat+0.5d0*dot_product(grid(s)%u2+grid(s)%s2(k),x_flux(:,k,1,s,nx,j)) &
                     +sum(x_flux(:,k,2,s,nx,j))
          if (enable_vibration) heat = heat+sum(x_flux(:,k,3,s,nx,j))
        end do
      end do
      wall(:,j) = wall(:,j)+interval*density_ref*c0**2*[force,c0*heat]
    end do
  end subroutine accumulate_mixture

  subroutine write_mixture(filename,values)
    character(*), intent(in) :: filename
    real(8), intent(in) :: values(6,ng,nx,nr)
    real(8) :: u,v,rho,number,pressure,tr,rot,rot_dof_sum,number_species(ng),vibrational(1)
    integer :: unit,i,j,s
    character(:), allocatable :: header

    header = '# x_m r_m n_m3 rho_kg_m3 ux_m_s ur_m_s Ttr_K Trot_K p_Pa'
    do s = 1,ng
      header = header//' rho_'//trim(species(s)%name)//'_kg_m3'
    end do
    if (enable_vibration) header = header//' Evib_J_m3'
    open(newunit=unit,file=filename,status='replace',action='write')
    write(unit,'(a)') header//' active'
    do j = 1,nr
      do i = 1,nx
        rho = sum(values(1,:,i,j))
        number_species = density_ref*values(1,:,i,j)/species%mass_kg
        number = sum(number_species)
        u = 0.d0
        v = 0.d0
        pressure = 0.d0
        tr = 0.d0
        rot = 0.d0
        if (rho > 1.d-30) then
          u = sum(values(2,:,i,j))/rho
          v = sum(values(3,:,i,j))/rho
          pressure = density_ref*c0*c0/3.d0*max(0.d0,sum(values(4,:,i,j))-rho*(u*u+v*v))
          tr = pressure/(number*boltzmann)
          rot_dof_sum = dot_product(number_species,dble(species%rot_dof))
          if (rot_dof_sum > 0.d0) rot = 2.d0*density_ref*c0*c0*sum(values(5,:,i,j))/(rot_dof_sum*boltzmann)
        end if
        vibrational(1) = density_ref*c0*c0*sum(values(6,:,i,j))
        write(unit,'(*(es24.15e3,1x))') axial_center(i),radial_center(j),number,density_ref*rho,u*c0,v*c0,tr,rot, &
          pressure,density_ref*values(1,:,i,j),vibrational(1:merge(1,0,enable_vibration)),merge(1.d0,0.d0,fluid(i,j))
      end do
    end do
    close(unit)
  end subroutine write_mixture
end submodule axisymmetric_mixture

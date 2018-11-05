program create_laser
  
  use specifications
  use read_parameters
  
  implicit none
  
  integer :: ilasers
  integer :: ixyz
  integer :: it
  integer :: NE
  integer :: max_polarizaton_index(1)
  
  real(kind=8) :: t
  real(kind=8) :: bandwidth_factor
  real(kind=8), allocatable :: envelope(:)
  real(kind=8), allocatable :: env(:,:)
  real(kind=8), allocatable :: momentary_frequency(:)
  real(kind=8), allocatable :: mom_freq(:,:)
  
  complex(kind=8),allocatable :: laser(:,:)
  complex(kind=8),allocatable :: laser_t(:)
  complex(kind=8),allocatable :: Wig(:,:)
  
  call read_params
  
  allocate (laser(Nt,3), STAT=allocatestatus)
  if (allocatestatus /= 0) stop "*** Not enough memory 1 ***"
  allocate (laser_t(Nt), STAT=allocatestatus)
  if (allocatestatus /= 0) stop "*** Not enough memory 2 ***"
  allocate (envelope(Nt), STAT=allocatestatus)
  if (allocatestatus /= 0) stop "*** Not enough memory 3 ***"
  allocate (env(Nt,Nlasers), STAT=allocatestatus)
  if (allocatestatus /= 0) stop "*** Not enough memory 4 ***"
  allocate (momentary_frequency(Nt), STAT=allocatestatus)
  if (allocatestatus /= 0) stop "*** Not enough memory 5 ***"
  allocate (mom_freq(Nt,Nlasers), STAT=allocatestatus)
  if (allocatestatus /= 0) stop "*** Not enough memory 6 ***"
  laser = 0.
  write(6,*) Nt
  
  do ilasers = 1, Nlasers
    call field_transform(laser_t,type_envelope(ilasers),field_strength(ilasers),fwhm(ilasers), &
                         pulse_begin(ilasers),pulse_center(ilasers),pulse_center2(ilasers), &
                         pulse_end(ilasers),omega_0(ilasers),phase(ilasers),b_1(ilasers), &
                         b_2(ilasers),b_3(ilasers),b_4(ilasers),dt,t0,Nt,envelope,momentary_frequency, &
                         E_flip(ilasers)) 
    env(:,ilasers) = sqrt(dble(laser_t(:))**2+imag(laser_t(:))**2)
    mom_freq(:,ilasers) = momentary_frequency(:)
    do ixyz = 1,3
      do it = 1,Nt
        if (abs(laser_t(it)) < threshold(ilasers)) then
          laser_t(it) = 0.
        endif
        laser(it,ixyz) = laser(it,ixyz) + polarization(ixyz,ilasers)*laser_t(it)
      enddo
    enddo
  enddo
  
  write(6,*) 'Writing out'
  
  open (10,file='laser.inp')
  do it = 1,Nt
    t = t0 + (it-1) * dt
    if (realvalued) then
      write(10,'(107(e16.8))') t*au2fs, &
                             dble(laser(it,1)), 0.d0, &
                             dble(laser(it,2)), 0.d0, &
                             dble(laser(it,3)), 0.d0, &
                             (mom_freq(it,ilasers),ilasers=1,Nlasers) 
    else
      write(10,'(107(e16.8))') t*au2fs, &
                             dble(laser(it,1)), imag(laser(it,1)), &
                             dble(laser(it,2)), imag(laser(it,2)), &
                             dble(laser(it,3)), imag(laser(it,3)), &
                             (mom_freq(it,ilasers),ilasers=1,Nlasers)
    endif 
  enddo
  close (10)
  
  open (10,file='pol.inp')
    write(10,*) Nlasers
    do ilasers = 1, Nlasers
      write(10,'(3(e16.8))') (polarization(ixyz,ilasers),ixyz=1,3)
    enddo
  close (10)

  open (10,file='env.inp')
    do it = 1,Nt
      t = t0 + (it-1) * dt
      write(10,'(100(e16.8))') t*au2fs, (env(it,ilasers),ilasers=1,Nlasers)
    enddo
  close (10)
  
  
  open (10,file='omega.dat')
    do ilasers = 1, Nlasers
      write(10,*) omega_0(ilasers)
    enddo
  close (10)
  
  if (Nlasers == 1 .and. l_wigner) then
    max_polarizaton_index = MAXLOC(polarization,DIM=1)   
    laser_t(:) = laser(:,max_polarizaton_index(1))
  
    NE = INT(Nt/2) 
    allocate(Wig(Nt,NE))
  
    call wigner(Nt,dt,t0,laser_t,NE,Wig,wigner_low_t_limit,wigner_high_t_limit,wigner_low_E_limit,wigner_high_E_limit)
  
    if (type_envelope(1)==1 .and. l_husimi) then
      bandwidth_factor = 8.*log(2.) / (fwhm(1)**2)
      call husimi(Nt,dt,t0,NE,Wig,bandwidth_factor,wigner_low_t_limit,wigner_high_t_limit,wigner_low_E_limit,wigner_high_E_limit)
    endif
  endif
  
  deallocate (polarization)
  deallocate (type_envelope)
  deallocate (field_strength)
  deallocate (fwhm)
  deallocate (pulse_begin)
  deallocate (pulse_center)
  deallocate (pulse_center2)
  deallocate (pulse_end)
  deallocate (omega_0)
  deallocate (phase)
  deallocate (b_2)
  deallocate (b_3)
  deallocate (b_4)
  deallocate (E_flip)
  deallocate (laser)
  deallocate (laser_t)
  
end



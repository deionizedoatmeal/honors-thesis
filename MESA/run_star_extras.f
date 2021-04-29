! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      !****************************************************
      ! written by ian k. bania (advised by cosmin ilie)
      ! colgate university
      ! email: ibania@colgate.edu, cilie@colgate.edu
      ! last updated dec. 1st 2020
      !****************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      
      implicit none
      ! global variables
      double precision :: DM_rate, DM_energy_rate, Mchi, frac
      double precision :: rho_chi, vbar, Hfrac
      double precision, dimension(100000) :: cell_DM_energy_rate
      double precision, dimension(100000) :: cell_DM_annihlation_rate
      logical :: sigma_xenon
      
      !****************************************************
      !
      ! array definitons:
      ! -> cell_DM_annihlation_rate(k) [1/s]
      ! array with annihlations per unit time for each cell
      ! -> cell_DM_energy_rate(k) [ergs/s]
      ! array with total energy rate for each cell
      ! -> extra_heat(k) [ergs/s/g]
      ! array with energy rate per unit mass for each cell
      !
      ! variable definitons:
      ! -> DM_energy_rate
      ! variable with total energy rate for star (L_DM) [erg/s]
      ! -> DM_rate
      ! capture rate of dark matter [1/s]
      !
      !****************************************************
      !
      ! quick flowchart/roadmap of the functions ian built
      ! extras_check_model:
      !  -> calls DM_energy
      !
      ! DM_energy:
      !  -> calls phi 
      !  -> calls phimid 
      !  -> calls phicent
      !  -> calls trap_int
      !  -> calls captureN_pureH
      !  -> assigns to extra_heat variables
      !
      ! captureN_pureH: (DM capture rate)
      !  -> calls gamma_incl 
      !
      ! gamma_incl: (lower incomplete gamma func.)
      !  -> calls trap_int 
      !
      ! phi: (gravitational potential)
      !  -> calls trap_int
      !
      ! phimid: (gravitational potential at COM of cell)
      !  -> calls trap_int
      !
      ! phicent: (gravitational potential at inside of innermost cell)
      !  -> calls trap_int
      !
      ! trap_int: (trapozoid rule numerical integral)
      !
      !****************************************************

      ! these routines are called by the standard run_star check_model
      contains
      
      subroutine extras_controls(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
               
            ! this is the place to set any procedure pointers you want to change
            ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
            ! the extras functions in this file will not be called
            ! unless you set their function pointers as done below.
            ! otherwise we use a null_ version which does nothing (except warn).

            s% extras_startup => extras_startup
            s% extras_start_step => extras_start_step
            s% extras_check_model => extras_check_model
            s% extras_finish_step => extras_finish_step
            s% extras_after_evolve => extras_after_evolve
            s% how_many_extra_history_columns => how_many_extra_history_columns
            s% data_for_extra_history_columns => data_for_extra_history_columns
            s% how_many_extra_profile_columns => how_many_extra_profile_columns
            s% data_for_extra_profile_columns => data_for_extra_profile_columns  

            s% how_many_extra_history_header_items => how_many_extra_history_header_items
            s% data_for_extra_history_header_items => data_for_extra_history_header_items
            s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
            s% data_for_extra_profile_header_items => data_for_extra_profile_header_items
            ! consider energy from captured DM annihilation
            s% other_energy_implicit => DM_energy
      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
            integer, intent(in) :: id
            logical, intent(in) :: restart
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
      end subroutine extras_startup
      

      integer function extras_start_step(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            ! TODO: is this necesary
            ! calculate and apply energy from DM
            ! call DM_energy(id, ierr)
            if (ierr /= 0) return
            extras_check_model = keep_going         
            if (.false. .and. s% star_mass_h1 < 0.35d0) then
               ! stop when star hydrogen mass drops to specified level
               extras_check_model = terminate
               write(*, *) 'have reached desired hydrogen mass'
               return
            end if
            ! if you want to check multiple conditions, it can be useful
            ! to set a different termination code depending on which
            ! condition was triggered.  MESA provides 9 customizeable
            ! termination codes, named t_xtra1 .. t_xtra9.  You can
            ! customize the messages that will be printed upon exit by
            ! setting the corresponding termination_code_str value.
            ! termination_code_str(t_xtra1) = 'my termination condition'

            ! by default, indicate where (in the code) MESA terminated
            if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            ! add 7 custom variables to history
            how_many_extra_history_columns = 7
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
            integer, intent(in) :: id, n
            character (len=maxlen_history_column_name) :: names(n)
            real(dp) :: vals(n)
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return

            ! adding my custom variables to history
            names(1) = 'DM_rate'
            vals(1) = DM_rate
            names(2) = 'DM_energy_rate'
            vals(2) = DM_energy_rate
            names(3) = 'Mchi'
            vals(3) = Mchi
            names(4) = 'frac'
            vals(4) = frac
            names(5) = 'Hfrac'
            vals(5) = Hfrac
            names(6) = 'rho_chi'
            vals(6) = rho_chi
            names(7) = 'vbar'
            vals(7) = vbar

            ! note: do NOT add the extras names to history_columns.list
            ! the history_columns.list is only for the built-in history column options.
            ! it must not include the new column names you are adding here.
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            ! add my 11 custom profile columns
            how_many_extra_profile_columns = 11
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
            integer, intent(in) :: id, n, nz
            character (len=maxlen_profile_column_name) :: names(n)
            real(dp) :: vals(nz,n)
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            integer :: k
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return

            names(1) = 'extra_heat'
            names(2) = 'cell_DM_energy_rate'
            names(3) = 'cell_DM_annihlation_rate'
            names(4) = 'd_extra_heat_dlnR00'
            names(5) = 'd_extra_heat_dlnRp1'
            names(6) = 'd_extra_heat_dlndm1'
            names(7) = 'd_extra_heat_dlnd00'
            names(8) = 'd_extra_heat_dlndp1'
            names(9) = 'd_extra_heat_dlnTm1'
            names(10) = 'd_extra_heat_dlnT00'
            names(11) = 'd_extra_heat_dlnTp1'
            do k = 1, nz
               vals(k,1) = s% extra_heat(k)
               vals(k,2) = cell_DM_energy_rate(k)
               vals(k,3) = cell_DM_annihlation_rate(k)
               vals(k,4) = s% d_extra_heat_dlnR00(k)
               vals(k,5) = s% d_extra_heat_dlnRp1(k)
               vals(k,6) = s% d_extra_heat_dlndm1(k)
               vals(k,7) = s% d_extra_heat_dlnd00(k)
               vals(k,8) = s% d_extra_heat_dlndp1(k)
               vals(k,9) = s% d_extra_heat_dlndm1(k)
               vals(k,10) = s% d_extra_heat_dlnd00(k)
               vals(k,11) = s% d_extra_heat_dlndp1(k)
            end do

            ! note: do NOT add the extra names to profile_columns.list
            ! the profile_columns.list is only for the built-in profile column options.
            ! it must not include the new column names you are adding here.
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
            integer, intent(in) :: id, n
            character (len=maxlen_history_column_name) :: names(n)
            real(dp) :: vals(n)
            type(star_info), pointer :: s
            integer, intent(out) :: ierr
            ierr = 0
            call star_ptr(id,s,ierr)
            if(ierr/=0) return

            ! here is an example for adding an extra history header item
            ! also set how_many_extra_history_header_items
            ! names(1) = 'mixing_length_alpha'
            ! vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
            integer, intent(in) :: id, n
            character (len=maxlen_profile_column_name) :: names(n)
            real(dp) :: vals(n)
            type(star_info), pointer :: s
            integer, intent(out) :: ierr
            ierr = 0
            call star_ptr(id,s,ierr)
            if(ierr/=0) return

            ! here is an example for adding an extra profile header item
            ! also set how_many_extra_profile_header_items
            ! names(1) = 'mixing_length_alpha'
            ! vals(1) = s% mixing_length_alpha
      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            extras_finish_step = keep_going

            ! to save a profile, 
               ! s% need_to_save_profiles_now = .true.
            ! to update the star log,
               ! s% need_to_update_history_now = .true.

            ! see extras_check_model for information about custom termination codes
            ! by default, indicate where (in the code) MESA terminated
            if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      

      ! calculates energy from DM annihilation
      subroutine DM_energy(id, ierr)
            use const_def
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            double precision :: k_GeV, k_cgs, G, denom, numer
            double precision :: Mchi_cgs
            double precision, dimension(:), allocatable :: array(:)
            double precision, dimension(:), allocatable :: n_capchi(:)
            double precision, dimension(3) :: rad
            double precision, dimension(3) :: numerator
            integer :: k, inc, j

            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return

            ! VARIABLES:
            ! -> DM_energy_rate
            ! variable with total energy rate for star (L_DM) [erg/s]
            ! -> DM_rate
            ! capture rate of dark matter [1/s]

            ! ARRAYS:
            ! -> cell_DM_annihlation_rate(k) [1/s]
            ! array with annihlations per unit time for each cell
            ! -> cell_DM_energy_rate(k) [ergs/s]
            ! array with total energy rate for each cell
            ! -> extra_heat(k) [ergs/s/g]
            ! array with energy rate per unit mass for each cell

            ! read DM variables from inlist
            Mchi = s% X_CTRL(1)
            rho_chi = s% X_CTRL(2)
            vbar = s% X_CTRL(3)
            sigma_xenon = s% X_LOGICAL_CTRL(1) 

            ! constants
            k_GeV = 8.617d-14             ! boltzman const in GeV/k  
            k_cgs = 1.381d-16             ! boltzman constant in ergs/K
            Mchi_cgs = Mchi * 1.783d-24   ! DM mass in g

            ! allocate our array for the integral
            allocate (array(s% nz + 1))
            allocate (n_capchi(s% nz + 1))

            ! integral for the denominator in EQ. 5.8
            ! since the factor of n_c * 2 * Mchi * <sigma*nu>_ann * 4pi
            ! cancles when we divide these two integrals, we can just leave it off
            ! integrate captured DM number density over all cells:
            ! k = 1, top of outer most cell
            ! k = nz, top of inner most cell
            ! k = 0, bottom of inner most cell i.e. the center
            inc = 0
            ! TODO: OMP
            do k = (s% nz + 1), 1, -1
               inc = inc + 1
               ! if we are at the center of the star
               if (k.eq.(s% nz + 1)) then
                  ! radius at center is 0
                  array(inc) = 0.0d0
                  ! at center of star e^0 = 1
                  n_capchi(inc) = 1.0d0

               else
               ! this array has the radius at the outer boundery of each cell
               ! plus the final element is 0
               array(inc) = s% r(k)

               ! this is function being evaled at the outer edge of cell k
               ! n_capchi(inc) = exp((-2.0d0* Mchi_cgs * phi(id,k))/(k_cgs*&
               ! &10d0**(s% log_center_temperature)))
               n_capchi(inc) = ((s% r(k))**2) * exp((-2.0d0* Mchi_cgs * phi(id,k))/(k_cgs*&
               &10d0**(s% log_center_temperature)))

               end if
            end do

            ! call the trapezoid integral approx function
            ! this is the denominator of ians thesis Eq. 5.8
            ! TODO: integrals should be volume and not radius 
            denom = trap_int(array, n_capchi)

            ! fraction of annihilation products that thermalize
            frac = 2d0 / 3d0

            ! calculate the hydrogen fraction from MESA
            ! Hfrac = (s% star_mass_h1)/(s% star_mass)
            Hfrac = 0.75d0 ! temporary fix

            ! TODO: code thermalization timescale in?
            ! tau_th = ?

            ! initialize DM capture rate
            ! this varibale gets assigned in captureN_pureH
            DM_rate = 0.0d0

            ! call the subroutine which calculates the DM capture rate
            call captureN_pureH(s% mstar,s% r(1),Mchi,rho_chi,vbar,&
            &sigma_xenon,Hfrac,DM_rate)

            ! assume 2/3rds of annihilation products thermalize
            ! TODO: why 624.2??? I should have written this down
            DM_energy_rate = (624.2d0)*frac*DM_rate*mchi ! [erg/s]

            ! eval Eq. 5.8 for all cells
            ! TODO: OMP
            do k = 1, s% nz
               ! we have previously calculated the integral in the denominator
               ! and we calculated the capture rate for this timestep
               ! now we need to calculate the integral in the numerator
               ! this produces a value of Eq 5.8 corresponding to each cell
               ! to do this we do a trapezoid inegral with only three
               ! points, the bottom edge of cell, the center of mass of the
               ! cell, then top edge of the cell

               ! edge conditon for innermost cell
               if (k.eq.(s% nz)) then
                  ! phi, phimid, and phicent eval phi at outer edge, center
                  ! of mass, and inner edge of the innermost cell
                  rad(1) = s% R_center 
                  rad(2) = s% rmid(k) 
                  rad(3) = s% r(k) 
                  numerator(1) = exp((-2d0*Mchi_cgs*phicent(id))/(k_cgs*10d0**&
                  &(s% log_center_temperature)))
                  numerator(2) = exp((-2d0*Mchi_cgs*phimid(id,k))/(k_cgs*10d0**&
                  &(s% log_center_temperature))) 
                  numerator(3) = exp((-2d0*Mchi_cgs*phi(id,k))/(k_cgs*10d0**&
                  &(s% log_center_temperature))) 

               ! general case we can use phi(k+1) and r(k+1)
               else
                  rad(1) = s% r(k + 1) 
                  rad(2) = s% rmid(k) 
                  rad(3) = s% r(k) 
                  numerator(1) = (rad(1)**2) * exp((-2d0*Mchi_cgs*phi(id,k+1))/(k_cgs*10d0**&
                  &(s% log_center_temperature)))
                  numerator(2) = (rad(2)**2) * exp((-2d0*Mchi_cgs*phimid(id,k))/(k_cgs*10d0**&
                  &(s% log_center_temperature))) 
                  numerator(3) = (rad(3)**2) * exp((-2d0*Mchi_cgs*phi(id,k))/(k_cgs*10d0**&
                  &(s% log_center_temperature))) 
               end if
               
               ! now do our trapezoid rule integral for 3 points
               numer = trap_int(rad, numerator)

               ! if the numerator for second innermost cell is 0
               ! then for greater accuracy we should set the innermost
               ! cells numerator to the denominator
               ! TODO: change this from exactly zero to very small
               if ((k.eq.(s% nz - 1)).and.(numer.eq.0)) then
                  ! DM profile is contained entirely in the first cell
                  ! so set it to the denominator
                  ! i.e. all DM heating is applied to inner cell
                  cell_DM_annihlation_rate(s% nz) = (DM_rate/2.0d0)
                  cell_DM_energy_rate(s% nz) = &
                  &2d0 * Mchi * cell_DM_annihlation_rate(s% nz)
                  s% extra_heat(s% nz) = (cell_DM_energy_rate(s% nz))/&
                  &(s% dm(s% nz))

                  ! now we have to set all other cells to zero
                  do j = 1, (s% nz - 1)
                     cell_DM_annihlation_rate(j) = 0d0
                     cell_DM_energy_rate(j)= 0d0
                     s% extra_heat(j) = 0d0
                  end do

                  ! break the do loop so we do not re-eval the other cells
                  exit
               end if

               ! specifically this is number of annihlations,
               ! NOT the number of particles that annihlate
               cell_DM_annihlation_rate(k) = (DM_rate/2.0d0)*(numer/denom)

               ! energy rate for each cell
               cell_DM_energy_rate(k) =2d0*Mchi*cell_DM_annihlation_rate(k)

               ! energy rate per unit mass for each cell:
               ! MESAs extra_heat array and partials derivatives
               s% extra_heat(k) = (cell_DM_energy_rate(k)) / (s% dm(k))
            end do

            ! partial derivatives of extra_heat(k)
            ! TODO: eval these analytically for added stability?
            ! dlnR00(k) is partial of Q w.r.t. ln R(k)
            ! dlnRp1(k) is partial of Q w.r.t. ln R(k+1)
            ! dlndm1(k) is partial of Q w.r.t. ln rho(k-1)
            ! dlnd00(k) is partial of Q w.r.t. ln rho(k)
            ! dlndp1(k) is partial of Q w.r.t. ln rho(k+1)

            ! do loop for assinging partials
            !$OMP PARALLEL DO PRIVATE (k)
            do k = 1, s% nz
               ! PARTIALS w.r.t. ln(radius)
               ! edge case for central cell
               if (k.eq.(s% nz)) then
                  ! w.r.t. ln R(k)
                  s% d_extra_heat_dlnR00(k) = (s% extra_heat(k) - &
                  &s% extra_heat(k-1))/(s% lnR(k) - s% lnR(k-1))

                  ! w.r.t. ln R(k+1)
                  s% d_extra_heat_dlnRp1(k) = (s% extra_heat(k-1) - &
                  &s% extra_heat(k-2))/(s% lnR(k) - s% lnR(k-1))

               ! edge case for second from central cell 
               else if (k.eq.(s% nz - 1)) then
                  ! w.r.t. ln R(k)
                  s% d_extra_heat_dlnR00(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnR(k+1) - s% lnR(k-1))

                  ! w.r.t. ln R(k+1)
                  s% d_extra_heat_dlnRp1(k) = (s% extra_heat(k) - &
                  &s% extra_heat(k-1))/(s% lnR(k+1) - s% lnR(k))

               ! edge case for surface cell
               else if (k.eq.1) then
                  ! w.r.t. ln R(k)
                  s% d_extra_heat_dlnR00(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k))/(s% lnR(k+1) - s% lnR(k))

                  ! w.r.t. ln R(k+1)
                  s% d_extra_heat_dlnRp1(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k))/(s% lnR(k+2) - s% lnR(k+1))

               ! general case
               else
                  ! dlnr00(k) = partial of Q w.r.t. ln R(k)
                  s% d_extra_heat_dlnR00(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnR(k+1) - s% lnR(k-1))

                  ! dlnrp1(k) = partial of Q w.r.t. ln R(k+1)
                  s% d_extra_heat_dlnRp1(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnR(k+2) - s% lnR(k))
               end if 

               ! PARTIALS w.r.t. ln density
               ! edge case for central cell
               if (k.eq.(s% nz)) then
                  ! w.r.t. ln rho(k-1)
                  s% d_extra_heat_dlndm1(k) = (s% extra_heat(k) - &
                  &s% extra_heat(k-1))/(s% lnd(k-1) - s% lnd(k-2))

                  ! w.r.t. ln rho(k)
                  s% d_extra_heat_dlnd00(k) = (s% extra_heat(k) - &
                  &s% extra_heat(k-1))/(s% lnd(k) - s% lnd(k-1))

                  ! w.r.t. ln rho(k+1)
                  s% d_extra_heat_dlndp1(k) = (s% extra_heat(k-1) - &
                  &s% extra_heat(k-2))/(s% lnd(k) - s% lnd(k-1))

               ! edge case for second from central cell 
               else if (k.eq.(s% nz - 1)) then
                  ! w.r.t. ln rho(k-1)
                  s% d_extra_heat_dlndm1(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnd(k) - s% lnd(k-2))

                  ! w.r.t. ln rho(k)
                  s% d_extra_heat_dlnd00(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnd(k+1) - s% lnd(k-1))

                  ! w.r.t. ln rho(k+1)
                  s% d_extra_heat_dlndp1(k) = (s% extra_heat(k) - &
                  &s% extra_heat(k-1))/(s% lnd(k+1) - s% lnd(k))

               ! edge case for surface cell
               else if (k.eq.1) then
                  ! w.r.t. ln rho(k-1)
                  s% d_extra_heat_dlndm1(k) = (s% extra_heat(k+2) - &
                  &s% extra_heat(k+1))/(s% lnd(k+1) - s% lnd(k))

                  ! w.r.t. ln rho(k)
                  s% d_extra_heat_dlnd00(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k))/(s% lnd(k+1) - s% lnd(k))

                  ! w.r.t. ln rho(k+1)
                  s% d_extra_heat_dlndp1(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k))/(s% lnd(k+2) - s% lnd(k+1))

               ! edge case for second from surface cell 
               else if (k.eq.(2)) then
                  ! w.r.t. ln rho(k-1)
                  s% d_extra_heat_dlndm1(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k))/(s% lnd(k) - s% lnd(k-1))

                  ! w.r.t. ln rho(k)
                  s% d_extra_heat_dlnd00(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnd(k+1) - s% lnd(k-1))

                  ! w.r.t. ln rho(k+1)
                  s% d_extra_heat_dlndp1(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnd(k+2) - s% lnd(k))

               ! general case
               else
                  ! dlndm1(k) = partial of Q w.r.t. ln rho(k-1)
                  s% d_extra_heat_dlndm1(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnd(k) - s% lnd(k-2))

                  ! dlnd00(k) = partial of Q w.r.t. ln rho(k)
                  s% d_extra_heat_dlnd00(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnd(k+1) - s% lnd(k-1))

                  ! dlndp1(k) = partial of Q w.r.t. ln rho(k+1)
                  s% d_extra_heat_dlndp1(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnd(k+2) - s% lnd(k))
               end if 

               ! PARTIALS w.r.t. ln Temp
               ! edge case for central cell
               if (k.eq.(s% nz)) then
                  ! w.r.t. ln T(k-1)
                  s% d_extra_heat_dlntm1(k) = (s% extra_heat(k) - &
                  &s% extra_heat(k-1))/(s% lnT(k-1) - s% lnT(k-2))

                  ! w.r.t. ln T(k)
                  s% d_extra_heat_dlnt00(k) = (s% extra_heat(k) - &
                  &s% extra_heat(k-1))/(s% lnT(k) - s% lnT(k-1))

                  ! w.r.t. ln T(k+1)
                  s% d_extra_heat_dlntp1(k) = (s% extra_heat(k-1) - &
                  &s% extra_heat(k-2))/(s% lnT(k) - s% lnT(k-1))

               ! edge case for second from central cell 
               else if (k.eq.(s% nz - 1)) then
                  ! w.r.t. ln T(k-1)
                  s% d_extra_heat_dlntm1(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnT(k) - s% lnT(k-2))

                  ! w.r.t. ln T(k)
                  s% d_extra_heat_dlnt00(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnT(k+1) - s% lnT(k-1))

                  ! w.r.t. ln T(k+1)
                  s% d_extra_heat_dlntp1(k) = (s% extra_heat(k) - &
                  &s% extra_heat(k-1))/(s% lnT(k+1) - s% lnT(k))

               ! edge case for surface cell
               else if (k.eq.1) then
                  ! w.r.t. ln T(k-1)
                  s% d_extra_heat_dlntm1(k) = (s% extra_heat(k+2) - &
                  &s% extra_heat(k+1))/(s% lnT(k+1) - s% lnT(k))

                  ! w.r.t. ln T(k)
                  s% d_extra_heat_dlnt00(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k))/(s% lnT(k+1) - s% lnT(k))

                  ! w.r.t. ln T(k+1)
                  s% d_extra_heat_dlntp1(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k))/(s% lnT(k+2) - s% lnT(k+1))

               ! edge case for second from surface cell 
               else if (k.eq.(2)) then
                  ! w.r.t. ln T(k-1)
                  s% d_extra_heat_dlntm1(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k))/(s% lnT(k) - s% lnT(k-1))

                  ! w.r.t. ln T(k)
                  s% d_extra_heat_dlnt00(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnT(k+1) - s% lnT(k-1))

                  ! w.r.t. ln T(k+1)
                  s% d_extra_heat_dlntp1(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnT(k+2) - s% lnT(k))

               ! general case
               else
                  ! dlndm1(k) = partial of Q w.r.t. ln T(k-1)
                  s% d_extra_heat_dlntm1(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnT(k) - s% lnT(k-2))

                  ! dlnd00(k) = partial of Q w.r.t. ln T(k)
                  s% d_extra_heat_dlnt00(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnT(k+1) - s% lnT(k-1))

                  ! dlndp1(k) = partial of Q w.r.t. ln T(k+1)
                  s% d_extra_heat_dlntp1(k) = (s% extra_heat(k+1) - &
                  &s% extra_heat(k-1))/(s% lnT(k+2) - s% lnT(k))
               end if 

               ! dlnTm1(k) is partial of Q w.r.t. ln T(k-1)
               ! dlnT00(k) is partial of Q w.r.t. ln T(k)
               ! dlnTp1(k) is partial of Q w.r.t. ln T(k+1)
               ! s% d_extra_heat_dlntm1(k) = 0d0
               ! s% d_extra_heat_dlnt00(k) = 0d0
               ! s% d_extra_heat_dlntp1(k) = 0d0
            end do
            !$OMP END PARALLEL DO

            ! clean up memory after from DM captured density integral
            deallocate (array)
            deallocate (n_capchi)
            return
      end subroutine DM_energy


      ! calucaltes gravitation potenital for cell k
      function phi(id, k) result (value)
            ! specifically the potenital at the outer edge cell k
            type (star_info), pointer :: s
            integer, intent(in) :: id
            integer :: ierr
            integer :: k, i, inc
            double precision :: G, integral, value
            double precision, dimension(:), allocatable :: array(:)
            double precision, dimension(:), allocatable :: func(:)
            ierr = 0
            call star_ptr(id, s, ierr)

            G = 6.6743d-8           ! gravitational const. in cgs

            ! allocate our array
            allocate (array(s% nz - k + 2))
            allocate (func(s% nz - k + 2))

            ! across all cells inteior of k
            ! generate an array from radius(k) -> 0
            ! then evaluate the integrand at each point in mesh
            inc = 0
            do i = (s% nz + 1), k, -1
               inc = inc + 1
               ! if we are at the center of the star
               if (i.eq.(s% nz + 1)) then
                  ! radius at center is 0
                  array(inc) = s% R_center
                  ! divide by zero check
                  if (s% R_center.eq.0) then 
                     ! GM(r)/r = 0 at center
                     func(inc) = 0d0
                  else
                     ! if there is mass or space interior to central cell
                     func(inc) = (G*(s% M_center))/(s% R_center) 
                  end if
               else
                  ! array is the radius at the outer edge of each cell
                  array(inc) = s% r(i)
                  ! func is the function being evaled at each radius
                  ! GM(r)/r
                  func(inc) = (G * (s% m(i))) / (s% r(i)) 
               end if
            end do

            ! call the trapozoid inegral approx function thing
            integral = trap_int(array, func)

            ! clean up memory
            deallocate (func)
            deallocate (array)

            value = integral 
      end function phi


      ! calucaltes gravitation potenital for inner edge of central cell
      function phicent(id) result (value)
            ! this function is almost always zero
            ! rarely MESA will have some mass and radius inteiror of the
            ! innermost cell, so in those cases we calculate grav poten.
            type (star_info), pointer :: s
            integer, intent(in) :: id
            integer :: ierr
            integer :: k, i, inc
            double precision :: G, integral, value
            double precision, dimension(:), allocatable :: array(:)
            double precision, dimension(:), allocatable :: func(:)
            ierr = 0
            call star_ptr(id, s, ierr)

            G = 6.6743d-8           ! gravitational const. in cgs

            ! allocate our array
            allocate (array(2))
            allocate (func(2))

            ! radius at center is ussually 0
            array(1) = 0d0
            array(2) = s% R_center
            func(1) = 0d0
            ! divide by zero check
            if (s% R_center.eq.0) then 
               ! GM(r)/r = 0 at center
               func(2) = 0d0
            else
               ! if there is mass or space interior to central cell
               func(2) = (G*(s% M_center))/(s% R_center) 
            end if

            ! call the trapozoid inegral approx function thing
            integral = trap_int(array, func)

            ! clean up memory
            deallocate (func)
            deallocate (array)
            value = integral 
      end function phicent


      ! calucaltes gravitation potenital at cener of mass of cell k
      function phimid(id, k) result (value)
            ! specifically the potenital at the COM of cell k
            type (star_info), pointer :: s
            integer, intent(in) :: id
            integer :: ierr
            integer :: k, inc, i
            double precision :: G, integral, value
            double precision, dimension(:), allocatable :: array(:)
            double precision, dimension(:), allocatable :: func(:)
            ierr = 0
            call star_ptr(id, s, ierr)

            G = 6.6743d-8           ! gravitational const. in cgs

            ! allocate our array
            allocate (array(s% nz - k + 2))
            allocate (func(s% nz - k + 2))

            ! across all cells inteior of k
            ! generate an array from radius(k) -> 0
            ! then evaluate the integrand at each point in mesh
            inc = 0
            do i = (s% nz + 1), k, -1
               inc = inc + 1
               ! if we are at the center of the star
               if (i.eq.(s% nz + 1)) then
                  ! radius at center is 0
                  array(inc) = s% R_center

                  ! divide by zero check
                  if (s% R_center.eq.0) then 
                     ! GM(r)/r = 0 at center
                     func(inc) = 0d0
                  else
                     ! if there is mass inside of innermost cell
                     func(inc) = (G*(s% M_center))/(s% R_center)
                  end if
               else
               ! array is the radius at the outer edge of each cell
               array(inc) = s% rmid(i)
               ! func is the function being evaled at each radius
               ! [G (M(r)-CellMass(r)/2)] / [(r center of mass)]
               ! (s% m(i) - s% dm(i)/2) is the mass enclosed at COM of cell k
               func(inc) = ( G * (s% m(i) - ((s% dm(i))/2d0)))/(s% rmid(i))
               end if
            end do

            ! call the trapozoid inegral approx function thing
            integral = trap_int(array, func)

            ! clean up memory and return
            deallocate (func)
            deallocate (array)
            value = integral 
      end function phimid


      ! calculates the lower incomplete gamma function
      function gamma_incl(s, x) result (value)
            ! incomplete plete gamma function is defined as:
            ! γ(s,x) = ∫ t^(s−1) e^(−t) dt evaluated 0 -> x 
            ! approximation. note that the mesh spacing of x does not have to be uniform.
            double precision :: value  ! gamma(s, x)
            double precision :: x, s, lg_func, integral
            double precision, dimension(:), allocatable :: array(:)
            double precision, dimension(:), allocatable :: func(:)
            integer :: i, mesh

            ! set mesh size and allocate
            mesh = 1000       !TODO: optimize mesh denisty and spacing
            allocate (array(mesh))
            allocate (func(mesh))

            ! generate arrays
            !$OMP PARALLEL DO PRIVATE(i, lg_func) SHARED(array, func, mesh, s, x)
            do i=1,mesh
               array(i) = (dble(i)/dble(mesh)) * x
               func(i) = (array(i)**(s-1d0)) * exp(-1d0*array(i))

               ! check for infinity overflow 
               if (func(i).gt.0.8d308) then
                  ! approximate in log form
                  lg_func = log10(array(i))*(S-1d0) + log10(exp(1d0))*(-1d0*array(i))
                  if (lg_func.gt.307d0) then
                     func(i) = 1d307
                  else
                     func(i) = 10d0**lg_func
                  end if
               end if

               ! check for NaN overflow
               if (func(i).ne.func(i)) then
                  ! approximate in log form
                  lg_func = log10(array(i))*(S-1d0) + log10(exp(1d0))*(-1d0*array(i))
                  if (lg_func.lt.-307d0) then
                     func(i) = 1d-307
                  else
                     func(i) = 10d0**lg_func
                  end if
               end if
               ! assign to arrays in thread safe manner
               !!$OMP ATOMIC WRITE
               !array(i) = arr
               !!$OMP ATOMIC WRITE
               !func(i) = fun
            end do
            !$OMP END PARALLEL DO

            ! take trapezoid rule integral 
            integral = trap_int(array, func) 
            ! check for infinity overflow 
            if (integral.gt.0.8d308) then
               integral = 1d307
            end if
            ! check for NaN overflow
            if (integral.lt.0.8d-307) then
               integral = 1d-307
            end if

            ! clean up memory and return
            deallocate (func)
            deallocate (array)
            value = integral
      end function gamma_incl


      ! heaviside step function, x=0 -> 1/2 
      function heavi(x) result (value)
            double precision :: value, x
            value = 0.5 * (sign(1d0,x) + 1)
      end function heavi


      ! calculates the integral of an array y with respect to x using the trapezoid
      function trap_int(x, y) result (value)
            ! approximation. note that the mesh spacing of x does not have to be uniform.
            double precision :: x(:)        ! variable x
            double precision :: y(size(x))  ! function y(x)
            double precision :: value       ! integral ∫y(x)·dx

            ! integrate using the trapezoidal rule
            associate(n => size(x))
                  value = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
            end associate
      end function trap_int


      ! caluclates the capture rate based on various DM and stellar
      ! parameters, assigns that rate to DM_rate 
      subroutine captureN_pureH(M,R,Mchi,rho_chi,vbar,sigma_xenon,&
      &Hfrac,DM_rate)
            integer :: N, less_count, Ncut 
            double precision :: G, mn, mn_grams, R, M, beta 
            double precision :: Mchi, rho_chi, vbar, nchi, A2N
            double precision :: Vesc, Vstar, n_H, tau, sigma
            double precision :: pi, Hfrac, VN
            double precision :: DM_rate, zavg, pn_tau
            double precision :: Cn, Ctot
            logical :: sigma_xenon

            ! define constants
            G = 6.6743d-8           ! gravitational const. in cgs
            mn = 0.93827d0          ! mass of nucleons (protons) [GeV]
            mn_grams = 1.6726d-24   ! mass of nucleons [g]
            pi = 3.14159d0          ! pi
            zavg = 2.0d0

            ! do we use the xenon bounds for cross section?
            if (sigma_xenon .eqv. .TRUE.) then
                sigma = (1.26d-40)*(Mchi / 1.0d8)
            end if

            ! calculate vesc, vstar, n_H, nchi, beta, tau
            Vesc = sqrt(2d0*G*M/R)            ! escape velocity [cm/s]
            Vstar = pi*(4d0/3d0)*(R**3d0)     ! volume of star [cm^3]
            n_H = ((Hfrac)*M)/(mn_grams*Vstar)   ! H # density [cm^-3]
            beta = (4d0*Mchi*mn)/((Mchi+mn)**2d0)   ! reduced mass
            nchi = rho_chi/mchi
            tau = 2d0*R*sigma*n_H
            ! TODO: optimize this?
            Ncut = ceiling(tau*2d0 + 10d0)
            ! if (Ncut.gt.170) then
            !    Ncut = 170
            ! end if

            ! initialize total cap rate
            Ctot = 0.0d0

            ! PARTIAL CAPTURE SUM
            !$OMP PARALLEL DO PRIVATE(N, VN, Cn, A2N) SHARED(Ctot, tau, Vesc, beta, mn, Mchi, vbar, G, nchi, M, R, zavg, pi)
            do N=1,Ncut
            ! calculate p_N(tau)
            ! check if tau >> 1 or tau << 1
            if (tau.gt.100) then
               ! from PRD 2020 Eq. A3b
               pn_tau = (2d0/(tau**2d0)) * (N + 1) * heavi(tau - N)
            else
               ! from PRD 2020 Eq. 2
               ! TODO: this returns 0 above N = 170, this is because 170! is the largest
               ! factorial fortran double precsion can store (10^308)
               ! so if there is a significant amount of capture past this, we are boned
               pn_tau = (2d0/(tau**2d0)) * (gamma_incl(dble(N+2),tau)/gamma(dble(N+1)))
            end if

                ! calculate velocity after scattering
                VN = Vesc*(1.0d0 - (beta/zavg))**((-1d0*N)/2d0)

                ! calculate the full partial capture rate and append
                ! check for loss of signficance error
                if (Mchi.gt.1d16) then
                  ! use eq. 23 from bramante 2017 
                  ! print *, "high mass"
                  A2N = (3d0 * (Vesc**2d0) * N * mn )/((vbar**2d0) * Mchi)
                  Cn = sqrt(24d0*pi)*pn_tau*G*nchi*M*R*(1d0/vbar)*(1d0-&
                  &(1d0+((2d0*A2N*(vbar**2d0))/(3d0*(Vesc**2d0))))*exp(-1d0*A2N))
                else
                  ! bramenate 2017 eq. 22 
                  Cn = pi*(R**2d0)*pn_tau*((sqrt(2.d0)*nchi)/&
                  &(sqrt(3d0*pi)*vbar))*((((2d0*(vbar**2d0))+(3d0*Vesc**2d0))-((2d0*vbar**2d0)+&
                  &(3d0*VN**2d0)) * exp(-1d0*((3d0*((VN**2d0)-(Vesc**2d0)))/(2d0*(vbar**2d0))))))
                end if

                ! add partial capture to total capture
                !$OMP ATOMIC UPDATE
                Ctot = Ctot + Cn
                !$OMP END ATOMIC
            end do
            !$OMP END PARALLEL DO

            ! assign Ctot to variable that gets passed back
            DM_rate = Ctot
            return
      end subroutine captureN_pureH

      
      subroutine extras_after_evolve(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
      end subroutine extras_after_evolve


      end module run_star_extras
      

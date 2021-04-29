! ******************************************************* !
!     writen in fall 2020 by ian k. bania, advised        !
!     by cosmin ilie at colgate university                !
! ******************************************************* !

program DM_heating

   ! global variables
   double precision :: DM_rate, DM_energy_rate, Mchi, frac
   double precision :: rho_chi, vbar, Hfrac
   double precision :: bass
   double precision, dimension(100000) :: cell_DM_energy_rate
   double precision, dimension(100000) :: cell_DM_annihlation_rate
   logical :: sigma_xenon
   CHARACTER(len=1) :: arg

   !*******************************************
   !!! \/ EDIT ME TO CHANGE DM VARIABLES \/ !!!
   !*******************************************
   ! make sure these are all double precision variables
   ! i.e. 2.12 * 10^6 is writen: 2.12d6
   Mchi = 1.0d12     ! [GeV] TODO:significant error above 10^13 GeV
   rho_chi = 1.0d9   ! []
   vbar = 1.0d6      ! []
   sigma_xenon = .TRUE. ! true or false 
   !*******************************************
   !!! /\ EDIT ME TO CHANGE DM VARIABLES /\ !!!
   !*******************************************

   ! initialize DM capture rate
   ! this varibale gets assigned in captureN_pureH
   DM_rate = 0.0d0
   CALL get_command_argument(1, arg)
   !$OMP PARALLEL DO PRIVATE(Bass, DM_rate)
   do N =1,943
   ! do N =3,4
   ! do N =-100,100
      bass = (1.0d0)*(1.05d0)**(dble(N))
      ! bass = 1.0d14

      if ( arg == "1" ) then
         CALL captureN_pureH(1.99d33, 6.96d10, bass, 1.0d10, 1.0d6, .TRUE., 0.75d0, DM_rate)
         print *, Bass, DM_rate
      end if

      if ( arg == "2" ) then
         CALL captureN_pureH(9.94d33, 8.37d10, bass, 1.0d10, 1.0d6, .TRUE., 0.75d0, DM_rate)
         print *, Bass, DM_rate
      end if

      if ( arg == "3" ) then
         CALL captureN_pureH(2.98d34, 1.03d11, bass, 1.0d10, 1.0d6, .TRUE., 0.75d0, DM_rate)
         print *, Bass, DM_rate
      end if

      if ( arg == "4" ) then
         CALL captureN_pureH(3.98d34, 1.15d11, bass, 1.0d10, 1.0d6, .TRUE., 0.75d0, DM_rate)
         print *, Bass, DM_rate
      end if

      if ( arg == "5" ) then
         CALL captureN_pureH(1.99d35, 2.97d11, bass, 1.0d10, 1.0d6, .TRUE., 0.75d0, DM_rate)
         print *, Bass, DM_rate
      end if

      if ( arg == "6" ) then
         CALL captureN_pureH(7.95d35, 6.28d11, bass, 1.0d10, 1.0d6, .TRUE., 0.75d0, DM_rate)
         print *, Bass, DM_rate
      end if

      if ( arg == "7" ) then
         CALL captureN_pureH(1.99d36, 1.03d12, bass, 1.0d10, 1.0d6, .TRUE., 0.75d0, DM_rate)
         print *, Bass, DM_rate
      end if

      ! gamma func?
      if ( arg == "8" ) then
         DM_rate = gamma_incl(3d0, Bass)
         print *, Bass, DM_rate
      end if
   end do
   !$OMP END PARALLEL DO
   contains


   ! calculates the lower incomplete gamma function
   function gamma_incl(s, x) result (value)
      ! incomplete plete gamma function is defined as:
      ! γ(s,x) = ∫ t^(s−1) e^(−t) dt evaluated 0 -> x 
      ! approximation. note that the mesh spacing of x does not have to be uniform.
      double precision :: value  ! gamma(s, x)
      double precision :: x, s, lg_func, integral
      double precision :: arr, fun
      double precision, dimension(:), allocatable :: array(:)
      double precision, dimension(:), allocatable :: func(:)
      integer :: i, mesh
   
      ! set mesh size and allocate
      mesh = 1000
      ! use  100000 at highest mass
      allocate (array(mesh))
      allocate (func(mesh))

      ! generate arrays
      !$OMP PARALLEL DO PRIVATE(i, arr, fun, lg_func) SHARED(array, func, mesh, s, x)
      do i=1,mesh
         ! print *, "made it here", i
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

         ! array(i) = arr
         ! func(i) = fun
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
      ! if (isnan(func(i))) then
      ! if (func(i).ne.func(i)) then !!!TODO: this breaks everything?
         ! integral = 1d-307
      ! end if

      ! check for NaN overflow
      if (integral.lt.0.8d-308) then
         integral = 1d-307
      end if

      ! clean up memory and return
      deallocate (func)
      deallocate (array)
      value = integral
   end function gamma_incl


   ! calculates the integral of an array y with respect to x using the trapezoid
   function trap_int(x, y) result (value)
      ! approximation. note that the mesh spacing of x does not have to be uniform.
      double precision :: x(:)            ! variable x
      double precision :: y(size(x))      ! function y(x)
      double precision :: value, tmpsum   ! integral ∫y(x)·dx

      ! integrate using the trapezoidal rule
      associate(n => size(x))
         value = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
      end associate
   end function trap_int

   ! heaviside step function, x=0 -> 1/2 
   function heavi(x) result (value)
      double precision :: value, x
      value = 0.5 * (sign(1d0,x) + 1)
   end function heavi

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
            ! from PRD'20 Eq. A3b
            pn_tau = (2d0/(tau**2d0)) * (N + 1) * heavi(tau - N)
         else
            ! from PRD'20 Eq. 2
            ! TODO: this returns 0 above N = 170, this is because 170! is the largest
            ! factorial fortran double precsion can store (10^308)
            ! so if there's a significant amount of capture past this, we're boned
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

         ! print *, N, Cn
         ! add partial capture to total capture
         !$OMP ATOMIC UPDATE
         Ctot = Ctot + Cn
         !$OMP END ATOMIC
      end do
      !$OMP END PARALLEL DO

      ! assign Ctot to variable that gets passed back
      DM_rate = Ctot
      Ctot = 0.0d0
      return
   end subroutine captureN_pureH


end program DM_heating

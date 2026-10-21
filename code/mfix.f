#include "error.inc"

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MFIX                                                    !
!  Author: M. Syamlal                                 Date: 29-JAN-92  !
!                                                                      !
!  Purpose: The main module in the MFIX program                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!
!> \mainpage Multiphase Flow with Interphase eXchanges
!!
!! MFIX is a general-purpose computer code developed at the National
!! Energy Technology Laboratory, NETL, for describing the hydrodynamics,
!! heat transfer, and chemical reactions in fluid-solid systems.
!!
!! It has been used for describing bubbling and circulating fluidized
!! beds and spouted beds. MFiX calculations give transient data on the
!! three-dimensional distribution of pressure, velocity, temperature,
!! and species mass fractions. MFiX code is based on a generally
!! accepted set of multiphase flow equations. The code is used as a
!! "test-stand" for testing and developing multiphase flow constitutive
!!  equations.
!!
!! \section Notice
!! Neither the United States Government nor any agency thereof, nor any
!! of their employees, makes any warranty, expressed or implied, or
!! assumes any legal liability or responsibility for the accuracy,
!! completeness, or usefulness of any information, apparatus, product,
!! or process disclosed or represents that its use would not infringe
!! privately owned rights.
!!
!! * MFIX is provided without any user support for applications in the
!!   user's immediate organization. It should not be redistributed in
!!   whole or in part.
!!
!! * The use of MFIX is to be acknowledged in any published paper based
!!   on computations using this software by citing the MFIX theory
!!   manual. Some of the submodels are being developed by researchers
!!   outside of NETL. The use of such submodels is to be acknowledged
!!   by citing the appropriate papers of the developers of the submodels.
!!
!! * The authors would appreciate receiving any reports of bugs or other
!!   difficulties with the software, enhancements to the software, and
!!   accounts of practical applications of this software.
!!
!! \section Disclaimer
!! This report was prepared as an account of work sponsored by an agency
!! of the United States Government. Neither the United States Government
!! nor any agency thereof, nor any of their employees, makes any
!! warranty, express or implied, or assumes any legal liability or
!! responsibility for the accuracy, completeness, or usefulness of any
!! information, apparatus, product, or process disclosed, or represents
!! that its use would not infringe privately owned rights. Reference
!! herein to any specific commercial product, process, or service by
!! trade name, trademark, manufacturer, or otherwise does not
!! necessarily constitute or imply its endorsement, recommendation, or
!! favoring by the United States Government or any agency thereof. The
!! views and opinions of authors expressed herein do not necessarily
!! state or reflect those of the United States Government or any
!! agency thereof.

    subroutine run_mfix(mfix_dat_filename, mfix_dat_size)
        use check_batch_queue_end_mod, only: check_batch_queue_end
        use compar, only:adjust_partition, mype, pe_io
        use des_time_march, only: des_time_init, des_time_step, des_time_end, factor, exit_loop
        use discretelement, only: des_continuum_coupled, discrete_element
        use error_manager, only: init_error_manager
        use exit, only: exit_flag, check_exit_flag
        use funits, only: unit_dat
        use iso_c_binding
        use iterate, only: converged, diverged, adjustdt
        use iterate, only: iterate_init, do_iteration, post_iterate
        use iterate, only: log_converged, log_diverged, nit, max_nit
        use main, only: get_data, initialize, finalize
        use main, only: print_flags
        use output_manager, only: print_granular_step_info, update_time_step_info_float
        use output_manager, only: write_outputs
        use pause, only: wait_while_paused
        use pic_time_march_mod, only: pic_time_init, pic_time_step, pic_time_end
        use pic_time_march_mod, only: s_time, tend_pic_loop
        use run, only:  dt, dem_solids, pic_solids, steady_state, time, tstop, mfix_input_filename
        use step, only: check_low_dt, chem_mass, time_step_init, time_step_end
        use time_mod, only: get_elapsed_time

        implicit none
! Path to input file
        character(len=*), intent(in) :: mfix_dat_filename
        integer, intent(in) :: mfix_dat_size

!-----------------------------------------------
! Local variables
!-----------------------------------------------

        double precision :: step_start_time
        integer :: ii

        character(mfix_dat_size) :: mfix_dat
        integer :: filename_len
        integer :: ios

        filename_len = len_trim(mfix_dat_filename)
        mfix_input_filename = mfix_dat_filename(1:filename_len)

        step_start_time = get_elapsed_time() ! first call sets timer to 0

! Dynamic load balance loop
        do

! Open input file. Report errors if the file is not located or
! there are difficulties opening it.

            open(unit=unit_dat, file=mfix_dat_filename, status='old', iostat=ios, &
                form="unformatted", access="stream", action="read")
            if(ios /= 0) then
                if(mype == pe_io) write (*,1001) mfix_dat_filename
#ifdef __GFORTRAN__
                call disable_backtrace
#endif
         error stop
1001     format(2/,1X,70('*')/' From: run_mfix',/' Error 1001: ',    &
              'Unable to open input data file: ',A/,'Aborting.',/1x,70('*'),2/)
      endif

      read(unit_dat) mfix_dat
      close(unit_dat)

! Read input data, check data, do computations for IC and BC locations
! and flows, and set geometry parameters such as X, X_E, DToDX, etc.
      call get_data(mfix_dat)

      if (check_exit_flag()) return

! Initialize the simulation
      call initialize(mfix_dat)

      if (check_exit_flag()) return

! Time march loop.
      dt_loop: DO WHILE (TIME + 0.1*DT < TSTOP)

          call wait_while_paused

          if(des_continuum_coupled .or. .not.discrete_element) then
              call run_fluid
          endif

          if (dem_solids) then
              call run_dem
              if(.not.des_continuum_coupled) then
! Pure DEM skips time_step_end; save a requested stop between outputs.
! Partition adjustment already writes its checkpoint in finalize.
                  if(exit_flag .and. .not.adjust_partition) then
                      call write_outputs(.true., .false.)
                  endif
                  exit
              endif
              if (check_exit_flag()) then
                  call time_step_end
                  exit
              endif
          endif

          if (pic_solids) then
              call run_pic
              if(.not.des_continuum_coupled) exit
              if (check_exit_flag()) then
                  call time_step_end
                  exit
              endif
          endif

! Terminate MFIX normally before batch queue terminates.
          call check_batch_queue_end(exit_flag)
          call time_step_end

! Transient or steady state simulation
          if (steady_state .or. adjust_partition) exit

          if (check_exit_flag()) exit

      enddo dt_loop

      call finalize
      if(.not.adjust_partition) exit
  enddo

contains

    subroutine run_fluid
        step_start_time = get_elapsed_time()
        call time_step_init(mfix_dat)
        do
            call iterate_init
            do while (nit<max_nit .and. .not.(converged.or.diverged))
                nit = nit + 1
                call do_iteration(mfix_dat)
                if (check_exit_flag()) return
                call wait_while_paused
            enddo
            call post_iterate
            if (steady_state) exit
            if (.not.adjustdt(mfix_dat)) exit
        enddo

! Exit if DT < DT_MIN
        call check_low_dt
        if (check_exit_flag()) return

! Stiff chemistry solver
        call chem_mass
        if (check_exit_flag()) return
        call print_walltime("Timestep walltime, fluid solver:", step_start_time)

    end subroutine run_fluid

    subroutine run_dem
        step_start_time = get_elapsed_time()
        call des_time_init
        do ii = 1, factor
            call des_time_step(ii)
            if (mod(factor-ii, 10) == 0) then
               call check_batch_queue_end(exit_flag)
               exit_loop = exit_loop .or. exit_flag
            endif
            if ( exit_loop ) exit
            call wait_while_paused
            call print_granular_step_info(ii)
         enddo
         call des_time_end
         call print_walltime("Timestep walltime, DEM solver:  ", step_start_time)
      end subroutine run_dem

      subroutine run_pic
! number of PIC time steps
         integer :: pic_iters

         step_start_time = get_elapsed_time()
         call pic_time_init

! If the current time in the discrete loop exceeds the current time in
! the continuum simulation, exit the Lagrangian loop
         pic_iters = 0
         do while(s_time.lt.tend_pic_loop)
            pic_iters  = pic_iters + 1
            call pic_time_step(pic_iters)
            if (mod(pic_iters, 10) == 0) then
               call check_batch_queue_end(exit_flag)
               exit_loop = exit_loop .or. exit_flag
            endif
            if ( exit_loop ) exit
            call wait_while_paused
         enddo

! 2026-03-18 cgw should pic_time_end and print_walltime be inside the loop?
! it seems that they are only called at the end of the run
         call pic_time_end(pic_iters)
         call print_walltime("Timestep walltime, PIC solver:  ", step_start_time)
      end subroutine run_pic

    end subroutine run_mfix

subroutine print_walltime(label, step_start_time)
   use error_manager, only: err_msg, log_message, loglevel_info
   use time_mod, only: get_elapsed_time
   use output_manager, only:  cgdem, superdem, update_time_step_info_float

   implicit none
   character(len=*), intent(in) :: label
   double precision, intent(in) :: step_start_time
   double precision :: step_time

   step_time = get_elapsed_time() - step_start_time
#ifndef QUIET
   write(err_msg, "(A, f9.3, A)" ) label, step_time, " s"
   call log_message(__FILE__, __LINE__, LOGLEVEL_INFO, header=.False., footer=.False.)
#endif

! Caution: The label has now a long form, the type of solver is determined
! from characters 20 to 22 of the label.
! This will need to be adjusted if the label changes
! XXX FIXME

   select case(trim(label(20:22)))
   case('flu')
       call update_time_step_info_float(10,step_time)

   case('DEM')
       call update_time_step_info_float(11,step_time)
       if(cgdem) call update_time_step_info_float(13,step_time)
       if(superdem) call update_time_step_info_float(14,step_time)

   case('PIC')
       call update_time_step_info_float(12,step_time)

   case('CGDEM')
       call update_time_step_info_float(13,step_time)

   end select

end subroutine print_walltime


program mfix

      use compar, only: mype, pe_io
      use error_manager, only: loglevel, loglevel_error, loglevel_info, loglevel_status, loglevel_warning
      use exit, only: mfix_exit
      use main, only: finalize
      use fs_util, only: file_exists
      use funits, only: file_size
      use iso_c_binding, only: c_bool, c_char
      use parallel_mpi, only: parallel_init

      implicit none
      character(len=4096) :: mfix_dat_filename = "mfix.dat"
      integer :: mfix_dat_size
      integer :: ii

      call parse_command_line_args

      if (.not. file_exists(mfix_dat_filename)) then
         if(mype == pe_io) write(*, 1000) trim(mfix_dat_filename)
1000     format(2/,1x,70('*')/' From: mfix.f',/' Error 1000: ',    &
              'Input data file does not exist: ',A/,'Aborting.',/1x,   &
              70('*'),2/)

         call mfix_exit(mype, .true.)
      endif

      mfix_dat_size = file_size(mfix_dat_filename)

      ! Invoke MPI/SMP initialization and get rank info.
      call parallel_init

      call run_mfix(mfix_dat_filename, mfix_dat_size)
      call finalize
   contains

   subroutine parse_command_line_args

      use main, only: add_command_line_keyword, print_flags
      use run, only: id_version

      implicit none

      logical :: reading_mfix_dat = .false.
      logical :: reading_log_level = .false.
      character(len=80) :: mfix_log_level = "INFO"
      character(len=1000) :: tmp

      do ii=1, command_argument_count()

         if (reading_mfix_dat) then
            call get_command_argument(ii, mfix_dat_filename)
            reading_mfix_dat = .false.
            cycle

         else if (reading_log_level) then
            call get_command_argument(ii, mfix_log_level)
            reading_log_level = .false.
            cycle

         else
            call get_command_argument(ii,tmp)
         endif

         if (tmp=='-h'.or.tmp=='--help') then
            call print_usage
            stop
         else if (tmp=='-v'.or.tmp=='--version') then
            print *, id_version
            stop
         else if (tmp=='-p'.or.tmp=='--print-flags') then
            call print_flags
            stop
         else if (tmp=='-f'.or.tmp=='--file') then
            reading_mfix_dat = .true.
            cycle
         else if (tmp=='-l'.or.tmp=='--log') then
            reading_log_level = .true.
            cycle
         else if (index(tmp, '-')==1 .or. index(tmp, '=')==0) then
            ! argument does not start with - or contain =
            print *, "Unknown option: ", tmp
            call print_usage
#ifdef __gfortran__
      call disable_backtrace
#endif
            error stop
         else
            call add_command_line_keyword(tmp)
         endif
      enddo

      ! TODO get rid of this nonsense
      if (mfix_log_level == "ERROR") then
         loglevel = loglevel_error
      else if (mfix_log_level == "WARNING") then
         loglevel = loglevel_warning
      else if (mfix_log_level == "STATUS") then
         loglevel = loglevel_status
      else if (mfix_log_level == "INFO") then
         loglevel = loglevel_info
      else
         print *, "Invalid log level: ", mfix_log_level
         print *, "Must be one of: ERROR, WARNING, STATUS, INFO"
         stop
      end if

   end subroutine parse_command_line_args

   subroutine print_usage

      character(len=1000) :: solver_executable

      call get_command_argument(0, solver_executable)

      print *, "Usage: ", trim(solver_executable), " [-h,--help] [-p,--print-flags] &
         &[-f,--file <filename>] [<KEYWORD>=<VALUE> ...]"

      print *, "       -h,--help: display this help message"
      print *, "       -p,--print-flags: print flags solver was &
         &built with, such as: dmp mkl python smp"
      print *, "       -f,--file <filename>: &
         &specify filename of input file (Default: mfix.dat)"
      print *, "       <KEYWORD>=<VALUE>: specify keyword on &
         &command line, overrides values in input file (mfix.dat or <RUN_NAME>.mfx)"
      print *, "       -l,--log <loglevel>: log message level, one of {ERROR, WARNING, STATUS, INFO}, default INFO"
      print *, "       -v,--version: print version info"
   end subroutine print_usage

end program mfix

module sub
  use func
  use,intrinsic :: iso_fortran_env
  use floating_point_kinds, only : zero, one
  implicit none
  public :: read_prm, read_mat, output, discard
contains
subroutine read_prm(nin, omg, omax, verbose, rmax, tol, at, directory)
  implicit none

  real(real64), intent(inout) :: omg, tol
  integer, intent(inout) :: at, nin, omax, rmax, verbose
  character(*), intent(out) :: directory
  integer :: i, ierror, length, narg, satus
  integer :: set_at = 0, set_fi = 0, set_nin = 0, set_omg = 0, set_tol = 0, set_omax = 0, set_rmax = 0
  character*100 argv

  narg = command_argument_count()

  do i = 1, narg
    call get_command_argument(i, argv, length, satus)
    if (satus > 0) then
      error stop "in get_command_argument"
    endif
    if (index(argv(1:5), "--at=") == 1) then        ! Flag for automatic parameter tuning
      read(argv(6:len_trim(argv)), '(i1)', iostat=ierror) at          
      if (ierror /= 0) then          
        error stop "option --at="
      endif
      set_at = 1
      if (verbose == 1) then
        if (at == 0) then
          write(*, '(a)') "  Disabled automatic paramete tuning"
        else
          write(*, '(a)') "  Enabled automatic paramete tuning"
        endif
      endif
    endif
  enddo

  do i = 1, narg
    call get_command_argument(i, argv, length, satus)
    if (satus > 0) then
      error stop "in get_command_argument"
    endif
    if (index(argv(1:5), "--at=") == 1) then        ! Flag for automatic parameter tuning
    elseif (index(argv(1:6), "--nin=") == 1) then  ! Maximum number of inner iterations
      read(argv(7:len_trim(argv)), '(i10)', iostat=ierror) nin
      if (ierror /= 0) then          
        error stop "option --nin="
      endif
      set_nin = 1
      if (verbose == 1) then
        if (at == 0) then
          write(*, '(a, i0)') "  # of inner iterations: ", nin
        else
          write(*, '(a, i0)') "  Maximum # of inner iterations: ", nin
        endif
      endif
    else if (index(argv(1:6), "--omg=") == 1) then  ! Acceleration (relaxation) parameter 
      read(argv(7:len_trim(argv)), '(f16.10)', iostat=ierror) omg
      if (ierror /= 0) then          
        error stop "option --omg="
      endif
      set_omg = 1
      if (verbose == 1) then
        if (at == 0) then
          write(*, '(a, f0.2)') "  Acceleration parameter (omega): ", omg
        else
          write(*, '(a, f0.2)') "  Initial acceleration parameter (omega): ", omg
          endif
      endif
    else if (index(argv(1:6), "--tol=") == 1) then  ! Stopping criterion
      read(argv(7:len_trim(argv)), '(e16.10e2)', iostat=ierror) tol
      if (ierror /= 0) then          
        error stop "option --tol="
      endif
      set_tol = 1
      if (verbose == 1) then
        write(*, '(a, es8.2e2)') "  Stopping criterion: ", tol
      endif
    else if (index(argv(1:7), "--omax=") == 1) then ! Maximum # of outer iterations
       read(argv(8:len_trim(argv)), '(i10)', iostat=ierror) omax
       if (ierror /= 0) then          
          error stop "option --omax="
       endif
       set_omax = 1
       if (verbose == 1) then
        write(*, '(a, i0)') "  Maximum # of outer iterations: ", omax
      endif
    else if (index(argv(1:7), "--rmax=") == 1) then ! Maximum # of restart cycles
      read(argv(8:len_trim(argv)), '(i10)', iostat=ierror) rmax
      if (ierror /= 0 ) then          
        error stop "option --rmax="
      endif
      set_rmax = 1
      if (verbose == 1) then
        write(*, '(a, i0)') "  Maximum # of restarts: ", rmax
      endif
    else if (index(argv(1:2), "-v") == 1) then          
    else if (index(argv(1:12), "--directory=") == 1) then      
      read(argv(13:len_trim(argv)), '(a)', iostat=ierror) directory
      if (ierror /= 0 ) then          
        error stop "option --directory="
      endif
      set_fi = 1
      if (verbose == 1) then
        write(*, '(a, a)') "  Directory: ", trim(directory)
      endif
    else
      error stop "Unknown argument"
    endif
  enddo 

! Set the default values of parameters (if not set above)
  if (set_at == 0) then
    at = 1
    if (verbose == 1) then
      write(*, '(a)') "  Enebled automatic paramete tuning (default)"
    endif
  endif

  if (set_nin == 0) then
    nin = 50
    if (verbose == 1) then
      write(*, '(a, i0)') "  Maximum # of inner iterations (default): ", nin
    endif
  endif

  if (set_omg == 0) then
    omg = 1.0d0
    if (verbose == 1) then
      write(*, '(a, f0.2)') "  Initial Value of omega (default): ", omg
    endif
  endif    

  if (set_tol == 0) then
    tol = 1.0d-8
    if (verbose == 1) then
      write(*, '(a, e16.10e2)') "  Stoppin  g criterion (default)", tol
    endif
  endif

  if (set_omax == 0) then
    omax = 800
    if (verbose == 1) then
      write(*, '(a, i0)') "  Maximum # of outer iterations (default): ", omax
    endif
  endif

  if (set_rmax == 0) then
    rmax = 0
    if (verbose == 1) then
      write(*, '(a, i0)') "  Maximum # of restarts (default): ", rmax
    endif
  endif  

  if (set_fi == 0) then
    directory = "RANDL7T"
    if (verbose == 1) then      
      write(*, '(a, a)') "  Directory: ", directory
    endif
  endif   

end subroutine read_prm
!---------------------------------------------------------------------------

subroutine read_mat(directory, AR, ja, ip, m, n, x0, b)  
  implicit none

  character(*), intent(in) :: directory
  real(real64), allocatable, intent(inout) :: AR(:), b(:), x0(:)

  integer, allocatable, intent(out) :: ja(:), ip(:)
  integer, intent(inout) :: m, n
  integer i, is, k1, k2, l, nnz
  integer :: fi_AR = 14, fi_ip = 15, fi_ja = 16  

  open(fi_AR, file = trim(directory)//'/AR.crs', action = 'read', iostat = is, status="old")
  if (is /= 0) error stop 'cannot open file AR.crs'
  read(fi_AR, *) m, n, nnz

  allocate(AR(nnz))
  allocate(ja(nnz))
  allocate(ip(m+1))
  allocate(b(m))
  allocate(x0(n))

! Load AC
  read(fi_AR, *) (AR(l), l = 1, nnz)
  close(fi_AR)

! Load jp
  open(fi_ip, file = trim(directory)//'/ip.crs', action = 'read', iostat = is, status="old")
  if (is /= 0) stop 'cannot open ip file'
  read(fi_ip, *) (ip(i), i = 1, m+1)
  close(fi_ip)

! Load ia
  open(fi_ja, file = trim(directory)//'/ja.crs', action = 'read', iostat = is, status="old")
  if (is /= 0) stop 'cannot open ja file'
  read(fi_ja, *) (ja(l), l = 1, nnz)
  close(fi_ja)

! Load b
  do i = 1, m
    k1 = ip(i)
    k2 = ip(i+1)-1
    b(i) = sum(AR(k1:k2))
  enddo

end subroutine read_mat
!-----------------------------------------------------------

subroutine output(AR, ja, ip, m, n, &
                  b, &
                  nin, omg, &
                  x, &
                  iter, omax, riter, relres, t_tot, &
                  verbose)
  use func
  implicit none

  real(real64), intent(in) :: AR(:), b(:)
  real(real64), intent(in) :: t_tot, x(n)
  real(real64), intent(inout) :: relres(iter)
  real(real64) r(m), ATr(n)
  real(real64), intent(in) :: omg
  real(real64) nrm_b, nrm_r, tmp, nrmATr, nrmATb

  integer, intent(in) :: ja(:), ip(:)
  integer, intent(in) :: iter, m, n, nin, omax, verbose, riter  
  integer :: info = 10, reshis = 11, sol = 12
  integer i, is, j, l, k, k1, k2

  nrm_b = nrm2(b(1:m), m)

  ATr(1:n) = zero
  do i = 1, m
    tmp = b(i)
    do l = ip(i), ip(i+1)-1
      j = ja(l)
      ATr(j) = ATr(j) + tmp*AR(l)
    enddo
  enddo

  nrmATb = nrm2(ATr(1:n), n)

   do i = 1, m
      k1 = ip(i)
      k2 = ip(i+1)-1
      r(i) = sum(AR(k1:k2)*x(ja(k1:k2)))
  enddo

  r(1:m) = b(1:m) - r(1:m)

  nrm_r = nrm2(r(1:m), m)

  ATr(1:n) = zero
  do i = 1, m
    tmp = r(i)
    do l = ip(i), ip(i+1)-1
      j = ja(l)
      ATr(j) = ATr(j) + tmp*AR(l)
    enddo
  enddo

  nrmATr = nrm2(ATr(1:n), n)

  tmp = one / nrm_b
  relres(1:iter) = tmp * relres(1:iter)

  if (verbose == 1) then

    write(*, '(a)') 'Results:'
    write(*, '(a, f6.3)') '  omega: ', omg
    write(*, '(a, i0)') '  # of outer iterations: ', iter
    write(*, '(a, i0)') '  # of inner iterations: ', nin
    write(*, '(a, i0)') '  # of restarts: ', riter
    write(*, '(a, i0)') '  Restart cycle: ', omax
    write(*, '(a, f0.2)') '   CPU time: ', t_tot
    write(*, '(a, es8.2e2)') '  Relative residual: ', relres(iter)
    write(*, '(a, es8.2e2)') ' Actual relative residual (ATr): ', nrmATr / nrmATb
    write(*, '(a, es8.2e2)') ' Actual relative residual (r): ', nrm_r / nrm_b    

  else

    write(*, '(a)') 'Results (omega time iter):'
    write(*, '(f0.2, a, f0.2, a, i0)') omg, ' ', t_tot, ' ', iter

  endif

  open(info, file='info.dat', action='write', iostat=is, status='replace')
  write(info, '(a, f6.3)') '  omega: ', omg
  write(info, '(a, i0)') '  # of outer iterations: ', iter
  write(info, '(a, i0)') '  # of inner iterations: ', nin
  write(info, '(a, i0)') '  # of restarts: ', riter
  write(info, '(a, i0)') '  Restart cycle: ', omax
  write(info, '(a, f0.2)') '   CPU time: ', t_tot
  write(info, '(a, es8.2e2)') '  Relative residual: ', relres(iter)
  write(info, '(a, es8.2e2)') ' Actual relative residual (ATr): ', nrmATr / nrmATb
  write(info, '(a, es8.2e2)') ' Actual relative residual (r): ', nrm_r / nrm_b 
  close(info)

  open(reshis, file='reshis.dat', action='write', iostat=is, status='replace')
  if (is /= 0) stop 'cannot open reshis.dat file'
  do k = 1, iter
    write(reshis, *) k, relres(k)
  enddo
  close(reshis)

  open(sol, file='solution.dat', action='write', iostat=is, status='replace')
  if (is /= 0) stop 'cannot open solution.dat file'
  do j = 1, n
    write(sol, *) x(j)
  enddo
  close(sol)

end subroutine output

subroutine discard(AC, ia, jp, b)
  implicit none
  real(real64), allocatable, intent(inout) :: AC(:), b(:)
  integer,  allocatable, intent(inout) :: ia(:), jp(:)

  deallocate(AC, jp, ia)
  deallocate(b)

end subroutine discard
end module sub
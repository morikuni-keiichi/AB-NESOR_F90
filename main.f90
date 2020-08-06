program main
  use, intrinsic :: iso_fortran_env
  use floating_point_kinds, only : zero
  use sub
  use solver
  implicit none

  real(real64), allocatable :: b(:), relres(:), x(:), x0(:), AR(:)
  integer,  allocatable :: ja(:), ip(:)
  real(real64) omg, t1, t2, t_tot, tol
  integer :: at, conv, err, i, iter, length, logcsv = 20, m, n, narg, nin, omax, verbose = 0, riter, rmax, status
  character*100 argv, directory    

  narg = command_argument_count()
  do i = 1, narg
    call get_command_argument(i, argv, length, status)
    if (index(argv(1:2), "-v") == 1) then
      verbose = 1  
    endif
  enddo

! Read parameters  
  if (verbose == 1) then
    write(*, '(a)') "Started setting the values of parameters"
  endif
  call read_prm(nin, omg, omax, verbose, rmax, tol, at, directory)
  if (verbose == 1) then
    write(*, '(a)') "Finished setting the values of parameters"
    write(*, *)
  endif 

! Read matrix data
  call read_mat(directory, AR, ja, ip, m, n, x0, b)

  allocate(x(n), stat=err)
  if (err /= 0) then
    error stop "allocate x"
  endif

  allocate(relres(omax*(rmax+1)))
  if (err /= 0) then
    error stop "allocate relres"
  endif

  open(logcsv, file='log.csv', position='append')

  t_tot = zero

  write(*, '(a)') 'AB-GMRES started' 

  call cpu_time(t1)

  call ABGMRES(AR, ja, ip, m, n, &
                x0, b, &
                at, nin, omg, &
                tol, omax, rmax, &
                x, &
                iter, riter, relres, conv, verbose)

  call cpu_time(t2)

  if (conv == 1) then
    write(*, '(a)') "AB-GMRES converged"
  else
    write(*, '(a)') 'AB-GMRES failed to converge'
  endif  
  write(*, *) 

  t_tot = t_tot + t2 - t1
  
  call output(AR, ja, ip, m, n, &
              b, &
              nin, omg, & 
              x, &
              iter, omax, riter, relres, t_tot, &
              verbose)

  write(logcsv, '(i9, a1, i2, a1, f4.1, a1, f16.4)') Iter, ',', nin, ',', omg, ',', t_tot
  close(logcsv)

  call discard(AR, ja, ip, b)

end program main
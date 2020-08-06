module solver
  use func
  use,intrinsic :: iso_fortran_env
  use floating_point_kinds
  implicit none
  private :: NESOR, atNESOR
  public :: ABGMRES
contains
  subroutine ABGMRES(AR, ja, ip, m, n, &
                      x0, b, &
                      at, nin, omg, &
                      tol, omax, rmax, &
                      x, &
                      iter, riter, relres, conv, verbose)

  real(real64) H(omax+1, omax), V(m, omax+1)
  real(real64), intent(in) :: b(:), AR(:)
  real(real64), intent(out) :: relres(:), x(:)  
  real(real64) c(omax), g(omax+1), r(m), s(omax), &
               w(m), x0(n), y(omax), ATei(m)  
  real(real64), intent(out) :: tol, omg  
  real(real64) :: beta, eps, tmp

  integer, intent(in) :: ja(:), ip(:)
  integer, intent(in) :: at, m, n, omax, rmax, verbose
  integer, intent(out) :: nin
  integer, intent(out) :: iter, riter, conv
  integer :: bd = 0, i, iter_tot = 1, j, k, k1, k2, l, p

  eps = epsilon(eps)
  conv = 0

!  The initial approximate solution is equal to zero
  x0(1:n) = zero

! The 2-norm of a_i, i=1,...,m, required for NE-SOR
  do i = 1, m
    k1 = ip(i)
    k2 = ip(i+1)-1
    tmp = sum(AR(k1:k2) * AR(k1:k2))
    if (tmp /= zero) then
      ATei(i) = one / tmp
    else
      write(*, *) 'warning: ||ai|| = 0.0 at i =', i
      stop
    endif
  enddo

! Tol = cri_res * ||b||_2
  tol = tol * nrm2(b(1:m), m)

  do p = 0, rmax

    do i = 1, m
      k1 = ip(i)
      k2 = ip(i+1)-1
      r(i) = sum(AR(k1:k2)*x0(ja(k1:k2)))
    enddo

  ! r0 = b - A*x0
    r(1:m) = b(1:m) - r(1:m)

    beta = nrm2(r(1:m), m)

    g(1) = beta

  ! Normalization of the vector v_1
    V(1:m, 1) = (one/beta) * r(1:m)

    do k = 1, omax

      if (k > 1) then        
      ! NE-SOR inner-iteration preconditioning without automatic parameter tuning              
        call NESOR(AR, ja, ip, m, n, V(1:m, k), ATei, nin, omg, x(1:n))
      else if (at > 0) then        
        if (verbose == 1) then
          write(*, '(a)') '  Automatic parameter tuning started'
        endif
        call atNESOR(AR, ja, ip, m, n, V(1:m, k), ATei, nin, omg, x)
        if (verbose == 1) then
          write(*, '(a)') '  Automatic parameter tuning finished'
        endif
      else
      ! NE-SOR inner-iteration preconditioning without automatic parameter tuning        
        call NESOR(AR, ja, ip, m, n, V(1:m, k), ATei, nin, omg, x(1:n))        
      endif
      
    ! w = A*zj
      do i = 1, m
        k1 = ip(i)
        k2 = ip(i+1)-1
        w(i) = sum(AR(k1:k2) * x(ja(k1:k2)))
      enddo    

    ! Modified Gram-Schmidt orthogonalization
      do l = 1, k
        tmp  = sum(w(1:m) * V(1:m, l))
        w(1:m) = w(1:m) - tmp*V(1:m, l)
        H(l, k) = tmp
      enddo      

      tmp = nrm2(w(1:m), m)

      if (tmp > eps) then
      ! Normalization of the vector v_{k+1}
        H(k+1, k) = tmp        
        V(1:m, k+1) = (one/tmp) * w(1:m)
      else
        write(*, *) 'BREAKDOW at step', k
        write(*, *) 'h_k+1, k =', H(k+1, k)
        bd = 1
      endif

    ! Application of Givens rotation
      do i = 1, k-1
        tmp =  c(i)*H(i, k) + s(i)*H(i+1, k)
        H(i+1, k) = -s(i)*H(i, k) + c(i)*H(i+1, k)
        H(i, k) = tmp
      enddo

    ! Construction of Givens rotations
      call rotg(H(k, k), H(k+1, k), c(k), s(k))

      H(k, k) = one / H(k, k)

      g(k+1) = -s(k) * g(k)
      g(k)   =  c(k) * g(k)

      relres(iter_tot) = abs(g(k+1))

    !  Convergence check
      if (relres(iter_tot)<=Tol .or. bd == 1) then

        !	Derivation of the approximate solution x_k
        !	Backward substitution
        y(k) = g(k) * H(k, k)
        do i = k-1, 1, -1
          y(i) = (g(i) - sum(H(i, i+1:k) * y(i+1:k))) * H(i, i)
        enddo

        w = matmul(V(:, 1:k), y(1:k))

      ! NE-SOR inner-iteration preconditioning        
        call NESOR(AR, ja, ip, m, n, w, ATei, nin, omg, x)    

        x(1:n) = x0(1:n) + x(1:n)

        if (bd == 1) then
          iter = iter_tot
          riter = p
          return
        endif

        do i = 1, m
          k1 = ip(i)
          k2 = ip(i+1)-1
          r(i) = sum(AR(k1:k2)*x(ja(k1:k2)))
        enddo

        ! r_k = b - A*x_k
        r(1:m) = b(1:m) - r(1:m)

        ! write(*, *) nrm2(r(1:m), m), tol

        if (nrm2(r(1:m), m) <= tol) then
          
          conv = 1
          iter = iter_tot
          Riter = p          

          return

        endif

      endif

      iter_tot = iter_tot + 1

    enddo

    y(omax) = g(omax) * H(omax, omax)
    do i = omax-1, 1, -1
      y(i) = (g(i) - sum(H(i, i+1:omax) * y(i+1:omax))) * H(i, i)
    enddo

    w = matmul(V(:, 1:omax), y(1:omax))

  ! NE-SOR inner-iteration preconditioning    
    call NESOR(AR, ja, ip, m, n, w, ATei, nin, omg, x)            

    x0(1:n) = x0(1:n) + x(1:n)

    iter = iter_tot
    Riter = p  

  enddo

  x = x0

  iter = iter_tot - 1
  Riter = p - 1

  end subroutine ABGMRES

             
  subroutine NESOR(AR, ja, ip, m, n, rhs, ATei, nin, omg, x)

  real(real64), intent(in) :: AR(:), ATei(:)
  real(real64), intent(out) :: x(:)
  real(real64), intent(in) :: rhs(:)
  integer, intent(in) :: ja(:), ip(:)
  integer, intent(in) :: m, n, nin
  real(real64), intent(in) :: omg
  real(real64) d
  integer i, j, k, l, k1, k2

  x(1:n) = zero

  do k = 1, nin
    do i = 1, m
      k1 = ip(i)
      k2 = ip(i+1)-1
      d = omg*(rhs(i) - sum(AR(k1:k2)*x(ja(k1:k2))))*ATei(i)
      do l = k1, k2
        j = ja(l)
        x(j) = x(j) + d*AR(l)
      enddo
    enddo
  enddo

  end subroutine NESOR


  subroutine atNESOR(AR, ja, ip, m, n, rhs, ATei, nin, omg, x)
  
  real(real64), intent(in) :: AR(:), ATei(:)  
  real(real64), intent(out) :: x(:)
  real(real64), intent(inout) :: rhs(:)
  integer, intent(in) :: ja(:), ip(:)
  integer, intent(in) :: m, n
  integer, intent(out) :: nin
  real(real64), intent(out) :: omg
  real(real64) :: d, e, r(m), norm_r1, norm_r2 = zero, y(n)
  integer i, j, k, k1, k2, kk, l

  omg = one

  x(1:n) = zero
  y(1:n) = zero

  do k = 1, 50

    do i = 1, m
      k1 = ip(i)
      k2 = ip(i+1)-1
      d = omg*(rhs(i) - sum(AR(k1:k2)*x(ja(k1:k2))))*ATei(i)
      do l = k1, k2
        j = ja(l)
        x(j) = x(j) + d*AR(l)
      enddo
    enddo

    d = maxval(abs(x(1:n)))
    e = maxval(abs(x(1:n) - y(1:n)))

    if (e < 1.0d-1 * d .or. i == 50) then
      nin = k
      exit
    endif

    y(1:n) = x(1:n)

  enddo

  do k = 19, 1, -1

    omg = 1.0d-1 * dble(k)

    x(1:n) = zero

    do kk = 1, nin

      do i = 1, m
        k1 = ip(i)
        k2 = ip(i+1)-1
        d = omg*(rhs(i) - sum(AR(k1:k2)*x(ja(k1:k2))))*ATei(i)
        do l = k1, k2
          j = ja(l)
          x(j) = x(j) + d*AR(l)
        enddo
      enddo
    enddo

    do i = 1, m
      k1 = ip(i)
      k2 = ip(i+1)-1
      r(i) = sum(AR(k1:k2)*x(ja(k1:k2)))
    enddo

    r(1:m) = rhs(1:m) - r(1:m)

    norm_r1 = nrm2(r(1:m), m)

    if (k < 19) then
      if (norm_r1 > norm_r2) then
        omg = omg + 1.0d-1
        x(1:n) = y(1:n)
        return
      elseif (k == 1) then
        omg = 0.1d0
        return
      endif
    endif

    norm_r2 = norm_r1

    y(1:n) = x(1:n)

  enddo

  end subroutine atNESOR
end module solver
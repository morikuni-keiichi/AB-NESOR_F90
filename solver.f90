module solver
	use globvar
	use func 
	implicit none
	private
	public :: ABGMRES
contains
	subroutine ABGMRES(Iter, Riter, RelRes, x)
	
	real(dp) H(omax+1, omax), V(m, omax+1)
	real(dp), intent(out) :: RelRes(:), x(:)
	real(dp) c(omax), g(omax+1), r(m), s(omax), w(m), x0(n), y(omax)
	real(dp) :: beta, deps = 1.0d-15, tmp, Tol 
	integer, intent(out) :: Iter, Riter
	integer :: bd = 0, i, iter_tot = 1, k, k1, k2, l, p

!	The initial approximate solution is equal to zero
	x0(1:n) = zero

! The 2-norm of a_i, i=1,...,m, required for NE-SOR
	do i = 1, m
		k1 = ip(i)
		k2 = ip(i+1)-1
		tmp = sum(AR(k1:k2)*AR(k1:k2))
		if (tmp /= zero) then
			ATei(i) = one / tmp 
		else
			write(*, *) 'warning: ||ai|| = 0.0 at i =', i
			stop
		endif
	enddo

! Tol = cri_res * ||b||_2
	Tol = cri_res * nrm2(b(1:m), m) 

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
		tmp = one / beta
		V(1:m, 1) = tmp * r(1:m)

		do k = 1, omax
	
			if (k > 1) then
			! NE-SOR inner-iteration preconditioning without automatic parameter tuning 
				call NESOR(V(1:m, k), x(1:n))
			elseif (at > 0) then
				write(*, *) 'Automatic NE-SOR inner-iteration parameter tuning'
				call atNESOR(V(1:m, k), x(1:n))
				write(*, *) 'Tuned'
			else 
			! NE-SOR inner-iteration preconditioning without automatic parameter tuning 
				call NESOR(V(1:m, k), x(1:n))
			endif
			
		! w = A*zj
			do i = 1, m
				k1 = ip(i)
				k2 = ip(i+1)-1
				w(i) = sum(AR(k1:k2)*x(ja(k1:k2)))
			enddo
	
		! Modified Gram-Schmidt orthogonalization
			do l = 1, k
				H(l, k)	= sum(w(1:m)*V(1:m, l)) 
				w(1:m) = w(1:m) - H(l, k)*V(1:m, l)
			enddo
	
			H(k+1, k) = nrm2(w(1:m), m) 
	
		!	write(*, *) k, H(k+1, k)
	
			if (H(k+1, k) > deps ) then
			! Normalization of the vector v_{k+1}
				tmp = one / H(k+1, k)
				V(1:m, k+1) = tmp * w(1:m)
			else
				write(*, *) 'BREAKDOW at step', k
				write(*, *) 'h_k+1, k =', H(k+1, k)
				bd = 1
			endif
	
		! Application of Given's rotation
			do l = 1, k-1
				tmp =  c(l)*H(l, k) + s(l)*H(l+1, k)
				H(l+1, k) = -s(l)*H(l, k) + c(l)*H(l+1, k)
				H(l, k) = tmp
			enddo
		
		! Construction of Givens rotations
			call rotg(H(k, k), H(k+1, k), c(k), s(k))
			
			g(k+1) = -s(k) * g(k)
			g(k)	 =  c(k) * g(k)
	
			RelRes(iter_tot) = abs(g(k+1)) 
	
		!	write(*, *) k, RelRes(k) / beta 
	
			H(k, k) = one / H(k, k)
	
		!	Convergence check
			if ((RelRes(iter_tot))<=Tol .or. (k == omax .and. p == rmax) .or. bd == 1) then
				y(k) = g(k) * H(k, k)
				do l = k-1, 1, -1
					y(l) = (g(l) - sum(H(l, l+1:k)*y(l+1:k))) * H(l, l)
				enddo
	
				w(1:m) = matmul(V(1:m, 1:k), y(1:k))
		
			! NE-SOR inner-iteration preconditioning	
				call NESOR(w(1:m), x(1:n))

				x(1:n) = x0(1:n) + x(1:n)
				
				Iter = iter_tot
				Riter = p
		
				return
			endif

			iter_tot = iter_tot + 1
			
		enddo

		y(omax) = g(omax) * H(omax, omax)
		do l = omax-1, 1, -1
			y(l) = (g(l) - sum(H(l, l+1:omax)*y(l+1:omax))) * H(l, l)
		enddo
	
		w(1:m) = matmul(V(1:m, 1:omax), y(1:omax))
	
	! NE-SOR inner-iteration preconditioning	
		call NESOR(w(1:m), x(1:n))
		
		x0(1:n) = x0(1:n) + x(1:n)

	enddo

	end subroutine ABGMRES


	subroutine NESOR(rhs, x)

	real(dp), intent(out) :: x(:)
	real(dp), intent(in) :: rhs(:)
	real(dp) d
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


	subroutine atNESOR(rhs, x)
	
	real(dp), intent(inout) :: rhs(:)
	real(dp), intent(out) :: x(:)
	real(dp) :: d, e, r(m), norm_r1, norm_r2 = zero, y(n)
	integer i, j, k, k1, k2, kk, l

	omg = one 
	
	x(1:n) = zero
	y(1:n) = zero

	do k = 1, 100

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

		if (e < 1.0d+1 ** (-1.0d0) * d .or. i == 100) then
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

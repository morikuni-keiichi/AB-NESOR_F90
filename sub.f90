subroutine read_prm()
	use globvar
	implicit none

	integer :: fi = 11, is
	
	open(fi, file='prm.dat', action='read', iostat = is)
	if (is /= 0) stop 'cannot open a parameter file'
	read(fi, *) cri_res
	read(fi, *) omg
	read(fi, *) omax
	read(fi, *) nin 
	read(fi, *) rmax
	read(fi, *) at
	read(fi, *) omode
	close(fi)

end subroutine read_prm
!---------------------------------------------------------------------------

subroutine read_mat()
	use globvar
	implicit none
	
	integer i, is, k1, k2, l
	integer :: fo_AR = 11, fo_ja = 12, fo_ip = 13

	open(fo_AR, file = 'RANDL7T/AR.crs', action = 'read', iostat = is)
	if (is /= 0) stop 'cannot open AR file'
	read(fo_AR, *) m, n, nnz
	
	allocate(AR(nnz))
	allocate(ja(nnz))
	allocate(ip(m+1))
	allocate(b(m))
	allocate(ATei(m))

! Load AR
	read(fo_AR, *) (AR(l), l = 1, nnz)
	close(fo_AR)

! Load ja
	open(fo_ja, file ='RANDL7T/ja.crs', action = 'read', iostat = is)
	if (is /= 0) stop 'cannot open ja file'
	read(fo_ja, *) (ja(l), l = 1, nnz)
	close(fo_ja)

! Load ip
	open(fo_ip, file = 'RANDL7T/ip.crs', action = 'read', iostat = is)
	if (is /= 0) stop 'cannot open ip file'
	read(fo_ip, *) (ip(i), i = 1, m+1)
	close(fo_ip)

! Load b
	do i = 1, m
		k1 = ip(i)
		k2 = ip(i+1)-1
		b(i) = sum(AR(k1:k2))
	enddo

end subroutine read_mat 
!-----------------------------------------------------------

subroutine output(Iter, Riter, t_tot, RelRes, x)
	use func
	use globvar
	implicit none

	real(dp), intent(in) :: t_tot, x(n)
	real(dp), intent(inout) :: RelRes(Iter)
	real(dp) r(m), ATr(n)
	real(dp) norm_b, norm_r, tmp, norm_ATr, norm_ATb
	integer, intent(in) :: Iter, Riter
	integer :: info = 10, reshis = 11, sol = 12 
	integer i, is, j, l, k, k1, k2

	norm_b = nrm2(b(1:m), m) 

	ATr(1:n) = zero
	do i = 1, m
		tmp = b(i)
		do l = ip(i), ip(i+1)-1
			j = ja(l)
			ATr(j) = ATr(j) + tmp*AR(l)
		enddo
	enddo

	norm_ATb = nrm2(ATr(1:n), n) 
	
	do i = 1, m
	 k1 = ip(i)
	 k2 = ip(i+1)-1
	 r(i) = sum(AR(k1:k2)*x(ja(k1:k2))) 
	enddo

	r(1:m) = b(1:m) - r(1:m)

	norm_r = nrm2(r(1:m), m) 

	ATr(1:n) = zero
	do i = 1, m
		tmp = r(i)
		do l = ip(i), ip(i+1)-1
			j = ja(l)
			ATr(j) = ATr(j) + tmp*AR(l)
		enddo
	enddo

	norm_ATr = nrm2(ATr(1:n), n) 

	tmp = one / norm_b
	RelRes(1:Iter) = tmp * RelRes(1:Iter)

	if (omode == 0) then

		write(*, '(a, f6.3)') '  omega: ', omg
		write(*, *) ' # of outer iterations: ', Iter
		write(*, *) ' # of inner iterations: ', nin
		write(*, *) ' # of restarts: ', Riter
		write(*, *) ' Restart cycle: ', omax
		write(*, '(a, f16.5)') '  CPU time: ', t_tot
		write(*, *) ' Relative residual: ', RelRes(Iter)
		write(*, *) ' Actual relative residual (ATr): ', norm_ATr / norm_ATb
		write(*, *) ' Actual relative residual (r): ', norm_r / norm_b
		write(*, *) " "

	else

		write(*, '(f6.3, f16.5, i9)') omg, t_tot, Iter
	
	endif

	open(info, file='info.dat', action='write', iostat=is, status='replace')
	write(info, '(a, f6.3)') '  omega: ', omg
	write(info, *) ' # of outer iterations: ', Iter
	write(info, *) ' # of inner iterations: ', nin
	write(info, *) ' # of restarts: ', Riter
	write(info, *) ' Restart cycle: ', omax		
	write(info, '(a, f16.5)') '  CPU time: ', t_tot
	write(info, *) ' Relative residual: ', RelRes(Iter)
	write(info, *) ' Actual relative residual (ATr): ', norm_ATr / norm_ATb	
	write(info, *) 'Actual relative residual (r): ', norm_r / norm_b
	close(info)
	
	open(reshis, file='reshis.dat', action='write', iostat=is, status='replace')
	if (is /= 0) stop 'cannot open reshis.dat file'	
	do k = 1, Iter
		write(reshis, *) k, RelRes(k) 
	enddo
	close(reshis)
	
	open(sol, file='solution.dat', action='write', iostat=is, status='replace')
	if (is /= 0) stop 'cannot open solution.dat file'	
	do j = 1, n
		write(sol, *) x(j)
	enddo
	close(sol)

end subroutine output

subroutine discard()
	use globvar
	implicit none
	
	deallocate(AR, ja, ip)
	deallocate(b)
	deallocate(ATei)

end subroutine discard


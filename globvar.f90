module globvar
	intrinsic                   ::      selected_real_kind
   integer,  parameter, public :: dp = selected_real_kind(15)
	real(dp), parameter, public :: zero = 0.0_dp, one = 1.0_dp
	real(dp), allocatable, save :: AR(:), b(:), ATei(:)
	integer, allocatable, save :: ja(:), ip(:)
	real(dp), save :: cri_res, omg
	integer, save :: at, m, n, nin, nnz, omax, omode, rmax
end module globvar

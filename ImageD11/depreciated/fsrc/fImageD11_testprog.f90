
 program fImageD11_testprog
    integer, parameter :: n = 10000000
    real,allocatable :: xlylzl(:,:),omega(:)
    real :: wvln,wedge,chi,t(3),dt,ubi(3,3)
    real(8), allocatable :: gv(:,:), drlv2(:)
    integer, allocatable :: labels(:), seed(:)
    integer :: omegasign,ns
    integer(8) :: start, end, rate
    allocate(xlylzl(3,n))
    allocate(omega(n))
    allocate(gv(3,n))
    call random_seed(size=ns)
    allocate(seed(ns))
    seed=42
    call random_seed(put=seed)
    call random_number(xlylzl)
    call random_number(omega)
    omega= omega*180
    t = (/ 0.1,0.2,0.3  /)
    wedge = 1.0
    chi = 2.0
    omegasign=1
    call system_clock(count_rate=rate)
    call system_clock(start)
    call compute_gv( xlylzl, omega, omegasign, wvln, wedge, chi, &
        t, gv, n)
    call system_clock(end)
    dt = real(end-start)/rate
    write(*,*)"Time: ", dt, n/dt/1e6,"Mgvecs computed per second"
    ig=1
    ubi = reshape( (/ 1., 2., 3., 3.,2.,1.,1.,-3.,1. /), (/ 3, 3 /))
    tol=0.01
    allocate(labels(n))
    labels=-1
    allocate(drlv2(n))
    drlv2 = 10.
    call system_clock(start)
    call assign( ubi, gv, tol, drlv2, labels, ig, n)
    call system_clock(end)
    dt = real(end-start)/rate
    write(*,*)"Time: ", dt, n/dt/1e6,"Mgvecs assigned per second"
 end program fImageD11_testprog

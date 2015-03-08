


! compute_gv for refinegrains.py

subroutine compute_gv( xlylzl, omega, wvln, wedge, chi, t, gv, n )
  use omp_lib
  implicit none
  real, intent(in) :: xlylzl(3,n), omega(n), wvln, wedge, chi, t(3)
  real, intent(out):: gv(3,n)
  integer, intent(in) ::  n
  real :: sc,cc,sw,cw,wmat(3,3),cmat(3,3), mat(3,3), u(3),d(3),v(3)
  real :: modyz, o(3), rtth, reta, co, so, ds, sth, cth, k(3)
  real, parameter :: PI=3.141592653589793,RAD=PI/180.0,DEG=180.0/PI
  integer :: i

  ! Fill in rotation matrix of wedge, chi
  sw = sin(wedge*RAD)
  cw = cos(wedge*RAD)
  sc = sin(chi*RAD)
  cc = cos(chi*RAD)
  wmat = RESHAPE ((/ cw,0.,-sw,0.,1.,0.,sw,0.,cw /), (/3,3/))
  cmat = RESHAPE ((/ 1.,0.,0.,0.,cc,sc,0.,-sc,cc /), (/3,3/))
  mat = matmul(cmat, wmat)
  write(*,*)'threads',omp_get_max_threads()
!  write(*,*)'mat',mat
!$omp parallel do private(so,co,u,o,d,modyz,reta,rtth,sth,cth,ds,v,k)
  do i=1,n
     ! Compute translation + rotation for grain origin
     so = sin(RAD*omega(i))
     co = cos(RAD*omega(i))
     u(1) =  co*t(1) - so*t(2)
     u(2) =  so*t(1) + co*t(2)
     u(3) = t(3)
     ! grain origin, difference vec, |yz| component
     o = matmul(transpose(mat),u)
     d = xlylzl(:,i) - o
     modyz  = sqrt(d(2)*d(2) + d(3)*d(3))
     ! tth/eta in radians
     reta = atan2(-d(2),d(3))
     rtth = atan2( modyz, d(1) )
     ! k-vector
     sth = sin(rtth/2.0)
     cth = cos(rtth/2.0)
     ds = 2*sth/wvln
     k(1) = -ds*sth
     k(2) = -ds*sin(reta)*cth
     k(3) =  ds*cos(reta)*cth
     v = matmul(mat,k)
     gv(1,i)= co*v(1) + so*v(2)
     gv(2,i)=-so*v(1) + co*v(2)
     gv(3,i)=v(3)
  enddo
 

end subroutine compute_gv
  
! f2py -m fImageD11 -c fImageD11.f90 --opt=-O3 --f90flags="-fopenmp" -lgomp -lpthread
! export OMP_NUM_THREADS=12
! python tst.py ../test/nac_demo/peaks.out_merge_t200 ../test/nac_demo/nac.prm 

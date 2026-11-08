program mc_box_sphere
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  integer :: ndim(18), n, nprob, boxtype, numeval, itrial, i
  real(dp) :: seed, f, fbest, nrm
  real(dp) :: x(50), xbest(50), l(50), u(50)
  character(len=64) :: arg
  real(dp), external :: drand
  external :: objfcn

  data ndim /3,6,3,2,3,10,6,4,4,2,4,3,10,10,12,2,4,10/

  if (command_argument_count() < 3) then
     write(*,*) 'Uso: ./mc_box_sphere NPROB BOXTYPE NUMEVAL'
     stop 1
  end if
  call get_command_argument(1, arg); read(arg, *) nprob
  call get_command_argument(2, arg); read(arg, *) boxtype
  call get_command_argument(3, arg); read(arg, *) numeval

  if (nprob < 1 .or. nprob > 18) stop 2
  n = ndim(nprob)
  call getbounds(boxtype, n, l, u)

  seed  = 123456.0_dp
  fbest = huge(1.0_dp)
  xbest = 0.0_dp

  do itrial = 1, numeval
     ! 1. sortear un punto en la caja
     do i = 1, n
        x(i) = l(i) + (u(i) - l(i)) * drand(seed)
     end do

     ! 2. proyectarlo a la esfera unidad
     nrm = sqrt(sum(x(1:n)**2))
     if (nrm < 1.0d-12) cycle
     x(1:n) = x(1:n) / nrm

     ! 3. evaluar f
     call objfcn(n, x, f, nprob)
     if (f < fbest) then
        fbest = f
        xbest(1:n) = x(1:n)
     end if
  end do

  write(*,'(a,es24.8e2)') 'fbest = ', fbest
  write(*,'(a,es24.8e2)') '|xbest| = ', sqrt(sum(xbest(1:n)**2))
  write(*,*) 'xbest ='
  write(*,'(es24.8e2)') xbest(1:n)
end program mc_box_sphere

subroutine getbounds(boxtype, n, l, u)
  implicit none
  integer, intent(in) :: boxtype, n
  real(kind=8), intent(out) :: l(50), u(50)

  if (boxtype == 1) then
     l(1:n) = -1.0d0;  u(1:n) = 1.0d0
  else if (boxtype == 2) then
     l(1:n) = -3.0d0;  u(1:n) = 3.0d0
  else if (boxtype == 3) then
     l(1:n) = -5.0d0;  u(1:n) = 5.0d0
  else
     write(*,*) 'Undefined box type!'
     stop 3
  end if

  ! Variante no simetrica (la caja si cambia las direcciones):
  ! l(1:n) = -1.0d0
  ! u(1:n) = 1.0d0 + 2.0d0*(boxtype-1)    ! 1, 3, 5
end subroutine getbounds

function drand(ix)
  implicit none
  real(kind=8) :: drand
  real(kind=8), intent(inout) :: ix
  real(kind=8) :: a,p,b15,b16,xhi,xalo,leftlo,fhi,k
  data a/16807.d0/,b15/32768.d0/,b16/65536.d0/,p/2147483647.d0/
  xhi= ix/b16
  xhi= xhi - dmod(xhi,1.d0)
  xalo= (ix-xhi*b16)*a
  leftlo= xalo/b16
  leftlo= leftlo - dmod(leftlo,1.d0)
  fhi= xhi*a + leftlo
  k= fhi/b15
  k= k - dmod(k,1.d0)
  ix= (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
  if (ix.lt.0) ix= ix + p
  drand= ix*4.656612875d-10
  return
end function drand

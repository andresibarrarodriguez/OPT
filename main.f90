program rsd
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: maxfc = 100000        ! presupuesto de evaluaciones de f
  real(dp), parameter :: eps = 1.0d-6

  integer :: n, nprob, k, fcnt, verbose, maxit
  integer :: ndim(18)
  real(dp) :: x(50), xn(50), g(50), rg(50), d(50)
  real(dp) :: f, fn, t, gn, f0, cost
  real(dp) :: c1, t0, beta
  logical :: conv, fail
  character(len=64) :: arg
  external :: initpt, projinit, objfcn, rgrad, retract

  data ndim /3,6,3,2,3,10,6,10,10,2,4,3,10,10,12,2,4,10/

  ! ---- argumentos: nprob c1 t0 beta maxit [verbose] ----
  if (command_argument_count() < 5) then
     write(*,*) 'Uso: ./rsd NPROB C1 T0 BETA MAXIT [VERBOSE]'
     stop 1
  end if
  call get_command_argument(1, arg); read(arg, *) nprob
  call get_command_argument(2, arg); read(arg, *) c1
  call get_command_argument(3, arg); read(arg, *) t0
  call get_command_argument(4, arg); read(arg, *) beta
  call get_command_argument(5, arg); read(arg, *) maxit

  verbose = 0
  if (command_argument_count() >= 6) then
     call get_command_argument(6, arg); read(arg, *) verbose
  end if

  if (nprob < 1 .or. nprob > 18) stop 2
  if (c1 <= 0.0_dp .or. c1 >= 1.0_dp .or. beta <= 0.0_dp .or. &
      beta >= 1.0_dp .or. t0 <= 0.0_dp .or. maxit < 1) then
     write(*,*) 'Parametros fuera de rango'
     stop 3
  end if
  n = ndim(nprob)

  call initpt(n, x, nprob, 1.0d0)
  call projinit(n, x)

  fcnt = 0
  conv = .false.
  fail = .false.
  f0   = 0.0_dp

  do k = 0, maxit
     call objfcn(n, x, f, nprob);  fcnt = fcnt + 1
     if (k == 0) f0 = f
     call rgrad(n, x, g, rg, nprob)
     gn = sqrt(sum(rg(1:n)**2))

     if (gn < eps) then
        conv = .true.
        exit
     end if
     if (k == maxit .or. fcnt >= maxfc) exit

     d(1:n) = -rg(1:n)

     t = t0
     do
        call retract(n, x, d, t, xn)
        call objfcn(n, xn, fn, nprob);  fcnt = fcnt + 1
        if (fn <= f - c1 * t * gn**2) exit
        t = beta * t
        if (t < 1.0d-16 .or. fcnt >= maxfc) then
           fail = .true.
           exit
        end if
     end do
     if (fail) exit

     x(1:n) = xn(1:n)
  end do

  ! recalcular f y gradiente en el punto final
  call objfcn(n, x, f, nprob)
  call rgrad(n, x, g, rg, nprob)
  gn = sqrt(sum(rg(1:n)**2))

  cost = (f - f0) / (abs(f0) + 1.0_dp)

  if (verbose == 1) then
     print '(a,i3)',      ' Problema    = ', nprob
     print '(a,i8)',      ' n           = ', n
     print '(a,i8)',      ' Iteraciones = ', k
     print '(a,i8)',      ' Eval. de f  = ', fcnt
     print '(a,es24.16)', ' f0          = ', f0
     print '(a,es24.16)', ' f           = ', f
     print '(a,es24.16)', ' |grad|      = ', gn
     print '(a,es24.16)', ' costo       = ', cost
     print '(a,es24.16)', ' |x|         = ', sqrt(sum(x(1:n)**2))
     print *, ' x ='
     print '(es24.16)', x(1:n)
     if (conv) then
        print *, ' Convergio'
     else
        print *, ' NO convergio (iteraciones, presupuesto o linea de busqueda)'
     end if
  else
     print '(es24.16)', cost
  end if
end program rsd

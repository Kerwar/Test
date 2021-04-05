program FEM
  use number_types_mod , only: DP 
  use point_mod        , only: Point, idx, swapNodes
  use quadrilateral_mod, only: Square
  use param_mod        , only: Param
  use matrix_mod       , only: GMRES, swapColumns
  use ogpf

  implicit none

! DIRICHLET BOUNDARY CONDITIONS RANT, READ IF IN NEED TO UNDERSTAND
! Let us recall that at the end we are solving a system of equation Ax=b.
! The nodes that are in the D.B already have they value set so there is no 
! need to calculate its solution. We just need to be careful because they 
! influence the other nodes in their elements. For that we are going to 
! rearrenge the A matrix in the following way:
! [A_int A_int_bc]
! [0       I     ]
! the block A_int contains the matrix of the coefficients of only the interior
! nodes. The block A_int_bc contains the columns of the matrix A for 
! the correspondent node that is in the boundary no change needed, this is 
! because they only have influence in other nodes therefor the rows can be 
! changed to the identity row and b to the Dirichlet values

  type(Point) , allocatable :: nodes(:)
  type(Square), allocatable :: elements(:)
  type(Param)               :: parametr
  type(gpf) :: gp 
  real(DP), allocatable :: xgrid(:), ygrid(:)
  real(DP), allocatable :: x(:,:), y(:,:), z(:,:)
  real(DP), allocatable :: A(:,:), b(:), AA(:)
  real(DP), allocatable :: x0(:), sol(:)
  integer , allocatable :: JA(:), IA(:), IDBC(:)

  real(DP) :: help, tStart, tEnd
  integer :: nEl, nN, nDBC
  integer :: i, j, k, l 

  call cpu_time(tStart) 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PARAMETER SET UP AND GRID SET UP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call parametr%read_param("datainput")

  xgrid = linspace(parametr%xMin, parametr%xMax, parametr%nx)
  ygrid = linspace(parametr%yMin, parametr%yMax, parametr%ny)

  call meshgrid(x, y, xgrid, ygrid)

  nN = parametr%nx * parametr%ny 
  allocate(nodes(nN), z(size(x,1), size(x,2)), IDBC(parametr%ny))
  call cpu_time(tEnd)
  print *, "Mesh dode", tEnd -tStart

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NODES SET UP AND BOUNDARY FLAGGING !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  l = 0
  do i = 1, size(x, 1)
    do j = 1, size(x, 2)
      nodes(idx(i,j)) = Point(x(i,j),y(i,j))
      if (nodes(idx(i,j))%x == parametr%xMin) then
        l = l + 1
        IDBC(l) = idx(i,j)
        nDBC = nDBC + 1
      end if 
    end do
  end do
  call cpu_time(tEnd)
  print *, "Nodes node", tEnd -tStart

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ELEMENTS SET UP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nEl = 0
  allocate(elements((size(x, 1)-1)*(size(x, 2)-1)))
  do i = 1, size(x,1) -1
    do j = 1, size(x,2) -1
      nEl = nEl + 1
      elements(nEl) = Square(nodes(idx(i,j)), nodes(idx(i,j+1)), &
      nodes(idx(i+1,j+1)), nodes(idx(i+1,j)))  
      call elements(nEl)%nodesIndex(idx(i,j), idx(i+1,j), &
      idx(i+1,j+1), idx(i,j+1))
      call elements(nEL)%set_K(parametr)
    end do
  end do
  call cpu_time(tEnd)
  print *, "Elements done", tEnd -tStart

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! MATRIX SET UP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(A(nN, nN), b(nN), x0(nN), sol(nN))

  A = 0.0_DP
  b = 0.0_DP
  ! do i = 1, size(x,1)
  !   do j = 1, size(x,2)
  !     if (nodes(idx(i,j))%x < -20) then
  !       x0 (idx(i,j))= 0.0_DP
  !     else if (nodes(idx(i,j))%x > 20) then
  !       x0 (idx(i,j))= 1.0_DP
  !     else
  !       x0 (idx(i,j))= (2/3*parametr%nx - i ) * 0.0_DP + &
  !       (i - 1/3 * parametr%nx ) * 3/parametr%nx
  !     end if
  !   end do
  ! end do
  x0 = 0.0_DP
  sol = 0.0_DP
  
  
  
  do k = 1, nEl 
    
    do i = 1, 4
      do j = 1, 4
        A(elements(k)%aI(i), elements(k)%aI(j)) = A(elements(k)%aI(i), elements(k)%aI(j)) + &
        elements(k)%K(i,j)
      end do
      b(elements(k)%aI(i)) = b(elements(k)%aI(i)) + elements(k)%b(i)
    end do
  end do
  
  call cpu_time(tEnd)
  print *, "Matrix done", tEnd -tStart

  l = nN
  do k = 1, nDBC 
    call swapColumns(A, b, iDBC(k), l)
    call swapNodes(nodes(IDBC(k)), nodes(l))
    
    A(l,:) = 0.0_DP
    A(l,l) = 1.0_DP
    b(IDBC(k)) = b(l)
    ! Set the right DBC
    b(l) = 0.0_DP
    
    l = l - 1
  end do 
  
  call cpu_time(tEnd)
  print *, "BC applied", tEnd - tStart
  
  k = 1
  do i = 1, nN - nDBC
    do j = 1, nN - nDBC
      if (abs(A(i,j)) >= 10D-8) k = k+1
    end do
  end do

  allocate(AA(k), JA(k), IA(nN-nDBC+1))

  
  k = 1
  do i = 1, nN - nDBC
    IA(i) = k
    do j = 1, nN - nDBC
      if (abs(A(i,j)) >= 10D-8) then
        AA(k) = A(i,j)
        JA(k) = j
        k = k+1
      end if 
    end do
  end do
  IA(nN -nDBC +1) = k
  call cpu_time(tEnd)
  print *, "Everything ready", tEnd - tStart
  sol = 0.0_DP
  sol(1:nN-nDBC) = GMRES(AA, JA, IA, b, x0(1:nN-nDBC), 10d-10, nN-nDBC, nN)

  l = nN
  do k = 1, nDBC 
    call swapNodes(nodes(IDBC(k)), nodes(l))
    help = sol(IDBC(k))
    sol(IDBC) = sol(l)
    sol(l) = help
    l = l - 1
  end do 

  do i = 1, size(x, 1)
    do j = 1, size(x, 2)
      z(i,j) = sol(idx(i,j))
    end do
  end do

  print *, "Everything done", tEnd - tStart
  call gp%title("Temperature in a channel")
  call gp%xlabel("x Position")
  call gp%ylabel("y Position")

  call gp%surf(x, y, z, lspec='t "default color spec"' ) ! color palette: gnuplot default
  
  deallocate(nodes, elements, A, b, x0, sol, AA, JA, IA, x, y, z)

end program FEM


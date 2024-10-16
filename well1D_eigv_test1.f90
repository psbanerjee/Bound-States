program schrodinger_fd
  use omp_lib   ! OpenMP library
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: N = 500  ! Number of grid points
  real(dp), parameter :: x_min = -2.0_dp  ! Minimum x value
  real(dp), parameter :: x_max = 2.0_dp   ! Maximum x value
  real(dp), parameter :: dx = (x_max - x_min) / (N - 1)  ! Grid spacing
  real(dp), allocatable :: x(:), V(:)  ! x values and potential
  real(dp), allocatable :: H(:,:), eigenvectors(:,:)  ! Hamiltonian and eigenvectors
  real(dp), allocatable :: eigenvalues(:)  ! Array to store eigenvalues
  integer :: i, j
  real(dp) :: norm_factor

  ! Allocate arrays
  allocate(x(N), V(N), H(N, N), eigenvalues(N), eigenvectors(N, N))

  ! Set the number of threads
  call omp_set_num_threads(5)   ! Use 10 cores

  ! Set up the grid points
  !$omp parallel do private(i)
  do i = 1, N
     x(i) = x_min + (i - 1) * dx
     !V(i) = 0.5_dp * x(i)**2  ! Potential V(x) = 0.5 * x^2
     if (abs(x(i)) < 1.0d0) then
     V(i)=0.0d0
     else
     V(i)=1000.0d0
     endif
  end do
  !$omp end parallel do

  ! Construct the Hamiltonian matrix using finite differences
  !$omp parallel do private(i)
  do i = 1, N
     ! Diagonal terms (kinetic + potential energy)
     H(i, i) = 1.0_dp / dx**2 + V(i)
     ! Off-diagonal terms (kinetic energy)
     if (i > 1) then
        H(i, i - 1) = -0.5_dp / dx**2
        H(i - 1, i) = -0.5_dp / dx**2
     end if
  end do
  !$omp end parallel do


  ! Diagonalize the Hamiltonian matrix to find eigenvalues and eigenvectors
  call diagonalize(H, eigenvalues, N)

  ! Print the first few eigenvalues (energy levels)
!  print *, 'First few energy levels:'
  print *, 'Square root of ratio of eigenvalues with eigenvalue of ground state:'
  do i = 1, 10
 !    print *, 'E(', i, ') = ', eigenvalues(i)
      print *, 'sqrt(E(', i, ')/E(1)) = ', sqrt(eigenvalues(i)/eigenvalues(1))
  end do

  ! Normalize and print the first few wavefunctions
!  print *, 'Wavefunctions for the first few states:'
!  do i = 1, 3
     ! Normalize the wavefunction
!     norm_factor = sqrt(sum(eigenvectors(:, i)**2) * dx)
!     eigenvectors(:, i) = eigenvectors(:, i) / norm_factor

     ! Print the normalized wavefunction
!     print *, 'Wavefunction for state ', i
!     do j = 1, N
!        print *, x(j), eigenvectors(j, i)
!     end do
!  end do

contains

  ! Subroutine to diagonalize a real symmetric matrix
  subroutine diagonalize(H, eigenvalues, N)
    implicit none
    integer, intent(in) :: N
    real(dp), intent(inout) :: H(N, N)
    real(dp), intent(out) :: eigenvalues(N)
!    real(dp), intent(out) :: eigenvectors(N, N)
    integer :: info, lwork
    real(dp), allocatable :: work(:)

    ! Query for the optimal workspace size
    lwork = -1
    allocate(work(1))
    call dsyev('V', 'U', N, H, N, eigenvalues, work, lwork, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))

    ! Diagonalize the matrix
    call dsyev('V', 'U', N, H, N, eigenvalues, work, lwork, info)
    deallocate(work)

    ! Copy the eigenvectors
    !eigenvectors = H
  end subroutine diagonalize

end program schrodinger_fd


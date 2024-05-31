program array
    implicit none

    ! 1D array of integers
    integer, dimension(4) :: array1D

    ! 2D Array
    integer, dimension(4,4) :: array2D

    ! Alternate Declarations
    integer :: alt1D(4)
    integer :: alt1Dcustom(-2:2)
    integer :: alt1Dzero(0:3)
    integer :: alt2D(4,4)

    print *, "Array 1D: ", array1D
    print *, "Array 2D: ", array2D
    print *, "Alt 1D: ", alt1D
    print *, "Alt 1D custom indices: ", alt1Dcustom
    print *, "Alt 1D custom indices: ", alt1Dzero
    print *, "Alt 2D: ", alt2D

end program array
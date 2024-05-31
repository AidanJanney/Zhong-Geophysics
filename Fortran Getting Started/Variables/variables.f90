program variables
    implicit none ! all variables will be explicitly declared, ALWAYS use

    integer :: amount
    real :: pi
    complex :: frequency
    character :: initial
    logical :: test
    ! variable names are case insensitive

    amount = 10
    pi = 3.141592653589
    frequency = (1.0, 2.0) ! (real, imaginary)
    initial = 'A' ! ' or " works
    test = .true. ! weird syntax for booleans

    ! can assign with declaraction, e.g. integer :: amount = 1, but this is bad practice
    ! this lets variables keep values between procedure calls (huh?), which can be confusing

    ! OUTPUT
    print *, "Amount: ", amount
    print *, "Pi: ", pi
    print *, "Frequency: ", frequency
    print *, "Initial: ", initial
    print *, "Test: ", test

end program variables 
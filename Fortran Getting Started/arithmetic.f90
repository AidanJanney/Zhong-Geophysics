program cylinderArithmetic
    ! Calculate the Volume and Surface Are of a cylinder when given the radius and height by the user
    implicit none

    real :: pi
    real :: radius
    real :: height
    real :: volume
    real :: surfaceArea

    ! Constant(s)
    pi = 3.141592653589793238462643383279502
    ! Get User Input
    print *, "Please enter radius of the cylinder:"
    read(*,*) radius

    print *, "Please enter the height of the cylinder:"
    read(*,*) height

    ! Calculate things
    volume = pi*(radius**2)*height
    surfaceArea = 2*pi*(radius**2) + 2*pi*radius*height

    ! Output
    print *, "Volume: ", volume
    print *, "Surface Area: ", surfaceArea
end program cylinderArithmetic
program float
    ! Use kind parameters, in this example the iso_fortran_env intrinsic module 
    ! provides 32-bit and 64-bit floating-point types
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none

    real(sp) :: EXfloat32
    real(dp) :: EXfloat64

    EXfloat32 = 1.0_sp
    EXfloat64 = 2.0_dp
    ! NEED TO USE: explicit suffix, i.e. _32 and _64
    print *, "Testing"
    print *, "Single Precision (32 bits): ", EXfloat32
    print *, "Double Precision (64 bits): ", EXfloat64
end program float

! Alternative compatiable with C
! program float
!     use, intrinsic :: iso_c_binding, only: sp=>c_float, dp=>c_double
!     implicit none
!   
!     real(sp) :: float32
!     real(dp) :: float64
  
! end program float
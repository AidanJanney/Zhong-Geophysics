module module_test
    implicit none
    integer :: num = 2

    contains
        subroutine curr_val()
            print *, "Num = ", num
        end subroutine
end module module_test

program module_example
    implicit none
    real :: test_num

    block
        use module_test !, only: n ! would only retrieve the n variable, leaving out the subroutines
        real :: local_num ! local scope var
        local_num = 4.0
        test_num = local_num**num
        print *, "Local Num: ", local_num
        print *, "Test Num: ", test_num
    end block
    print *, "Local Num no longer accessible, but Test Num: ", test_num
end program module_example
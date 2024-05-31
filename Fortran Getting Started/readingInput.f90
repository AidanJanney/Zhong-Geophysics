program read_input
    implicit none
    integer :: grad_year

    print *, "What year will you graduate? "
    read(*,*) grad_year ! standard input or stdin
    
    print *, "You will graduate in ", grad_year

end program read_input
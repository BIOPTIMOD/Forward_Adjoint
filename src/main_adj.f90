program  bio_optical_adjoint
use bioptimod_memory,  only: nw, wavelength, Rrs0p_to_Eu0m, Ed0m, Es0m, Eu0m, sunz, a_init, b_init, bb_init
use adj_3stream, only: compute_3stream_adjoint
!local
integer           :: i
CHARACTER(len=32) :: arg
double precision  :: Q
double precision  :: Ed0mOASIM(nw), Es0mOASIM(nw)
double precision  :: Rrs0p_sat(nw), Eu0m_sat(nw)

! Init

Q = 4.0D0

if (nw .EQ. 1) then

    if (iargc() .NE. 5) then 
            write(*,*) 'STOP :: Arguments number must be 5 '
            STOP
    endif

    CALL getarg(1, arg)
    READ(arg,fmt=*) wavelength(1)
    write(*,*) 'Wavelenght ', wavelength(1)

    CALL getarg(2, arg)
    READ(arg,fmt=*) Rrs0p_sat(1)
    write(*,*) 'Rrs',  Rrs0p_sat(1)

    CALL getarg(3, arg)
    READ(arg,fmt=*) Ed0mOASIM(1)
    write(*,*) 'Ed0m', Ed0mOASIM(1)

    CALL getarg(4, arg)
    READ(arg,fmt=*) Es0mOASIM(1)
    write(*,*) 'Es0m', Es0mOASIM(1)
    write(*,*) "________________"

    CALL getarg(5, arg)
    READ(arg,fmt=*) sunz
    write(*,*) 'sunz', sunz
    write(*,*) "________________"

else

    open (unit=15, file="surfdata.txt", status='old',    &
          access='sequential', form='formatted', action='read' )

    do i=1, nw

       read(15,*) wavelength(i), Rrs0p_sat(i), Ed0mOASIM(i), Es0mOASIM(i)

       write(*,*) 'Wavelenght ', wavelength(i)
       write(*,*) 'Rrs',  Rrs0p_sat(i)
       write(*,*) 'Ed0m', Ed0mOASIM(i)
       write(*,*) 'Es0m', Es0mOASIM(i)
       write(*,*) "________________"

    end do

    close(unit=15)

endif

open (unit=16, file="init_params.txt", status='old',    &
          access='sequential', form='formatted', action='read' )

   read(16,*) a_init, b_init, bb_init

   write(*,*) 'init a ', a_init
   write(*,*) 'init b',  b_init
   write(*,*) 'init bb', bb_init

close(unit=16)

Ed0m(:) = Ed0mOASIM(:)
Es0m(:) = Es0mOASIM(:)

call Rrs0p_to_Eu0m(Rrs0p_sat, Ed0mOASIM, Es0mOASIM, Q, Eu0m_sat)

write(*,*) "Eu0m_sat", Eu0m_sat

Eu0m(:) = Eu0m_sat(:)

! compute
call compute_3stream_adjoint()

!finalize
end program bio_optical_adjoint


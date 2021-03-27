program  bio_optical_adjoint
use bioptimod_memory,  only: nw, Rrs0p_to_Eu0m, Ed0m, Es0m, Eu0m
use adj_3stream, only: compute_3stream_adjoint
!local
integer :: i
double precision :: Q
double precision :: Ed0mOASIM(nw), Es0mOASIM(nw)
double precision :: Rrs0p_sat(nw), Eu0m_sat(nw)

! Init

Q = 4.0D0

do i=1, nw
   Rrs0p_sat(i) = 0.05D0
   Ed0mOASIM(i) = 0.35D0
   Es0mOASIM(i) = 0.35D0
end do

Ed0m(:) = Ed0mOASIM(:)
Es0m(:) = Es0mOASIM(:)

call Rrs0p_to_Eu0m(Rrs0p_sat, Ed0mOASIM, Es0mOASIM, Q, Eu0m_sat)

write(*,*) "Eu0m_sat", Eu0m_sat

Eu0m(:) = Eu0m_sat(:)

! compute
call compute_3stream_adjoint()

!finalize
end program bio_optical_adjoint


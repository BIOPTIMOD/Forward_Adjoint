program  bio_optical_forward
use bioptimod_memory,  only: nw, wavelength, read_command_line, parse_command_line, nlev, nphy, Ed0m, Es0m, Eu0m, &
                             allocate_bio_optical_parameters
use adj_3stream, only: solve_direct
!local
integer :: num_args, ix
character(len=12), dimension(:), allocatable :: args

integer :: i,m,col
double precision :: Q
double precision :: Ed0mOASIM(nw), Es0mOASIM(nw)
double precision :: Rrs0p_sat(nw), Eu0m_sat(nw)
double precision, allocatable :: z(:) !layer boundaries (depth levels); z(1)=0 (must be), z(n+1) = bottom
double precision, allocatable :: chl(:,:),C(:,:) 

! read arguments (nlev and nphy)
call read_command_line
call parse_command_line

! Init
allocate(z(nlev+1))
allocate(chl(nlev,nphy),C(nlev,nphy))
call allocate_bio_optical_parameters(nw, nphy)


Q = 4.0D0

! Read vertical grid mesh
open (unit=15, file="z.txt", status='old',    &
      access='sequential', form='formatted', action='read' )

do i=1, nlev+1

   read(15,*) z(i)
   write(*,*) z(i)
   write(*,*) "________________"

end do

close(unit=15)

! Read chlorophyll concentration
open (unit=16, file="chl.txt", status='old',    &
      access='sequential', form='formatted', action='read' )

do i = 1,nlev
   read(16,*) (chl(i,col),col=1,nphy)
   write(*,*) (chl(i,col),col=1,nphy)
   write(*,*) "________________"
end do

close(unit=16)


! compute
!call solve_direct(nlev+1, z, nlev, z, nlt, a, b, bb, rd, rs, ru, vd, vs, vu, EdOASIM, EsOASIM, E, E_ave)

!finalize
end program bio_optical_forward


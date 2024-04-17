program  bio_optical_forward
use bioptimod_memory,  only: nlt, wavelength, read_command_line, parse_command_line, nlev, nphy, Ed0m, Es0m, Eu0m, &
                             allocate_bio_optical_parameters, read_coefficients, read_1d_ascii, read_2d_ascii, &
                             compute_total_a_b_bb, &
                             rd, rs, ru, vs, vu
use adj_3stream, only: solve_direct
!local
double precision :: Ed_0m(nlt), Es_0m(nlt)
double precision :: Rrs0p_sat(nlt), Eu0m_sat(nlt)
double precision, allocatable :: z(:) !layer boundaries (depth levels); z(1)=0 (must be), z(n+1) = bottom
double precision, allocatable :: chl(:,:),C(:,:),nap(:),cdom(:) 
double precision, allocatable :: a(:,:), b(:,:), bb(:,:)
double precision, allocatable :: vd(:,:)
double precision, allocatable :: E(:,:,:), E_ave(:,:,:)

! read arguments (nlev and nphy)
call read_command_line
call parse_command_line

! Allocate
allocate(z(nlev+1))
allocate(chl(nlev,nphy),C(nlev,nphy))
allocate(nap(nlev),cdom(nlev))
allocate(a(nlev,nlt),b(nlev,nlt),bb(nlev,nlt))
allocate(vd(nlev,nlt))
allocate(E(3,nlev+1,nlt),E_ave(3,nlev,nlt))

call allocate_bio_optical_parameters(nlt, nphy, nlev)

!Init
call read_coefficients() ! read values for rd, rs, ru, vs, vu
call lidata_test(nlt,nchl)

vd(:,:)=1.0

! read vertical mesh
call read_1d_ascii("z.txt", nlev+1, z)

! Read chlorophyll and carbon C concentration
call read_2d_ascii("chl.txt", nlev, nphy, chl)
call read_2d_ascii("C.txt", nlev, nphy, C)

! Read nap, cdom concentration
call read_1d_ascii("nap.txt", nlev, nap)
call read_1d_ascii("cdom.txt", nlev, cdom)

! Compute total absorption (a), total scattering (b), total back scattering (bb)
call compute_total_a_b_bb(nlt, nphy, nlev, chl, C, cdom, nap, a, b, bb)
write(*,*) 'a', a(:,:)
write(*,*) 'b', b(:,:)
write(*,*) 'bb',bb(:,:)

! Retrieve boundary values at the surface
Ed_0m(:)=1.0
Es_0m(:)=1.0

! compute
call solve_direct(nlev+1, z, nlev, z, nlt, a, b, bb, rd, rs, ru, vd, vs, vu, Ed_0m, Es_0m, E, E_ave)

! print output
write(*,*) 'Ed', E(1,:,:)
write(*,*) 'Es', E(2,:,:)
write(*,*) 'Eu', E(3,:,:)
write(*,*) 'Ed_ave', E_ave(1,:,:)
write(*,*) 'Es_ave', E_ave(2,:,:)
write(*,*) 'Eu_ave', E_ave(3,:,:)

!finalize
end program bio_optical_forward


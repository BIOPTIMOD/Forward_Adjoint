MODULE bioptimod_memory

!   PUBLIC:: test_3stream_adjoint, get_derivs

!   PRIVATE
    integer         , parameter :: nw  = 1 ! Number of walength solved 
! downward planar irradiance [direct and diffuse] just below the sea surface
! provided by 
    double precision            :: Ed0m(nw), Es0m(nw), Eu0m(nw)
! Remote Sensing reflectance just above the Sea surface
    double precision            :: Rrs0p(nw)

    CONTAINS

subroutine  Rrs0p_to_Eu0m(inRrs0p, inEd0m, inEs0m, inQ, outEu0m)
 implicit none
  double precision, intent(IN)  :: inRrs0p(nw), inEd0m(nw), inEs0m(nw), inQ  
  double precision, intent(OUT) :: outEu0m(nw)  

! local
  double precision              :: Rrs0m(nw)
! this subroutine applies two corrections 
! 1) derive Rrs0m using the correction by Lee et al. 2002
  Rrs0m(:) = inRrs0p(:)
! 2) correct for the Raman Scattering at surface

  outEu0m(:)  = Rrs0m(:) * (inEd0m(:) + inEs0m(:)) * inQ

end subroutine Rrs0p_to_Eu0m

END MODULE bioptimod_memory



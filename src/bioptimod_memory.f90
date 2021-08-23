MODULE bioptimod_memory

!   PUBLIC:: test_3stream_adjoint, get_derivs

!   PRIVATE
    integer         , parameter :: nw  = 1 ! Number of walength solved 
! downward planar irradiance [direct and diffuse] just below the sea surface
! provided by 
    double precision            :: wavelength(nw)
    double precision            :: Ed0m(nw), Es0m(nw), Eu0m(nw)
! Remote Sensing reflectance just above the Sea surface
    double precision            :: Rrs0p(nw)
    double precision            :: sunz ! solar zenith angle
    double precision            :: a_init, b_init, bb_init 
    double precision,parameter  :: rd_init = 1.0D0
    double precision,parameter  :: rs_init = 1.5D0
    double precision,parameter  :: ru_init = 3.0D0
    double precision,parameter  :: vs_init = 0.83D0
    double precision,parameter  :: vu_init = 0.4D0

    CONTAINS

subroutine  Rrs0p_to_Eu0m(inRrs0p, inEd0m, inEs0m, inQ, outEu0m)
 implicit none
  double precision, intent(IN)  :: inRrs0p(nw), inEd0m(nw), inEs0m(nw), inQ  
  double precision, intent(OUT) :: outEu0m(nw)  
  double precision, parameter   :: T=0.52D0,GammaQ=1.7D0

! local
  double precision              :: Rrs0m(nw)
! this subroutine applies two corrections 
! 1) derive Rrs0m using the correction by Lee et al. 2002
  Rrs0m(:) = inRrs0p(:)/(T+GammaQ*inRrs0p(:))
! 2) correct for the Raman Scattering at surface

  outEu0m(:)  = Rrs0m(:) * (inEd0m(:) + inEs0m(:)) * inQ

end subroutine Rrs0p_to_Eu0m

END MODULE bioptimod_memory



MODULE bioptimod_memory

!   PUBLIC:: test_3stream_adjoint, get_derivs

!   PRIVATE
    integer         , parameter :: nw  = 1 ! Number of walength solved 
    character(1024) :: command_line ! parser argument variable
    integer :: nlev
    integer :: nphy
! downward planar irradiance [direct and diffuse] just below the sea surface
! provided by 
    integer                       :: wavelength(nw)
    double precision              :: Ed0m(nw), Es0m(nw), Eu0m(nw)
! Remote Sensing reflectance just above the Sea surface
    double precision              :: Rrs0p(nw)
    double precision, allocatable :: lam(:), lam1(:), lam2(:), aw(:), bw(:), bbw(:), acdom(:), apoc(:), bpoc(:), bbpoc(:) 
    double precision, allocatable :: ac(:,:), ac_ps(:,:), bc(:,:), bbc(:,:)

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


subroutine  allocate_bio_optical_parameters(nlt, nphy)
  implicit none
  integer, intent(IN)  :: nlt,nphy

  allocate(lam(nlt),lam1(nlt),lam2(nlt),aw(nlt),bw(nlt),bbw(nlt))
  allocate(acdom(nlt),apoc(nlt),bpoc(nlt),bbpoc(nlt))
  allocate(ac(nphy,nlt),ac_ps(nphy,nlt),bc(nphy,nlt),bbc(nphy,nlt))

end subroutine allocate_bio_optical_parameters

     subroutine read_command_line
       integer :: exenamelength
       integer :: io, io2

       command_line = ""
       call get_command(command = command_line,status = io)
       if (io==0) then
         call get_command_argument(0,length = exenamelength,status = io2)
         if (io2==0) then
           command_line = "&cmd "//adjustl(trim(command_line(exenamelength+1:)))//" /"
         else
           command_line = "&cmd "//adjustl(trim(command_line))//" /"
         end if
       else
         write(*,*) io,"Error getting command line."
       end if
     end subroutine

     subroutine parse_command_line
       character(256) :: msg
       namelist /cmd/  nlev, nphy
       integer :: io

       if (len_trim(command_line)>0) then
         msg = ''
         read(command_line,nml = cmd,iostat = io,iomsg = msg)
         if (io/=0) then
           error stop "Error parsing the command line or cmd.conf " // msg
         end if
       end if
     end subroutine

END MODULE bioptimod_memory



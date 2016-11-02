
  module wmpi
#if defined (MPI)
     include 'mpif.h'
#endif
  end module wmpi


  module para
     !> Some global parameters 
     !
     !> Copyright (c) 2010 QuanSheng Wu. All rights reserved.
     !
     !> add namelist for convenience  June 5th 2016 by QuanSheng Wu

     use wmpi
     implicit none

     integer,parameter :: stdout= 8

     !> define the file index to void the same index in different subroutines
     integer, public, save :: outfileindex= 10

     character(80) :: Hrfile
     namelist / TB_FILE / Hrfile

     !> control parameters
     logical :: BulkBand_calc    ! Flag for bulk energy band calculation
     logical :: BulkFS_calc      ! Flag for bulk 3D fermi surface in 3D BZ calculation
     logical :: BulkFS_plane_calc ! Flag for bulk fermi surface for a fix k plane calculation
     logical :: BulkGap_cube_calc  ! Flag for Gap_cube calculation
     logical :: BulkGap_plane_calc ! Flag for Gap_plane calculation
     logical :: SlabBand_calc  ! Flag for 2D slab energy band calculation
     logical :: WireBand_calc  ! Flag for 1D wire energy band calculation
     logical :: SlabSS_calc    ! Flag for surface state ARPES spectrum calculation
     logical :: Dos_calc       ! Flag for density of state calculation
     logical :: JDos_calc      ! Flag for joint density of state calculation
     logical :: SlabArc_calc   ! Flag for surface state fermi-arc calculation
     logical :: SlabQPI_calc   ! Flag for surface state QPI spectrum calculation
     logical :: SlabSpintexture_calc ! Flag for surface state spin-texture calculation
     logical :: WannierCenter_calc  ! Flag for Wilson loop calculation
     logical :: BerryPhase_calc   ! Flag for Berry phase calculation
     logical :: BerryCurvature_calc ! Flag for Berry curvature calculation
     logical :: EffectiveMass_calc  ! Flag for effective mass calculation
     
     namelist / Control / BulkBand_calc,BulkFS_calc,  BulkFS_Plane_calc, BulkGap_plane_calc, &
                          BulkGap_cube_calc, SlabBand_calc, WireBand_calc, &
                          SlabSS_calc, SlabArc_calc, SlabSpintexture_calc, &
                          WannierCenter_calc,BerryPhase_calc,BerryCurvature_calc, &
                          Dos_calc, JDos_calc, EffectiveMass_calc, SlabQPI_calc


     
     integer,parameter :: Dp=kind(1.0d0) ! double precision  

     
     integer :: Nslab  ! Number of slabs for 2d Slab system
     integer :: Nslab1 ! Number of slabs for 1D wire system
     integer :: Nslab2 ! Number of slabs for 1D wire system

     integer :: Np !> Number of princple layers for surface green's function
     
     integer, public, save :: ijmax=6

     integer :: Ndim !> Leading dimension of surface green's function 

     integer :: Numoccupied !> Number of occupied bands for bulk unit cell
    

     integer :: Ntotch !> Number of electrons
    
     integer :: Num_wann  ! Number of Wannier functions

     integer :: Nrpts ! Number of R points
    

   
     integer :: Nk   ! number of k points for different use
     integer :: Nk1  ! number of k points for different use
     integer :: Nk2  ! number of k points for different use
     integer :: Nk3  ! number of k points for different use

     integer, public, save :: Nr1=5
     integer, public, save :: Nr2=5
     integer, public, save :: Nr3=2

   
     integer,parameter :: kmesh(2)=(/200 , 200/)  ! kmesh used for spintexture
     integer,parameter :: knv=kmesh(1)*kmesh(2)  ! number of k points used for spintexture


     integer :: Soc  ! A parameter to control soc;  Soc=0 means no spin-orbit coupling; Soc>0 means spin-orbit coupling

    
     real(Dp) :: eta     ! used to calculate dos epsilon+i eta
     real(Dp) :: Eta_Arc ! used to calculate dos epsilon+i eta

   
     integer :: OmegaNum   ! The number of energy slices between OmegaMin and OmegaMax
     
     real(dp) :: OmegaMin, OmegaMax ! omega interval 

  
     real(Dp) :: E_arc ! Fermi energy for arc calculation

   
     real(Dp) :: Gap_threshold  ! threshold value for output the the k points data for Gap3D

     !> namelist parameters
     namelist /PARAMETERS/ Eta_Arc, OmegaNum, OmegaMin, OmegaMax, &
        E_arc, Nk1, Nk2, Nk3, NP, Gap_threshold

    
     real(Dp) :: E_fermi  ! Fermi energy, search E-fermi in OUTCAR for VASP, set to zero for Wien2k

     real(dp) :: surf_onsite  !> surface onsite energy shift
    
     real(dp) :: Bx, By, Bz !> magnetic field (in Tesla)

     !> system parameters namelist
     namelist / SYSTEM / Soc, E_fermi, Bx, By, Bz, surf_onsite, &
        Nslab, Nslab1, Nslab2, Numoccupied, Ntotch

     real(dp),parameter :: alpha= 1.20736d0*1D-6  !> e/2/h*a*a   a=1d-10m, h is the planck constant then the flux equals alpha*B*s

     !> some parameters related to atomic units
     real(dp),parameter :: bohr2atomic=0.529177211d0  ! Bohr to atomic length unit
     real(dp),parameter :: eV2Hartree= 1d0/27.211385d0 ! electron Voltage to Hartree unit

    
     real(dp),parameter :: Pi= 3.14159265359d0  ! circumference ratio pi  
     real(dp),parameter :: half= 0.5d0  ! 1/2
     real(dp),parameter :: zero= 0.0d0  ! 0
     real(dp),parameter :: one = 1.0d0  ! 1
     real(dp),parameter :: eps3= 1e-3   ! 0.001
     real(dp),parameter :: eps6= 1e-6   ! 0.000001
     real(dp),parameter :: eps9= 1e-9   ! 0.000000001

     real(Dp),parameter :: Ka(2)=(/1.0d0,0.0d0/)  
     real(Dp),parameter :: Kb(2)=(/0.0d0,1.0d0/)

     real(Dp),public, save :: Ra2(2)
     real(Dp),public, save :: Rb2(2)

     real(Dp),public, save :: Ka2(2)
     real(Dp),public, save :: Kb2(2)

     
     real(dp),public, save :: Rua(3) ! three  primitive vectors in Cartsien coordinatec, a
     real(dp),public, save :: Rub(3) ! three  primitive vectors in Cartsien coordinatec, b
     real(dp),public, save :: Ruc(3) ! three  primitive vectors in Cartsien coordinatec, c

    
     real(dp),public, save :: Rua_new(3) !> three primitive vectors in new coordinate system, see slab part
     real(dp),public, save :: Rub_new(3) !> three primitive vectors in new coordinate system, see slab part
     real(dp),public, save :: Ruc_new(3) !> three primitive vectors in new coordinate system, see slab part

    
     real(dp),public, save :: Kua(3)   ! three reciprocal primitive vectors, a
     real(dp),public, save :: Kub(3)   ! three reciprocal primitive vectors, b
     real(dp),public, save :: Kuc(3)   ! three reciprocal primitive vectors, c

     real(dp),public, save :: Urot(3, 3)  ! Rotate matrix for the new primitve cell 

     ! k list for 3D case band
     integer :: nk3lines  ! Howmany k lines for bulk band calculation
     integer :: nk3_band  ! Howmany k points for each k line
     character(4), allocatable :: k3line_name(:) ! Name of the K points
     real(dp),allocatable :: k3line_stop(:)  ! Connet points
     real(dp),allocatable :: k3line_start(:, :) ! Start point for each k line
     real(dp),allocatable :: k3line_end(:, :) ! End point for each k line
     real(dp),allocatable :: K3list_band(:, :) ! coordinate of k points for bulk band calculation in kpath mode
     real(dp),allocatable :: K3len(:)  ! put all k points in a line in order to plot the bands 
     real(dp),allocatable :: K3points(:, :) ! coordinate of k points for bulk band calculation in cube mode

     namelist / KPATH_BULK / nk3lines, k3line_name, k3line_start


     !>> top surface atoms
     integer :: NtopAtoms, NtopOrbitals  ! Select atoms on the top surface for surface state output
     integer, allocatable :: TopAtoms(:)  ! Select atoms on the top surface for surface state output
     integer, allocatable :: TopOrbitals(:)  ! Orbitals on the top surface for surface state output

     !>> bottom surface atoms
     integer :: NBottomAtoms, NBottomOrbitals ! Select atoms on the bottom  surface for surface state output
     integer, allocatable :: BottomAtoms(:) ! Select atoms on the bottom  surface for surface state output
     integer, allocatable :: BottomOrbitals(:) ! Orbitals on the bottom  surface for surface state output

     !>> effective mass

   
     real(dp), public, save :: dk_mass  ! k step for effective mass calculation
     integer , public, save :: iband_mass   ! the i'th band for effective mass calculation
     real(dp), public, save :: k_mass(3) ! the k point for effective mass calculation

     !>  klist for 2D case include all 2D system
     integer :: nk2lines  ! Number of k lines for 2D slab band calculation
     integer :: knv2 ! Number of k points for each k line
     real(dp) :: kp(2, 32) ! start k point for each k line
     real(dp) :: ke(2, 32) ! end k point for each k line
     real(dp) :: k2line_stop(32)
     character(4) :: k2line_name(32)
     real(dp),allocatable :: k2len(:)
     real(dp),allocatable :: k2_path(:, :)

     !> kpoints plane for 2D system--> arcs  
     real(dp) :: K2D_start(2)   ! start k point for 2D system calculation, like arcs
     real(dp) :: K2D_vec1(2) ! the 1st k vector for the k plane
     real(dp) :: K2D_vec2(2) ! the 2nd k vector for the k plane

     !> kpoints plane for 3D system --> gapshape
     real(dp) :: K3D_start(3) ! the start k point for the 3D k plane
     real(dp) :: K3D_vec1(3) ! the 1st k vector for the 3D k plane
     real(dp) :: K3D_vec2(3) ! the 2nd k vector for the 3D k plane
     real(dp) :: K3D_vec3(3)

     !> kpoints plane for 3D system --> gapshape3D
     real(dp) :: K3D_start_cube(3) ! the start k point for the k cube
     real(dp) :: K3D_vec1_cube(3) ! the 1st k vector for the k cube
     real(dp) :: K3D_vec2_cube(3) ! the 2nd k vector for the k cube
     real(dp) :: K3D_vec3_cube(3) ! the 3rd k vector for the k cube

     integer, allocatable     :: irvec(:,:)   ! R coordinates

     complex(dp), allocatable :: HmnR(:,:,:)   ! Hamiltonian m,n are band indexes

     integer, allocatable     :: ndegen(:)  ! degree of degeneracy of R point
 
     complex(dp),parameter    :: zi=(0.0d0, 1.0d0)    ! complex constant 0+1*i
     complex(dp),parameter    :: pi2zi=(0.0d0, 6.283185307179586d0) ! 2*pi*zi

     integer :: cpuid  ! CPU id for mpi
     integer :: num_cpu  ! Number of processors for mpi

# if defined (MPI)
     integer, parameter :: mpi_in= mpi_integer
     integer, parameter :: mpi_dp= mpi_double_precision
     integer, parameter :: mpi_dc= mpi_double_complex
     integer, parameter :: mpi_cmw= mpi_comm_world
# endif 

     real(dp), public, save :: Umatrix(3, 3)    ! a matrix change old primitive cell to new primitive cell which can define new surface, it is a 3*3 matrix

     integer :: Num_atoms  ! number of atoms in one primitive cell
     character(10) :: AngOrBohr  ! Angstrom unit to Bohr unit
     character(10) :: DirectOrCart ! Whether direct coordinates or Cartisen coordinates
     character(10), allocatable :: Atom_name(:)  ! Atom's name
     real(dp) :: CellVolume ! Cell volume
     real(dp) :: kCubeVolume
     real(dp) :: PrimitiveCellVolume
     real(dp), allocatable :: Atom_position(:, :)  ! Atom's position, only the atoms which have Wannier orbitals
     real(dp), allocatable :: Atom_position_direct(:, :)
     real(dp), allocatable :: wannier_centers_cart(:, :)
     real(dp), allocatable :: wannier_centers_direct(:, :)

     integer :: max_projs
     integer, allocatable :: nprojs(:) ! Number of projectors for each atoms
     character(10), allocatable :: proj_name(:, :) 

     !> the start index for each atoms, only consider the spinless component
     integer, allocatable :: orbitals_start(:)

     !> symmetry operator apply on function basis
     complex(dp), allocatable :: inversion(:, :)
     complex(dp), allocatable :: mirror_x(:, :)
     complex(dp), allocatable :: mirror_z(:, :)
     complex(dp), allocatable :: glide(:, :)
     
     !> symmetry operator apply on coordinate system
     real(dp), allocatable :: inv_op(:, :)
     real(dp), allocatable :: mirror_z_op(:, :)
     real(dp), allocatable :: mirror_x_op(:, :)
     real(dp), allocatable :: mirror_y_op(:, :)
     real(dp), allocatable :: glide_y_op(:, :)

 end module para

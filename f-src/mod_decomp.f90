!==============================================================================
!> \file mod_decomp.f90
!> Modules for domain decomposition
!>
!> The file holds modules which are used to perform domain decomposition within 
!> the material structure process chain
!>
!> \author Ralf Schneider
!> \date 22.01.2010

!==============================================================================
!> Module with routines for domain decomposition
module decomp

  USE global_std
  use puredat

  implicit none

  !> Type which holds informations about the global domain decomposition and 
  !> selected sub-domain
  Type tDecomp
     
     !> Number of domain in global numbering scheme
     Integer(kind=ik) :: nn

     !> Maximum number of domains in global numbering sheme
     Integer(kind=ik) :: nn_D_max
     !> Number of points in complete decomposed volume
     Integer(kind=ik), Dimension(3) :: x_VD

     !> Position of Domain on each axis in global numbering sheme
     Integer(kind=ik) :: nn_1, nn_2, nn_3
     !> Number of points in domain
     Integer(kind=ik) :: no_dat_D

     !> Points on each domain edge
     Integer(kind=ik), Dimension(3) :: x_D
     !> Physical dimension of domain
     Real(kind=rk), Dimension(3) :: x_D_phy
     !> Number of Domains on each axis
     Integer(kind=ik), Dimension(3) :: nn_D
     !> Lower left and upper right corner of domain
     Integer(kind=ik), Dimension(3) :: xa_n, xe_n
     !> Lower left and upper right corner of extended domain
     Integer(kind=ik), Dimension(3) :: xa_n_ext, xe_n_ext

     !> Grid spacing / Voxel size
     Real(kind=rk), Dimension(3) :: delta

     !> Boundary points (Domain extension)
     Integer(kind=ik), Dimension(3) :: bpoints = (/1_ik, 1_ik, 1_ik/)

  End type tDecomp
  
  !> Type to hold the description of a scalar field on a regular rectangular 
  !> grid
  Type tScalar_field

     !> Field content description
     Character, Dimension(:), Allocatable :: desc
     !> Field origin (in voxel coordinates)
     Integer(Kind=4), Dimension(3) :: orig
     !> Field dimension (in voxel coordinates) 
     !> => No. of voxels on each domain axis
     Integer(Kind=4), Dimension(3) :: vdim
     Real(kind=rk)  , Dimension(3) :: delta, shift

  End Type tScalar_field

  !============================================================================
  !== Interfaces
  !> Interface: Overloading of domain decomposition routines
  !> \author Ralf Schneider
  !> \date 22.01.2010
  Interface calc_decomp

     Module Procedure calc_decomp_general
     Module Procedure calc_decomp_domain
     !Module Procedure calc_decomp_read_general_calc_rest

  End Interface calc_decomp

Contains

  !============================================================================
  !> Function to calculate a domain decomposition for a scalar field
  !>
  !> The Function calculates all attributes of the type tDecomp from the
  !> description of a scalar field by tScalar_Field. The selection and size of
  !> the subdomain is done by the input paramaters nn and x_D
  Function calc_decomp_domain(nn, x_D, phi_desc, un) Result(dc)

    Integer(Kind=ik)              , Intent(In) :: nn
    Integer(Kind=ik), Dimension(3), Intent(In) :: x_D
    Type(tScalar_Field)           , Intent(in) :: phi_desc
    Integer                       , Intent(In) :: un

    Type(tDecomp) :: dc

    !--------------------------------------------------------------------------
    
    dc%x_D = x_D
    dc%nn  = nn

    !** Calc physical doman size **********************************************
    dc%x_D_phy = (x_D)*phi_desc%delta
    dc%delta   = phi_desc%delta

    !** Calc number of domains on each axis ***********************************
    dc%nn_D = INT((phi_desc%vdim - 2*dc%bpoints) / dc%x_D)

    dc%nn_D_max = dc%nn_D(1)*dc%nn_D(2)*dc%nn_D(3)

    !** Calc number of points on each axis in decomposed volume **************
    dc%x_VD     = dc%nn_D*x_D

    dc%nn_3 = INT( dc%nn / ( dc%nn_D(1)*dc%nn_D(2) ))
    dc%nn_2 = INT(( dc%nn - dc%nn_D(1)*dc%nn_D(2)*dc%nn_3 ) / ( dc%nn_D(1) ))
    dc%nn_1 = INT( dc%nn - dc%nn_D(1)*dc%nn_D(2)*dc%nn_3 - dc%nn_D(1)*dc%nn_2 )

    dc%xa_n(1) = dc%x_D(1) * dc%nn_1 + dc%bpoints(1) + 1
    dc%xa_n(2) = dc%x_D(2) * dc%nn_2 + dc%bpoints(2) + 1
    dc%xa_n(3) = dc%x_D(3) * dc%nn_3 + dc%bpoints(3) + 1
    
    dc%xe_n(1) = dc%x_D(1) * (dc%nn_1 + 1) + dc%bpoints(1)
    dc%xe_n(2) = dc%x_D(2) * (dc%nn_2 + 1) + dc%bpoints(2)
    dc%xe_n(3) = dc%x_D(3) * (dc%nn_3 + 1) + dc%bpoints(3)
    
    dc%xa_n_ext = dc%xa_n - dc%bpoints
    dc%xe_n_ext = dc%xe_n + dc%bpoints

    dc%no_dat_D = (dc%xe_n_ext(1)-dc%xa_n_ext(1)+1) * &
                  (dc%xe_n_ext(2)-dc%xa_n_ext(2)+1) * &
                  (dc%xe_n_ext(3)-dc%xa_n_ext(3)+1)

    Write(un,"(A,3(',',I0))"   )'nn       ',dc%nn
    Write(un,"(A,3(',',I0))"   )'x_D      ',dc%x_D
    Write(un,"(A,3(',',E18.9))")'x_D_phy  ',dc%x_D_phy
    Write(un,"(A,3(',',E18.9))")'delta    ',dc%delta
    Write(un,"(A,3(',',I0))"   )'x_VD     ',dc%x_VD
    Write(un,"(A,1(',',I0))"   )'no_dat_D ',dc%no_dat_D
    Write(un,"(A,3(',',I0))"   )'nn_D     ',dc%nn_D
    Write(un,"(A,1(',',I0))"   )'nn_D_max ',dc%nn_D_max
    Write(un,"(A,3(',',I0))"   )'nn_i     ',dc%nn_1, dc%nn_2, dc%nn_3
    Write(un,"(A,3(',',I0))"   )'xa_n     ',dc%xa_n
    Write(un,"(A,3(',',I0))"   )'xe_n     ',dc%xe_n
    Write(un,"(A,3(',',I0))"   )'xa_n_ext ',dc%xa_n_ext
    Write(un,"(A,3(',',I0))"   )'xe_n_ext ',dc%xe_n_ext

  End Function calc_decomp_domain

  !============================================================================
  !> Function to calculate a general domain decomposition for a scalar field
  !> 
  !> The function returns the global decomposition parameters of a scalar field
  !> The decomposition is determined by the physical size of the subdomains
  !> The description of phi is passed as a puredat tBranch structure
  Function calc_decomp_general(x_D_phy, phi_desc) Result(dc)

    Real(Kind=rk)   , Dimension(3), Intent(In)    :: x_D_phy
    Type(tBranch)                 , Intent(inOut) :: phi_desc

    Type(tBranch)                                 :: dc

    Real(Kind=rk)   , Dimension(:), Allocatable   :: delta
    Integer(kind=4) , Dimension(:), Allocatable   :: vdim

    Integer(kind=ik), Dimension(3)                :: x_D, nn_D
    Integer(kind=ik), Dimension(1)                :: nn_D_max, no_dat_D

    Integer(kind=ik), Dimension(3), parameter     :: bpoints=[1_ik,1_ik,1_ik]

    !--------------------------------------------------------------------------

    !** Get phi description ***************************************************
    call open_stream_files(phi_desc, "read" , "old")

    call pd_load_leaf(phi_desc%streams,phi_desc,"Grid spacings",                  delta)  
    call pd_load_leaf(phi_desc%streams,phi_desc,"Number of voxels per direction", vdim)

    call close_stream_files(phi_desc)

    !** Generate global domain decomposition **********************************

    call raise_tree("",dc)

    call raise_branch("Global domain decomposition",0,8,dc)

    call raise_leaves(no_leaves=8, &
         desc   = ["nn_D_max", "x_VD    ", "no_dat_D", "x_D     "  , &
                   "x_D_phy ", "nn_D    ", "delta   ", "bpoints "] , &
         dat_ty = [4_1 ,4_1 ,4_1 ,4_1 ,5_1 ,4_1 ,5_1 ,4_1 ] , &
         dat_no = [1_ik,3_ik,1_ik,3_ik,3_ik,3_ik,3_ik,3_ik] , &
         branch = dc)

    call open_stream_files(dc, "write", "replace")

    !** Calc voxel number on each domain edge *********************************
    x_D = Nint(x_D_phy/delta)
    call pd_store(dc%streams,dc,"x_D",x_D)

    !** Calc physical doman size **********************************************
    call pd_store(dc%streams,dc,"x_D_phy",Real(x_D,rk)*delta)
    call pd_store(dc%streams,dc,"delta", delta)

    !** Calc number of domains on each axis ***********************************
    nn_D = INT((vdim - 2*bpoints) / x_D)
    call pd_store(dc%streams,dc,"nn_D",nn_D)
    
    nn_D_max = nn_D(1)*nn_D(2)*nn_D(3) - 1
    call pd_store(dc%streams,dc,"nn_D_max",nn_D_max)

    !** Calc number of points on each axis in decomposed volume **************
    call pd_store(dc%streams,dc,"x_VD",nn_D*x_D)
    
    no_dat_D = (x_D(1) + 2*bpoints(1)) * (x_D(2) + 2*bpoints(2)) * &
               (x_D(3) + 2*bpoints(3))
    call pd_store(dc%streams,dc,"no_dat_D",no_dat_D)

    !** Boundary points on each axis in decomposed volume *********************
    call pd_store(dc%streams,dc,"bpoints",bpoints)

    call close_stream_files(dc,.TRUE.)

    call write_tree(dc)

  End Function calc_decomp_general

  !============================================================================
  !> Function to calculate a general domain decomposition for a scalar field
  !> 
  !> The function returns the global decomposition parameters of a scalar field
  !> The decomposition is determined by the physical size of the subdomains
  !> The description of phi is passed as a puredat tBranch structure
  Function calc_general_ddc_params(x_D_phy_in, phi_desc) Result(dc)

    Real(Kind=rk), Dimension(3), Intent(In) :: x_D_phy_in
    Type(tBranch), Intent(inOut) :: phi_desc

    Type(tBranch):: dc

    Real(Kind=rk) , Dimension(:), Allocatable :: delta
    Integer(kind=4) , Dimension(:), Allocatable :: vdim

    Integer(kind=ik), Dimension(3) :: x_D, nn_D
    Real(kind=rk)   , Dimension(3) :: x_D_phy
    Integer(kind=ik), Dimension(1) :: nn_D_max, no_dat_D

    Integer(kind=ik), Dimension(3), parameter :: bpoints=[1_ik,1_ik,1_ik]

    !--------------------------------------------------------------------------

    !** Get phi description ***************************************************
    call open_stream_files(phi_desc, "read" , "old")

    call pd_load_leaf(phi_desc%streams,phi_desc,"Grid spacings", delta)  
    call pd_load_leaf(phi_desc%streams,phi_desc,"Number of voxels per direction", vdim)

    call close_stream_files(phi_desc)

    !** Generate global domain decomposition **********************************
    call raise_tree("",dc)

    call raise_branch("Global domain decomposition",0,0,dc)

!!$    call raise_leaves(no_leaves=8, &
!!$         desc   = ["nn_D_max", "x_VD    ", "no_dat_D", "x_D     "  , &
!!$                   "x_D_phy ", "nn_D    ", "delta   ", "bpoints "] , &
!!$         dat_ty = [4_1 ,4_1 ,4_1 ,4_1 ,5_1 ,4_1 ,5_1 ,4_1 ] , &
!!$         dat_no = [1_ik,3_ik,1_ik,3_ik,3_ik,3_ik,3_ik,3_ik] , &
!!$         branch = dc)

    
!!$    call open_stream_files(dc,       "write", "replace")

    !** Calc voxel number on each domain edge *********************************
    x_D = Nint(x_D_phy_in/delta)
    Call add_leaf_to_branch(dc,"x_D",3_ik,x_D)
    !Call pd_store(dc%streams,dc,"x_D",x_D)

    !** Calc physical doman size **********************************************
    !call pd_store(dc%streams,dc,"x_D_phy",Real(x_D,rk)*delta)
    x_D_phy = Real(x_D,rk)*delta
    Call add_leaf_to_branch(dc, "x_D_phy", 3_ik, x_D_phy)
    !call pd_store(dc%streams,dc,"delta", delta)
    Call add_leaf_to_branch(dc,"delta",3_ik,delta)
    
    !** Calc number of domains on each axis ***********************************
    nn_D = INT((vdim - 2*bpoints) / x_D)
    !call pd_store(dc%streams,dc,"nn_D",nn_D)
    Call add_leaf_to_branch(dc,"nn_D",3_ik,nn_D)
    
    nn_D_max = nn_D(1)*nn_D(2)*nn_D(3) - 1
    !call pd_store(dc%streams,dc,"nn_D_max",nn_D_max)
    Call add_leaf_to_branch(dc,"nn_D_max",1_ik,nn_D_max)
    
    !** Calc number of points on each axis in decomposed volume **************
    !call pd_store(dc%streams,dc,"x_VD",nn_D*x_D)
    Call add_leaf_to_branch(dc,"x_VD",3_ik,nn_D*x_D)
    
    no_dat_D = (x_D(1) + 2*bpoints(1)) * (x_D(2) + 2*bpoints(2)) * &
               (x_D(3) + 2*bpoints(3))
    !call pd_store(dc%streams,dc,"no_dat_D",no_dat_D)
    Call add_leaf_to_branch(dc,"no_dat_D",1_ik,no_dat_D)

    !** Boundary points on each axis in decomposed volume *********************
    !call pd_store(dc%streams,dc,"bpoints",bpoints)
    Call add_leaf_to_branch(dc,"bpoints",3_ik,bpoints)
    
    !call close_stream_files(dc,.TRUE.)

    !call write_tree(dc)

  End Function calc_general_ddc_params

!!$  !============================================================================
!!$  !> Function to calculate the domain decomposition parameters for a subdomain
!!$  !> 
!!$  !> The function reads the  global decomposition parameters of a scalar field
!!$  !> from an input file and calculates the parameters for the spcific subdomain
!!$  !> selected by the input parameter nn
!!$  Function calc_decomp_domain_from_general(nn, g_ddc) Result(ddc)
!!$
!!$    Integer(Kind=ik)              , Intent(In) :: nn
!!$    Type(tBranch)                 , Intent(In) :: g_ddc
!!$
!!$    Type(tDecomp)                              :: ddc
!!$
!!$    Type(tLeaf), Pointer                       :: 
!!$
!!$    !--------------------------------------------------------------------------
!!$    Read(un_if,*)tmp_char,dc%nn
!!$
!!$    Read(un_if,*)tmp_char,dc%nn_1, dc%nn_2, dc%nn_3
!!$    Read(un_if,*)tmp_char,dc%xa_n
!!$    Read(un_if,*)tmp_char,dc%xe_n
!!$    Read(un_if,*)tmp_char,dc%xa_n_ext
!!$    Read(un_if,*)tmp_char,dc%xe_n_ext
!!$
!!$    dc%nn  = nn
!!$
!!$    dc%nn_3 = INT( dc%nn / ( dc%nn_D(1)*dc%nn_D(2) ))
!!$    dc%nn_2 = INT(( dc%nn - dc%nn_D(1)*dc%nn_D(2)*dc%nn_3 ) / ( dc%nn_D(1) ))
!!$    dc%nn_1 = INT( dc%nn - dc%nn_D(1)*dc%nn_D(2)*dc%nn_3 - dc%nn_D(1)*dc%nn_2 )
!!$
!!$    dc%xa_n(1) = dc%x_D(1) * dc%nn_1 + dc%bpoints(1) + 1
!!$    dc%xa_n(2) = dc%x_D(2) * dc%nn_2 + dc%bpoints(2) + 1
!!$    dc%xa_n(3) = dc%x_D(3) * dc%nn_3 + dc%bpoints(3) + 1
!!$    
!!$    dc%xe_n(1) = dc%x_D(1) * (dc%nn_1 + 1) + dc%bpoints(1)
!!$    dc%xe_n(2) = dc%x_D(2) * (dc%nn_2 + 1) + dc%bpoints(2)
!!$    dc%xe_n(3) = dc%x_D(3) * (dc%nn_3 + 1) + dc%bpoints(3)
!!$    
!!$    dc%xa_n_ext = dc%xa_n - dc%bpoints
!!$    dc%xe_n_ext = dc%xe_n + dc%bpoints
!!$
!!$    Write(un_lf,"(A,1(',',I0))"   )'nn       ',dc%nn
!!$    Write(un_lf,"(A,3(',',I0))"   )'x_D      ',dc%x_D
!!$    Write(un_lf,"(A,3(',',E18.9))")'x_D_phy  ',dc%x_D_phy
!!$    Write(un_lf,"(A,3(',',E18.9))")'delta    ',dc%delta
!!$    Write(un_lf,"(A,3(',',I0))"   )'x_VD     ',dc%x_VD
!!$    Write(un_lf,"(A,1(',',I0))"   )'no_dat_D ',dc%no_dat_D
!!$    Write(un_lf,"(A,3(',',I0))"   )'nn_D     ',dc%nn_D
!!$    Write(un_lf,"(A,1(',',I0))"   )'nn_D_max ',dc%nn_D_max
!!$    Write(un_lf,"(A,3(',',I0))"   )'nn_i     ',dc%nn_1, dc%nn_2, dc%nn_3
!!$    Write(un_lf,"(A,3(',',I0))"   )'xa_n     ',dc%xa_n
!!$    Write(un_lf,"(A,3(',',I0))"   )'xe_n     ',dc%xe_n
!!$    Write(un_lf,"(A,3(',',I0))"   )'xa_n_ext ',dc%xa_n_ext
!!$    Write(un_lf,"(A,3(',',I0))"   )'xe_n_ext ',dc%xe_n_ext
!!$
!!$  End Function calc_decomp_read_general_calc_rest

  !============================================================================
  !> Subroutine for reading a selected subdomain from a scalar field
  Subroutine domain_I2_from_pd_file_I2(unit, ddc, field, phi)

    Integer,             Intent(In) :: unit
    Type(tDecomp)      , Intent(In) :: ddc
    Type(tScalar_Field), Intent(In) :: field

    Integer(Kind=2), Dimension(*), intent(out) :: phi

    integer(Kind=ik)                :: ii, jj, nn, ma

    nn = 1
    
    Do jj = ddc%xa_n_ext(3), ddc%xe_n_ext(3)

       ma = jj * field%vdim(1) * field%vdim(2) + &
            ddc%xa_n_ext(2) * field%vdim(1)    + &
            ddc%xa_n_ext(1)

       write(*,*)ma

       Do ii = 0, ddc%x_D(2) + 2*ddc%bpoints(2) - 1

          Read(unit, pos=((ma+field%vdim(1)*ii)-1)*2+1) &
               phi(nn : nn + ddc%x_D(1) + 2 * ddc%bpoints(1) - 1)

          nn = nn + ddc%x_D(1) + 2 * ddc%bpoints(1) 

       End Do
    End Do

  End Subroutine domain_I2_from_pd_file_I2

  !============================================================================
  !> Subroutine for reading a selected subdomain from a scalar field
  Subroutine domain_I4_from_pd_file_I4(unit, ddc, field, phi, phi_branch)

    Integer,             Intent(In)               :: unit
    Type(tDecomp)      , Intent(In)               :: ddc
    Type(tScalar_Field), Intent(In)               :: field
    type(tBranch), pointer , Intent(In),optional  :: phi_branch

    Integer(Kind=4), Dimension(*), intent(out)    :: phi

    integer(Kind=ik)                :: ii, jj, nn, ma, lb = 0,kk

    character(len=256)              :: fname

    if(present(phi_branch)) then
       lb = phi_branch%leaves(6)%lbound
       !Write(*,*)"In domain_I4_from_pd_file_I4 Argument phi_branch = ",lb
    Else
       lb = 0
       !Write(*,*)"In domain_I4_from_pd_file_I4 Argument phi_branch not present"
    End if

    nn = 1

    inquire(unit=unit,name=fname)
    close(unit)
    !** intel
    !open(unit=unit,file=trim(fname),access="direct",action="read",status="old",recl=1)    
    !** gfortran
    open(unit=unit,file=trim(fname),access="direct",action="read",status="old",recl=4)
!!$
!!$    write(*,*)field%vdim
!!$    write(*,*)ddc%xa_n_ext
!!$    write(*,*)ddc%xe_n_ext
!!$    write(*,*)ddc%x_D
!!$    write(*,*)ddc%bpoints

    Do jj = ddc%xa_n_ext(3), ddc%xe_n_ext(3)

       ma = (jj-1) * field%vdim(1) * field%vdim(2) + &
            (ddc%xa_n_ext(2)-1) * field%vdim(1)    + &
            ddc%xa_n_ext(1) + lb - 1     

       !write(*,*)jj,ma,ddc%x_D(2) + 2*ddc%bpoints(2) - 1

       Do ii = 0, ddc%x_D(2) + 2*ddc%bpoints(2) - 1
     
          !write(*,*)ii,((ma-1 + (field%vdim(1))*ii)*4)+1

!!$          Read(unit, pos=((ma-1 + (field%vdim(1))*ii)*4)+1) &
!!$               phi(nn : nn + ddc%x_D(1) + 2 * ddc%bpoints(1) - 1)
          
          Do kk = 0, ddc%x_D(1) + 2 * ddc%bpoints(1) - 1
             Read(unit, rec=((ma-1 + (field%vdim(1))*ii))+1+kk) &
                  phi(nn+kk)
          end Do

          nn = nn + ddc%x_D(1) + 2 * ddc%bpoints(1) 

       End Do
    End Do

  End Subroutine domain_I4_from_pd_file_I4

  !============================================================================
  !> The Function ocnverts a puredat branch to tScalar_field
  Function branch_to_type(phi_branch, tree) Result(field)

    !** Prototype of branch "Scalar data structure"
    !** +-- Branch : Scalar data structure
    !**      +-- Leaf : Field content description
    !**      |    +-- Type of data         : 6
    !**      |
    !**      +-- Leaf : Number of voxels per direction
    !**      |    +-- Type of data         : 3
    !**      |
    !**      +-- Leaf : Grid spacings
    !**      |    +-- Type of data         : 5
    !**      |
    !**      +-- Leaf : Origin
    !**      |    +-- Type of data         : 3
    !**      |
    !**      +-- Leaf : Origin shift against global coordinate system
    !**      |    +-- Type of data         : 5
    !**      |
    !**      +-- Leaf : Scalar data
    !**           +-- Type of data         : 2

    Type(tBranch), Intent(InOut)   :: tree
    Type(tBranch), Intent(In)      :: phi_branch

    Type(tScalar_field)            :: field

    Call open_stream_files(tree, 'read', 'old')

    Allocate(field%desc(phi_branch%leaves(1)%dat_no))
    Read(tree%streams%units(6), pos=phi_branch%leaves(1)%lbound)field%desc

    Read(tree%streams%units(3), pos=(phi_branch%leaves(2)%lbound-1)*4+1)field%vdim
    Read(tree%streams%units(5), pos=(phi_branch%leaves(3)%lbound-1)*8+1)field%delta
    Read(tree%streams%units(3), pos=(phi_branch%leaves(4)%lbound-1)*4+1)field%orig
    Read(tree%streams%units(5), pos=(phi_branch%leaves(5)%lbound-1)*8+1)field%shift

  End Function branch_to_type

end module decomp

!==============================================================================
!> \file mod_vtkio.f90
!> Module for vtk output in Fortran
!>
!> The file holds a module with routines for vtk output
!>
!> \author Ralf Schneider
!> \date 11.06.2012

!==============================================================================
!> VTK output module
!> 
!> Module for VTK output of different data structures. 
!> The routines perform binary output of datasets in BIG ENDIAN see 
!>  http://vtk.1045678.n5.nabble.com/VTK-0012819-Wrong-byte-order-for-long-
!>         64bit-in-vtkDataReader-vtkDataWriter-td5092339.html
!> for details
Module vtkio

  Implicit None

  Private give_new_unit
  Private file_err

  Integer, Parameter, Private :: ik=8, rk=8, mcl=256

  Integer, Dimension(20), Parameter :: topo_FE_to_VTK = &
       [1,3,5,7, 13,15,17,19, 2,4,6,8, 14,16,18,20, 9,10,11,12]
  
contains

  !============================================================================
  !> Subroutine which writes a vtk unstructured grid
  !> 
  !> Subroutine which writes a vtk unstructured grid with mixed cell list     
  !> The Cell type is determined by the cell_type field
  !> Currently the following vtk cell types are supported:
  !> 3  => 2 Node VTK_LINE
  !> 5  => 3 Node VTK_TRIANGLE
  !> 9  => 4 Node VTK_QUAD
  !> 10 => 4 Node VTK_TETRA
  !> 12 => 8 Node VTK_HEXAHEDRON
  !> 13 => 6 Node VTK_WEDGE
  subroutine write_vtk_unstructured_grid_mixed_element_lists(filename, nodes, el_types, topo, &
       ext_nn, ext_en)

    Real(Kind=rk)    , Dimension(:,:), intent(in) :: nodes
    Integer(Kind=ik) , Dimension(:)  , intent(in) :: el_types
    Integer(Kind=ik) , Dimension(:)  , intent(in) :: topo

    Integer(Kind=ik) , Dimension(:)  , intent(in), Optional :: ext_nn
    Integer(Kind=ik) , Dimension(:)  , intent(in), Optional :: ext_en

    Character(Len=*)                 , Intent(in) :: filename

    Integer(Kind=ik)                  :: no_nodes, no_elems
    character(len=mcl)                :: tmp_line
    integer(kind=4)                   :: un_out
    integer(kind=ik)                  :: ii, ii_topo

    Integer(Kind=ik) , Dimension(:)  , Allocatable :: cref_nodes

    un_out = give_new_unit()

    no_nodes = size(ext_nn)
    no_elems = size(el_types)

    call Open_as_big_endian_stream(unit=un_out, file=trim(filename), &
         action='write', status='replace')

    call write_vtk_head(un_out)

    write(un_out) 'DATASET UNSTRUCTURED_GRID',achar(10)

    call write_vtk_nodelist(nodes, un_out)

    tmp_line=""
    write(tmp_line,"(A,1X,I0,1X,I0,A)")'CELLS',no_elems,size(topo)+size(ext_en),achar(10)
    write(un_out)trim(tmp_line)

    If (present(ext_nn)) then

       Allocate(cref_nodes(minval(ext_nn):maxval(ext_nn)))
       
       cref_nodes = -1
       
       Do ii = 0, no_nodes-1
          cref_nodes(ext_nn(ii+1)) = ii
       End Do

       ii_topo = 1
       Do ii = 1, no_elems

         Select Case (el_types(ii))

             !** Tets *******************************
          Case( 3)
             Write(un_out) 2_4 , Int(cref_nodes(topo(ii_topo:ii_topo+1)),4)
             ii_topo = ii_topo + 2
          Case( 5)
             Write(un_out) 3_4 , Int(cref_nodes(topo(ii_topo:ii_topo+2)),4)
             ii_topo = ii_topo + 3
          Case( 9)
             Write(un_out) 4_4 , Int(cref_nodes(topo(ii_topo:ii_topo+3)),4)
             ii_topo = ii_topo + 4
          Case(10)
             Write(un_out) 4_4 , Int(cref_nodes(topo(ii_topo:ii_topo+3)),4)
             ii_topo = ii_topo + 4
          Case(12)
             Write(un_out) 8_4 , Int(cref_nodes(topo(ii_topo:ii_topo+7)),4)
             ii_topo = ii_topo + 8
          Case(13)
             Write(un_out) 6_4 , Int(cref_nodes(topo(ii_topo:ii_topo+5)),4)
             ii_topo = ii_topo + 6
          case default
             Write(*,*)"Found element type ",el_types(ii)," not supported by "
             Write(*,*)"write_vtk_unstructured_grid_mixed_element_lists"
             Write(*,*)
             stop

          End Select

       End Do

    Else

       ii_topo = 1
       Do ii = 1, no_elems

          Select Case (el_types(ii))

             !** Tets *******************************
          Case( 3)
             Write(un_out) 2_4 , Int(topo(ii_topo:ii_topo+1)-1,4)
             ii_topo = ii_topo + 2
          Case( 5)
             Write(un_out) 3_4 , Int(topo(ii_topo:ii_topo+2)-1,4)
             ii_topo = ii_topo + 3
          Case( 9)
             Write(un_out) 4_4 , Int(topo(ii_topo:ii_topo+3)-1,4)
             ii_topo = ii_topo + 4
          Case(10)
             Write(un_out) 4_4 , Int(topo(ii_topo:ii_topo+3)-1,4)
             ii_topo = ii_topo + 4
          Case(12)
             Write(un_out) 8_4 , Int(topo(ii_topo:ii_topo+7)-1,4)
             ii_topo = ii_topo + 8
          Case(13)
             Write(un_out) 6_4 , Int(topo(ii_topo:ii_topo+5)-1,4)
             ii_topo = ii_topo + 6
          case default
             Write(*,*)"Found element type ",el_types(ii)," not supported by "
             Write(*,*)"write_vtk_unstructured_grid_mixed_element_lists"
             Write(*,*)
          End Select

       End Do

    End If

    tmp_line=""
    write(tmp_line,"(A,1X,I0,A)")'CELL_TYPES',no_elems,achar(10,1)
    write(un_out)trim(tmp_line)

    Write(un_out)Int(el_types,4)
    
    call close_big_endian_stream(un_out)

  end subroutine write_vtk_unstructured_grid_mixed_element_lists

  !============================================================================
  !> Subroutine which writes a vtk unstructured grid
  !>
  !> Subroutine which writes a vtk unstructured grid. No mixed cell lists are 
  !> supported. The Cell type is determined by the length of the first 
  !> dimension of elems ==> size(elems(:,1)). Currently vtk cell types
  !> 12 => 8 Node Hexaedra and 
  !>  5 => 3 Nodes Triangle are supported.
  subroutine write_vtk_unstructured_grid(filename, nodes, elems)

    Real(Kind=rk)    , Dimension(:,:), intent(in) :: nodes
    Integer(Kind=ik) , Dimension(:,:), intent(in) :: elems

    Character(Len=*)                 , Intent(in) :: filename

    Integer(Kind=ik)                  :: no_nodes, no_elems
    character(len=mcl)                :: tmp_line
    integer(kind=4)                   :: un_out, nnpe, cell_type
    integer(kind=ik)                  :: ii
    integer(kind=4) , Dimension(:), allocatable :: cell_type_array

    un_out = give_new_unit()

    if (size(elems(:,1)) == 8) then
       nnpe  = 8
       cell_type = 12
    Else if (size(elems(:,1)) == 3) then
       nnpe  = 3
       cell_type = 5
    Else if (size(elems(:,1)) == 20) then
       nnpe  = 20
       cell_type = 25
    Else
       Write(*,*)"Output of elements with ",size(elems(:,1))," nodes is not supported"
       write(*,*)"Program halted"
       Stop
    End if

    no_nodes = size(nodes(1,:))
    no_elems = size(elems(1,:))

    Allocate(cell_type_array(no_elems))
    cell_type_array = cell_type

    call Open_as_big_endian_stream(unit=un_out, file=trim(filename), &
         action='write', status='replace')

    call write_vtk_head(un_out)

    write(un_out) 'DATASET UNSTRUCTURED_GRID',achar(10)

    call write_vtk_nodelist(nodes, un_out)

    tmp_line=""
    write(tmp_line,"(A,1X,I0,1X,I0,A)")'CELLS',no_elems,no_elems*(nnpe+1),achar(10)
    write(un_out)trim(tmp_line)

    if (cell_type == 25) then
       Do ii = 1, no_elems
          Write(un_out)nnpe,Int(elems(:,ii)-1,4)
       End Do
    Else
       Do ii = 1, no_elems
          Write(un_out)nnpe,Int(elems(:,ii)-1,4)
       End Do
    End if

    tmp_line=""
    write(tmp_line,"(A,1X,I0,A)")'CELL_TYPES',no_elems,achar(10,1)
    write(un_out)trim(tmp_line)

    Write(un_out)cell_type_array

    deallocate(cell_type_array)

    call close_big_endian_stream(un_out)

  end subroutine write_vtk_unstructured_grid

  !============================================================================
  !> Subroutine which writes a vtk structured_points dataset
  subroutine write_vtk_polydata(filename, grids)

    Real(kind=rk)   , Dimension(:,:)          , intent(in) :: grids
    character(len=*)                          , intent(in) :: filename

    integer(kind=4)                  :: un_out

    integer(kind=ik)                 :: no_points

    no_points = size(grids(1,:))

    un_out = give_new_unit()
    
    call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
                                   action='write', status='replace')

    call write_vtk_head(un_out)

    write(un_out) 'DATASET POLYDATA',achar(10)

    call  write_vtk_nodelist(grids, un_out)
    
    call close_big_endian_stream(un_out)

  end subroutine write_vtk_polydata

  !============================================================================
  !> Subroutine which writes a vtk structured_points grid                    
  subroutine write_vtk_structured_points(filename, &
       extend, spacing, origin)

    character(len=*)                  , intent(in) :: filename
    Real(kind=rk)   , Dimension(3)    , intent(in) :: spacing,origin
    integer(kind=ik), Dimension(3)    , intent(in) :: extend

    integer(kind=4)                  :: un_out

    un_out = give_new_unit()
    call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
         action='write', status='replace')

    call write_vtk_head(un_out)

    call write_vtk_structured_points_head(un_out,  extend, &
                                          spacing, origin )    

    call close_big_endian_stream(un_out)

  end subroutine write_vtk_structured_points

  !============================================================================
  !> Subroutine which writes a head for a vtk binary file
  subroutine write_vtk_head(un_out)

    integer(kind=4)                           , intent(in) :: un_out

    write(un_out) '# vtk DataFile Version 3.0',achar(10)
    write(un_out) 'vtk output',achar(10)
    write(un_out) 'BINARY',achar(10)

  End subroutine write_vtk_head

  !============================================================================
  !> Subroutine which writes a head for vtk binary data output
  subroutine write_data_head(filename, size_matrix, location)

    character(len=*)          , intent(in) :: filename
    integer(kind=ik)          , intent(in) :: size_matrix
    character(len=*), optional, intent(in) :: location
        
    integer(kind=4)                        :: un_out
    character(len=mcl)                     :: tmp_pointdata
    character(len=10)                      :: loc_location
    
    if (present(location)) then
       loc_location = trim(location)
    Else
       loc_location = "POINT_DATA"
    End if

    if( (loc_location /= "POINT_DATA") .AND. (trim(loc_location) /= "CELL_DATA")) then
       write(*,'(A)')"In write_vtk_data_int4_scalar_1D"
       write(*,'(A)')"Given actual argument "//trim(loc_location)
       write(*,'(A)')"for dummy parameter location is not supported"
       write(*,'(A)')"As location parameter only POINT_DATA or CELL_DATA are accepted"
       write(*,'(A)')"Output data header is skipped."

    Else
       
       un_out = give_new_unit()
       call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
            action='write', status='old', &
            position='append')
       
       write(tmp_pointdata,'(A,1X,I0,A)') trim(loc_location), size_matrix,achar(10)
       write(un_out)trim(tmp_pointdata)
       
       call close_big_endian_stream(un_out)
    End if

  End subroutine write_data_head

  !============================================================================
  !> Subroutine which writes a head for a vtk structured_points dataset
  subroutine write_vtk_structured_points_head(un_out,  extend, &
       spacing, origin   ) 

    integer(kind=4)                   , intent(in) :: un_out
    Real(kind=rk)   , Dimension(3)    , intent(in) :: spacing,origin
    integer(kind=ik), Dimension(3)    , intent(in) :: extend

    character(len=mcl)               :: tmp_line, tmp_origin,tmp_real

    write(un_out) 'DATASET STRUCTURED_POINTS',achar(10)

    !# ---------
    !# Dimension
    tmp_real=""
    write(tmp_real,'(A,I0,1X,I0,1X,I0,A)')'DIMENSIONS ',&
         extend(1),extend(2),extend(3),achar(10)
    write(un_out) trim(tmp_real)

    !# ---------
    !# spacing
    tmp_line=""
    Write(tmp_line,'(A,F0.9,1X,F0.9,1X,F0.9,A)')'SPACING ', &
         spacing(1),spacing(2),spacing(3),achar(10)
    write(un_out)trim(tmp_line)

    !# ---------
    !# origin
    tmp_origin=""
    write(tmp_origin,'(A,F0.9,1X,F0.9,1X,F0.9,A)') 'ORIGIN ',&
         origin(1),origin(2),origin(3),achar(10)
    write(un_out)trim(tmp_origin)

  End subroutine write_vtk_structured_points_head

  !============================================================================
  !> Subroutine which writes a vtk unstructured grid nodelist
  subroutine write_vtk_nodelist(nodes, un_out)

    Real(Kind=rk)    , Dimension(:,:), intent(in) :: nodes
    Integer(Kind=4)                  , intent(in) :: un_out

    Integer(Kind=ik)                  :: no_nodes
    character(len=mcl)                :: tmp_line
   
    no_nodes = size(nodes(1,:))

    tmp_line=""
    write(tmp_line,"(A,1X,I0,1X,A,A)")'POINTS',no_nodes,'double',achar(10)
    write(un_out)trim(tmp_line)

    Write(un_out)nodes


  end subroutine write_vtk_nodelist

  !============================================================================
  !> Subroutine which writes vtk scalar data from a 1D int4 matrix
  subroutine write_vtk_data_int4_scalar_1D (matrix, filename, desc, head, location)

    Integer(kind=4) , Dimension(:)            , intent(in) :: matrix
    character(len=*)                          , intent(in) :: filename
    character(len=*), optional                , intent(in) :: desc
    Logical         , optional                , intent(in) :: head
    character(len=*), optional                , intent(in) :: location

    character(len=mcl)               :: tmp_pointdata
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_head
    character(len=10)                :: loc_location

    if (present(head)) then
       loc_head = head
    Else
       loc_head = .FALSE.
    End if

    if (present(desc)) then
       loc_desc = desc(1:min(12,len_trim(desc)))
    Else
       loc_desc = "scalars_i4"
    End if

    if (present(location)) then
       loc_location = trim(location)
    Else
       loc_location = "POINT_DATA"
    End if

    if( (loc_location /= "POINT_DATA") .AND. (trim(loc_location) /= "CELL_DATA")) then
       write(*,'(A)')"In write_vtk_data_int4_scalar_1D"
       write(*,'(A)')"Given actual argument "//trim(loc_location)
       write(*,'(A)')"for dummy parameter location is not supported"
       write(*,'(A)')"As location parameter only POINT_DATA or CELL_DATA are accepted"
       write(*,'(A)')"Output of point data is skipped."
       GOTO 1000
    End if

    un_out = give_new_unit()
    call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
                                   action='write', status='old', &
                                   position='append')
    If (loc_head) then
       !# -------
       !# Header 
       write(tmp_pointdata,'(A,1X,I0,A)') trim(loc_location), size(matrix),achar(10)
       write(un_out)trim(tmp_pointdata)

    End If

    !# Scalars
    write(un_out) 'SCALARS '//trim(loc_desc)//' int',achar(10)
    !# -------
    !# Lockup_table
    write(un_out) 'LOOKUP_TABLE default',achar(10)

    write(un_out) matrix
    write(un_out)achar(10)

1000 continue

    call close_big_endian_stream(un_out)

  end subroutine write_vtk_data_int4_scalar_1D

  !============================================================================
  !> Subroutine which writes vtk scalar data from a 1D int4 matrix
  subroutine write_vtk_data_int8_scalar_1D (matrix, filename, desc, head,location)

    Integer(kind=8) , Dimension(:)            , intent(in) :: matrix
    character(len=*)                          , intent(in) :: filename
    character(len=*), optional                , intent(in) :: desc
    Logical         , optional                , intent(in) :: head
    character(len=*), optional                , intent(in) :: location

    character(len=mcl)               :: tmp_pointdata
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_head
    character(len=10)                :: loc_location

    if (present(head)) then
       loc_head = head
    Else
       loc_head = .FALSE.
    End if

    if (present(desc)) then
       loc_desc = desc(1:min(12,len_trim(desc)))
    Else
       loc_desc = "scalars_i8"
    End if

    if (present(location)) then
       loc_location = trim(location)
    Else
       loc_location = "POINT_DATA"
    End if

    if( (loc_location /= "POINT_DATA") .AND. (trim(loc_location) /= "CELL_DATA")) then
       write(*,'(A)')"In write_vtk_data_int8_scalar_1D"
       write(*,'(A)')"Given actual argument "//trim(loc_location)
       write(*,'(A)')"for dummy parameter location is not supported"
       write(*,'(A)')"As location parameter only POINT_DATA or CELL_DATA are accepted"
       write(*,'(A)')"Output of point data is skipped."
       GOTO 1000
    End if

    un_out = give_new_unit()
    call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
                                   action='write', status='old', &
                                   position='append')
    If (loc_head) then
       !# -------
       !# Header 
       write(tmp_pointdata,'(A,1X,I0,A)') trim(loc_location), size(matrix),achar(10)
       write(un_out)trim(tmp_pointdata)

    End If

    !# Scalars
    write(un_out) 'SCALARS '//trim(loc_desc)//' long',achar(10)
    !# -------
    !# Lockup_table
    write(un_out) 'LOOKUP_TABLE default',achar(10)

    write(un_out) matrix
    write(un_out)achar(10)

1000 continue

    call close_big_endian_stream(un_out)

  end subroutine write_vtk_data_int8_scalar_1D

  !============================================================================
  !> Subroutine which writes vtk scalar data from a 1D int4 matrix
  subroutine write_vtk_data_real8_scalar_1D (matrix, filename, desc, head,location)

    Real(kind=8)    , Dimension(:)            , intent(in) :: matrix
    character(len=*)                          , intent(in) :: filename
    character(len=*), optional                , intent(in) :: desc
    Logical         , optional                , intent(in) :: head
    character(len=*), optional                , intent(in) :: location

    character(len=mcl)               :: tmp_pointdata
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical(kind=8)                  :: loc_head
    character(len=10)                :: loc_location

    if (present(head)) then
       loc_head = head
    Else
       loc_head = .FALSE.
    End if

    if (present(desc)) then
       loc_desc = desc(1:min(12,len_trim(desc)))
    Else
       loc_desc = "scalars_r8"
    End if

    if (present(location)) then
       loc_location = trim(location)
    Else
       loc_location = "POINT_DATA"
    End if

    if( (loc_location /= "POINT_DATA") .AND. (trim(loc_location) /= "CELL_DATA")) then
       write(*,'(A)')"In write_vtk_data_real8_scalar_1D"
       write(*,'(A)')"Given actual argument "//trim(loc_location)
       write(*,'(A)')"for dummy parameter location is not supported"
       write(*,'(A)')"As location parameter only POINT_DATA or CELL_DATA are accepted"
       write(*,'(A)')"Output of point data is skipped."
       GOTO 1000
    End if

    un_out = give_new_unit()
    call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
         action='write', status='old', &
         position='append')
    If (loc_head) then
       !# -------
       !# Header 
       write(tmp_pointdata,'(A,1X,I0,A)') trim(loc_location), size(matrix),achar(10)
       write(un_out)trim(tmp_pointdata)

    End If

    !# Scalars
    write(un_out) 'SCALARS '//trim(loc_desc)//' double',achar(10)
    !# -------
    !# Lockup_table
    write(un_out) 'LOOKUP_TABLE default',achar(10)

    write(un_out) matrix
    !    write(un_out)achar(10)

1000 continue

    call close_big_endian_stream(un_out)

  end subroutine write_vtk_data_real8_scalar_1D
  
  !============================================================================
  !> Subroutine which writes vtk scalar data from a 1D int4 matrix
  subroutine write_vtk_data_real8_vector_1D (matrix, filename, desc, head,location)

    Real(kind=8)    , Dimension(:,:)          , intent(in) :: matrix
    character(len=*)                          , intent(in) :: filename
    character(len=*), optional                , intent(in) :: desc
    Logical         , optional                , intent(in) :: head
    character(len=*), optional                , intent(in) :: location

    character(len=mcl)               :: tmp_pointdata
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_head
    character(len=10)                :: loc_location

    if (present(head)) then
       loc_head = head
    Else
       loc_head = .FALSE.
    End if

    if (present(desc)) then
       loc_desc = desc(1:min(12,len_trim(desc)))
    Else
       loc_desc = "vectors_r8"
    End if

    if (present(location)) then
       loc_location = trim(location)
    Else
       loc_location = "POINT_DATA"
    End if

    if( (loc_location /= "POINT_DATA") .AND. (trim(loc_location) /= "CELL_DATA")) then
       write(*,'(A)')"In write_vtk_data_real8_vector_1D"
       write(*,'(A)')"Given actual argument "//trim(loc_location)
       write(*,'(A)')"for dummy parameter location is not supported"
       write(*,'(A)')"As location parameter only POINT_DATA or CELL_DATA are accepted"
       write(*,'(A)')"Output of point data is skipped."
       GOTO 1000
    End if

    un_out = give_new_unit()
    call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
                                   action='write', status='old', &
                                   position='append')
    If (loc_head) then
       !# -------
       !# Header 
       write(tmp_pointdata,'(A,1X,I0,A)') trim(loc_location), size(matrix(1,:)),achar(10)
       write(un_out)trim(tmp_pointdata)

    End If

    !# Vectors ---
    write(un_out) 'VECTORS '//trim(loc_desc)//' double',achar(10)
    write(un_out) matrix
    !   write(un_out)achar(10)
    
1000 continue

    call close_big_endian_stream(un_out)

  end subroutine write_vtk_data_real8_vector_1D

  !============================================================================
  !> Subroutine which writes vtk scalar data from a 1D int4 matrix
  subroutine write_vtk_data_real8_tensor_1D (matrix, filename, desc, head,location)

    Real(kind=8)    , Dimension(:,:,:)        , intent(in) :: matrix
    character(len=*)                          , intent(in) :: filename
    character(len=*), optional                , intent(in) :: desc
    Logical         , optional                , intent(in) :: head
    character(len=*), optional                , intent(in) :: location

    character(len=mcl)               :: tmp_pointdata
    character(len=12)                :: loc_desc
    integer(kind=4)                  :: un_out
    logical                          :: loc_head
    character(len=10)                :: loc_location

    if (present(head)) then
       loc_head = head
    Else
       loc_head = .FALSE.
    End if

    if (present(desc)) then
       loc_desc = desc(1:min(12,len_trim(desc)))
    Else
       loc_desc = "tensors_r8"
    End if

    if (present(location)) then
       loc_location = trim(location)
    Else
       loc_location = "POINT_DATA"
    End if

    if( (loc_location /= "POINT_DATA") .AND. (trim(loc_location) /= "CELL_DATA")) then
       write(*,'(A)')"In write_vtk_data_real8_tensor_1D"
       write(*,'(A)')"Given actual argument "//trim(loc_location)
       write(*,'(A)')"for dummy parameter location is not supported"
       write(*,'(A)')"As location parameter only POINT_DATA or CELL_DATA are accepted"
       write(*,'(A)')"Output of point data is skipped."
       GOTO 1000
    End if

    un_out = give_new_unit()
    call open_as_big_endian_stream(unit=un_out, file=trim(filename), &
                                   action='write', status='old', &
                                   position='append')
    If (loc_head) then
       !# -------
       !# Header 
       write(tmp_pointdata,'(A,1X,I0,A)') trim(loc_location), size(matrix),achar(10)
       write(un_out)trim(tmp_pointdata)

    End If

    !# Scalars
    write(un_out) 'TENSORS '//trim(loc_desc)//' double',achar(10)
    write(un_out) matrix
 !   write(un_out)achar(10)

1000 continue

    call close_big_endian_stream(un_out)

  end subroutine write_vtk_data_real8_tensor_1D

  !############################################################################
  !> \name vtkio private auxilliary routines
  !> @{
  !> These routines exist in most codes. So they are declared private to aviod
  !> interference with other implementations

  !============================================================================
  !> Function which returns new free unit
  function give_new_unit() result(new_unit)

    Integer(kind=4) :: new_unit

    Integer(kind=4) :: ii

    Logical :: unit_is_open

    Do ii = 3000, huge(new_unit)-1

       inquire(unit=ii, opened=unit_is_open)

       if( .not.unit_is_open ) then
          new_unit = ii
          Exit
       end if

    End Do

    if ( unit_is_open ) then

       WRITE(*,*)
       WRITE(*,*)'Something bad and unexpected happened during search ',&
            'for free unit'
       WRITE(*,*)'Could not find a new unit between 3000 and huge(Int(kind=4))'
       WRITE(*,*)' '
       WRITE(*,*)'PROGRAM STOPPED'
       STOP
    END IF

  End function give_new_unit

  !============================================================================
  !> Subroutine for I/O error handling while operating on files
  SUBROUTINE file_err(in_file,io_stat)

    INTEGER             :: io_stat
    CHARACTER (LEN=*)   :: in_file

    IF (io_stat /= 0) Then
       WRITE(*,*)
       WRITE(*,"(80('='))")
       WRITE(*,"('EE ',A,T77,' EE')")   'Operation on file :'       
       WRITE(*,"('EE ',A          )")   in_file
       WRITE(*,"('EE ',A,T77,' EE')")   'faild !!'
       WRITE(*,"('EE ',A,I0,T77,' EE')")'With I/O Status ',io_stat
       WRITE(*,"('EE PROGRAM STOPPED ..... ',T77,' EE',/,'<',77('='),'>')")
       STOP
    End IF

  END SUBROUTINE file_err

  !============================================================================
  !> Subroutine for opening files with big endian encoding
  !> 
  Subroutine open_as_big_endian_stream(unit,file,action,status,position)

    integer(kind=4) ,intent(in)           :: unit
    character(len=*),intent(in)           :: file,action,status
    character(len=*),intent(in), optional :: position

    !integer(kind=4)             :: ier
    character(len=mcl)          :: loc_pos

    If (present(position)) then
       loc_pos = position
    Else
       loc_pos='rewind'
    End If

    !**************************************************************************
    !** GFortran, Intel implementation
    Open(unit=unit, file=trim(file), action=trim(action), &
         status=trim(status), &
         access="stream", convert="big_endian", position=trim(loc_pos))

!!$    !**************************************************************************
!!$    !** CRAY-Fortran implementation
!!$    call asnunit (unit,"-N swap_endian",ier)
!!$    IF (ier /= 0) Then
!!$       WRITE(*,*)
!!$       WRITE(*,"(80('='))")
!!$       WRITE(*,"('EE ',A,T77,' EE')")  'Asign operation on file :'       
!!$       WRITE(*,"('EE ',A          )")  file
!!$       WRITE(*,"('EE ',A,I0,T77,' EE')")   'Connected to unit :',unit
!!$       WRITE(*,"('EE ',A,T77,' EE')")   'faild !!'
!!$       WRITE(*,"('EE ',A,I0,T77,' EE')")'With error flag ',ier
!!$       WRITE(*,"('EE PROGRAM STOPPED ..... ',T77,' EE',/,'<',77('='),'>')")
!!$       STOP
!!$    End IF
!!$
!!$    Open(unit=unit, file=trim(file), action=trim(action), &
!!$         status=trim(status), access="stream")

  End Subroutine open_as_big_endian_stream

  !============================================================================
  !> Subroutine for closing files with big endian encoding
  !> 
  Subroutine close_big_endian_stream(unit)

    integer(kind=4) ,intent(in) :: unit
    !integer(kind=4)             :: ier

    !**************************************************************************
    !** GFortran, Intel implementation
    CLOSE(unit=unit)

!!$
!!$    !**************************************************************************
!!$    !** CRAY-Fortran implementation
!!$    call asnunit (unit,"-R",ier)
!!$    IF (ier /= 0) Then
!!$       WRITE(*,*)
!!$       WRITE(*,"(80('='))")
!!$       WRITE(*,"('EE ',A,T77,' EE')")  'Asign release on unit :',unit
!!$       WRITE(*,"('EE ',A,T77,' EE')")  'faild !!'
!!$       WRITE(*,"('EE ',A,I0,T77,' EE')")'With error flag ',ier
!!$       WRITE(*,"('EE PROGRAM STOPPED ..... ',T77,' EE',/,'<',77('='),'>')")
!!$       STOP
!!$    End IF
!!$
!!$    CLOSE(unit=unit)

  End Subroutine close_big_endian_stream
  !> @}
  !# End of memeber group "vtkio private auxilliary routines" #################

End Module vtkio

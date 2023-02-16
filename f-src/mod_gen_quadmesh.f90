Module gen_quadmesh

Use decomp
Use chain_routines
USE user_interaction

implicit none

!------------------------------------------------------------------------------
! Types
!------------------------------------------------------------------------------
Type tconreg
    Integer(kind=ik) :: color
    Logical :: tb   = .FALSE.

    Type(tconreg), Pointer    :: next => null()

End type tconreg

contains

    Subroutine gen_quadmesh_from_phi(delta, ddc, loc_ddc, llimit, elt_micro, &
        nodes, elems, HU_magnitudes, node_col , elem_col , no_nodes, &
        no_elems, typeraw, phi_ik1, phi_ik2, phi_ik4)

    !------------------------------------------------------------------------------
    ! Parameters
    !------------------------------------------------------------------------------                             
    Integer(1), Intent(In), optional, Dimension(1:,1:,1:) :: phi_ik1
    Integer(2), Intent(In), optional, Dimension(1:,1:,1:) :: phi_ik2
    Integer(4), Intent(In), optional, Dimension(1:,1:,1:) :: phi_ik4
    Character(Len=scl) :: typeraw
    Real(rk)  , Intent(in), Dimension(3) :: delta

    !------------------------------------------------------------------------------
    ! Domain decomposition
    !------------------------------------------------------------------------------
    Type(tBranch), Intent(In) :: ddc
    Type(tBranch), Intent(In) :: loc_ddc
    
    !------------------------------------------------------------------------------
    ! Iso Value limit
    !------------------------------------------------------------------------------
    Integer(Kind=ik)   , Intent(In) :: llimit

    !------------------------------------------------------------------------------
    ! Type of elements in micro mesh
    !------------------------------------------------------------------------------
    Character(len=*)   , Intent(In) :: elt_micro

    !------------------------------------------------------------------------------
    ! Fiels which are generated and passed out
    !------------------------------------------------------------------------------
    Real(Kind=rk)   , Dimension(:,:), Allocatable, Intent(Out) :: nodes
    Integer(Kind=ik), Dimension(:,:), Allocatable, Intent(Out) :: elems
    Integer(Kind=ik), DIMENSION(:)  , ALLOCATABLE, INTENT(OUT) :: HU_magnitudes
    Integer(Kind=ik), Dimension(:)  , Allocatable, Intent(Out) :: elem_col
    Integer(Kind=ik), Dimension(:)  , Allocatable, Intent(Out) :: node_col

    !------------------------------------------------------------------------------
    ! Sizes 
    !------------------------------------------------------------------------------
    Integer(Kind=ik), Intent(Out) :: no_nodes
    Integer(Kind=ik), Intent(Out) :: no_elems

    Integer(Kind=ik), Dimension(:), Allocatable :: nodes_no, node_cref

    Type(tconreg), Pointer :: cregs, start_cregs, tmp_creg

    Integer(kind=ik) :: ii, jj, kk, ll
    Integer(kind=ik) :: lb_nodes_no, ub_nodes_no
    Integer(kind=ik) :: min_col, max_col

    Integer(ik) :: min_val_phi, max_val_phi, val_phi
    Integer(Kind=ik), Dimension(3)  :: x_D_nodes
    Real(Kind=rk)   , Dimension(3)  :: x_min, x_max
    Integer         , Dimension(27) :: el_nn

    integer :: no_elem_nodes, alloc_stat, phi_stat

    Logical :: Change, next_exists

    Integer(kind=ik), Allocatable, Dimension(:) :: xa_n, xe_n
    Integer(kind=ik), Allocatable, Dimension(:) :: x_VD, x_D
    Character(len=9) :: nn_char
    Integer(ik) :: ddc_nn

    phi_stat=0_ik
    SELECT CASE(TRIM(ADJUSTL(typeraw)))
        CASE('ik1'); IF(.NOT. PRESENT(phi_ik1)) phi_stat=1_ik
        CASE('ik2'); IF(.NOT. PRESENT(phi_ik2)) phi_stat=1_ik
        CASE('ik4'); IF(.NOT. PRESENT(phi_ik4)) phi_stat=1_ik
    END SELECT

    IF (phi_stat == 1_ik)  THEN
        CALL print_err_stop(std_out, "Phi missing. Check your implementation!", 1)
    END IF 


    !--------------------------------------------------------------------------
    call pd_get(loc_ddc,"nn",ddc_nn)
    write(nn_char,'(I0)')ddc_nn
    
    !------------------------------------------------------------------------------
    ! Generate quadmesh
    !------------------------------------------------------------------------------
    Select Case (timer_level)
        Case (3)
            call start_timer("  +-- Generating quadmesh "//trim(nn_char))
        Case (2)
            call start_timer("  +-- Generating quadmesh "//trim(nn_char))
        Case default
            continue
    End Select

    If (out_amount /= "PRODUCTION" ) then
       Write(un_lf, FMT_MSG_SEP)
       Write(un_lf, FMT_MSG ) 'Generating Quadmesh'
       Write(un_lf,*)
    End If
    
    call pd_get(loc_ddc,"xa_n",xa_n)
    call pd_get(loc_ddc,"xe_n",xe_n)
    call pd_get(ddc,"x_VD",x_VD)
    call pd_get(ddc,"x_D",x_D)
    
    SELECT CASE(TRIM(ADJUSTL(typeraw)))
        CASE('ik1')
            min_val_phi = INT(minval(phi_ik1), ik)
            max_val_phi = INT(maxval(phi_ik1), ik)            
        CASE('ik2')
            min_val_phi = INT(minval(phi_ik2), ik)
            max_val_phi = INT(maxval(phi_ik2), ik)
        CASE('ik4')
            min_val_phi = INT(minval(phi_ik4), ik)
            max_val_phi = INT(maxval(phi_ik4), ik)
    END SELECT

    !------------------------------------------------------------------------------
    ! Check the domain for iso value inclusion
    !------------------------------------------------------------------------------
    !If ( (min_val_phi > llimit) .Or. (max_val_phi < llimit) ) Then
    If ( (max_val_phi < llimit) ) Then

        If (out_amount /= "PRODUCTION" ) Write(un_lf,FMT_MSG_SEP)
        Write(un_lf,"('EE ',A,I0,A,T77,' EE')")'Isovalue = ',llimit,' not enclosed in field'
        Write(un_lf,"('EE ',A,I0,T77,  ' EE')")'Minimum value in PHI = ',min_val_phi
        Write(un_lf,"('EE ',A,I0,T77,  ' EE')")'Maximum value in PHI = ',max_val_phi
        If (out_amount /= "PRODUCTION" ) Write(un_lf,FMT_MSG_SEP)

        no_nodes = 0
        no_elems = 0

        Select Case (timer_level)
        Case (3)
            call end_timer("  +-- Generating quadmesh "//trim(nn_char))
        Case (2)
            call end_timer("  +-- Generating quadmesh "//trim(nn_char))
        Case default
            continue
        End Select
      
        GOTO 1000

    Else

        Write(un_lf,FMT_MSG_xAI0)'Minimal value in phi = ',min_val_phi
        Write(un_lf,FMT_MSG_xAI0)'Maximal value in phi = ',max_val_phi

    End If

    !------------------------------------------------------------------------------
    ! Allocate required fields
    !------------------------------------------------------------------------------

    ! For hexaedral elements with linear trial functions
    if (elt_micro == "HEX08") then

        no_elem_nodes = 8

        x_D_nodes = x_D + 1

        Allocate(elems( 8,     x_D(1)*x_D(2)*x_D(3)), stat=alloc_stat)
        call alloc_err("elems", alloc_stat)

    ! For hexaedral elements with quadratic trial functions
    else if (elt_micro == "HEX20") then

        no_elem_nodes = 20

        x_D_nodes = x_D * 2 + 1

        Allocate(elems(20, x_D(1)*x_D(2)*x_D(3)),stat=alloc_stat)
        call alloc_err("elems",alloc_stat)

    Else

        write(std_out,FMT_ERR)"El_Type not supported in gen_quadmesh"
        write(std_out,FMT_ERR)"Program stopped "
        stop

    End if

    !------------------------------------------------------------------------------
    ! Allocate list of element gradients
    !------------------------------------------------------------------------------
    Allocate(HU_magnitudes(x_D(1)*x_D(2)*x_D(3)), stat=alloc_stat)
    call alloc_err("HU_magnitudes", alloc_stat)
    HU_magnitudes = 0_ik

    lb_nodes_no = 0
    ub_nodes_no = x_D_nodes(1) * x_D_nodes(2) * x_D_nodes(3)
    elems = 0

    Allocate(nodes_no(lb_nodes_no:ub_nodes_no),stat=alloc_stat)
    call alloc_err("nodes_no",alloc_stat)

    Write(un_lf, FMT_MSG_xAI0)'Allocated nodes numbers with index range : ',lb_nodes_no,'-',ub_nodes_no

    nodes_no  = 0
    elems     = 0
    no_elems  = 1 

    !------------------------------------------------------------------------------
    ! Generate Elements 
    !------------------------------------------------------------------------------
    ! 20 Node hexaeral elements
    If (elt_micro == "HEX20") then

        Do kk = 1,  x_D(3)
            Do jj =  1, x_D(2)
                Do ii =  1, x_D(1)
                
                    ! Transform voxel at position ii,jj,kk to Hexa element
                    SELECT CASE(TRIM(ADJUSTL(typeraw)))
                        CASE('ik1')
                            val_phi = INT(phi_ik1(ii,jj,kk), ik)
                        CASE('ik2')
                            val_phi = INT(phi_ik2(ii,jj,kk), ik)
                        CASE('ik4')
                            val_phi = INT(phi_ik4(ii,jj,kk), ik)
                    END SELECT

                    If (val_phi >= llimit) Then
                    
                    el_nn( 1) = (2*ii-2) + (2*jj-2) * (x_D_nodes(1)) + (2*kk-2) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn( 2) = (2*ii-1) + (2*jj-2) * (x_D_nodes(1)) + (2*kk-2) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn( 3) = (2*ii  ) + (2*jj-2) * (x_D_nodes(1)) + (2*kk-2) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn( 4) = (2*ii  ) + (2*jj-1) * (x_D_nodes(1)) + (2*kk-2) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn( 5) = (2*ii  ) + (2*jj  ) * (x_D_nodes(1)) + (2*kk-2) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn( 6) = (2*ii-1) + (2*jj  ) * (x_D_nodes(1)) + (2*kk-2) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn( 7) = (2*ii-2) + (2*jj  ) * (x_D_nodes(1)) + (2*kk-2) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn( 8) = (2*ii-2) + (2*jj-1) * (x_D_nodes(1)) + (2*kk-2) * (x_D_nodes(1)) * (x_D_nodes(2))
                    
                    el_nn( 9) = (2*ii-2) + (2*jj-2) * (x_D_nodes(1)) + (2*kk-1) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn(10) = (2*ii  ) + (2*jj-2) * (x_D_nodes(1)) + (2*kk-1) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn(11) = (2*ii  ) + (2*jj  ) * (x_D_nodes(1)) + (2*kk-1) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn(12) = (2*ii-2) + (2*jj  ) * (x_D_nodes(1)) + (2*kk-1) * (x_D_nodes(1)) * (x_D_nodes(2))
                    
                    el_nn(13) = (2*ii-2) + (2*jj-2) * (x_D_nodes(1)) + (2*kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn(14) = (2*ii-1) + (2*jj-2) * (x_D_nodes(1)) + (2*kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn(15) = (2*ii  ) + (2*jj-2) * (x_D_nodes(1)) + (2*kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn(16) = (2*ii  ) + (2*jj-1) * (x_D_nodes(1)) + (2*kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn(17) = (2*ii  ) + (2*jj  ) * (x_D_nodes(1)) + (2*kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn(18) = (2*ii-1) + (2*jj  ) * (x_D_nodes(1)) + (2*kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn(19) = (2*ii-2) + (2*jj  ) * (x_D_nodes(1)) + (2*kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn(20) = (2*ii-2) + (2*jj-1) * (x_D_nodes(1)) + (2*kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))

                    nodes_no(el_nn( 1:20))  = -1
                    elems(:,no_elems)       = el_nn(1:20)

                    HU_magnitudes(no_elems) = val_phi
                    no_elems = no_elems + 1
 
                    End If
                    
                End Do
            End Do
        End Do
    
    !------------------------------------------------------------------------------
    ! 8 Node hexaedral elements
    !------------------------------------------------------------------------------
    else if (elt_micro == "HEX08") then

        Do kk = 1,  x_D(3)
            Do jj =  1, x_D(2)
                Do ii =  1, x_D(1)
                
                    ! Transform voxel at position ii,jj,kk to Hexa element
                    SELECT CASE(TRIM(ADJUSTL(typeraw)))
                    CASE('ik1')
                        val_phi = INT(phi_ik1(ii,jj,kk), ik)
                    CASE('ik2')
                        val_phi = INT(phi_ik2(ii,jj,kk), ik)
                    CASE('ik4')
                        val_phi = INT(phi_ik4(ii,jj,kk), ik)
                    END SELECT

                    If (val_phi >= llimit) Then

                    el_nn(1) =  ii-1 + (jj-1) * (x_D_nodes(1)) + (kk-1) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn(2) =  ii   + (jj-1) * (x_D_nodes(1)) + (kk-1) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn(3) =  ii   + (jj  ) * (x_D_nodes(1)) + (kk-1) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn(4) =  ii-1 + (jj  ) * (x_D_nodes(1)) + (kk-1) * (x_D_nodes(1)) * (x_D_nodes(2))
                    
                    el_nn(5) =  ii-1 + (jj-1) * (x_D_nodes(1)) + (kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn(6) =  ii   + (jj-1) * (x_D_nodes(1)) + (kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn(7) =  ii   + (jj  ) * (x_D_nodes(1)) + (kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
                    el_nn(8) =  ii-1 + (jj  ) * (x_D_nodes(1)) + (kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
                    
                    nodes_no(el_nn(1:8)) = -1
                    elems(:,no_elems)    = el_nn(1:8)
                    
                    HU_magnitudes(no_elems) = val_phi
                    no_elems = no_elems + 1

                    End If
                    
                End Do
            End Do
        End Do

    End If

    no_elems  = no_elems  - 1
    no_nodes  = 0

    !------------------------------------------------------------------------------
    ! Renumber nodes
    !------------------------------------------------------------------------------
    Do ii = 1, no_elems
        Do jj = 1, no_elem_nodes
            If (nodes_no(elems(jj,ii)) == -1) then
                no_nodes = no_nodes + 1
                nodes_no(elems(jj,ii)) = no_nodes
            End If
        End Do
    End Do

    write(un_lf,FMT_MSG_xAI0)'No nodes found in domain : ',no_nodes

    Do ii = 1, no_elems
        elems(1:no_elem_nodes,ii) = nodes_no(elems( 1:no_elem_nodes,ii))       
    End Do

    !------------------------------------------------------------------------------
    ! Calculate Physical positions of nodes
    !------------------------------------------------------------------------------
    Allocate(nodes(3,no_nodes))
    
    If (elt_micro == "HEX20") then

        Do ll = lb_nodes_no, ub_nodes_no
            
            If (nodes_no(ll) /= 0) then
                
                kk =  ll / (x_D_nodes(1) * x_D_nodes(2))                        !+ xa_n(3)
                jj = (ll - kk * x_D_nodes(1) * x_D_nodes(2)) / x_D_nodes(1)     !+ xa_n(2)
                ii =  ll - kk * x_D_nodes(1) * x_D_nodes(2) - jj * x_D_nodes(1) !+ xa_n(1)

                nodes(3,nodes_no(ll)) = (Real(kk,rk)/2._rk + Real(xa_n(3)-1,rk))*delta(3) 
                nodes(2,nodes_no(ll)) = (Real(jj,rk)/2._rk + Real(xa_n(2)-1,rk))*delta(2) 
                nodes(1,nodes_no(ll)) = (Real(ii,rk)/2._rk + Real(xa_n(1)-1,rk))*delta(1)
    
            End If
        
        End Do

        Else if (elt_micro == "HEX08") then

        Do ll = lb_nodes_no, ub_nodes_no
            
            If (nodes_no(ll) /= 0) then

                kk =  ll / (x_D_nodes(1) * x_D_nodes(2))                        !+ xa_n(3)
                jj = (ll - kk * x_D_nodes(1) * x_D_nodes(2)) / x_D_nodes(1)     !+ xa_n(2)
                ii =  ll - kk * x_D_nodes(1) * x_D_nodes(2) - jj * x_D_nodes(1) !+ xa_n(1)

                nodes(3,nodes_no(ll)) = (Real(kk+ xa_n(3)-1,rk))*delta(3) 
                nodes(2,nodes_no(ll)) = (Real(jj+ xa_n(2)-1,rk))*delta(2) 
                nodes(1,nodes_no(ll)) = (Real(ii+ xa_n(1)-1,rk))*delta(1) 

            End If

        End Do

    End If

    deallocate(nodes_no)

    Select Case (timer_level)
    Case (3)
        call end_timer  ("  +-- Generating quadmesh "//trim(nn_char))
        call Start_timer("  +-- Coloring connected domains "//trim(nn_char))
    Case (2)
        call end_timer  ("  +-- Generating quadmesh "//trim(nn_char))
        call Start_timer("  +-- Coloring connected domains "//trim(nn_char))
    Case default
        continue
    End Select

    Write(un_lf,FMT_MSG_SEP)
    Write(un_lf,FMT_MSG)    'Coloring connected domains'


    !------------------------------------------------------------------------------
    ! Only colorize if more than 0 elements.
    !------------------------------------------------------------------------------
    IF (no_nodes .LE. 8_ik) GOTO 1000

    !------------------------------------------------------------------------------
    ! Color nodes
    !------------------------------------------------------------------------------
    Allocate(node_col(no_nodes))

    Do ii = 1, no_nodes
        node_col(ii) = ii
    End Do

    Change = .TRUE.

    DO WHILE (CHANGE)

        Change = .FALSE.

        Do ii = 1, no_elems

            min_col = minval(node_col(elems(:,ii)))

            max_col = maxval(node_col(elems(:,ii)))

            if (min_col /= max_col) then
                node_col(elems(:,ii)) = max_col
                change = .TRUE.
            End if

        End Do

    End DO


    !------------------------------------------------------------------------------
    ! Color elements
    !------------------------------------------------------------------------------
    Allocate(elem_col(no_elems))

    Do ii = 1, no_elems
        elem_col(ii) = node_col(elems(1,ii))
    End Do

    !------------------------------------------------------------------------------
    ! Setup colorlist
    !------------------------------------------------------------------------------
    Allocate(cregs)
    cregs%next => null()
    start_cregs => cregs

    cregs%color = node_col(1)

    Do ii = 1, no_nodes

        !------------------------------------------------------------------------------
        ! Check for different colors
        !------------------------------------------------------------------------------
        if (node_col(ii) /= cregs%color) then

            ! Check whether the different color exists in list
            tmp_creg => start_cregs
            Do While ( Associated(tmp_creg%next) )
                if ( tmp_creg%color == node_col(ii) ) then
                    exit
                Else
                    tmp_creg => tmp_creg%next
                End if
            End Do

            !------------------------------------------------------------------------------
            ! If color not catched yet
            !------------------------------------------------------------------------------
            If (tmp_creg%color /= node_col(ii)) then
                allocate(cregs%next)
                cregs => cregs%next
                cregs%color = node_col(ii)
            End If

       End if

    End Do

    !------------------------------------------------------------------------------
    ! Check which domains are connected to the cube faces
    !------------------------------------------------------------------------------
    x_min(1) = minval(nodes(1,:))
    x_max(1) = maxval(nodes(1,:))
    x_min(2) = minval(nodes(2,:))
    x_max(2) = maxval(nodes(2,:))
    x_min(3) = minval(nodes(3,:))
    x_max(3) = maxval(nodes(3,:))

    tmp_creg => start_cregs
    Do ii = 1, no_nodes
        if (node_col(ii) == tmp_creg%color) then
            if ((nodes(1,ii) -x_min(1)) < (0.1_rk* delta(1))) then
                tmp_creg%tb = .TRUE.
                Exit
            End if
            if ((nodes(1,ii) -x_max(1)) < (0.1_rk* delta(1))) then
                tmp_creg%tb = .TRUE.
                Exit
            End if
            if ((nodes(2,ii) -x_min(2)) < (0.1_rk* delta(2))) then
                tmp_creg%tb = .TRUE.
                Exit
            End if
            if ((nodes(2,ii) -x_max(2)) < (0.1_rk* delta(2))) then
                tmp_creg%tb = .TRUE.
                Exit
            End if
            if ((nodes(3,ii) -x_min(3)) < (0.1_rk* delta(3))) then
                tmp_creg%tb = .TRUE.
                Exit
            End if
            if ((nodes(3,ii) -x_max(3)) < (0.1_rk* delta(3))) then
                tmp_creg%tb = .TRUE.
                Exit
            End if
        End if
    End Do

    Do While ( Associated(tmp_creg%next) )
        tmp_creg => tmp_creg%next
        Do ii = 1, no_nodes
            if (node_col(ii) == tmp_creg%color) then
                if (abs(nodes(1,ii) -x_min(1)) < (0.1_rk* delta(1))) then
                    tmp_creg%tb = .TRUE.
                    Exit
                End if
                if (abs(nodes(1,ii) -x_max(1)) < (0.1_rk* delta(1))) then
                    tmp_creg%tb = .TRUE.
                    Exit
                End if
                if (abs(nodes(2,ii) -x_min(2)) < (0.1_rk* delta(2))) then
                    tmp_creg%tb = .TRUE.
                    Exit
                End if
                if (abs(nodes(2,ii) -x_max(2)) < (0.1_rk* delta(2))) then
                    tmp_creg%tb = .TRUE.
                    Exit
                End if
                if (abs(nodes(3,ii) -x_min(3)) < (0.1_rk* delta(3))) then
                    tmp_creg%tb = .TRUE.
                    Exit
                End if
                if (abs(nodes(3,ii) -x_max(3)) < (0.1_rk* delta(3))) then
                    tmp_creg%tb = .TRUE.
                    Exit
                End if
            End if
        End Do
    End Do

    !------------------------------------------------------------------------------
    ! Debug ** Output color list
    !------------------------------------------------------------------------------
    tmp_creg => start_cregs
    Do While ( Associated(tmp_creg%next) )
       Write(un_lf,"('MM ',I0,L1)")tmp_creg%color, tmp_creg%tb
       tmp_creg => tmp_creg%next
    End Do
    Write(un_lf,"('MM ',I0,L1)")tmp_creg%color, tmp_creg%tb
    !------------------------------------------------------------------------------
    ! Debug
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Make color lists binary for output
    !------------------------------------------------------------------------------
    tmp_creg    => start_cregs
    next_exists =  .TRUE.
    Do While (next_exists)

        if (tmp_creg%tb) then
            Do ii = 1, no_nodes
                if (node_col(ii) == tmp_creg%color) node_col(ii) = -1
            End Do
            Do ii = 1, no_elems
                if (elem_col(ii) == tmp_creg%color) elem_col(ii) = -1
            End Do
        end if

        If (Associated(tmp_creg%next)) then
            tmp_creg => tmp_creg%next
        Else
            next_exists = .FALSE.
        End If

    End Do

    !------------------------------------------------------------------------------
    ! Delete unconnected elements
    !------------------------------------------------------------------------------
    jj = 1
    kk = 1

    Do while (jj < no_elems)

        if (elem_col(jj) /= -1) then
            jj = jj + 1
        else      
            jj = jj + 1
            kk = kk + 1  
        end if

        elem_col(kk) = elem_col(jj)
        elems(:,kk)  = elems(:,jj)

        HU_magnitudes(kk)  = HU_magnitudes(jj)

    End Do

    no_elems = kk
    
    !------------------------------------------------------------------------------
    ! Delete unconnected nodes
    !------------------------------------------------------------------------------
    allocate(node_cref(no_nodes))

    jj = 1
    kk = 1
    if (node_col(jj) == -1) then
        node_col(kk)  = node_col(jj)
        nodes(:,kk)   = nodes(:,jj)
        node_cref(jj) = kk
    End if

    Do 

        if (node_col(jj) /= -1) then
            node_cref(jj) = 0
            node_col(jj)  = 0
            jj = jj + 1
            else      
            jj = jj + 1
            kk = kk + 1  
        end if

        node_col(kk)  = node_col(jj)
        nodes(:,kk)   = nodes(:,jj)
        node_cref(jj) = kk

        if (jj == no_nodes) exit

    End Do

    no_nodes = kk

    !------------------------------------------------------------------------------
    ! renumber nodes in elem list
    !------------------------------------------------------------------------------
    do ii = 1, no_elems
        elems(1:no_elem_nodes,ii) = node_cref(elems(1:no_elem_nodes,ii))
    end do

    deallocate(node_cref)

    Select Case (timer_level)
    Case (3)
    
        call end_timer("  +-- Coloring connected domains "//trim(nn_char))
    Case (2)
        call end_timer("  +-- Coloring connected domains "//trim(nn_char))
    Case default
        continue
    End Select
   
    !------------------------------------------------------------------------------
    ! Error with Iso-Value continue
    !------------------------------------------------------------------------------
    1000 Continue

    If (out_amount /= "PRODUCTION" ) then
       Write(un_lf,FMT_MSG_xAI0)"Remaining nodes   : ",no_nodes
       Write(un_lf,FMT_MSG_xAI0)"Remaining elements: ",no_elems
    End If
    
  End Subroutine gen_quadmesh_from_phi
  
!   Subroutine quadmesh_from_phi(phi, delta, ddc, llimit, elt_micro, &
!                                nodes, elems, node_col , elem_col , &
!                                no_nodes, no_elems)

!     !-- Parameters ------------------------------------------------------------

!     ! Iso Field *************************************************************
!     Integer(4)    , Intent(In), Dimension(1:,1:,1:) :: Phi
!     Real(rk)      , Intent(in), Dimension(3)        :: delta

!     ! Domain decomposition **************************************************
!     Type(tBranch)      , Intent(In)                   :: ddc
    
!     ! Iso Value limit *******************************************************
!     Integer(ik)   , Intent(In)                 :: llimit

!     ! Type of elements in micro mesh ****************************************
!     Character(len=*)   , Intent(In)                 :: elt_micro

!     ! Fiels which are generated and passed out ******************************
!     Real(rk)   , Dimension(:,:), Allocatable, Intent(Out) :: nodes
!     Integer(ik), Dimension(:,:), Allocatable, Intent(Out) :: elems
!     Integer(ik), Dimension(:)  , Allocatable, Intent(Out) :: elem_col
!     Integer(ik), Dimension(:)  , Allocatable, Intent(Out) :: node_col

!     ! Sizes *****************************************************************
!     Integer(ik)                             , Intent(Out) :: no_nodes
!     Integer(ik)                             , Intent(Out) :: no_elems

!     !--------------------------------------------------------------------------
!     Integer(ik), Dimension(:), Allocatable :: nodes_no, node_cref

!     Type(tconreg), Pointer                      :: cregs, start_cregs, tmp_creg

!     Integer(ik)                            :: ii, jj, kk, ll
!     Integer(ik)                            :: lb_nodes_no, ub_nodes_no
!     Integer(ik)                            :: min_col, max_col

!     Integer                                     :: min_val_phi, max_val_phi
!     Integer(ik), Dimension(3)              :: x_D_nodes
!     Real(rk)   , Dimension(3)              :: x_min, x_max
!     Integer         , Dimension(27)             :: el_nn

!     integer                                     :: no_elem_nodes, alloc_stat

!     Logical                                     :: Change, next_exists

!     Integer(ik), Allocatable, Dimension(:) :: xa_n, xe_n
!     Integer(ik), Allocatable, Dimension(:) :: x_VD, x_D

!     !--------------------------------------------------------------------------

!     ! Generate quadmesh *****************************************************
!     call start_timer("+-- Generating quadmesh")

!     Write(un_lf,FMT_MSG_SEP)
!     Write(un_lf,FMT_MSG)    'Generating Quadmesh'
!     Write(un_lf,*)
!     GOTO 1000
!     call pd_get(ddc,"xa_n",xa_n)
!     call pd_get(ddc,"xe_n",xe_n)
!     call pd_get(ddc,"x_VD",x_VD)
!     call pd_get(ddc,"x_D",x_D)

!     min_val_phi = minval(phi)
!     max_val_phi = maxval(phi)

!     ! Check the domain for iso value inclusion ******************************
!     !If ( (min_val_phi > llimit) .Or. (max_val_phi < llimit) ) Then
!     If ( (max_val_phi < llimit) ) Then

!        Write(un_lf,FMT_ERR_SEP)
!        Write(un_lf,FMT_ERR_xAI0) 'Isovalue = ',llimit,' not enclosed in field'
!        Write(un_lf,FMT_ERR_xAI0) 'Minimum value in PHI = ',min_val_phi
!        Write(un_lf,FMT_ERR_xAI0) 'Maximum value in PHI = ',max_val_phi
!        Write(un_lf,FMT_ERR_SEP)

!        no_nodes = 0
!        no_elems = 0

!        GOTO 1000

!     Else

!        Write(un_lf,FMT_MSG_xAI0)'Minimal value in phi = ',min_val_phi
!        Write(un_lf,FMT_MSG_xAI0)'Maximal value in phi = ',max_val_phi

!     End If

!     !************************************************************************
!     ! Allocate needed fields ************************************************

!     ! For hexaedral elements with linear trial functions ********************
!     if (elt_micro == "HEX08") then

!        no_elem_nodes = 8

!        x_D_nodes = x_D + 1

!        Allocate(elems( 8, x_D(1)*x_D(2)*x_D(3)),stat=alloc_stat)
!        call alloc_err("elems",alloc_stat)

!     ! For hexaedral elements with quadratic trial functions *****************
!     else if (elt_micro == "HEX20") then

!        no_elem_nodes = 20

!        x_D_nodes = x_D * 2 + 1

!        Allocate(elems(20, x_D(1)*x_D(2)*x_D(3)),stat=alloc_stat)
!        call alloc_err("elems",alloc_stat)

!     Else

!        write(*,'(A)')"El_Type not supported in gen_quadmesh !!!"
!        write(*,'(A)')"Program stopped                       !!!"
!        stop

!     End if

!     lb_nodes_no = 0
!     ub_nodes_no = x_D_nodes(1) * x_D_nodes(2) * x_D_nodes(3)
!     elems = 0

!     Allocate(nodes_no(lb_nodes_no:ub_nodes_no),stat=alloc_stat)
!     call alloc_err("nodes_no",alloc_stat)

!     If (out_amount /= "PRODUCTION" ) then
!        Write(un_lf, FMT_MSG_xAI0)'Allocated nodes numbers with index range : ',lb_nodes_no,'-',ub_nodes_no
!     End If
    
!     nodes_no  = 0
!     elems     = 0
!     no_elems  = 1 

!     !************************************************************************
!     ! Generate Elements 

!     ! 20 Node hexaeral elements *********************************************
!     If (elt_micro == "HEX20") then

!        Do kk = 1,  x_D(3)
!           Do jj =  1, x_D(2)
!              Do ii =  1, x_D(1)
                
!                 ! Transform voxel at position ii,jj,kk to Hexa element ***********
!                 If (Phi(ii,jj,kk) >= llimit) Then

!                    el_nn( 1) = (2*ii-2) + (2*jj-2) * (x_D_nodes(1)) + (2*kk-2) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn( 2) = (2*ii-1) + (2*jj-2) * (x_D_nodes(1)) + (2*kk-2) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn( 3) = (2*ii  ) + (2*jj-2) * (x_D_nodes(1)) + (2*kk-2) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn( 4) = (2*ii  ) + (2*jj-1) * (x_D_nodes(1)) + (2*kk-2) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn( 5) = (2*ii  ) + (2*jj  ) * (x_D_nodes(1)) + (2*kk-2) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn( 6) = (2*ii-1) + (2*jj  ) * (x_D_nodes(1)) + (2*kk-2) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn( 7) = (2*ii-2) + (2*jj  ) * (x_D_nodes(1)) + (2*kk-2) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn( 8) = (2*ii-2) + (2*jj-1) * (x_D_nodes(1)) + (2*kk-2) * (x_D_nodes(1)) * (x_D_nodes(2))
                   
!                    el_nn( 9) = (2*ii-2) + (2*jj-2) * (x_D_nodes(1)) + (2*kk-1) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn(10) = (2*ii  ) + (2*jj-2) * (x_D_nodes(1)) + (2*kk-1) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn(11) = (2*ii  ) + (2*jj  ) * (x_D_nodes(1)) + (2*kk-1) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn(12) = (2*ii-2) + (2*jj  ) * (x_D_nodes(1)) + (2*kk-1) * (x_D_nodes(1)) * (x_D_nodes(2))
                  
!                    el_nn(13) = (2*ii-2) + (2*jj-2) * (x_D_nodes(1)) + (2*kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn(14) = (2*ii-1) + (2*jj-2) * (x_D_nodes(1)) + (2*kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn(15) = (2*ii  ) + (2*jj-2) * (x_D_nodes(1)) + (2*kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn(16) = (2*ii  ) + (2*jj-1) * (x_D_nodes(1)) + (2*kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn(17) = (2*ii  ) + (2*jj  ) * (x_D_nodes(1)) + (2*kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn(18) = (2*ii-1) + (2*jj  ) * (x_D_nodes(1)) + (2*kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn(19) = (2*ii-2) + (2*jj  ) * (x_D_nodes(1)) + (2*kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn(20) = (2*ii-2) + (2*jj-1) * (x_D_nodes(1)) + (2*kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))

!                    nodes_no(el_nn( 1:20)) = -1
!                    elems(:,no_elems)      = el_nn(1:20)
!                    no_elems               = no_elems + 1

!                 End If
                
!              End Do
!           End Do
!        End Do
    
!     ! 8 Node hexaedral elements *********************************************
!     else if (elt_micro == "HEX08") then

!        Do kk = 1,  x_D(3)
!           Do jj =  1, x_D(2)
!              Do ii =  1, x_D(1)
                
!                 ! Transform voxel at position ii,jj,kk to Hexa element ***********
!                 If (Phi(ii,jj,kk) >= llimit) Then

!                    el_nn(1) =  ii-1 + (jj-1) * (x_D_nodes(1)) + (kk-1) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn(2) =  ii   + (jj-1) * (x_D_nodes(1)) + (kk-1) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn(3) =  ii   + (jj  ) * (x_D_nodes(1)) + (kk-1) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn(4) =  ii-1 + (jj  ) * (x_D_nodes(1)) + (kk-1) * (x_D_nodes(1)) * (x_D_nodes(2))
                   
!                    el_nn(5) =  ii-1 + (jj-1) * (x_D_nodes(1)) + (kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn(6) =  ii   + (jj-1) * (x_D_nodes(1)) + (kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn(7) =  ii   + (jj  ) * (x_D_nodes(1)) + (kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
!                    el_nn(8) =  ii-1 + (jj  ) * (x_D_nodes(1)) + (kk  ) * (x_D_nodes(1)) * (x_D_nodes(2))
                   
!                    nodes_no(el_nn(1:8)) = -1
!                    elems(:,no_elems)    = el_nn(1:8)
!                    no_elems             = no_elems + 1

!                 End If
                
!              End Do
!           End Do
!        End Do

!     End If

!     no_elems  = no_elems  - 1
!     no_nodes  = 0

!     ! Renumber nodes **********************************************************
!     Do ii = 1, no_elems
!        Do jj = 1, no_elem_nodes
!           If (nodes_no(elems(jj,ii)) == -1) then
!              no_nodes = no_nodes + 1
!              nodes_no(elems(jj,ii)) = no_nodes
!           End If
!        End Do
!     End Do

!     If (out_amount /= "PRODUCTION" ) then
!        write(un_lf,FMT_MSG_xAI0)'No nodes found in domain : ',no_nodes
!     End If
    
!     Do ii = 1, no_elems
!        elems(1:no_elem_nodes,ii) = nodes_no(elems( 1:no_elem_nodes,ii))       
!     End Do

!     ! Calculate Physical positions of nodes ***********************************
!     Allocate(nodes(3,no_nodes))
    
!     If (elt_micro == "HEX20") then

!        Do ll = lb_nodes_no, ub_nodes_no
          
!           If (nodes_no(ll) /= 0) then
             
!              kk =  ll / (x_D_nodes(1) * x_D_nodes(2))                        !+ xa_n(3)
!              jj = (ll - kk * x_D_nodes(1) * x_D_nodes(2)) / x_D_nodes(1)     !+ xa_n(2)
!              ii =  ll - kk * x_D_nodes(1) * x_D_nodes(2) - jj * x_D_nodes(1) !+ xa_n(1)

!              nodes(3,nodes_no(ll)) = (Real(kk,rk)/2._rk + Real(xa_n(3)-1,rk))*delta(3) 
!              nodes(2,nodes_no(ll)) = (Real(jj,rk)/2._rk + Real(xa_n(2)-1,rk))*delta(2) 
!              nodes(1,nodes_no(ll)) = (Real(ii,rk)/2._rk + Real(xa_n(1)-1,rk))*delta(1)
 
!           End If
       
!        End Do

!     Else if (elt_micro == "HEX08") then

!        Do ll = lb_nodes_no, ub_nodes_no
          
!           If (nodes_no(ll) /= 0) then

!              kk =  ll / (x_D_nodes(1) * x_D_nodes(2))                        !+ xa_n(3)
!              jj = (ll - kk * x_D_nodes(1) * x_D_nodes(2)) / x_D_nodes(1)     !+ xa_n(2)
!              ii =  ll - kk * x_D_nodes(1) * x_D_nodes(2) - jj * x_D_nodes(1) !+ xa_n(1)

!              nodes(3,nodes_no(ll)) = (Real(kk+ xa_n(3)-1,rk))*delta(3) 
!              nodes(2,nodes_no(ll)) = (Real(jj+ xa_n(2)-1,rk))*delta(2) 
!              nodes(1,nodes_no(ll)) = (Real(ii+ xa_n(1)-1,rk))*delta(1) 

!           End If

!        End Do

!     End If

    

!     deallocate(nodes_no)
!     call end_timer("+-- Generating quadmesh")

!     ! Color connected domains *************************************************
!     Call start_timer("+-- Coloring connected domains")

!     If (out_amount /= "PRODUCTION" ) then
!        Write(un_lf,FMT_MSG_SEP)
!        Write(un_lf,FMT_MSG)'Coloring connected domains'
!     End If
    
!     ! Color nodes *************************************************************
!     Allocate(node_col(no_nodes))

!     Do ii = 1, no_nodes
!        node_col(ii) = ii
!     End Do

!     Change = .TRUE.

!     DO WHILE (CHANGE)

!        Change = .FALSE.

!        Do ii = 1, no_elems

!           min_col = minval(node_col(elems(:,ii)))

!           max_col = maxval(node_col(elems(:,ii)))

!           if (min_col /= max_col) then
!              node_col(elems(:,ii)) = max_col
!              change = .TRUE.
!           End if

!        End Do

!     End DO

!     ! Color elements **********************************************************
!     Allocate(elem_col(no_elems))

!     Do ii = 1, no_elems
!        elem_col(ii) = node_col(elems(1,ii))
!     End Do

!     ! Setup colorlist *********************************************************
!     Allocate(cregs)
!     cregs%next => null()
!     start_cregs => cregs

!     cregs%color = node_col(1)

!     Do ii = 1, no_nodes

!        ! Check for different colors *******************************************
!        if (node_col(ii) /= cregs%color) then

!           ! Check whether the different color exists in list ******************
!           tmp_creg => start_cregs
!           Do While ( Associated(tmp_creg%next) )
!              if ( tmp_creg%color == node_col(ii) ) then
!                 exit
!              Else
!                 tmp_creg => tmp_creg%next
!              End if
!           End Do

!           ! If color not catched yet ******************************************
!           If (tmp_creg%color /= node_col(ii)) then
!              allocate(cregs%next)
!              cregs => cregs%next
!              cregs%color = node_col(ii)
!           End If

!        End if

!     End Do

!     ! Check which domains are connected to the cube faces *********************

!     ! Rechnen und nicht suchen !! *********************************************
!     x_min(1) = minval(nodes(1,:))
!     x_max(1) = maxval(nodes(1,:))
!     x_min(2) = minval(nodes(2,:))
!     x_max(2) = maxval(nodes(2,:))
!     x_min(3) = minval(nodes(3,:))
!     x_max(3) = maxval(nodes(3,:))

!     tmp_creg => start_cregs
!     Do ii = 1, no_nodes
!        if (node_col(ii) == tmp_creg%color) then
!           if ((nodes(1,ii) -x_min(1)) < (0.1_rk* delta(1))) then
!              tmp_creg%tb = .TRUE.
!              Exit
!           End if
!           if ((nodes(1,ii) -x_max(1)) < (0.1_rk* delta(1))) then
!              tmp_creg%tb = .TRUE.
!              Exit
!           End if
!           if ((nodes(2,ii) -x_min(2)) < (0.1_rk* delta(2))) then
!              tmp_creg%tb = .TRUE.
!              Exit
!           End if
!           if ((nodes(2,ii) -x_max(2)) < (0.1_rk* delta(2))) then
!              tmp_creg%tb = .TRUE.
!              Exit
!           End if
!           if ((nodes(3,ii) -x_min(3)) < (0.1_rk* delta(3))) then
!              tmp_creg%tb = .TRUE.
!              Exit
!           End if
!           if ((nodes(3,ii) -x_max(3)) < (0.1_rk* delta(3))) then
!              tmp_creg%tb = .TRUE.
!              Exit
!           End if
!        End if
!     End Do

!     Do While ( Associated(tmp_creg%next) )
!        tmp_creg => tmp_creg%next
!        Do ii = 1, no_nodes
!           if (node_col(ii) == tmp_creg%color) then
!              if (abs(nodes(1,ii) -x_min(1)) < (0.1_rk* delta(1))) then
!                 tmp_creg%tb = .TRUE.
!                 Exit
!              End if
!              if (abs(nodes(1,ii) -x_max(1)) < (0.1_rk* delta(1))) then
!                 tmp_creg%tb = .TRUE.
!                 Exit
!              End if
!              if (abs(nodes(2,ii) -x_min(2)) < (0.1_rk* delta(2))) then
!                 tmp_creg%tb = .TRUE.
!                 Exit
!              End if
!              if (abs(nodes(2,ii) -x_max(2)) < (0.1_rk* delta(2))) then
!                 tmp_creg%tb = .TRUE.
!                 Exit
!              End if
!              if (abs(nodes(3,ii) -x_min(3)) < (0.1_rk* delta(3))) then
!                 tmp_creg%tb = .TRUE.
!                 Exit
!              End if
!              if (abs(nodes(3,ii) -x_max(3)) < (0.1_rk* delta(3))) then
!                 tmp_creg%tb = .TRUE.
!                 Exit
!              End if
!           End if
!        End Do
!     End Do

!     ! Debug ** Output color list **********************************************
!     If (out_amount /= "PRODUCTION" ) then
!        tmp_creg => start_cregs
!        Do While ( Associated(tmp_creg%next) )
!           Write(un_lf,"('MM ',I0,L1)")tmp_creg%color, tmp_creg%tb
!           tmp_creg => tmp_creg%next
!        End Do
!        Write(un_lf,"('MM ',I0,L1)")tmp_creg%color, tmp_creg%tb
!     End If
!     ! Debug *****************************************************************

!     ! Make color lists binary for output ************************************
!     tmp_creg    => start_cregs
!     next_exists =  .TRUE.
!     Do While (next_exists)

!        if (tmp_creg%tb) then
!           Do ii = 1, no_nodes
!              if (node_col(ii) == tmp_creg%color) node_col(ii) = -1
!           End Do
!           Do ii = 1, no_elems
!              if (elem_col(ii) == tmp_creg%color) elem_col(ii) = -1
!           End Do
!        end if

!        If (Associated(tmp_creg%next)) then
!           tmp_creg => tmp_creg%next
!        Else
!           next_exists = .FALSE.
!        End If

!     End Do

!     ! Delete unconnected elements *******************************************
!     jj = 1
!     kk = 1

!     Do while (jj < no_elems)

!        if (elem_col(jj) /= -1) then
!           jj = jj + 1
!        else      
!           jj = jj + 1
!           kk = kk + 1  
!        end if

!        elem_col(kk) = elem_col(jj)
!        elems(:,kk)  = elems(:,jj)

!     End Do

!     no_elems = kk
    
!     ! Delete unconnected nodes **********************************************
!     allocate(node_cref(no_nodes))

!     jj = 1
!     kk = 1
!     if (node_col(jj) == -1) then
!        node_col(kk)  = node_col(jj)
!        nodes(:,kk)   = nodes(:,jj)
!        node_cref(jj) = kk
!     End if

!     Do 

!        if (node_col(jj) /= -1) then
!           node_cref(jj) = 0
!           node_col(jj)  = 0
!           jj = jj + 1
!         else      
!           jj = jj + 1
!           kk = kk + 1  
!        end if

!        node_col(kk)  = node_col(jj)
!        nodes(:,kk)   = nodes(:,jj)
!        node_cref(jj) = kk

!        if (jj == no_nodes) exit

!     End Do

!     no_nodes = kk

!     ! renumber nodes in elem list *******************************************
!     do ii = 1, no_elems
!        elems(1:no_elem_nodes,ii) = node_cref(elems(1:no_elem_nodes,ii))
!     end do

!     deallocate(node_cref)

!     ! Error with Iso-Value continue *****************************************
! 1000 Continue

!     Write(un_lf,FMT_MSG_xAI0)"Remaining nodes   : ",no_nodes
!     Write(un_lf,FMT_MSG_xAI0)"Remaining elements: ",no_elems

!     Call end_timer("+-- Coloring connected domains")

!   End Subroutine quadmesh_from_phi

End Module gen_quadmesh

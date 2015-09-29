module hi_SIEPPEM

    use hi_const
    use hi_mod_funcs
    use hi_target_func
    
    implicit none 
    
    !include 'hi_const.f90'    
    !include 'hi_funcs.f90'    
    integer,protected    ::  num_dim,num_node,num_nrml,num_elem
    integer,protected    ::  elem_nd_count,num_target_func
    real(8),protected    ::  hi_beta 


    real(8),allocatable,private ::  node_matrix(:,:),normal_matrix(:,:),src_local_list(:,:)
    !node_matrix(1:3,node_id),normal_matrix(1:3,nrml_id)
    integer,allocatable,private ::  elem_matrix(:,:),src_flag(:)
    !elem_matrix(1:8,elem_id)
    real(8),allocatable,private ::  full_mesh_matrix(:,:,:)
    !full_mesh_matrix(1:3,1:8,elem_id)

    
    integer,private :: model_readed_flag = 0 ! 0 for not readed
    integer,private :: external_src_ctr_flag = 0! if or not use external src ctr input

    integer,private,parameter :: NPW = 4
    real(8),private,allocatable :: value_list(:,:)
    integer,private :: n_pwr_g = -1
    real(8),allocatable,private :: cnr_glb_mtx(:,:) !corner_global_matrix
    real(8),private :: src_lcl_preset(2),src_glb_preset(3)
   
    real(8),private ::  src_glb(3),src_ctr_glb(3)
    !src coordinate in global,src center point in global
contains
    include 'test4.f90'
    include 'test3.f90'
    include 'test2.f90'
    include 'test.f90'
        
    subroutine read_model_from_WAVDUT()

        
        !use MVAR_MOD
        use mesh
        implicit none

        integer :: i,j,ie,id,tmp,tmp1

        print *,"------------Start Reading Model from WAVDUT-------------"

        num_dim = 3
        num_node = NNODE
        num_elem = NELEM
        num_nrml = NNODED
        elem_nd_count = 8 !NCN(IELEM)!! to be changed
        hi_beta = 3.
        num_target_func = 8

        allocate(node_matrix(num_dim,num_node))
        allocate(normal_matrix(num_dim,num_nrml))

        allocate(elem_matrix(elem_nd_count,num_elem))
        allocate(value_list(num_elem,num_target_func))
        allocate(cnr_glb_mtx(num_dim,elem_nd_count))

        if (num_dim == 2) ngl = 1

        forall (i = 1:num_dim,j=1:num_node)
            node_matrix(i,j) = XYZ(i,j)
        end forall

        forall (i = 1:num_dim,j=1:num_nrml)
            normal_matrix(i,j) = DXYZ(i,j)
        end forall

        forall(i = 1:num_elem,j = 1:elem_nd_count)
            elem_matrix(j,i) = NCON(i,j)!triangle
            ! there is index order problem
!              elem_mtx_nrml(i,j) = NCOND(i,j)
        end forall

        ! for first and third edge, vertical component unchanged,1st comp change
                ! for second and forth edge, horizontal component unchanged,2rd comp change
                !       4----7----3
                !       |         |
                !       8    9    6
                !       |         |  new position
                !       1----5----2

                !         7     6     5

                !         8           4 old position

                !         1     2     3
        !switch order
        print *,"Transfering elem matrix from old order to new order"
            do i = 1,num_elem
                if (elem_nd_count .eq. 8) then
!                     tmp = elem_matrix(2,i)
!                     elem_matrix(2,i) = elem_matrix(5,i)
!                     tmp1 = elem_matrix(3,i)
!                     elem_matrix(3,i) = tmp
!                     tmp = elem_matrix(4,i)
!                     elem_matrix(4,i) = elem_matrix(6,i)
!                     elem_matrix(5,i) = tmp1
!                     elem_matrix(6,i) = elem_matrix(7,i)
!                     elem_matrix(7,i) = tmp



                endif
                if (NCN(i) .eq. 8) then
                tmp = elem_matrix(2,i)
                elem_matrix(2,i) = elem_matrix(3,i)
                elem_matrix(3,i) = elem_matrix(5,i)
                elem_matrix(5,i) = tmp
                tmp = elem_matrix(4,i)
                elem_matrix(4,i) = elem_matrix(7,i)
                elem_matrix(7,i) = elem_matrix(6,i)
                elem_matrix(6,i) = tmp
            endif
            end do
        !-------------------- Data Manipulation------------------

        allocate(full_mesh_matrix(num_dim,elem_nd_count,num_elem))
            
        forall (ie = 1:num_elem,id = 1:elem_nd_count)

                full_mesh_matrix(1:num_dim,id,ie)=node_matrix(1:num_dim,elem_matrix(id,ie))
                ! reorganize nodes coordinate by element node order
        end forall

        !------------- Initialization-------------------------
        model_readed_flag = 0
        value_list = 0

        print *,"------------Stop Reading Model from WAVDUT-------------"
        !print *, node_matrix
!         print *,"------------elem matrix-------------"

        !print *, elem_matrix



    end subroutine read_model_from_WAVDUT

    subroutine swap_result(result)
        implicit none
        real(8) :: result(*),tmp,tmp1

                !       4----7----3
                !       |         |
                !       8    9    6
                !       |         |
                !       1----5----2

                !         7     6     5

                !         8           4

                !         1     2     3

        tmp = result(2)
        result(2) = result(5)
        tmp1 = result(3)
        result(3) = tmp
        tmp = result(4)
        result(4) = result(6)
        result(5) = tmp1
        result(6) = result(7)
        result(7) = tmp
    end subroutine


    subroutine read_model_from_DAT()

        implicit none
        integer :: ip,ie,tmp,i,id        

        if (model_readed_flag == 0) then
            print *,"------------Start Reading Model-------------"
            OPEN(5,FILE='SIEPPEM.DAT',STATUS='OLD')

            read (5,*) num_dim,num_node,num_elem,elem_nd_count,hi_beta,num_target_func
            ! dimesnion
            ! number of node
            ! number of element
            ! number of node per element
            ! beta is the power of r in target equation
            ! number of target func components

            allocate(node_matrix(num_dim,num_node))
            allocate(elem_matrix(elem_nd_count,num_elem))
            allocate(src_flag(num_elem))
            allocate(src_local_list(2,num_elem))
            allocate(value_list(num_elem,num_target_func))

            if (num_dim == 2) ngl = 1
    
         !    Input nodal coordinates and element connectivity
          
            DO IP = 1,num_node
                READ(5,*) tmp,(node_matrix(I,tmp),I=1,num_dim)                 ! Card set 2
            end do  


            DO IE = 1,num_elem
                READ(5,*) tmp,(elem_matrix(ID,tmp),ID=1,elem_nd_count),src_flag(tmp)    ! Card set 3
            end do
            !==========================
            !====src_flag
            ! if = 0 src not on elem
            ! if > 0 src is given in global coordinate, use node with id (src_flag)
            ! if < 0 src is given in local coordinate, use local src list given in card set 4

          
            READ(5,*) (src_glb(i),i=1,num_dim)                       ! Card set 4  
             ! read src x,y,z coord
             ! there seems a error, src_glb should be an array of global coordinate
             ! since src_glb cannot remain unchanged for different element
        
            DO IE=1,num_elem
                if (src_flag(ie) < 0) then
                    read(5,*) (src_local_list(i,ie),i=1,num_dim-1)
                end if
                ! src local position given in elements input order
            end do

            close(5)

            print *,"------------Finish Reading Model-------------"
            !------------Some initialisation of data

            allocate(full_mesh_matrix(num_dim,elem_nd_count,num_elem))
            
            forall (ie = 1:num_elem,id = 1:elem_nd_count)

                    full_mesh_matrix(1:num_dim,id,ie)=node_matrix(1:num_dim,elem_matrix(id,ie))
                    ! reorganize nodes coordinate by element node order
            end forall

            model_readed_flag = 0
            value_list = 0

        else
            print *,"--Attention! Reading process skipped,model already loaded---"
        end if

    end subroutine read_model_from_DAT

     subroutine set_npwg(number)
        
        implicit none

        integer,intent(in) :: number

        n_pwr_g = number
        print *,"n_power_g is set to ",number

    end subroutine

    subroutine set_src_preset(ksi,eta,glb,ctr_glb)
        
        implicit none

        real(8),intent(in) :: ksi,eta,glb(3),ctr_glb(3)

        src_lcl_preset(1) = ksi
        src_lcl_preset(2) = eta
        src_glb_preset = glb
        src_ctr_glb = ctr_glb
        print *,"src preseted as",src_lcl_preset
        print *,"src preseted as",src_glb_preset
        print *,"src_ctr_glb preseted as",src_ctr_glb


    end subroutine

    subroutine debug_test()
        
        implicit none

        print *,"elem matrix",elem_matrix

    end subroutine
    


   
!     subroutine RIM_ELEMS()

!         implicit none 

!         integer  :: num_edge,ie,id

! !         EXTERNAL INT_ELEM

!         IF(num_dim == 2) n_pwr_g = (elem_nd_count/3)*2               ! 0,  2
!         IF(num_dim == 3) n_pwr_g = elem_nd_count/2+(elem_nd_count/9)*2        ! 2,  4,  6 
!         ! refer to  Equ. 3-6-56 for parameter m
     
!         num_edge = 2 * (num_dim - 1 ) ! 4 -----how many edges

!         ! NDSID=2+elem_nd_count/8 !3 !removed not used=== July 25th====

!         src_ctr_glb = 0. !src center global, define the center of src for calculating r

!         do ie = 1,num_elem ! can introduce parallel here!!!!!!!!!!!!!!!!!
!             If (src_flag(ie) == 0) THEN    
!                 ! EVALUATE INTEGRAL OVER REGULAR ELEMENT
!                 !CALL ADAPTINT_ELEM(ie,src_ctr_glb,cnr_lcl_mtx,value_list(ie),GPR,GWR,GPL,GWL,INT_ELEM)
!                 print *," need evaluate integral over element,src not on element"            
!             else     
!                 CALL eval_SINGULAR_ELEM(ie,num_edge,value_list(ie,:),0)
!             end if 
!         end do

!         print *,"==========final result========="
!         print *,sum(value_list)
!         print *,"==============================="
!         print *,"end of RIM_ELEMS"
!     end subroutine

!     subroutine initialise_ELEMS(x0,y0,z0)

!         implicit none 

!         real,intent(in) :: x0,y0,z0
!         integer  :: num_edge,ie,id

!         !         EXTERNAL INT_ELEM

!         !-------------------------------------------

!         IF(num_dim == 2) n_pwr_g = (elem_nd_count/3)*2               ! 0,  2
!         IF(num_dim == 3) n_pwr_g = elem_nd_count/2+(elem_nd_count/9)*2        ! 2,  4,  6 
!         ! refer to  Equ. 3-6-56 for parameter m
     
!         num_edge = 2 * (num_dim - 1 ) ! 4 -----how many edges

!         external_src_ctr_flag = 1
        
!         if (external_src_ctr_flag .eq. 0) then
!             src_ctr_glb = 0. !src center global, define the center of src for calculating r
!         else
!             src_ctr_glb(1) = x0
!             src_ctr_glb(2) = y0
!             src_ctr_glb(3) = z0
!         end if

!         !         do ie = 1,num_elem ! can introduce parallel here!!!!!!!!!!!!!!!!!
!         !             If (src_flag(ie) == 0) THEN    
!         !                 ! EVALUATE INTEGRAL OVER REGULAR ELEMENT
!         !                 !CALL ADAPTINT_ELEM(ie,src_ctr_glb,cnr_lcl_mtx,value_list(ie),GPR,GWR,GPL,GWL,INT_ELEM)
!         !                 print *," need evaluate integral over element,src not on element"            
!         !             else     
!         !                 CALL eval_SINGULAR_ELEM(ie,num_edge,value_list(ie,:))
!         !             end if 
!         !         end do

!         !         print *,"==========final result========="
!         !         print *,sum(value_list)
!         !         print *,"==============================="
!         !         print *,"end of RIM_ELEMS"
!     end subroutine
end module
      



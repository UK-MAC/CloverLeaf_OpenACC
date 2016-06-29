!Crown Copyright 2012 AWE.
!
! This file is part of CloverLeaf.
!
! CloverLeaf is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! CloverLeaf is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! CloverLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Driver for the PdV update.
!>  @author Wayne Gaudin
!>  @details Invokes the user specified kernel for the PdV update.

MODULE loadbalance_module
  USE MPI
  USE load_balance_kernel_module
  IMPLICIT NONE

CONTAINS

  SUBROUTINE loadbalance(my_timer, my_step)


    USE clover_module
    USE data_module

    IMPLICIT NONE

    REAL(KIND=8) :: my_timer
    INTEGER      :: my_step


    INTEGER         :: neigh,tsk,tile
    REAL(KIND=8)    :: left_timer, right_timer

    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: array, neigh_array
    INTEGER, ALLOCATABLE,DIMENSION(:) :: my_int_array, neigh_int_arr
    INTEGER :: arr_size, neigh_arr_size
    INTEGER :: send(2), recv(2)
    INTEGER      :: status(MPI_STATUS_SIZE,2), rcv_status(MPI_STATUS_SIZE), err, right_task, left_task
    INTEGER      :: req_send(2), req_recv(2)

    send=0
    recv=0

    left_timer=0.0_8
    right_timer=0.0_8



    ALLOCATE(my_int_array(9), neigh_int_arr(9))


    CALL clover_swap_timers(my_timer, left_timer, right_timer)

    write(*,*) parallel%task, "Timers" ,my_timer, left_timer, right_timer, tiles_per_chunk
 ! Establish exchange pattern

    IF(MOD(my_step,2).EQ.1) THEN
      ! LB SEND  left -> right
      IF(chunk%chunk_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE) THEN
       IF(my_timer.GT.(right_timer*(1.1_8)))THEN
         send(2)=1
         !write(*,*) parallel%task, "Send right",my_timer, right_timer
       ENDIF
      END IF

      ! LB RECV  left -> right
      IF(chunk%chunk_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE) THEN
       IF(left_timer.GT.(my_timer*(1.1_8)))THEN
          recv(1)=1
         !write(*,*) parallel%task, "recv left",my_timer, left_timer
       ENDIF


      END IF

    ELSE
      ! LB SEND  right -> left
      IF(chunk%chunk_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE) THEN
       IF(my_timer.GT.(left_timer*(1.1_8)))THEN
         send(1)=1
         !write(*,*) parallel%task, "Send left",my_timer, left_timer
       ENDIF


      END IF

      ! LB RECV  right -> left
      IF(chunk%chunk_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE) THEN
       IF(right_timer.GT.(my_timer*(1.1_8)))THEN
         recv(2)=1
         !write(*,*) parallel%task, "recv right",my_timer, right_timer
       ENDIF
      END IF

    END IF

! Left to right




!IF(parallel%task.eq.0) send(2) =1
!IF(parallel%task.eq.1) recv(1) =1

IF(send(2).EQ.1)THEN

  CALL pack_tile(tiles_per_chunk,array, arr_size)
  CALL pack_tile_meta(tiles_per_chunk,my_int_array, arr_size)
  ! send ints to right
  right_task=chunk%chunk_neighbours(chunk_right) - 1
  CALL MPI_ISEND(my_int_array,9,MPI_INT,right_task,1, &
      MPI_COMM_WORLD,req_send(1),err)

  ! Send Data to right
  CALL MPI_ISEND(array,arr_size,MPI_DOUBLE_PRECISION,right_task,2, &
    MPI_COMM_WORLD,req_send(2),err)
END IF

IF(recv(1).EQ.1)THEN
  ! reorder tiles

  left_task=chunk%chunk_neighbours(chunk_left) - 1
  CALL MPI_RECV(neigh_int_arr,9,MPI_INT,left_task,1, &
      MPI_COMM_WORLD,rcv_status,err)

  CALL add_tile_start()
  ! recv size from left
  CALL unpack_tile_meta(1,neigh_int_arr, neigh_arr_size)
  CALL build_field_tile(1)
  ALLOCATE(neigh_array(neigh_arr_size))
  ! recv tile from left
  ! Send Data to right
  CALL MPI_RECV(neigh_array,neigh_arr_size,MPI_DOUBLE_PRECISION,left_task,2, &
      MPI_COMM_WORLD,rcv_status,err)
  CALL unpack_tile(1,neigh_array, neigh_arr_size)

  DEALLOCATE(neigh_array)


END IF


IF(send(2).EQ.1)THEN
  ! Wait

  CALL MPI_WAITALL(2,req_send,status,err)
  ! Kill patch
  CALL kill_field_tile(tiles_per_chunk)
  CALL remove_tile_end()
  ! send
  DEALLOCATE(array)
END IF




! Right to left



IF(send(1).EQ.1)THEN
  CALL pack_tile(1,array, arr_size)
  CALL pack_tile_meta(1,my_int_array, arr_size)
  ! send ints to right
  left_task=chunk%chunk_neighbours(chunk_left) - 1
  CALL MPI_ISEND(my_int_array,9,MPI_INT,left_task,1, &
      MPI_COMM_WORLD,req_send(1),err)
  ! Send Data to right
  CALL MPI_ISEND(array,arr_size,MPI_DOUBLE_PRECISION,left_task,2, &
      MPI_COMM_WORLD,req_send(2),err)
END IF

IF(recv(2).EQ.1)THEN
  ! reorder tiles

  right_task=chunk%chunk_neighbours(chunk_right) - 1
  CALL MPI_RECV(neigh_int_arr,9,MPI_INT,right_task,1, &
      MPI_COMM_WORLD,rcv_status,err)

  ! recv size from left
  CALL add_tile_end()
  CALL unpack_tile_meta(tiles_per_chunk,neigh_int_arr, neigh_arr_size)
  CALL build_field_tile(tiles_per_chunk)
  ALLOCATE(neigh_array(neigh_arr_size))
  ! recv tile from left
  ! Send Data to right
  CALL MPI_RECV(neigh_array,neigh_arr_size,MPI_DOUBLE_PRECISION,right_task,2, &
      MPI_COMM_WORLD,rcv_status,err)
  CALL unpack_tile(tiles_per_chunk,neigh_array, neigh_arr_size)

  DEALLOCATE(neigh_array)


END IF


IF(send(1).EQ.1)THEN
  ! Wait

  CALL MPI_WAITALL(2,req_send,status,err)
  ! Kill patch
  CALL kill_field_tile(1)
  CALL remove_tile_start()
  ! send
  DEALLOCATE(array)
END IF
!
!do tsk=0,parallel%max_task-1
!  IF(tsk.eq.parallel%task) THEN
!  do neigh=1,tiles_per_chunk
!   write(*,*) parallel%task, "Tiles" , neigh, chunk%tiles(neigh)%tp%t_left, chunk%tiles(neigh)%tp%t_right
!  END DO
!    CALL FLUSH(0)
!    CALL FLUSH(6)
!  END IF
!  CALL MPI_BARRIER(MPI_COMM_WORLD,err)
!END DO


    DEALLOCATE(my_int_array, neigh_int_arr)


!do tsk=0,parallel%max_task-1
!  IF(tsk.eq.parallel%task) THEN
!  do tile=1,tiles_per_chunk
!   CALL updatehost_tile(tile)
!   write(*,*) parallel%task, "Tiles" , tile, chunk%tiles(tile)%tp%t_left, chunk%tiles(tile)%tp%t_right, &
!                                       chunk%tiles(tile)%tp%tile_neighbours(TILE_LEFT), chunk%tiles(tile)%tp%tile_neighbours(TILE_RIGHT),&
!                                       chunk%tiles(tile)%tp%tile_neighbours(TILE_TOP), chunk%tiles(tile)%tp%tile_neighbours(TILE_BOTTOM)
!   write(*,*) SUM(chunk%tiles(tile)%tp%field%cellx), SUM(chunk%tiles(tile)%tp%field%celly), SUM(chunk%tiles(tile)%tp%field%energy0), &
!              SUM(chunk%tiles(tile)%tp%field%density0)
!   !write(*,*) chunk%tiles(tile)%tp%field%celly(:)
!chunk%tiles(tile)%tp%field%cellx(chunk%tiles(tile)%tp%t_xmin-2), &
!       chunk%tiles(tile)%tp%field%cellx(chunk%tiles(tile)%tp%t_xmax+2), chunk%tiles(tile)%tp%field%celly(chunk%tiles(tile)%tp%t_ymin-2),chunk%tiles(tile)%tp%field%celly(chunk%tiles(tile)%tp%t_ymax+2)

!  END DO
!    CALL FLUSH(0)
!    CALL FLUSH(6)
!  END IF
!  CALL MPI_BARRIER(MPI_COMM_WORLD,err)
!END DO




  END SUBROUTINE loadbalance


  SUBROUTINE pack_tile(tile, array, arr_size)
    USE clover_module
    USE data_module
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: tile
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(OUT) :: arr_size
    INTEGER :: xmin, xmax, ymin, ymax
    INTEGER :: xwidth_small, xwidth_large, yheight_small, yheight_large

    INTEGER :: index, x, y

    !call updatehost_tile(tile)

    xmin = chunk%tiles(tile)%tp%t_xmin
    xmax = chunk%tiles(tile)%tp%t_xmax
    ymin = chunk%tiles(tile)%tp%t_ymin
    ymax = chunk%tiles(tile)%tp%t_ymax

    xwidth_small = 5+(xmax-xmin)
    xwidth_large = 6+(xmax-xmin)
    yheight_small= 5+(ymax-ymin)
    yheight_large= 6+(ymax-ymin)

    arr_size = (xwidth_small * yheight_small) + & !field%density0
               (xwidth_small * yheight_small) + & !field%density1
               (xwidth_small * yheight_small) + & !field%energy0
               (xwidth_small * yheight_small) + & !field%energy1
               (xwidth_small * yheight_small) + & !field%pressure
               (xwidth_small * yheight_small) + & !field%viscosity
               (xwidth_small * yheight_small) + & !field%soundspeed
               (xwidth_large * yheight_large) + & !field%xvel0
               (xwidth_large * yheight_large) + & !field%xvel1
               (xwidth_large * yheight_large) + & !field%yvel0
               (xwidth_large * yheight_large) + & !field%yvel1
               (xwidth_large * yheight_small) + & !field%vol_flux_x
               (xwidth_large * yheight_small) + & !field%mass_flux_x
               (xwidth_small * yheight_large) + & !field%voll_flux_y
               (xwidth_small * yheight_large) + & !field%mass_flux_y
               (xwidth_small)                 + & !field%cellx
               (yheight_small)                + & !field%celly
               (xwidth_large)                 + & !field%vertexx
               (yheight_large)                + & !field%vertexy
               (xwidth_small)                 + & !field%celldx
               (yheight_small)                + & !field%celldy
               (xwidth_large)                 + & !field%vertexdx
               (yheight_large)                + & !field%vertexdy
               (xwidth_small * yheight_small) + & !field%volume
               (xwidth_large * yheight_small) + & !field%xarea
               (xwidth_small * yheight_large)     !field%yarea

    allocate(array(arr_size))

    index = 1

!$ACC DATA COPYOUT(array) 

    CALL device_pack(xmin,xmax,ymin,ymax,arr_size,array, &
      chunk%tiles(tile)%tp%field%density0   , &
      chunk%tiles(tile)%tp%field%density1   , &
      chunk%tiles(tile)%tp%field%energy0    , &
      chunk%tiles(tile)%tp%field%energy1    , &
      chunk%tiles(tile)%tp%field%pressure   , &
      chunk%tiles(tile)%tp%field%soundspeed , &
      chunk%tiles(tile)%tp%field%viscosity  , &
      chunk%tiles(tile)%tp%field%xvel0      , &
      chunk%tiles(tile)%tp%field%yvel0      , &
      chunk%tiles(tile)%tp%field%xvel1      , &
      chunk%tiles(tile)%tp%field%yvel1      , &
      chunk%tiles(tile)%tp%field%vol_flux_x , &
      chunk%tiles(tile)%tp%field%vol_flux_y , &
      chunk%tiles(tile)%tp%field%mass_flux_x, &
      chunk%tiles(tile)%tp%field%mass_flux_y, &
      chunk%tiles(tile)%tp%field%volume     , &
      chunk%tiles(tile)%tp%field%work_array1, &
      chunk%tiles(tile)%tp%field%work_array2, &
      chunk%tiles(tile)%tp%field%work_array3, &
      chunk%tiles(tile)%tp%field%work_array4, &
      chunk%tiles(tile)%tp%field%work_array5, &
      chunk%tiles(tile)%tp%field%work_array6, &
      chunk%tiles(tile)%tp%field%work_array7, &
      chunk%tiles(tile)%tp%field%cellx      , &
      chunk%tiles(tile)%tp%field%celly      , &
      chunk%tiles(tile)%tp%field%celldx     , &
      chunk%tiles(tile)%tp%field%celldy     , &
      chunk%tiles(tile)%tp%field%vertexx    , &
      chunk%tiles(tile)%tp%field%vertexdx   , &
      chunk%tiles(tile)%tp%field%vertexy    , &
      chunk%tiles(tile)%tp%field%vertexdy   , &
      chunk%tiles(tile)%tp%field%xarea      , &
      chunk%tiles(tile)%tp%field%yarea)




!$ACC END DATA






  END SUBROUTINE pack_tile





  SUBROUTINE pack_tile_meta(tile, int_array, array_size)
    USE clover_module
    USE data_module

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: tile, array_size
    INTEGER, INTENT(OUT) :: int_array(9)

    int_array(1)=array_size
    int_array(2)=chunk%tiles(tile)%tp%t_xmin
    int_array(3)=chunk%tiles(tile)%tp%t_xmax
    int_array(4)=chunk%tiles(tile)%tp%t_ymin
    int_array(5)=chunk%tiles(tile)%tp%t_ymax
    int_array(6)=chunk%tiles(tile)%tp%t_left
    int_array(7)=chunk%tiles(tile)%tp%t_right
    int_array(8)=chunk%tiles(tile)%tp%t_bottom
    int_array(9)=chunk%tiles(tile)%tp%t_top

  END SUBROUTINE pack_tile_meta


SUBROUTINE unpack_tile(tile, array, arr_size)
    USE clover_module
    USE data_module
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: tile
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(IN) :: arr_size
    INTEGER :: xmin, xmax, ymin, ymax

    INTEGER :: index, x, y

    xmin = chunk%tiles(tile)%tp%t_xmin
    xmax = chunk%tiles(tile)%tp%t_xmax
    ymin = chunk%tiles(tile)%tp%t_ymin
    ymax = chunk%tiles(tile)%tp%t_ymax


    CALL deepcopy_tile(tile)


!$ACC DATA COPYIN(array) 

   CALL device_unpack(xmin,xmax,ymin,ymax,arr_size,array, &
      chunk%tiles(tile)%tp%field%density0   , &
      chunk%tiles(tile)%tp%field%density1   , &
      chunk%tiles(tile)%tp%field%energy0    , &
      chunk%tiles(tile)%tp%field%energy1    , &
      chunk%tiles(tile)%tp%field%pressure   , &
      chunk%tiles(tile)%tp%field%soundspeed , &
      chunk%tiles(tile)%tp%field%viscosity  , &
      chunk%tiles(tile)%tp%field%xvel0      , &
      chunk%tiles(tile)%tp%field%yvel0      , &
      chunk%tiles(tile)%tp%field%xvel1      , &
      chunk%tiles(tile)%tp%field%yvel1      , &
      chunk%tiles(tile)%tp%field%vol_flux_x , &
      chunk%tiles(tile)%tp%field%vol_flux_y , &
      chunk%tiles(tile)%tp%field%mass_flux_x, &
      chunk%tiles(tile)%tp%field%mass_flux_y, &
      chunk%tiles(tile)%tp%field%volume     , &
      chunk%tiles(tile)%tp%field%work_array1, &
      chunk%tiles(tile)%tp%field%work_array2, &
      chunk%tiles(tile)%tp%field%work_array3, &
      chunk%tiles(tile)%tp%field%work_array4, &
      chunk%tiles(tile)%tp%field%work_array5, &
      chunk%tiles(tile)%tp%field%work_array6, &
      chunk%tiles(tile)%tp%field%work_array7, &
      chunk%tiles(tile)%tp%field%cellx      , &
      chunk%tiles(tile)%tp%field%celly      , &
      chunk%tiles(tile)%tp%field%celldx     , &
      chunk%tiles(tile)%tp%field%celldy     , &
      chunk%tiles(tile)%tp%field%vertexx    , &
      chunk%tiles(tile)%tp%field%vertexdx   , &
      chunk%tiles(tile)%tp%field%vertexy    , &
      chunk%tiles(tile)%tp%field%vertexdy   , &
      chunk%tiles(tile)%tp%field%xarea      , &
      chunk%tiles(tile)%tp%field%yarea)




!$ACC END DATA





  !CALL deepcopy_tile(tile)

END SUBROUTINE unpack_tile

  SUBROUTINE unpack_tile_meta(tile, int_array, array_size)

    USE clover_module
    USE data_module
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: tile, int_array(9)
    INTEGER, INTENT(OUT) :: array_size

    array_size = int_array(1)
    chunk%tiles(tile)%tp%t_xmin = int_array(2)
    chunk%tiles(tile)%tp%t_xmax = int_array(3)
    chunk%tiles(tile)%tp%t_ymin = int_array(4)
    chunk%tiles(tile)%tp%t_ymax = int_array(5)
    chunk%tiles(tile)%tp%t_left = int_array(6)
    chunk%tiles(tile)%tp%t_right = int_array(7)
    chunk%tiles(tile)%tp%t_bottom = int_array(8)
    chunk%tiles(tile)%tp%t_top = int_array(9)

  END SUBROUTINE unpack_tile_meta

  SUBROUTINE add_tile_start()

    USE clover_module
    USE data_module
  IMPLICIT NONE

  INTEGER :: tile, tmp_left, tmp_right

  DO tile=tiles_per_chunk, 1, -1
    chunk%tiles(tile+1)%tp=>chunk%tiles(tile)%tp

  END DO
  tiles_per_chunk=tiles_per_chunk+1

  CALL add_tile(1)
  call reset_neigh


  END SUBROUTINE add_tile_start

  SUBROUTINE add_tile_end()

    USE clover_module
    USE data_module
  IMPLICIT NONE

  tiles_per_chunk=tiles_per_chunk+1

  CALL add_tile(tiles_per_chunk)


    call reset_neigh


  END SUBROUTINE add_tile_end

  SUBROUTINE remove_tile_start()

    USE clover_module
    USE data_module
  IMPLICIT NONE


  INTEGER :: tile, tmp_left, tmp_right
  TYPE(tile_type), POINTER :: ltp

  ltp=>chunk%tiles(1)%tp

  tiles_per_chunk=tiles_per_chunk-1

  DO tile=1, tiles_per_chunk
    chunk%tiles(tile)%tp=>chunk%tiles(tile+1)%tp
  END DO

    call reset_neigh

  DEALLOCATE(ltp)

  END SUBROUTINE remove_tile_start

  SUBROUTINE remove_tile_end()

    USE clover_module
    USE data_module
  IMPLICIT NONE

  DEALLOCATE(chunk%tiles(tiles_per_chunk)%tp)
  tiles_per_chunk=tiles_per_chunk-1


    call reset_neigh

  END SUBROUTINE remove_tile_end

  SUBROUTINE add_tile(pos)

    USE clover_module
    USE data_module
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: pos
  TYPE(tile_type), POINTER :: ltp
   ALLOCATE(ltp)
   chunk%tiles(pos)%tp=>ltp

  END SUBROUTINE add_tile

  SUBROUTINE reset_neigh
  USE clover_module
    USE data_module
  IMPLICIT NONE

  INTEGER :: tile, tmp_left, tmp_right

  DO tile=1,tiles_per_chunk
    chunk%tiles(tile)%tp%tile_neighbours(TILE_RIGHT)=tile+1
    chunk%tiles(tile)%tp%tile_neighbours(TILE_LEFT)=tile-1
    chunk%tiles(tile)%tp%tile_neighbours(TILE_TOP)=EXTERNAL_TILE
    chunk%tiles(tile)%tp%tile_neighbours(TILE_BOTTOM)=EXTERNAL_TILE

  END DO

  chunk%tiles(tiles_per_chunk)%tp%tile_neighbours(TILE_RIGHT)=EXTERNAL_TILE
  chunk%tiles(1)%tp%tile_neighbours(TILE_LEFT)=EXTERNAL_TILE

  END SUBROUTINE reset_neigh


  SUBROUTINE deepcopy_tile(tile)

   USE openacc
    USE clover_module
    USE data_module
   
   IMPLICIT NONE

   INTEGER :: tile
    INTEGER :: xmin, xmax, ymin, ymax


    xmin = chunk%tiles(tile)%tp%t_xmin
    xmax = chunk%tiles(tile)%tp%t_xmax
    ymin = chunk%tiles(tile)%tp%t_ymin
    ymax = chunk%tiles(tile)%tp%t_ymax

    !CALL acc_copyin(chunk%tiles(tile))

    !CALL acc_attach(chunk%tiles(tile))

    !CALL acc_copyin(chunk%tiles(tile)%tp%field)

    !CALL acc_attach(chunk%tiles(tile)%tp%field)

    
    
    
    CALL acc_create(chunk%tiles(tile)%tp%field%density0)   
    CALL acc_create(chunk%tiles(tile)%tp%field%density1)   
    CALL acc_create(chunk%tiles(tile)%tp%field%energy0)    
    CALL acc_create(chunk%tiles(tile)%tp%field%energy1)    
    CALL acc_create(chunk%tiles(tile)%tp%field%pressure)   
    CALL acc_create(chunk%tiles(tile)%tp%field%soundspeed) 
    CALL acc_create(chunk%tiles(tile)%tp%field%viscosity)  
    CALL acc_create(chunk%tiles(tile)%tp%field%xvel0)      
    CALL acc_create(chunk%tiles(tile)%tp%field%yvel0)      
    CALL acc_create(chunk%tiles(tile)%tp%field%xvel1)      
    CALL acc_create(chunk%tiles(tile)%tp%field%yvel1)      
    CALL acc_create(chunk%tiles(tile)%tp%field%vol_flux_x) 
    CALL acc_create(chunk%tiles(tile)%tp%field%vol_flux_y) 
    CALL acc_create(chunk%tiles(tile)%tp%field%mass_flux_x)
    CALL acc_create(chunk%tiles(tile)%tp%field%mass_flux_y)
    CALL acc_create(chunk%tiles(tile)%tp%field%volume)     
    CALL acc_create(chunk%tiles(tile)%tp%field%work_array1)
    CALL acc_create(chunk%tiles(tile)%tp%field%work_array2)
    CALL acc_create(chunk%tiles(tile)%tp%field%work_array3)
    CALL acc_create(chunk%tiles(tile)%tp%field%work_array4)
    CALL acc_create(chunk%tiles(tile)%tp%field%work_array5)
    CALL acc_create(chunk%tiles(tile)%tp%field%work_array6)
    CALL acc_create(chunk%tiles(tile)%tp%field%work_array7)
    CALL acc_create(chunk%tiles(tile)%tp%field%cellx)      
    CALL acc_create(chunk%tiles(tile)%tp%field%celly)      
    CALL acc_create(chunk%tiles(tile)%tp%field%celldx)     
    CALL acc_create(chunk%tiles(tile)%tp%field%celldy)     
    CALL acc_create(chunk%tiles(tile)%tp%field%vertexx)    
    CALL acc_create(chunk%tiles(tile)%tp%field%vertexdx)   
    CALL acc_create(chunk%tiles(tile)%tp%field%vertexy)    
    CALL acc_create(chunk%tiles(tile)%tp%field%vertexdy)   
    CALL acc_create(chunk%tiles(tile)%tp%field%xarea)      
    CALL acc_create(chunk%tiles(tile)%tp%field%yarea)    
  
  END SUBROUTINE deepcopy_tile


  SUBROUTINE updatehost_tile(tile)

   USE openacc
    USE clover_module
    USE data_module

   IMPLICIT NONE

   INTEGER :: tile

    CALL acc_update_self(chunk%tiles(tile)%tp%field%density0)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%density1)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%energy0)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%energy1)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%pressure)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%soundspeed)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%viscosity)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%xvel0)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%yvel0)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%xvel1)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%yvel1)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%vol_flux_x)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%vol_flux_y)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%mass_flux_x)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%mass_flux_y)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%volume)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%cellx)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%celly)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%celldx)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%celldy)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%vertexx)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%vertexdx)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%vertexy)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%vertexdy)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%xarea)
    CALL acc_update_self(chunk%tiles(tile)%tp%field%yarea)


  END SUBROUTINE updatehost_tile

  END MODULE loadbalance_module

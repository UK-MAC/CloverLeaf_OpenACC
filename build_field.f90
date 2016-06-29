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

!>  @brief  Allocates the data for each mesh chunk
!>  @author Wayne Gaudin
!>  @details The data fields for the mesh chunk are allocated based on the mesh
!>  size.

SUBROUTINE build_field()

  USE clover_module

  IMPLICIT NONE

  INTEGER :: tile,j,k

  DO tile=1, tiles_per_chunk

    ALLOCATE(chunk%tiles(tile)%tp%field%density0  (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%density1  (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%energy0   (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%energy1   (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%pressure  (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%viscosity (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%soundspeed(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))

    ALLOCATE(chunk%tiles(tile)%tp%field%xvel0(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%xvel1(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%yvel0(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%yvel1(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))

    ALLOCATE(chunk%tiles(tile)%tp%field%vol_flux_x (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%mass_flux_x(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%vol_flux_y (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%mass_flux_y(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))

    ALLOCATE(chunk%tiles(tile)%tp%field%work_array1(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%work_array2(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%work_array3(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%work_array4(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%work_array5(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%work_array6(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%work_array7(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))

    ALLOCATE(chunk%tiles(tile)%tp%field%cellx   (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%celly   (chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%vertexx (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%vertexy (chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%celldx  (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%celldy  (chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%vertexdx(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%vertexdy(chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%volume  (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%xarea   (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%yarea   (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))

    ! Zeroing isn't strictly neccessary but it ensures physical pages
    ! are allocated. This prevents first touch overheads in the main code
    ! cycle which can skew timings in the first step

    !$OMP PARALLEL
    !$OMP DO
    DO k=chunk%tiles(tile)%tp%t_ymin-2,chunk%tiles(tile)%tp%t_ymax+3
      DO j=chunk%tiles(tile)%tp%t_xmin-2,chunk%tiles(tile)%tp%t_xmax+3
        chunk%tiles(tile)%tp%field%work_array1(j,k)=0.0
        chunk%tiles(tile)%tp%field%work_array2(j,k)=0.0
        chunk%tiles(tile)%tp%field%work_array3(j,k)=0.0
        chunk%tiles(tile)%tp%field%work_array4(j,k)=0.0
        chunk%tiles(tile)%tp%field%work_array5(j,k)=0.0
        chunk%tiles(tile)%tp%field%work_array6(j,k)=0.0
        chunk%tiles(tile)%tp%field%work_array7(j,k)=0.0

        chunk%tiles(tile)%tp%field%xvel0(j,k)=0.0
        chunk%tiles(tile)%tp%field%xvel1(j,k)=0.0
        chunk%tiles(tile)%tp%field%yvel0(j,k)=0.0
        chunk%tiles(tile)%tp%field%yvel1(j,k)=0.0
      ENDDO
    ENDDO
    !$OMP END DO

    !$OMP DO
    DO k=chunk%tiles(tile)%tp%t_ymin-2,chunk%tiles(tile)%tp%t_ymax+2
      DO j=chunk%tiles(tile)%tp%t_xmin-2,chunk%tiles(tile)%tp%t_xmax+2
        chunk%tiles(tile)%tp%field%density0(j,k)=0.0
        chunk%tiles(tile)%tp%field%density1(j,k)=0.0
        chunk%tiles(tile)%tp%field%energy0(j,k)=0.0
        chunk%tiles(tile)%tp%field%energy1(j,k)=0.0
        chunk%tiles(tile)%tp%field%pressure(j,k)=0.0
        chunk%tiles(tile)%tp%field%viscosity(j,k)=0.0
        chunk%tiles(tile)%tp%field%soundspeed(j,k)=0.0
        chunk%tiles(tile)%tp%field%volume(j,k)=0.0
      ENDDO
    ENDDO
    !$OMP END DO

    !$OMP DO
    DO k=chunk%tiles(tile)%tp%t_ymin-2,chunk%tiles(tile)%tp%t_ymax+2
      DO j=chunk%tiles(tile)%tp%t_xmin-2,chunk%tiles(tile)%tp%t_xmax+3
        chunk%tiles(tile)%tp%field%vol_flux_x(j,k)=0.0
        chunk%tiles(tile)%tp%field%mass_flux_x(j,k)=0.0
        chunk%tiles(tile)%tp%field%xarea(j,k)=0.0
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP DO
    DO k=chunk%tiles(tile)%tp%t_ymin-2,chunk%tiles(tile)%tp%t_ymax+3
      DO j=chunk%tiles(tile)%tp%t_xmin-2,chunk%tiles(tile)%tp%t_xmax+2
        chunk%tiles(tile)%tp%field%vol_flux_y(j,k)=0.0
        chunk%tiles(tile)%tp%field%mass_flux_y(j,k)=0.0
        chunk%tiles(tile)%tp%field%yarea(j,k)=0.0
      ENDDO
    ENDDO
    !$OMP END DO


    !$OMP DO
    DO j=chunk%tiles(tile)%tp%t_xmin-2,chunk%tiles(tile)%tp%t_xmax+2
      chunk%tiles(tile)%tp%field%cellx(j)=0.0
      chunk%tiles(tile)%tp%field%celldx(j)=0.0
    ENDDO
    !$OMP END DO
    !$OMP DO
    DO k=chunk%tiles(tile)%tp%t_ymin-2,chunk%tiles(tile)%tp%t_ymax+2
      chunk%tiles(tile)%tp%field%celly(k)=0.0
      chunk%tiles(tile)%tp%field%celldy(k)=0.0
    ENDDO
    !$OMP END DO

    !$OMP DO
    DO j=chunk%tiles(tile)%tp%t_xmin-2,chunk%tiles(tile)%tp%t_xmax+3
      chunk%tiles(tile)%tp%field%vertexx(j)=0.0
      chunk%tiles(tile)%tp%field%vertexdx(j)=0.0
    ENDDO
    !$OMP END DO
    !$OMP DO
    DO k=chunk%tiles(tile)%tp%t_ymin-2,chunk%tiles(tile)%tp%t_ymax+3
      chunk%tiles(tile)%tp%field%vertexy(k)=0.0
      chunk%tiles(tile)%tp%field%vertexdy(k)=0.0
    ENDDO
  !$OMP END DO

  !$OMP END PARALLEL  
 
  END DO
 
END SUBROUTINE build_field

SUBROUTINE build_field_tile(tile)
  USE clover_module

    IMPLICIT NONE

    INTEGER :: tile

    ALLOCATE(chunk%tiles(tile)%tp%field%density0  (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%density1  (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%energy0   (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%energy1   (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%pressure  (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%viscosity (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%soundspeed(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))

    ALLOCATE(chunk%tiles(tile)%tp%field%xvel0(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%xvel1(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%yvel0(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%yvel1(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))

    ALLOCATE(chunk%tiles(tile)%tp%field%vol_flux_x (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%mass_flux_x(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%vol_flux_y (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%mass_flux_y(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))

    ALLOCATE(chunk%tiles(tile)%tp%field%work_array1(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%work_array2(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%work_array3(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%work_array4(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%work_array5(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%work_array6(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%work_array7(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))

    ALLOCATE(chunk%tiles(tile)%tp%field%cellx   (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%celly   (chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%vertexx (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%vertexy (chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%celldx  (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%celldy  (chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%vertexdx(chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%vertexdy(chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))
    ALLOCATE(chunk%tiles(tile)%tp%field%volume  (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%xarea   (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+3, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+2))
    ALLOCATE(chunk%tiles(tile)%tp%field%yarea   (chunk%tiles(tile)%tp%t_xmin-2:chunk%tiles(tile)%tp%t_xmax+2, &
      chunk%tiles(tile)%tp%t_ymin-2:chunk%tiles(tile)%tp%t_ymax+3))

END SUBROUTINE build_field_tile

SUBROUTINE kill_field_tile(tile)
  USE clover_module
  USE OPENACC

    IMPLICIT NONE

    INTEGER :: tile

    CALL acc_delete(chunk%tiles(tile)%tp%field%density0)
    CALL acc_delete(chunk%tiles(tile)%tp%field%density1)
    CALL acc_delete(chunk%tiles(tile)%tp%field%energy0)
    CALL acc_delete(chunk%tiles(tile)%tp%field%energy1)
    CALL acc_delete(chunk%tiles(tile)%tp%field%pressure)
    CALL acc_delete(chunk%tiles(tile)%tp%field%viscosity)
    CALL acc_delete(chunk%tiles(tile)%tp%field%soundspeed)

    CALL acc_delete(chunk%tiles(tile)%tp%field%xvel0)
    CALL acc_delete(chunk%tiles(tile)%tp%field%xvel1)
    CALL acc_delete(chunk%tiles(tile)%tp%field%yvel0)
    CALL acc_delete(chunk%tiles(tile)%tp%field%yvel1)

    CALL acc_delete(chunk%tiles(tile)%tp%field%vol_flux_x)
    CALL acc_delete(chunk%tiles(tile)%tp%field%mass_flux_x)
    CALL acc_delete(chunk%tiles(tile)%tp%field%vol_flux_y)
    CALL acc_delete(chunk%tiles(tile)%tp%field%mass_flux_y)

    CALL acc_delete(chunk%tiles(tile)%tp%field%work_array1)
    CALL acc_delete(chunk%tiles(tile)%tp%field%work_array2)
    CALL acc_delete(chunk%tiles(tile)%tp%field%work_array3)
    CALL acc_delete(chunk%tiles(tile)%tp%field%work_array4)
    CALL acc_delete(chunk%tiles(tile)%tp%field%work_array5)
    CALL acc_delete(chunk%tiles(tile)%tp%field%work_array6)
    CALL acc_delete(chunk%tiles(tile)%tp%field%work_array7)

    CALL acc_delete(chunk%tiles(tile)%tp%field%cellx)
    CALL acc_delete(chunk%tiles(tile)%tp%field%celly)
    CALL acc_delete(chunk%tiles(tile)%tp%field%vertexx)
    CALL acc_delete(chunk%tiles(tile)%tp%field%vertexy )
    CALL acc_delete(chunk%tiles(tile)%tp%field%celldx )
    CALL acc_delete(chunk%tiles(tile)%tp%field%celldy )
    CALL acc_delete(chunk%tiles(tile)%tp%field%vertexdx)
    CALL acc_delete(chunk%tiles(tile)%tp%field%vertexdy)
    CALL acc_delete(chunk%tiles(tile)%tp%field%volume)
    CALL acc_delete(chunk%tiles(tile)%tp%field%xarea)
    CALL acc_delete(chunk%tiles(tile)%tp%field%yarea)

    DEALLOCATE(chunk%tiles(tile)%tp%field%density0)
    DEALLOCATE(chunk%tiles(tile)%tp%field%density1)
    DEALLOCATE(chunk%tiles(tile)%tp%field%energy0)
    DEALLOCATE(chunk%tiles(tile)%tp%field%energy1)
    DEALLOCATE(chunk%tiles(tile)%tp%field%pressure)
    DEALLOCATE(chunk%tiles(tile)%tp%field%viscosity)
    DEALLOCATE(chunk%tiles(tile)%tp%field%soundspeed)

    DEALLOCATE(chunk%tiles(tile)%tp%field%xvel0)
    DEALLOCATE(chunk%tiles(tile)%tp%field%xvel1)
    DEALLOCATE(chunk%tiles(tile)%tp%field%yvel0)
    DEALLOCATE(chunk%tiles(tile)%tp%field%yvel1)

    DEALLOCATE(chunk%tiles(tile)%tp%field%vol_flux_x)
    DEALLOCATE(chunk%tiles(tile)%tp%field%mass_flux_x)
    DEALLOCATE(chunk%tiles(tile)%tp%field%vol_flux_y)
    DEALLOCATE(chunk%tiles(tile)%tp%field%mass_flux_y)

    DEALLOCATE(chunk%tiles(tile)%tp%field%work_array1)
    DEALLOCATE(chunk%tiles(tile)%tp%field%work_array2)
    DEALLOCATE(chunk%tiles(tile)%tp%field%work_array3)
    DEALLOCATE(chunk%tiles(tile)%tp%field%work_array4)
    DEALLOCATE(chunk%tiles(tile)%tp%field%work_array5)
    DEALLOCATE(chunk%tiles(tile)%tp%field%work_array6)
    DEALLOCATE(chunk%tiles(tile)%tp%field%work_array7)

    DEALLOCATE(chunk%tiles(tile)%tp%field%cellx)
    DEALLOCATE(chunk%tiles(tile)%tp%field%celly)
    DEALLOCATE(chunk%tiles(tile)%tp%field%vertexx)
    DEALLOCATE(chunk%tiles(tile)%tp%field%vertexy )
    DEALLOCATE(chunk%tiles(tile)%tp%field%celldx )
    DEALLOCATE(chunk%tiles(tile)%tp%field%celldy )
    DEALLOCATE(chunk%tiles(tile)%tp%field%vertexdx)
    DEALLOCATE(chunk%tiles(tile)%tp%field%vertexdy)
    DEALLOCATE(chunk%tiles(tile)%tp%field%volume)
    DEALLOCATE(chunk%tiles(tile)%tp%field%xarea)
    DEALLOCATE(chunk%tiles(tile)%tp%field%yarea)

    

END SUBROUTINE kill_field_tile

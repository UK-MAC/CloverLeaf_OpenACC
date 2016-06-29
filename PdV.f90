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

MODULE PdV_module

CONTAINS

  SUBROUTINE PdV(predict)

    USE clover_module
    USE report_module
    USE PdV_kernel_module
    USE revert_module
    USE update_halo_module
    USE ideal_gas_module

    IMPLICIT NONE

    LOGICAL :: predict

    INTEGER :: prdct

    INTEGER :: tile
    INTEGER :: fields(NUM_FIELDS)

    REAL(KIND=8) :: kernel_time,timer

    error_condition=0 ! Not used yet due to issue with OpenA reduction

    IF(profiler_on) kernel_time=timer()

    IF(predict) THEN
      prdct=0
    ELSE
      prdct=1
    ENDIF

      DO tile=1,tiles_per_chunk


        CALL PdV_kernel(predict,                  &
          chunk%tiles(tile)%tp%t_xmin,      &
          chunk%tiles(tile)%tp%t_xmax,      &
          chunk%tiles(tile)%tp%t_ymin,      &
          chunk%tiles(tile)%tp%t_ymax,      &
          dt,                         &
          chunk%tiles(tile)%tp%field%xarea,      &
          chunk%tiles(tile)%tp%field%yarea,      &
          chunk%tiles(tile)%tp%field%volume ,    &
          chunk%tiles(tile)%tp%field%density0,   &
          chunk%tiles(tile)%tp%field%density1,   &
          chunk%tiles(tile)%tp%field%energy0,    &
          chunk%tiles(tile)%tp%field%energy1,    &
          chunk%tiles(tile)%tp%field%pressure,   &
          chunk%tiles(tile)%tp%field%viscosity,  &
          chunk%tiles(tile)%tp%field%xvel0,      &
          chunk%tiles(tile)%tp%field%xvel1,      &
          chunk%tiles(tile)%tp%field%yvel0,      &
          chunk%tiles(tile)%tp%field%yvel1,      &
          chunk%tiles(tile)%tp%field%work_array1 )


      ENDDO
  

    

    CALL clover_check_error(error_condition)
    IF(profiler_on) profiler%PdV=profiler%PdV+(timer()-kernel_time)

    IF(error_condition.EQ.1) THEN
      CALL report_error('PdV','error in PdV')
    ENDIF

    IF(predict)THEN
      IF(profiler_on) kernel_time=timer()
      DO tile=1,tiles_per_chunk
        CALL ideal_gas(tile,.TRUE.)
      ENDDO

      IF(profiler_on) profiler%ideal_gas=profiler%ideal_gas+(timer()-kernel_time)
      fields=0
      fields(FIELD_PRESSURE)=1
      CALL update_halo(fields,1)
    ENDIF

    IF ( predict ) THEN
      IF(profiler_on) kernel_time=timer()
      CALL revert()
      IF(profiler_on) profiler%revert=profiler%revert+(timer()-kernel_time)
    ENDIF

  END SUBROUTINE PdV

END MODULE PdV_module

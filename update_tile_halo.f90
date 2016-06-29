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

!>  @brief Driver for the halo updates
!>  @author Wayne Gaudin
!>  @details Invokes the kernels for the internal and external halo cells for
!>  the fields specified.

MODULE update_tile_halo_module

CONTAINS

  SUBROUTINE update_tile_halo(fields,depth)

    USE clover_module
    USE update_tile_halo_kernel_module

    IMPLICIT NONE

    INTEGER :: tile,fields(NUM_FIELDS), depth
    INTEGER :: t_left, t_right, t_up, t_down

        !TODO: fix the chunk comms phase
      !CALL clover_exchange(fields,depth)



    ! Update Top Bottom - Real to Real


    DO tile=1,tiles_per_chunk
      t_up   =chunk%tiles(tile)%tp%tile_neighbours(TILE_TOP)
      t_down =chunk%tiles(tile)%tp%tile_neighbours(TILE_BOTTOM)

      IF(t_up.NE.EXTERNAL_TILE) THEN
        call update_tile_halo_t_kernel(                                        &
          chunk%tiles(tile)%tp%t_xmin,                              &
          chunk%tiles(tile)%tp%t_xmax,                              &
          chunk%tiles(tile)%tp%t_ymin,                              &
          chunk%tiles(tile)%tp%t_ymax,                              &
          chunk%tiles(tile)%tp%field%density0,                      &
          chunk%tiles(tile)%tp%field%energy0,                       &
          chunk%tiles(tile)%tp%field%pressure,                      &
          chunk%tiles(tile)%tp%field%viscosity,                     &
          chunk%tiles(tile)%tp%field%soundspeed,                    &
          chunk%tiles(tile)%tp%field%density1,                      &
          chunk%tiles(tile)%tp%field%energy1,                       &
          chunk%tiles(tile)%tp%field%xvel0,                         &
          chunk%tiles(tile)%tp%field%yvel0,                         &
          chunk%tiles(tile)%tp%field%xvel1,                         &
          chunk%tiles(tile)%tp%field%yvel1,                         &
          chunk%tiles(tile)%tp%field%vol_flux_x,                    &
          chunk%tiles(tile)%tp%field%vol_flux_y,                    &
          chunk%tiles(tile)%tp%field%mass_flux_x,                   &
          chunk%tiles(tile)%tp%field%mass_flux_y,                   &
          chunk%tiles(t_up)%tp%t_xmin,                           &
          chunk%tiles(t_up)%tp%t_xmax,                           &
          chunk%tiles(t_up)%tp%t_ymin,                           &
          chunk%tiles(t_up)%tp%t_ymax,                           &
          chunk%tiles(t_up)%tp%field%density0,                   &
          chunk%tiles(t_up)%tp%field%energy0,                    &
          chunk%tiles(t_up)%tp%field%pressure,                   &
          chunk%tiles(t_up)%tp%field%viscosity,                  &
          chunk%tiles(t_up)%tp%field%soundspeed,                 &
          chunk%tiles(t_up)%tp%field%density1,                   &
          chunk%tiles(t_up)%tp%field%energy1,                    &
          chunk%tiles(t_up)%tp%field%xvel0,                      &
          chunk%tiles(t_up)%tp%field%yvel0,                      &
          chunk%tiles(t_up)%tp%field%xvel1,                      &
          chunk%tiles(t_up)%tp%field%yvel1,                      &
          chunk%tiles(t_up)%tp%field%vol_flux_x,                 &
          chunk%tiles(t_up)%tp%field%vol_flux_y,                 &
          chunk%tiles(t_up)%tp%field%mass_flux_x,                &
          chunk%tiles(t_up)%tp%field%mass_flux_y,                &
          fields,                                                   &
          depth)
      !ELSE
     
      END IF

      IF(t_down.NE.EXTERNAL_TILE) THEN
        call update_tile_halo_b_kernel(                                        &
          chunk%tiles(tile)%tp%t_xmin,                              &
          chunk%tiles(tile)%tp%t_xmax,                              &
          chunk%tiles(tile)%tp%t_ymin,                              &
          chunk%tiles(tile)%tp%t_ymax,                              &
          chunk%tiles(tile)%tp%field%density0,                      &
          chunk%tiles(tile)%tp%field%energy0,                       &
          chunk%tiles(tile)%tp%field%pressure,                      &
          chunk%tiles(tile)%tp%field%viscosity,                     &
          chunk%tiles(tile)%tp%field%soundspeed,                    &
          chunk%tiles(tile)%tp%field%density1,                      &
          chunk%tiles(tile)%tp%field%energy1,                       &
          chunk%tiles(tile)%tp%field%xvel0,                         &
          chunk%tiles(tile)%tp%field%yvel0,                         &
          chunk%tiles(tile)%tp%field%xvel1,                         &
          chunk%tiles(tile)%tp%field%yvel1,                         &
          chunk%tiles(tile)%tp%field%vol_flux_x,                    &
          chunk%tiles(tile)%tp%field%vol_flux_y,                    &
          chunk%tiles(tile)%tp%field%mass_flux_x,                   &
          chunk%tiles(tile)%tp%field%mass_flux_y,                   &
          chunk%tiles(t_down)%tp%t_xmin,                           &
          chunk%tiles(t_down)%tp%t_xmax,                           &
          chunk%tiles(t_down)%tp%t_ymin,                           &
          chunk%tiles(t_down)%tp%t_ymax,                           &
          chunk%tiles(t_down)%tp%field%density0,                   &
          chunk%tiles(t_down)%tp%field%energy0,                    &
          chunk%tiles(t_down)%tp%field%pressure,                   &
          chunk%tiles(t_down)%tp%field%viscosity,                  &
          chunk%tiles(t_down)%tp%field%soundspeed,                 &
          chunk%tiles(t_down)%tp%field%density1,                   &
          chunk%tiles(t_down)%tp%field%energy1,                    &
          chunk%tiles(t_down)%tp%field%xvel0,                      &
          chunk%tiles(t_down)%tp%field%yvel0,                      &
          chunk%tiles(t_down)%tp%field%xvel1,                      &
          chunk%tiles(t_down)%tp%field%yvel1,                      &
          chunk%tiles(t_down)%tp%field%vol_flux_x,                 &
          chunk%tiles(t_down)%tp%field%vol_flux_y,                 &
          chunk%tiles(t_down)%tp%field%mass_flux_x,                &
          chunk%tiles(t_down)%tp%field%mass_flux_y,                &
          fields,                                                   &
          depth)
      !ELSE

      END IF

  
    END DO

    ! Update Left Right - Ghost, Real, Ghost - > Real

    DO tile=1,tiles_per_chunk
      t_left   =chunk%tiles(tile)%tp%tile_neighbours(TILE_LEFT)
      t_right  =chunk%tiles(tile)%tp%tile_neighbours(TILE_RIGHT)
  
      IF(t_left.NE.EXTERNAL_TILE) THEN
        call update_tile_halo_l_kernel(                                        &
          chunk%tiles(tile)%tp%t_xmin,                              &
          chunk%tiles(tile)%tp%t_xmax,                              &
          chunk%tiles(tile)%tp%t_ymin,                              &
          chunk%tiles(tile)%tp%t_ymax,                              &
          chunk%tiles(tile)%tp%field%density0,                      &
          chunk%tiles(tile)%tp%field%energy0,                       &
          chunk%tiles(tile)%tp%field%pressure,                      &
          chunk%tiles(tile)%tp%field%viscosity,                     &
          chunk%tiles(tile)%tp%field%soundspeed,                    &
          chunk%tiles(tile)%tp%field%density1,                      &
          chunk%tiles(tile)%tp%field%energy1,                       &
          chunk%tiles(tile)%tp%field%xvel0,                         &
          chunk%tiles(tile)%tp%field%yvel0,                         &
          chunk%tiles(tile)%tp%field%xvel1,                         &
          chunk%tiles(tile)%tp%field%yvel1,                         &
          chunk%tiles(tile)%tp%field%vol_flux_x,                    &
          chunk%tiles(tile)%tp%field%vol_flux_y,                    &
          chunk%tiles(tile)%tp%field%mass_flux_x,                   &
          chunk%tiles(tile)%tp%field%mass_flux_y,                   &
          chunk%tiles(t_left)%tp%t_xmin,                           &
          chunk%tiles(t_left)%tp%t_xmax,                           &
          chunk%tiles(t_left)%tp%t_ymin,                           &
          chunk%tiles(t_left)%tp%t_ymax,                           &
          chunk%tiles(t_left)%tp%field%density0,                   &
          chunk%tiles(t_left)%tp%field%energy0,                    &
          chunk%tiles(t_left)%tp%field%pressure,                   &
          chunk%tiles(t_left)%tp%field%viscosity,                  &
          chunk%tiles(t_left)%tp%field%soundspeed,                 &
          chunk%tiles(t_left)%tp%field%density1,                   &
          chunk%tiles(t_left)%tp%field%energy1,                    &
          chunk%tiles(t_left)%tp%field%xvel0,                      &
          chunk%tiles(t_left)%tp%field%yvel0,                      &
          chunk%tiles(t_left)%tp%field%xvel1,                      &
          chunk%tiles(t_left)%tp%field%yvel1,                      &
          chunk%tiles(t_left)%tp%field%vol_flux_x,                 &
          chunk%tiles(t_left)%tp%field%vol_flux_y,                 &
          chunk%tiles(t_left)%tp%field%mass_flux_x,                &
          chunk%tiles(t_left)%tp%field%mass_flux_y,                &
          fields,                                                   &
          depth)
      !ELSE

      END IF

      IF(t_right.NE.EXTERNAL_TILE) THEN
        call update_tile_halo_r_kernel(                                        &
          chunk%tiles(tile)%tp%t_xmin,                              &
          chunk%tiles(tile)%tp%t_xmax,                              &
          chunk%tiles(tile)%tp%t_ymin,                              &
          chunk%tiles(tile)%tp%t_ymax,                              &
          chunk%tiles(tile)%tp%field%density0,                      &
          chunk%tiles(tile)%tp%field%energy0,                       &
          chunk%tiles(tile)%tp%field%pressure,                      &
          chunk%tiles(tile)%tp%field%viscosity,                     &
          chunk%tiles(tile)%tp%field%soundspeed,                    &
          chunk%tiles(tile)%tp%field%density1,                      &
          chunk%tiles(tile)%tp%field%energy1,                       &
          chunk%tiles(tile)%tp%field%xvel0,                         &
          chunk%tiles(tile)%tp%field%yvel0,                         &
          chunk%tiles(tile)%tp%field%xvel1,                         &
          chunk%tiles(tile)%tp%field%yvel1,                         &
          chunk%tiles(tile)%tp%field%vol_flux_x,                    &
          chunk%tiles(tile)%tp%field%vol_flux_y,                    &
          chunk%tiles(tile)%tp%field%mass_flux_x,                   &
          chunk%tiles(tile)%tp%field%mass_flux_y,                   &
          chunk%tiles(t_right)%tp%t_xmin,                           &
          chunk%tiles(t_right)%tp%t_xmax,                           &
          chunk%tiles(t_right)%tp%t_ymin,                           &
          chunk%tiles(t_right)%tp%t_ymax,                           &
          chunk%tiles(t_right)%tp%field%density0,                   &
          chunk%tiles(t_right)%tp%field%energy0,                    &
          chunk%tiles(t_right)%tp%field%pressure,                   &
          chunk%tiles(t_right)%tp%field%viscosity,                  &
          chunk%tiles(t_right)%tp%field%soundspeed,                 &
          chunk%tiles(t_right)%tp%field%density1,                   &
          chunk%tiles(t_right)%tp%field%energy1,                    &
          chunk%tiles(t_right)%tp%field%xvel0,                      &
          chunk%tiles(t_right)%tp%field%yvel0,                      &
          chunk%tiles(t_right)%tp%field%xvel1,                      &
          chunk%tiles(t_right)%tp%field%yvel1,                      &
          chunk%tiles(t_right)%tp%field%vol_flux_x,                 &
          chunk%tiles(t_right)%tp%field%vol_flux_y,                 &
          chunk%tiles(t_right)%tp%field%mass_flux_x,                &
          chunk%tiles(t_right)%tp%field%mass_flux_y,                &
          fields,                                                   &
          depth)
      !ELSE

      END IF


    END DO





  END SUBROUTINE update_tile_halo

END MODULE update_tile_halo_module

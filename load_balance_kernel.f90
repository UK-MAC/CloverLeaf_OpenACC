MODULE load_balance_kernel_module

  IMPLICIT NONE

  CONTAINS

 SUBROUTINE device_pack(xmin,xmax,ymin,ymax,asize,array,density0,density1,energy0,energy1,&
                        pressure,soundspeed,viscosity,xvel0,yvel0,xvel1,yvel1,vol_flux_x,vol_flux_y,mass_flux_x,&
                        mass_flux_y,volume,work_array1,work_array2,work_array3,work_array4,work_array5,work_array6,work_array7,&
                        cellx,celly,celldx,celldy,vertexx,vertexdx,vertexy,vertexdy,xarea,yarea)
   
   IMPLICIT NONE

    INTEGER               :: xmin,xmax,ymin,ymax, asize
    REAL(KIND=8), DIMENSION(asize) :: array
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+2) :: density0   
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+2) :: density1   
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+2) :: energy0    
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+2) :: energy1    
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+2) :: pressure  
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+2) :: soundspeed 
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+2) :: viscosity  
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: xvel0      
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: yvel0      
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: xvel1      
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: yvel1      
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+2) :: vol_flux_x 
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+3) :: vol_flux_y 
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+2) :: mass_flux_x
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+3) :: mass_flux_y
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+2) :: volume     
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: work_array1
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: work_array2
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: work_array3
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: work_array4
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: work_array5
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: work_array6
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: work_array7
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2) :: cellx      
    REAL(KIND=8), DIMENSION(ymin-2:ymax+2) :: celly      
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2) :: celldx     
    REAL(KIND=8), DIMENSION(ymin-2:ymax+2) :: celldy     
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3) :: vertexx    
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3) :: vertexdx   
    REAL(KIND=8), DIMENSION(ymin-2:ymax+3) :: vertexy    
    REAL(KIND=8), DIMENSION(ymin-2:ymax+3) :: vertexdy   
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+2) :: xarea      
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+3) :: yarea

    INTEGER :: index, offset, x ,y
    INTEGER :: xwidth_small, xwidth_large, yheight_small,yheight_large 


    xwidth_small = 5+(xmax-xmin)
    xwidth_large = 6+(xmax-xmin)
    yheight_small= 5+(ymax-ymin)
    yheight_large= 6+(ymax-ymin)

!$ACC DATA &
!$ACC PRESENT(array      ) &
!$ACC PRESENT(density0   ) &
!$ACC PRESENT(density1   ) &
!$ACC PRESENT(energy0    ) &
!$ACC PRESENT(energy1    ) &
!$ACC PRESENT(pressure   ) &
!$ACC PRESENT(soundspeed ) &
!$ACC PRESENT(viscosity  ) &
!$ACC PRESENT(xvel0      ) &
!$ACC PRESENT(yvel0      ) &
!$ACC PRESENT(xvel1      ) &
!$ACC PRESENT(yvel1      ) &
!$ACC PRESENT(vol_flux_x ) &
!$ACC PRESENT(vol_flux_y ) &
!$ACC PRESENT(mass_flux_x) &
!$ACC PRESENT(mass_flux_y) &
!$ACC PRESENT(volume     ) &
!$ACC PRESENT(work_array1) &
!$ACC PRESENT(work_array2) &
!$ACC PRESENT(work_array3) &
!$ACC PRESENT(work_array4) &
!$ACC PRESENT(work_array5) &
!$ACC PRESENT(work_array6) &
!$ACC PRESENT(work_array7) &
!$ACC PRESENT(cellx      ) &
!$ACC PRESENT(celly      ) &
!$ACC PRESENT(celldx     ) &
!$ACC PRESENT(celldy     ) &
!$ACC PRESENT(vertexx    ) &
!$ACC PRESENT(vertexdx   ) &
!$ACC PRESENT(vertexy    ) &
!$ACC PRESENT(vertexdy   ) &
!$ACC PRESENT(xarea      ) &
!$ACC PRESENT(yarea)

    offset=1
!$ACC KERNELS


!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        array(index) = density0(x,y)
      END DO
    END DO
    
!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_small)
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        array(index) = density1(x,y)
      END DO
    END DO
    
!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_small)
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        array(index) = energy0(x,y)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_small)
    
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        array(index) = energy1(x,y)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_small)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        array(index) = pressure(x,y)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_small)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        array(index) = viscosity(x,y)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_small)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        array(index) = soundspeed(x,y)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_small)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+3
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+3
        index=offset+((y-(ymin-2))*(xwidth_large))+(x-(xmin-2))
        array(index) = xvel0(x,y)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large * yheight_large)
    
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+3
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+3
        index=offset+((y-(ymin-2))*(xwidth_large))+(x-(xmin-2))
        array(index) = xvel1(x,y)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large * yheight_large)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+3
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+3
        index=offset+((y-(ymin-2))*(xwidth_large))+(x-(xmin-2))
        array(index) = yvel0(x,y)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large * yheight_large)
    
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+3
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+3
        index=offset+((y-(ymin-2))*(xwidth_large))+(x-(xmin-2))
        array(index) = yvel1(x,y)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large * yheight_large)
    
!$ACC KERNELS 
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+3
        index=offset+((y-(ymin-2))*(xwidth_large))+(x-(xmin-2))
        array(index) = vol_flux_x(x,y)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large * yheight_small)
    
!$ACC KERNELS 
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+3
        index=offset+((y-(ymin-2))*(xwidth_large))+(x-(xmin-2))
        array(index) = mass_flux_x(x,y)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large * yheight_small)
    
!$ACC KERNELS 
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+3
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        array(index) = vol_flux_y(x,y)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_large)
    
!$ACC KERNELS 
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+3
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        array(index) = mass_flux_y(x,y)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_large)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT PRIVATE(index)
    DO x=xmin-2, xmax+2
        index=offset+(x-(xmin-2))
      array(index) = cellx(x)
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT PRIVATE(index)
    DO y=ymin-2, ymax+2
        index=offset+(y-(ymin-2))
      array(index) = celly(y)
    END DO

!$ACC END KERNELS
    offset = offset + (yheight_small)
    
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT PRIVATE(index)
    DO x=xmin-2, xmax+3
        index=offset+(x-(xmin-2))
      array(index) = vertexx(x)
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT PRIVATE(index)
    DO y=ymin-2, ymax+3
        index=offset+(y-(ymin-2))
      array(index) = vertexy(y)
    END DO

!$ACC END KERNELS
    offset = offset + (yheight_large)
    
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT PRIVATE(index)
    DO x=xmin-2, xmax+2
        index=offset+(x-(xmin-2))
    
      array(index) = celldx(x)
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT PRIVATE(index)
    DO y=ymin-2, ymax+2
        index=offset+(y-(ymin-2))
      array(index) = celldy(y)
    END DO

!$ACC END KERNELS
    offset = offset + (yheight_small)
    
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT PRIVATE(index)
    DO x=xmin-2, xmax+3
    
        index=offset+(x-(xmin-2))
      array(index) = vertexdx(x)
      index=index+1
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large)
    
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT PRIVATE(index)
    DO y=ymin-2, ymax+3
    
        index=offset+(y-(ymin-2))
      array(index) = vertexdy(y)
    END DO

!$ACC END KERNELS
    offset = offset + (yheight_large)
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        array(index) = volume(x,y)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_small)
    
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+3
        index=offset+((y-(ymin-2))*(xwidth_large))+(x-(xmin-2))
        array(index) = xarea(x,y)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large * yheight_small)
    
    
!$ACC KERNELS
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+3
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        array(index) = yarea(x,y)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_large)
    
    
!$ACC END DATA




 END SUBROUTINE device_pack


 SUBROUTINE device_unpack(xmin,xmax,ymin,ymax,asize,array,density0,density1,energy0,energy1,&
                        pressure,soundspeed,viscosity,xvel0,yvel0,xvel1,yvel1,vol_flux_x,vol_flux_y,mass_flux_x,&
                        mass_flux_y,volume,work_array1,work_array2,work_array3,work_array4,work_array5,work_array6,work_array7,&
                        cellx,celly,celldx,celldy,vertexx,vertexdx,vertexy,vertexdy,xarea,yarea)
   
   IMPLICIT NONE

    INTEGER               :: xmin,xmax,ymin,ymax, asize
    REAL(KIND=8), DIMENSION(asize) :: array
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+2) :: density0   
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+2) :: density1   
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+2) :: energy0    
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+2) :: energy1    
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+2) :: pressure  
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+2) :: soundspeed 
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+2) :: viscosity  
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: xvel0      
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: yvel0      
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: xvel1      
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: yvel1      
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+2) :: vol_flux_x 
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+3) :: vol_flux_y 
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+2) :: mass_flux_x
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+3) :: mass_flux_y
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+2) :: volume     
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: work_array1
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: work_array2
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: work_array3
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: work_array4
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: work_array5
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: work_array6
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+3) :: work_array7
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2) :: cellx      
    REAL(KIND=8), DIMENSION(ymin-2:ymax+2) :: celly      
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2) :: celldx     
    REAL(KIND=8), DIMENSION(ymin-2:ymax+2) :: celldy     
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3) :: vertexx    
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3) :: vertexdx   
    REAL(KIND=8), DIMENSION(ymin-2:ymax+3) :: vertexy    
    REAL(KIND=8), DIMENSION(ymin-2:ymax+3) :: vertexdy   
    REAL(KIND=8), DIMENSION(xmin-2:xmax+3,ymin-2:ymax+2) :: xarea      
    REAL(KIND=8), DIMENSION(xmin-2:xmax+2,ymin-2:ymax+3) :: yarea

    INTEGER :: index, offset, x ,y
    INTEGER :: xwidth_small, xwidth_large, yheight_small,yheight_large 


    xwidth_small = 5+(xmax-xmin)
    xwidth_large = 6+(xmax-xmin)
    yheight_small= 5+(ymax-ymin)
    yheight_large= 6+(ymax-ymin)

!$ACC DATA &
!$ACC PRESENT(array) &
!$ACC PRESENT(density0   ) &
!$ACC PRESENT(density1   ) &
!$ACC PRESENT(energy0    ) &
!$ACC PRESENT(energy1    ) &
!$ACC PRESENT(pressure  ) &
!$ACC PRESENT(soundspeed ) &
!$ACC PRESENT(viscosity  ) &
!$ACC PRESENT(xvel0      ) &
!$ACC PRESENT(yvel0      ) &
!$ACC PRESENT(xvel1      ) &
!$ACC PRESENT(yvel1      ) &
!$ACC PRESENT(vol_flux_x ) &
!$ACC PRESENT(vol_flux_y ) &
!$ACC PRESENT(mass_flux_x) &
!$ACC PRESENT(mass_flux_y) &
!$ACC PRESENT(volume     ) &
!$ACC PRESENT(work_array1) &
!$ACC PRESENT(work_array2) &
!$ACC PRESENT(work_array3) &
!$ACC PRESENT(work_array4) &
!$ACC PRESENT(work_array5) &
!$ACC PRESENT(work_array6) &
!$ACC PRESENT(work_array7) &
!$ACC PRESENT(cellx      ) &
!$ACC PRESENT(celly      ) &
!$ACC PRESENT(celldx     ) &
!$ACC PRESENT(celldy     ) &
!$ACC PRESENT(vertexx    ) &
!$ACC PRESENT(vertexdx   ) &
!$ACC PRESENT(vertexy    ) &
!$ACC PRESENT(vertexdy   ) &
!$ACC PRESENT(xarea      ) &
!$ACC PRESENT(yarea)

    offset=1
!$ACC KERNELS


!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        density0(x,y) = array(index)
      END DO
    END DO
    
    offset = offset + (xwidth_small * yheight_small)
    
!$ACC END KERNELS
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        density1(x,y) = array(index)
      END DO
    END DO
    
!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_small)
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        energy0(x,y) = array(index)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_small)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        energy1(x,y) = array(index)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_small)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        pressure(x,y) = array(index)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_small)
    
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        viscosity(x,y) = array(index)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_small)
    
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        soundspeed(x,y) = array(index)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_small)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+3
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+3
        index=offset+((y-(ymin-2))*(xwidth_large))+(x-(xmin-2))
        xvel0(x,y) = array(index)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large * yheight_large)
    
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+3
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+3
        index=offset+((y-(ymin-2))*(xwidth_large))+(x-(xmin-2))
        xvel1(x,y) = array(index)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large * yheight_large)
    
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+3
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+3
        index=offset+((y-(ymin-2))*(xwidth_large))+(x-(xmin-2))
        yvel0(x,y) = array(index)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large * yheight_large)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+3
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+3
        index=offset+((y-(ymin-2))*(xwidth_large))+(x-(xmin-2))
        yvel1(x,y) = array(index)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large * yheight_large)
    
!$ACC KERNELS 
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+3
        index=offset+((y-(ymin-2))*(xwidth_large))+(x-(xmin-2))
        vol_flux_x(x,y) = array(index)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large * yheight_small)
    
!$ACC KERNELS 
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+3
        index=offset+((y-(ymin-2))*(xwidth_large))+(x-(xmin-2))
        mass_flux_x(x,y) = array(index)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large * yheight_small)
    
!$ACC KERNELS 
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+3
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        vol_flux_y(x,y) = array(index)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_large)
    
!$ACC KERNELS 
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+3
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        mass_flux_y(x,y) = array(index)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_large)
     
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT PRIVATE(index)
    DO x=xmin-2, xmax+2
        index=offset+(x-(xmin-2))
      cellx(x) = array(index)
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT PRIVATE(index)
    DO y=ymin-2, ymax+2
        index=offset+(y-(ymin-2))
      celly(y) = array(index)
    END DO

!$ACC END KERNELS
    offset = offset + (yheight_small)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT PRIVATE(index)
    DO x=xmin-2, xmax+3
        index=offset+(x-(xmin-2))
      vertexx(x) = array(index)
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT PRIVATE(index)
    DO y=ymin-2, ymax+3
        index=offset+(y-(ymin-2))
      vertexy(y) = array(index)
    END DO

!$ACC END KERNELS
    offset = offset + (yheight_large)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT PRIVATE(index)
    DO x=xmin-2, xmax+2
        index=offset+(x-(xmin-2))
    
      celldx(x) = array(index)
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT PRIVATE(index)
    DO y=ymin-2, ymax+2
        index=offset+(y-(ymin-2))
      celldy(y) = array(index)
    END DO

!$ACC END KERNELS
    offset = offset + (yheight_small)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT PRIVATE(index)
    DO x=xmin-2, xmax+3
    
        index=offset+(x-(xmin-2))
      vertexdx(x) = array(index)
      index=index+1
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT PRIVATE(index)
    DO y=ymin-2, ymax+3
    
        index=offset+(y-(ymin-2))
      vertexdy(y) = array(index)
    END DO

!$ACC END KERNELS
    offset = offset + (yheight_large)
    
!$ACC KERNELS 
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        volume(x,y) = array(index)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_small)
    
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+2
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+3
        index=offset+((y-(ymin-2))*(xwidth_large))+(x-(xmin-2))
        xarea(x,y) = array(index)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_large * yheight_small)
   
!$ACC KERNELS 
    
    
!$ACC LOOP INDEPENDENT
    DO y=ymin-2, ymax+3
!$ACC LOOP INDEPENDENT PRIVATE(index)
      DO x=xmin-2, xmax+2
        index=offset+((y-(ymin-2))*(xwidth_small))+(x-(xmin-2))
        yarea(x,y) = array(index)
      END DO
    END DO

!$ACC END KERNELS
    offset = offset + (xwidth_small * yheight_large)
    
    
!$ACC END DATA




 END SUBROUTINE device_unpack

END MODULE load_balance_kernel_module



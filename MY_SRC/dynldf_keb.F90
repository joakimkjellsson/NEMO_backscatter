MODULE dynldf_keb
   !!======================================================================
   !!                   ***  MODULE  dynldf_lap_blp  ***
   !! Ocean dynamics:  lateral viscosity trend (laplacian and bilaplacian)
   !!======================================================================
   !! History : 3.7  ! 2014-01  (G. Madec, S. Masson)  Original code, re-entrant laplacian
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_ldf_keb   : update the momentum trend with backscatter KE
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE ldfdyn         ! lateral diffusion: eddy viscosity coef.
   USE ldfslp         ! iso-neutral slopes 
   USE zdf_oce        ! ocean vertical physics
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC dyn_ldf_keb  ! called by dynldf.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dynldf_lap_blp.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_ldf_keb( kt, pub, pvb, pua, pva, psgs_ke )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dyn_ldf_lap  ***
      !!                       
      !! ** Purpose :   Compute the before horizontal momentum diffusive 
      !!      trend and add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The Laplacian operator apply on horizontal velocity is 
      !!      writen as :   grad_h( ahmt div_h(U )) - curl_h( ahmf curl_z(U) ) 
      !!
      !! ** Action : - pua, pva increased by the harmonic operator applied on pub, pvb.
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   kpass      ! =1/2 first or second passage
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pub, pvb   ! before velocity   [m/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   psgs_ke    ! sub-grid scale KE [m2/s2]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pua, pva   ! velocity trend    [m/s2]
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zsign        ! local scalars
      REAL(wp) ::   zua, zva     ! local scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   zcur, zdiv
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_ldf : KE backscatter with iso-level harmonic (laplacian) operator' 
         WRITE(numout,*) '~~~~~~~ '
      ENDIF
      !
      ! Calculate the Laplacian
      ! 
!!jk do this very easy. Spatially uniform anti-viscosity (replace ahm by zahm, local variable)
!!jk think we should have another set of ahm for backscatter, i.e. separate from lateral viscosity            
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         DO jj = 2, jpj
            DO ji = fs_2, jpi   ! vector opt.
               !                                      ! ahm * e3 * curl  (computed from 1 to jpim1/jpjm1)
!!gm open question here : e3f  at before or now ?    probably now...
!!gm note that ahmf has already been multiplied by fmask
               zcur(ji-1,jj-1) = ahmf(ji-1,jj-1,jk) * e3f_n(ji-1,jj-1,jk) * r1_e1e2f(ji-1,jj-1)       &
                  &     * (  e2v(ji  ,jj-1) * pvb(ji  ,jj-1,jk) - e2v(ji-1,jj-1) * pvb(ji-1,jj-1,jk)  &
                  &        - e1u(ji-1,jj  ) * pub(ji-1,jj  ,jk) + e1u(ji-1,jj-1) * pub(ji-1,jj-1,jk)  )
               !                                      ! ahm * div        (computed from 2 to jpi/jpj)
!!gm note that ahmt has already been multiplied by tmask
               zdiv(ji,jj)     = ahmt(ji,jj,jk) * r1_e1e2t(ji,jj) / e3t_b(ji,jj,jk)                                         &
                  &     * (  e2u(ji,jj)*e3u_b(ji,jj,jk) * pub(ji,jj,jk) - e2u(ji-1,jj)*e3u_b(ji-1,jj,jk) * pub(ji-1,jj,jk)  &
                  &        + e1v(ji,jj)*e3v_b(ji,jj,jk) * pvb(ji,jj,jk) - e1v(ji,jj-1)*e3v_b(ji,jj-1,jk) * pvb(ji,jj-1,jk)  )
            END DO  
         END DO  
         !
         DO jj = 2, jpjm1                             ! - curl( curl) + grad( div )
            DO ji = fs_2, fs_jpim1   ! vector opt.
               pua(ji,jj,jk) = pua(ji,jj,jk) + zsign * (                                                 &
                  &              - ( zcur(ji  ,jj) - zcur(ji,jj-1) ) * r1_e2u(ji,jj) / e3u_n(ji,jj,jk)   &
                  &              + ( zdiv(ji+1,jj) - zdiv(ji,jj  ) ) * r1_e1u(ji,jj)                     )
                  !
               pva(ji,jj,jk) = pva(ji,jj,jk) + zsign * (                                                 &
                  &                ( zcur(ji,jj  ) - zcur(ji-1,jj) ) * r1_e1v(ji,jj) / e3v_n(ji,jj,jk)   &
                  &              + ( zdiv(ji,jj+1) - zdiv(ji  ,jj) ) * r1_e2v(ji,jj)                     )
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      !
   END SUBROUTINE dyn_ldf_keb
   !!======================================================================
END MODULE dynldf_keb

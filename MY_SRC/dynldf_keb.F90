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

   SUBROUTINE dyn_ldf_keb( kt, pub, pvb, pua, pva, kpass )
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
      INTEGER                         , INTENT(in   ) ::   kpass      ! number of times we smooth
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pub, pvb   ! before velocity   [m/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pua, pva   ! velocity trend    [m/s2]
      REAL(wp), DIMENSION(jpi,jpj,jpk)                ::   zua, zva   
      !
      INTEGER  ::   ji, jj, jk,jp ! dummy loop indices
      REAL(wp) ::   zsign        ! local scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   zcur, zdiv
      REAL(wp) ::   rn_wgt, rn_wgt_c, rn_wgt_n ! local weights 
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
      zsign = 1.0_wp 
      !
      zua(:,:,:) = 0.
      zva(:,:,:) = 0.
      !
!!jk do this very easy. bhmt, bhmf set by ldf_dyn each time step using sgs_ke from last time step
!!jk bhmt, bhmf are negative here.
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         DO jj = 2, jpj
            DO ji = fs_2, jpi   ! vector opt.
               !                                      ! ahm * e3 * curl  (computed from 1 to jpim1/jpjm1)
               zcur(ji-1,jj-1) = bhmf(ji-1,jj-1,jk) * e3f_n(ji-1,jj-1,jk) * r1_e1e2f(ji-1,jj-1)       &
                  &     * (  e2v(ji  ,jj-1) * pvb(ji  ,jj-1,jk) - e2v(ji-1,jj-1) * pvb(ji-1,jj-1,jk)  &
                  &        - e1u(ji-1,jj  ) * pub(ji-1,jj  ,jk) + e1u(ji-1,jj-1) * pub(ji-1,jj-1,jk)  )               
               !PRINT*," JTK: bhmf, zcur: ",bhmf(ji-1,jj-1,jk),zcur(ji-1,jj-1)
               !                                      ! ahm * div        (computed from 2 to jpi/jpj)
!!gm note that ahmt has already been multiplied by tmask
               zdiv(ji,jj)     = bhmt(ji,jj,jk) * r1_e1e2t(ji,jj) / e3t_b(ji,jj,jk)                                         &
                  &     * (  e2u(ji,jj)*e3u_b(ji,jj,jk) * pub(ji,jj,jk) - e2u(ji-1,jj)*e3u_b(ji-1,jj,jk) * pub(ji-1,jj,jk)  &
                  &        + e1v(ji,jj)*e3v_b(ji,jj,jk) * pvb(ji,jj,jk) - e1v(ji,jj-1)*e3v_b(ji,jj-1,jk) * pvb(ji,jj-1,jk)  )
               !PRINT*," JTK: bhmt, zdiv: ",bhmt(ji-1,jj-1,jk), zdiv(ji,jj) 
            END DO  
         END DO  
         !
         DO jj = 2, jpjm1                             ! - curl( curl) + grad( div )
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zua(ji,jj,jk) = zsign * (                                                                 &
                  &              - ( zcur(ji  ,jj) - zcur(ji,jj-1) ) * r1_e2u(ji,jj) / e3u_n(ji,jj,jk)   &
                  &              + ( zdiv(ji+1,jj) - zdiv(ji,jj  ) ) * r1_e1u(ji,jj)                     )
                  !
               zva(ji,jj,jk) = zsign * (                                                                 &
                  &                ( zcur(ji,jj  ) - zcur(ji-1,jj) ) * r1_e1v(ji,jj) / e3v_n(ji,jj,jk)   &
                  &              + ( zdiv(ji,jj+1) - zdiv(ji  ,jj) ) * r1_e2v(ji,jj)                     )
               
               !IF ( abs(zua(ji,jj,jk)) > 1e-3 .or. isnan(zua(ji,jj,jk)) ) THEN
               IF (zua(ji,jj,jk) /= zua(ji,jj,jk)) THEN
                   PRINT*," Warning: zua is not valid: ",zua(ji,jj,jk),zcur(ji,jj),zdiv(ji,jj)
               ENDIF
               !IF ( abs(zva(ji,jj,jk)) > 1e-3 .or. isnan(zva(ji,jj,jk)) ) THEN
               IF (zva(ji,jj,jk) /= zva(ji,jj,jk)) THEN
      	           PRINT*," Warning: zva is not valid: ",zva(ji,jj,jk),zcur(ji,jj),zdiv(ji,jj)
               ENDIF
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      !
      ! Smooth the backscatter tendency
!!jk: this needs to be energy conservative      
      IF ( kpass > 0 ) THEN
         DO jp = 1, kpass 
            !
            rn_wgt_c = 1._wp
            rn_wgt_n = 1._wp
            rn_wgt   = 1._wp / ( 4._wp * rn_wgt_n + rn_wgt_c )
            !
            zua(2:jpim1, 2:jpjm1, 1:jpkm1) = rn_wgt * ( rn_wgt_n * zua(2:jpim1,   1:jpjm1-1, 1:jpkm1) & 
                                                      + rn_wgt_n * zua(2:jpim1,   3:jpj,     1:jpkm1) & 
                                                      + rn_wgt_n * zua(1:jpim1-1, 2:jpjm1,   1:jpkm1) & 
                                                      + rn_wgt_n * zua(3:jpi,     2:jpjm1,   1:jpkm1) & 
                                                      + rn_wgt_c * zua(2:jpim1,   2:jpjm1,   1:jpkm1) )
            !
            zva(2:jpim1, 2:jpjm1, 1:jpkm1) = rn_wgt * ( rn_wgt_n * zva(2:jpim1,   1:jpjm1-1, 1:jpkm1) &
                                                      + rn_wgt_n * zva(2:jpim1,   3:jpj,     1:jpkm1) &
                                                      + rn_wgt_n * zva(1:jpim1-1, 2:jpjm1,   1:jpkm1) & 
                                                      + rn_wgt_n * zva(3:jpi,     2:jpjm1,   1:jpkm1) & 
                                                      + rn_wgt_c * zva(2:jpim1,   2:jpjm1,   1:jpkm1) ) 
            !
            DO jk = 1, jpkm1 
               DO jj = 2, jpj
                  DO ji = fs_2, jpi
                     !IF ( abs(zua(ji,jj,jk)) > 1e-3 .or. isnan(zua(ji,jj,jk)) ) THEN
                     !    PRINT*," JTK: zua, pass: ",zua(ji,jj,jk), jp
                     !ENDIF
                     !IF ( abs(zva(ji,jj,jk)) > 1e-3 .or. isnan(zva(ji,jj,jk)) ) THEN
                     !    PRINT*," JTK: zva, pass: ",zva(ji,jj,jk), jp
                     !ENDIF
                  END DO
               END DO
            END DO
            !
            !jk: not sure about north fold sign
            CALL lbc_lnk_multi( 'dynldf_keb', zua, 'U', -1.,  zva, 'V', -1. )
         END DO
         !      
      ENDIF
      !
      ! Add tendency to pua, pva
      pua(:,:,:) = pua(:,:,:) + zua(:,:,:)
      pva(:,:,:) = pva(:,:,:) + zva(:,:,:)
      PRINT*," JTK: max zua, zva ",MAXVAL(zua),MAXVAL(zva)
      !
   END SUBROUTINE dyn_ldf_keb
   !!======================================================================
END MODULE dynldf_keb

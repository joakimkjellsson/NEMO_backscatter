MODULE dynldf
   !!======================================================================
   !!                       ***  MODULE  dynldf  ***
   !! Ocean physics:  lateral diffusivity trends 
   !!=====================================================================
   !! History :  2.0  ! 2005-11  (G. Madec)  Original code (new step architecture)
   !!            3.7  ! 2014-01  (F. Lemarie, G. Madec)  restructuration/simplification of ahm specification,
   !!                 !                                  add velocity dependent coefficient and optional read in file
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_ldf      : update the dynamics trend with the lateral diffusion
   !!   dyn_ldf_init : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE ldfdyn         ! lateral diffusion: eddy viscosity coef.
   USE dynldf_lap_blp ! lateral mixing   (dyn_ldf_lap & dyn_ldf_blp routines)
   USE dynldf_iso     ! lateral mixing                 (dyn_ldf_iso routine )
   USE dynldf_keb     ! KE backscatter  (dyn_ldf_keb routine)
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics   (trd_dyn      routine)
   !
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE iom 
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_ldf       ! called by step module 
   PUBLIC   dyn_ldf_init  ! called by opa  module 
   PUBLIC   dyn_ldf_alloc ! called by nemogcm.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dynldf.F90 10068 2018-08-28 14:09:04Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
   
   INTEGER FUNCTION dyn_ldf_alloc()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_ldf_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( sgs_ke(jpi,jpj,jpk), STAT=dyn_ldf_alloc )
         !
      IF( dyn_ldf_alloc /= 0 )   CALL ctl_warn('dyn_ldf_alloc: array allocate failed.')
   END FUNCTION dyn_ldf_alloc
   
   SUBROUTINE dyn_ldf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_ldf  ***
      !! 
      !! ** Purpose :   compute the lateral ocean dynamics physics.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER             ::   ji,jj,jk ! loop indices (Joakim)
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ztrdu, ztrdv
      ! Joakim added these
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   bu, bv   ! volume of u- and v-boxes
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   r1_bt    ! inverse of t-box volume
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zke      ! temporary KE array
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_ldf')
      !
      ! Allocate sub-grid scale energy if first step (Joakim)
      !IF( kt == nit000 .AND. ln_sgske ) THEN
      !   IF(lwp) WRITE(numout,*)
      !   IF(lwp) WRITE(numout,*) 'dyn_ldf : use sub-grid scale kinetic energy reservoir '
      !   IF( dyn_ldf_alloc() /= 0 )   CALL ctl_stop('STOP', 'dyn_ldf: failed to allocate arrays')
      !ENDIF
      !
      IF( l_trddyn .OR. ln_sgske )   THEN                      ! temporary save of momentum trends
         ALLOCATE( ztrdu(jpi,jpj,jpk) , ztrdv(jpi,jpj,jpk) )
         ztrdu(:,:,:) = ua(:,:,:) 
         ztrdv(:,:,:) = va(:,:,:) 
         IF ( ln_sgske ) THEN
            ALLOCATE( bu(jpi,jpj,jpk), bv(jpi,jpj,jpk), r1_bt(jpi,jpj,jpk), zke(jpi,jpj,jpk) )    ! inverse of t-box volume
         END IF
      ENDIF
      !
      SELECT CASE ( nldf_dyn )                   ! compute lateral mixing trend and add it to the general trend
      !
      CASE ( np_lap   )    ;   CALL dyn_ldf_lap( kt, ub, vb, ua, va, 1 )      ! iso-level    laplacian
      CASE ( np_lap_i )    ;   CALL dyn_ldf_iso( kt )                         ! rotated      laplacian
      CASE ( np_blp   )    ;   CALL dyn_ldf_blp( kt, ub, vb, ua, va    )      ! iso-level bi-laplacian
      !
      END SELECT     
      !
      !IF ( ln_kebs ) CALL dyn_ldf_keb( kt, ub, vb, ua, va ) ! KE backscatter 
      !
      IF( l_trddyn .OR. ln_sgske ) THEN                        ! save the horizontal diffusive trends for further diagnostics
!!jk: this is the contribution from viscosity
         ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
         ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
         CALL lbc_lnk_multi( 'dynldf', ztrdu, 'U', -1. , ztrdv, 'V', -1. )
         IF ( l_trddyn ) THEN
            CALL trd_dyn( ztrdu, ztrdv, jpdyn_ldf, kt )
         END IF
         ! 
	 IF ( ln_sgske ) THEN 	    
!!jk: Calculate dissipated KE. Taken from trdken routine
            DO jk = 1, jpkm1
               bu   (:,:,jk) =    e1e2u(:,:) * e3u_n(:,:,jk)
               bv   (:,:,jk) =    e1e2v(:,:) * e3v_n(:,:,jk)
               r1_bt(:,:,jk) = r1_e1e2t(:,:) / e3t_n(:,:,jk) * tmask(:,:,jk)
            END DO
            !
            zke(:,:,:) = 0._wp
            DO jk = 1, jpkm1
               DO jj = 2, jpj
                  DO ji = 2, jpi
!!jk: KE tendency from ldf dissipation [m2/s3]
                     zke(ji,jj,jk) = 0.5_wp * ( un(ji  ,jj,jk) * ztrdu(ji  ,jj,jk) * bu(ji  ,jj,jk)  &
                     &                        + un(ji-1,jj,jk) * ztrdu(ji-1,jj,jk) * bu(ji-1,jj,jk)  &
                     &                        + vn(ji,jj  ,jk) * ztrdv(ji,jj  ,jk) * bv(ji,jj  ,jk)  &
                     &                        + vn(ji,jj-1,jk) * ztrdv(ji,jj-1,jk) * bv(ji,jj-1,jk)  ) * r1_bt(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
!!jk: multiply KE tendency by rn_rdt. Not 2*rn_rdt
            sgs_ke(:,:,:) = sgs_ke(:,:,:) - zke(:,:,:) * rn_rdt  !! [m2/s2]
            CALL lbc_lnk_multi( 'dynldf', sgs_ke, 'T', -1.) 
            !
            IF ( ln_kebs ) THEN               
               ztrdu(:,:,:) = ua(:,:,:) 
               ztrdv(:,:,:) = va(:,:,:) 
               !
               CALL dyn_ldf_keb( kt, ub, vb, ua, va, 4 ) ! KE backscatter
               !               
!!jk: this is the contribution from backscatter
               ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
               ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
               CALL lbc_lnk_multi( 'dynldf', ztrdu, 'U', -1. , ztrdv, 'V', -1. )              
               !
               !DO jk = 1, jpkm1
               !   DO jj = 2, jpj
               !       DO ji = 2, jpi
               !          IF ( abs(ua(ji,jj,jk)) > 1e-3 .or. isnan(ua(ji,jj,jk)) ) THEN
               !             PRINT*," JTK: ua: ",ua(ji,jj,jk)
               !          ENDIF
               !          IF ( abs(va(ji,jj,jk)) > 1e-3 .or. isnan(va(ji,jj,jk)) ) THEN
               !             PRINT*," JTK: va: ",va(ji,jj,jk)
               !          ENDIF
               !       END DO
               !   END DO
               !END DO
               !
               IF ( l_trddyn ) THEN
                  CALL trd_dyn( ztrdu, ztrdv, jpdyn_keb, kt )
               END IF
               !
               zke(:,:,:) = 0._wp
               DO jk = 1, jpkm1
                  DO jj = 2, jpj
                     DO ji = 2, jpi
                        zke(ji,jj,jk) = 0.5_wp * ( un(ji  ,jj,jk) * ztrdu(ji  ,jj,jk) * bu(ji  ,jj,jk)  &
                        &                        + un(ji-1,jj,jk) * ztrdu(ji-1,jj,jk) * bu(ji-1,jj,jk)  &
                        &                        + vn(ji,jj  ,jk) * ztrdv(ji,jj  ,jk) * bv(ji,jj  ,jk)  &
                        &                        + vn(ji,jj-1,jk) * ztrdv(ji,jj-1,jk) * bv(ji,jj-1,jk)  ) * r1_bt(ji,jj,jk)
                     END DO
                  END DO
               END DO
!!jk: Now we take the KE from the sub-grid scale budget
!!jk: Again multiply by rn_rdt
               sgs_ke(:,:,:) = sgs_ke(:,:,:) - zke(:,:,:) * rn_rdt
               CALL lbc_lnk_multi( 'dynldf', sgs_ke, 'T', -1.)
            ENDIF
	    !  
         END IF
         !IF ( l_trddyn ) THEN
         !   CALL trd_dyn( ztrdu, ztrdv, jpdyn_ldf, kt )
         !END IF
         DEALLOCATE ( ztrdu , ztrdv )         
      ENDIF
      !                                          ! print sum trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' ldf  - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( ln_timing )   CALL timing_stop('dyn_ldf')
      !
   END SUBROUTINE dyn_ldf


   SUBROUTINE dyn_ldf_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_ldf_init  ***
      !! 
      !! ** Purpose :   initializations of the horizontal ocean dynamics physics
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                     !==  Namelist print  ==!
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_ldf_init : Choice of the lateral diffusive operator on dynamics'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_ldf: already read in ldfdyn module'
         WRITE(numout,*) '      see ldf_dyn_init report for lateral mixing parameters'
         WRITE(numout,*)
         !
         SELECT CASE( nldf_dyn )             ! print the choice of operator
         CASE( np_no_ldf )   ;   WRITE(numout,*) '   ==>>>   NO lateral viscosity'
         CASE( np_lap    )   ;   WRITE(numout,*) '   ==>>>   iso-level laplacian operator'
         CASE( np_lap_i  )   ;   WRITE(numout,*) '   ==>>>   rotated laplacian operator with iso-level background'
         CASE( np_blp    )   ;   WRITE(numout,*) '   ==>>>   iso-level bi-laplacian operator'
         END SELECT
      ENDIF
      !
   END SUBROUTINE dyn_ldf_init

   !!======================================================================
END MODULE dynldf

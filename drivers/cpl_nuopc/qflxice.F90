!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module qflxice

#if defined CCSMCOUPLED
!BOP
! !MODULE: ice
!
! !DESCRIPTION:
!  This module currently contains routines for computing sea ice 
!  formation / potential melting and the associated heat flux.
!  This heat flux is sent to the ice model via the flux coupler.
!
! !REVISION HISTORY:

! !USES:

   use par_kind
   use par_oce
   use dom_oce
   use oce
   use phycst
   use in_out_manager
   use iom
   use restart
   use lib_fortran
   use shr_frz_mod,        only: shr_frz_freezetemp_init, shr_frz_freezetemp

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_qflxice,          &
             increment_tlast_ice,   &
             ice_formation,         &
             ice_flx_to_coupler
             

! !PUBLIC DATA MEMBERS:

   logical, public :: &
      liceform,       &   ! flag to turn on/off ice formation
      lice_form_ts,   &   ! T ==> ice formation/melting computation time step
      lice_cpl_ts         ! T ==> ice flux to coupler time step

   logical :: &
      lactive_ice         ! T ==> ocn is coupled to an active ice model
                          ! F ==> ocn is coupled to a dummy ice model

   logical, parameter :: &
      lfw_as_salt_flx = .false.   ! treat fw flux as virtual salt flux
                                  ! even with var.thickness sfc layer

!   integer, parameter, public :: &
!      nn_nits = 2         ! ice formation/melting computed starting at nn_nits-1
!                          ! time steps before the coupling time step
!                          ! 1 ==> computed at the coupling time step only
!                          ! 2 ==> comp. at the coup. ts and 1 ts before
!                          ! n ==> comp. at the coup. ts and n-1 ts before

   real (wp), dimension(:,:), allocatable, public :: &
      QFLUX               ! internal ocn heat flux due to ice formation

   real (wp), dimension(:,:), allocatable, public :: &
      AQICE,             &! sum of accumulated ice heat flux since tlast
      QICE!,              &! tot column cooling from ice form (in C*m)
!      SALT_FREEZE         ! salt flux at T points due to frazil ice formation

!   real (wp), dimension(:,:), allocatable, public ::  &
!      FW_FREEZE,    &! water flux at T points due to frazil ice formation
!      SFLUX          ! salt flux at T points due to frazil ice formation

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

   real (wp), public ::  &
      tlast_ice           ! time since last ice flux computed

!   real (wp) ::  &
!      cp_over_lhfusion    ! rcp/lfus

!   real (wp) ::          &
!      hflux_factor

   ! FIXME: should be a namelist variable !
   integer :: &
      kmxice = 1          ! lowest level from which to integrate 
                          ! ice formation (1 => surface)

#include "domzgr_substitute.h90"

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE:
! !INTERFACE:

 subroutine init_qflxice(tfrz_option)

  character(len=lc), INTENT(IN)   :: tfrz_option     ! tfrz_option from driver
! !DESCRIPTION:
!  This routine initializes ice formation/melting related variables.
!  It must be called before initializing restarts because this module
!  add the accumulated ice heat flux to the restart file.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!   integer ::   &
!      nml_error             ! namelist i/o error flag

!   namelist /ice_nml/ kmxice, lactive_ice

!-----------------------------------------------------------------------
!
!  read input namelists
!
!-----------------------------------------------------------------------

!   cp_over_lhfusion = rho0_rcp/(rLfus*raufw)
!   hflux_factor     = r1_rho0_rcp

   CALL shr_frz_freezetemp_init(tfrz_option, lwp)

   kmxice           = 1
   lactive_ice      = .true.

   liceform     = .true.

   lice_form_ts = .false.
   lice_cpl_ts  = .false.

   if (liceform .and. lwp) then
      write(numout,'(a20,1pe10.3)') 'Ice salinity(PSU) = ', sice
      write(numout,'(a30,i3,a13)') 'Ice formation computed in top ', &
                kmxice, ' levels only.'
      write(numout,'(a)') 'tfreeze_option from driver = '//TRIM(tfrz_option)
   endif

   tlast_ice = 0.0_wp

   !***
   !*** allocate and initialize ice flux arrays
   !***

   allocate( QICE(jpi,jpj), & 
             AQICE(jpi,jpj), &
             QFLUX(jpi,jpj))!, &
!             FW_FREEZE(jpi,jpj),   &
!             SALT_FREEZE(jpi,jpj), &
!             SFLUX(jpi,jpj))

   QICE(:,:)  = 0.0_wp
   AQICE(:,:) = 0.0_wp
   QFLUX(:,:) = 0.0_wp

!   FW_FREEZE(:,:)   = 0.0_wp
!   SALT_FREEZE(:,:) = 0.0_wp
!   SFLUX(:,:) = 0.0_wp

   if (ln_rstart) then
      IF(lwp) WRITE(numout,*) ' nit000-1 sea ice freezing/melting potential fields read in restart file'
      call iom_get( numror, jpdom_auto, 'AQICE', AQICE(:,:) )
      call iom_get( numror, jpdom_auto, 'QFLUX', QFLUX(:,:) )
!!      IF (iom_varid( numror, 'SFLUX', ldstop = .FALSE. ) > 0) &
!        call iom_get( numror, jpdom_auto, 'SFLUX', SFLUX(:,:) )
!!      IF (iom_varid( numror, 'SALT_FREEZE', ldstop = .FALSE. ) > 0) &
!        call iom_get( numror, jpdom_auto, 'SALT_FREEZE', SALT_FREEZE(:,:) )
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine init_qflxice

!***********************************************************************
!BOP
! !IROUTINE:
! !INTERFACE:
   subroutine increment_tlast_ice(kt)

   integer, intent(in) :: kt

! !DESCRIPTION:
!  This subroutine increments tlast_ice in a nonthreaded region.

! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  increment time since last evaluation
!
!-----------------------------------------------------------------------

   tlast_ice = tlast_ice + rdt

!-----------------------------------------------------------------------
!EOC

 end subroutine increment_tlast_ice

!***********************************************************************
!BOP
! !IROUTINE:
! !INTERFACE:
   subroutine ice_formation( kt, Ktt )

! !DESCRIPTION:
!  This subroutine computes ocean heat flux to the sea-ice. it forms
!  the necessary ice in the ocean and adjusts the potential 
!  temperature and salinity fields accordingly. the logic of this 
!  subroutine is based on William Large''s 1-d model and is based
!  on a version from the NCOM model written by Gokhan Danabasoglu. 

! !REVISION HISTORY:
!  same as module
   implicit none
 
! !INPUT/OUTPUT PARAMETERS:
   integer, intent(in), optional :: kt   ! Time step
   integer, intent(in) :: Ktt            ! Time level
!EOP
!BOC


!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   integer :: &
      k,                 &! vertical level index
      n                   ! tracer index


   real (wp), dimension(jpi,jpj) :: &
     POTICE,           & ! potential amt of ice formation
     WORK1,            & ! work array
     WORK2,            & ! work array for salt flux
     TFRZ                ! freezing temp of water

   real (wp) ::        &
     r_rdttra

   real (wp) ::        &
     ref_val             ! tracer reference value


!========================
   if ( lice_form_ts ) then
!========================

!-----------------------------------------------------------------------
!
!  initialize flux to zero
!
!-----------------------------------------------------------------------
 
     QICE(:,:)   = 0.0_wp
     POTICE(:,:) = 0.0_wp
 
!-----------------------------------------------------------------------
!
!  compute frazil ice formation for sub-surface layers. if ice
!  forms in lower layers but layers above are warm - the heat is
!  used to melt the ice. the ice formation occurs at salinity, Si.
!  this volume is replaced with an equal volume at the salinity of
!  the layer above. the total ice heat flux is accumulated. 
!
!  WARNING: unless a monotone advection scheme is in place, 
!  advective errors could lead to temps that are far below freezing
!  in some locations and this scheme will form lots of ice.
!  ice formation should be limited to the top layer (kmxice=1)
!  if the advection scheme is not monotone.
!
!-----------------------------------------------------------------------
 
     do k=kmxice,2,-1
!
!     !***
!     !*** potice is the potential amount of ice formation 
!     !*** (potice>0) or melting (potice<0) in layer k
!     !***
!
       call tfreez(TFRZ(:,:),ts(:,:,k,jp_sal,Ktt))
       POTICE(:,:) = (TFRZ(:,:) - ts(:,:,k,jp_tem,Ktt))*e3t(:,:,k,Ktt)*tmask(:,:,k)
!
!     !***
!     !*** if potice < 0, use the heat to melt any ice
!     !*** from lower layers
!     !*** if potice > 0, keep on freezing (QICE < 0)
!     !***
!
       POTICE(:,:) = max(POTICE(:,:),QICE(:,:))
!
!     !***
!     !*** adjust tracer values based on freeze/melt
!     !***
!
       where (POTICE(:,:)>0.0_dp)
         ts(:,:,k,jp_tem,Ktt) = TFRZ(:,:)
       endwhere
!
!       if (lk_vvl .and. .not. lfw_as_salt_flx) then
!            tsa(:,:,k,jp_sal) = ( tsa(:,:,k,jp_sal)                               &
!            * (fse3t(:,:,k) + cp_over_lhfusion * QICE(:,:))         &
!            + cp_over_lhfusion * (tsa(:,:,k-1,jp_sal)                         &
!            * (POTICE(:,:) - QICE(:,:)) - sice * POTICE(:,:)) )/fse3t(:,:,k)
!       else
!            ref_val = soce - sice
!            if (ref_val /= 0.0_wp)  then
!               tsa(:,:,k,jp_sal) = tsa(:,:,k,jp_sal) &
!               + ref_val*POTICE(:,:)*cp_over_lhfusion/fse3t(:,:,k)
!            endif
!       endif
!
!       !*** accumulate freezing potential
       QICE(:,:) = QICE(:,:) - POTICE(:,:)
!
     enddo ! k loop
!
!-----------------------------------------------------------------------
!
!  now repeat the above algorithm for the surface layer. when fresh
!  water flux formulation is used, the surface layer does not get
!  any salt from other layers. instead, its volume changes. 
!
!-----------------------------------------------------------------------

     k = 1
     
     call tfreez(TFRZ(:,:),ts(:,:,k,jp_sal,Ktt))

     WORK1(:,:) = e3t(:,:,k,Ktt)

!     if (.not. lk_vvl)  &
!       WORK1 = WORK1 + ssha(:,:)

     POTICE(:,:) = (TFRZ(:,:) - ts(:,:,k,jp_tem,Ktt))*WORK1(:,:)*tmask(:,:,k)

     POTICE(:,:) = max(POTICE(:,:), QICE(:,:))

     where (POTICE(:,:)>0.0_dp)
       ts(:,:,k,jp_tem,Ktt) = TFRZ(:,:)
     endwhere

!     if (lk_vvl .and. .not. lfw_as_salt_flx) then
!        tsa(:,:,k,jp_sal) =  &
!           (tsa(:,:,k,jp_sal)*(WORK1(:,:) + cp_over_lhfusion*QICE(:,:)) - &
!           sice*QICE(:,:)*cp_over_lhfusion )/WORK1(:,:)
!     else
!        ref_val = soce - sice
!        if (ref_val /= 0.0_wp)  &
!           WORK2(:,:) = ref_val*POTICE(:,:)*cp_over_lhfusion/WORK1(:,:)
!!           tsa(:,:,k,jp_sal) = tsa(:,:,k,jp_sal) + WORK2(:,:)
!     endif

     QICE(:,:) = QICE(:,:) - POTICE(:,:)
!     SALT_FREEZE(:,:) = SALT_FREEZE(:,:) + WORK2(:,:)

!-----------------------------------------------------------------------
!
!  let any residual heat in the upper layer melt previously formed ice
!
!-----------------------------------------------------------------------
 
     AQICE(:,:) = AQICE(:,:) + QICE(:,:)

!-----------------------------------------------------------------------
!
!  recalculate freezing potential based on adjusted T.
!  only interested in melt potential now (POTICE < 0) - use this 
!  melt to offset any accumulated freezing (AQICE < 0) and
!  adjust T and S to reflect this melting. when freshwater flux
!  formulation, compute the associated freshwater flux instead of
!  adjusting S.
!
!-----------------------------------------------------------------------

!     call tfreez(TFRZ(:,:),tsa(:,:,k,jp_sal))
!
!     where (k <= mbkt(:,:))
!       POTICE(:,:) = (TFRZ(:,:) - tsa(:,:,k,jp_tem)) * WORK1(:,:) * tmask(:,:,k)
!     endwhere
!
!     POTICE(:,:) = max(POTICE(:,:), AQICE(:,:))
!
!     tsa(:,:,k,jp_tem) = tsa(:,:,k,jp_tem) + POTICE(:,:)/WORK1(:,:)
!
!!     if (lk_vvl .and. .not. lfw_as_salt_flx) then
!!        FW_FREEZE(:,:) = min(POTICE(:,:),QICE(:,:)) &
!!                  * cp_over_lhfusion / rdttra(k)
!!     else
!!        ref_val = soce - sice
!!        if (ref_val /= 0.0_wp) &
!!           WORK2(:,:) = ref_val*POTICE(:,:)*cp_over_lhfusion/WORK1(:,:)
!!!           tsa(:,:,k,jp_sal) = tsa(:,:,k,jp_sal) + WORK2(:,:)
!!     endif
!!
!!     AQICE(:,:) = AQICE(:,:) - POTICE(:,:)
!!!     SALT_FREEZE(:,:) = SALT_FREEZE(:,:) + WORK2(:,:)
!!

   endif ! time to do ice

   if (lrst_oce) then
      IF(lwp) WRITE(numout,*) 'qflxice : AQICE written in ocean restart file ',   &
          &                    'at it= ', kt,' date= ', ndastp
      IF(lwp) WRITE(numout,*) '~~~~'
      CALL iom_rstput( kt, nitrst, numrow, 'AQICE', AQICE(:,:) )
!      CALL iom_rstput( kt, nitrst, numrow, 'SALT_FREEZE', SALT_FREEZE(:,:) )
   endif


!-----------------------------------------------------------------------
!EOC

   end subroutine ice_formation

!***********************************************************************

   subroutine ice_flx_to_coupler( kt, Ktt )

!-----------------------------------------------------------------------
!
!  This subroutine sets up the ice formation / melting potential
!  heat fluxes to be sent to the coupler. ice formation heat flux
!  is accumulated for time averaging.
!
!-----------------------------------------------------------------------

   integer, intent(in) :: kt   ! Time step
   integer, intent(in) :: Ktt  ! Time level

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (wp), dimension(jpi,jpj) :: & 
      WORK1, WORK2         ! work arrays

   real (wp), dimension(jpi,jpj) :: &
      TFRZ                 ! freezing temp of water

!-----------------------------------------------------------------------

   if (lice_form_ts .AND. lice_cpl_ts) then
!-----------------------------------------------------------------------
!
!  compute the first layer thickness
!
!-----------------------------------------------------------------------

   call tfreez(TFRZ(:,:),ts(:,:,1,jp_sal,Ktt))
!   call tfreez(TFRZ(:,:),sn(:,:,1))

   WORK1(:,:) = e3t(:,:,1,Ktt)

!   if ( .not. lk_vvl ) &
!     WORK1 = WORK1 + sshn(:,:)

!-----------------------------------------------------------------------
!
!  first compute the melt potential
!
!-----------------------------------------------------------------------

   WORK2(:,:) = 0.0_wp
   WORK2(:,:) = (TFRZ(:,:) - ts(:,:,1,jp_tem,Ktt)) * WORK1(:,:) * tmask(:,:,1)

!-----------------------------------------------------------------------
!
!  adjust ice formation amount
!
!-----------------------------------------------------------------------
   AQICE(:,:) = AQICE(:,:)*0.5_wp  ! Adjust for leap-frog time stepping

!-----------------------------------------------------------------------
!
!  merge the ice formation and melt potential fluxes
!
!-----------------------------------------------------------------------
!   if (tlast_ice == 0.0_wp) then
!     SFLUX(:,:) = 0.0_wp
!   else
!     SFLUX(:,:) = SALT_FREEZE(:,:)*rho0*WORK1(:,:)*tmask(:,:,1)/tlast_ice
!   endif

   WORK1(:,:) = 0.0_wp
   where ( AQICE(:,:) < 0.0_wp ) 
     WORK1(:,:) = -AQICE(:,:)
   elsewhere
     WORK1(:,:) = WORK2(:,:)
   endwhere

   if (tlast_ice == 0.0_wp) then
     QFLUX(:,:) = 0.0_wp
   else
     QFLUX(:,:) = WORK1(:,:)*tmask(:,:,1)*rho0_rcp/tlast_ice
   endif

   lice_form_ts = .false.
   lice_cpl_ts  = .false.

   endif

   if (lrst_oce) then
      IF(lwp) WRITE(numout,*) 'qflxice : QFLUX written in ocean restart file ',   &
          &                    'at it= ', kt,' date= ', ndastp
      IF(lwp) WRITE(numout,*) '~~~~'
      CALL iom_rstput( kt, nitrst, numrow, 'QFLUX' , QFLUX(:,:) )
!      CALL iom_rstput( kt, nitrst, numrow, 'SFLUX' , SFLUX(:,:) )
   endif

   end subroutine ice_flx_to_coupler


!***********************************************************************
!BOP
! !IROUTINE:
! !INTERFACE:

 subroutine tfreez(TFRZ,SALT)

! !DESCRIPTION:
!  This function computes the freezing point of salt water.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (wp), dimension(jpi,jpj), intent(in) :: &
      SALT                ! salinity in model units (g/g)

! !OUTPUT PARAMETERS:

   real (wp), dimension(jpi,jpj), intent(out) :: &
      TFRZ                ! freezing temperature of water in deg C

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  call shr function to return freezing temp based on drv namelist
!  choice of minus1p8, linear_salt, or mushy algorithms.
!
!-----------------------------------------------------------------------

   TFRZ(:,:) = shr_frz_freezetemp(SALT(:,:))

!-----------------------------------------------------------------------
!EOC

 end subroutine tfreez


!BOP
! !IROUTINE:
! !INTERFACE:

   subroutine tmelt (TMLT,SALT)

! !DESCRIPTION:
!  This subroutine sets the melting point temperature of ice.
!  For now, TMLT is a separate routine than TFRZ.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (wp), dimension(jpi,jpj), intent(in) ::  &
     SALT                ! salinity in model units (PSU)

! !OUTPUT PARAMETERS:

   real (wp), dimension(jpi,jpj), intent(out) :: &
     TMLT                ! melting temperature in deg C
!EOP
!BOC

   if ( lactive_ice ) then
     TMLT(:,:) = 0.0_wp
   else
     call tfreez(TMLT(:,:),SALT(:,:))
   endif

!-----------------------------------------------------------------------
!EOC

   end subroutine tmelt


!***********************************************************************
#endif

end module qflxice

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

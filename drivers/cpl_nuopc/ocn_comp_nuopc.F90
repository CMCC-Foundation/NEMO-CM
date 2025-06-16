module ocn_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for POP
  !----------------------------------------------------------------------------
!$ use omp_lib, only : omp_set_num_threads
  use ESMF
  use NUOPC                 , only : NUOPC_CompDerive,NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC                 , only : NUOPC_CompFilterPhaseMap, NUOPC_IsUpdated, NUOPC_IsAtTime
  use NUOPC                 , only : NUOPC_CompAttributeGet, NUOPC_Advertise, NUOPC_CompSetClock
  use NUOPC                 , only : NUOPC_SetAttribute, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
  use NUOPC_Model           , only : model_routine_SS           => SetServices
  use NUOPC_Model           , only : SetVM
  use NUOPC_Model           , only : model_label_Advance        => label_Advance
  use NUOPC_Model           , only : model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model           , only : model_label_SetRunClock    => label_SetRunClock
  use NUOPC_Model           , only : model_label_CheckImport    => label_CheckImport
  use NUOPC_Model           , only : model_label_SetClock       => label_SetClock
  use NUOPC_Model           , only : model_label_Finalize       => label_Finalize
  use NUOPC_Model           , only : NUOPC_ModelGet
  use shr_file_mod          , only : shr_file_getLogUnit, shr_file_setLogUnit
  use shr_cal_mod           , only : shr_cal_date2ymd, shr_cal_ymd2date
  use shr_sys_mod           , only : shr_sys_flush, shr_sys_abort
  use shr_kind_mod          , only : cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_const_mod         , only : shr_const_spval
  use shr_string_mod        , only : shr_string_listGetNum, shr_string_listGetName
  use shr_mpi_mod           , only : shr_mpi_min, shr_mpi_max
  use ocn_communicator      , only : mpi_communicator_ocn
!  use MCT_vars_mod          , only : nemo_mct_init
  use perf_mod              , only : t_startf, t_stopf
!  use ocn_import_export     , only : ocn_advertise_fields, ocn_realize_fields
!  use ocn_import_export     , only : ocn_import, ocn_export, ocn_sum_buffer, tlast_coupled
!  use ocn_import_export     
  use nuopc_shr_methods     , only : chkerr, state_setscalar, state_getscalar, state_diagnose, alarmInit
  use nuopc_shr_methods     , only : set_component_logging, get_component_instance, log_clock_advance
  use nemogcm               , only : cform_aaa, nemo_init, nemo_closefile
  use domzgr                , only : dom_zgr
  use par_kind
  use par_oce
  use dom_oce
  use oce
  use lib_mpp
  use iom
  use lbclnk
  use in_out_manager
  use qflxice
  use geo2ocean
  use phycst
  use lib_fortran
  use sbccpl_cesm
  use eosbn2
  use timing
#if defined key_qco   ||   defined key_linssh
   USE stpmlf         ! NEMO time-stepping               (stp_MLF   routine)
#else
   USE step           ! NEMO time-stepping                 (stp     routine)
#endif


  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private                              ! By default make data private

  public  :: SetServices
  public  :: SetVM
  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelSetRunClock
  private :: ModelFinalize
  private :: lice_form

  private  :: ocn_advertise_fields
  private  :: ocn_realize_fields
  private  :: ocn_import
  private  :: ocn_export
  private  :: ocn_sum_buffer

  private :: fldlist_add
  private :: fldlist_realize
  private :: state_FldChk

  character(len=CL)   :: flds_scalar_name = ''
  integer             :: flds_scalar_num = 0
  integer             :: flds_scalar_index_nx = 0
  integer             :: flds_scalar_index_ny = 0
  integer             :: flds_scalar_index_precip_factor = 0

  character(len=lc)   :: tfrz_option_driver    ! tfrz_option from driver

  integer, parameter  :: dbug = 1
  logical             :: ldiag_cpl = .false.
  character(len=lc) :: runtype
  character(len=lc) :: message

  character(*), parameter :: u_FILE_u = &
       __FILE__

  integer (i4)  ::   &
      istp,           & ! time step counter (from nit000 to nitend)
      nn_ncpl,        & ! # of model time steps in 1 coupling time step
      ndt05             ! half-time step
  integer :: nproc
  real(wp), parameter ::  c0 = 0.0_wp
  real(wp), parameter ::  c1 = 1.0_wp

  type fld_list_type
    character(len=128) :: stdname
    integer :: ungridded_lbound = 0
    integer :: ungridded_ubound = 0
  end type fld_list_type

  integer, parameter       :: fldsMax = 100
  integer                  :: fldsToOcn_num = 0
  integer                  :: fldsFrOcn_num = 0
  type (fld_list_type)     :: fldsToOcn(fldsMax)
  type (fld_list_type)     :: fldsFrOcn(fldsMax)

  ! area correction factors for fluxes send and received from mediator
  real(wp), allocatable :: mod2med_areacor(:) ! ratios of model areas to input mesh areas
  real(wp), allocatable :: med2mod_areacor(:) ! ratios of input mesh areas to model areas

  interface state_getfldptr
     module procedure state_getfldptr_1d
     module procedure state_getfldptr_2d
  end interface state_getfldptr

  logical :: ocn2glc_coupling
  integer :: num_ocn2glc_levels
  integer, allocatable :: ocn2glc_levels(:)

  ! accumulated sum of send buffer quantities for averaging before being sent
  real(wp), allocatable :: sbuff_sum_u    (:,:) !(nx_block,ny_block,max_blocks_clinic)
  real(wp), allocatable :: sbuff_sum_v    (:,:) !(nx_block,ny_block,max_blocks_clinic)
  real(wp), allocatable :: sbuff_sum_t    (:,:) !(nx_block,ny_block,max_blocks_clinic)
  real(wp), allocatable :: sbuff_sum_s    (:,:) !(nx_block,ny_block,max_blocks_clinic)
  real(wp), allocatable :: sbuff_sum_dhdx (:,:) !(nx_block,ny_block,max_blocks_clinic)
  real(wp), allocatable :: sbuff_sum_dhdy (:,:) !(nx_block,ny_block,max_blocks_clinic)
!  real(r8) :: sbuff_sum_bld  (jpi,jpj) !(nx_block,ny_block,max_blocks_clinic)
  real(wp), allocatable :: sbuff_sum_co2  (:,:) !(nx_block,ny_block,max_blocks_clinic)
  !real(wp), allocatable :: sbuff_sum_t_depth (:,:,:)
  !real(wp), allocatable :: sbuff_sum_s_depth (:,:,:)

  ! tlast_coupled is incremented by delt every time pop_sum_buffer is called
  ! tlast_coupled is reset to 0 when ocn_export is called
  real (wp) :: tlast_coupled

  real (wp), save ::  area

!=======================================================================
contains
!=======================================================================

  subroutine SetServices(gcomp, rc)

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    character(len=*),parameter  :: subname='ocn_comp_nuopc:(SetServices) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    nproc = narea - 1

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_DataInitialize, &
         specRoutine=DataInitialize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_CheckImport, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_CheckImport, &
         specRoutine=ModelCheckImport, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)

    ! Arguments
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !--------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
         acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    use NUOPC, only : NUOPC_isConnected

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)     :: vm
    integer(i4)       :: iam
    integer(i4)       :: lmpicom
    integer(i4)       :: shrlogunit
    character(len=CL) :: logmsg
    character(len=CS) :: cvalue
    logical           :: isPresent, isSet
    character(len=*), parameter :: subname='ocn_comp_nuopc:(InitializeAdvertise) '
    !--------------------------------

    call ESMF_VMLogMemInfo("Entering "//trim(subname))

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, localPet=iam, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    mpi_communicator_ocn = lmpicom

    ! reset shr logging to my log file
    call set_component_logging(gcomp, iam==0, numout, shrlogunit, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    flds_scalar_name = trim(cvalue)
    call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue, *) flds_scalar_num
    write(logmsg,*) flds_scalar_num
    call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_scalar_index_nx
    write(logmsg,*) flds_scalar_index_nx
    call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_scalar_index_ny
    write(logmsg,*) flds_scalar_index_ny
    call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

!    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxPrecipFactor", value=cvalue, rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    read(cvalue,*) flds_scalar_index_precip_factor
!    if ( .not. flds_scalar_index_precip_factor > 0 ) then
!       call shr_sys_abort('flds_scalar_index_precip_factor must be > 0 for pop')
!    else
!       write(logmsg,*) flds_scalar_index_precip_factor
!       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_precip_factor = '//trim(logmsg), ESMF_LOGMSG_INFO)
!    end if

    ! Form of ocean freezing temperature
    ! 'minus1p8' = -1.8 C
    ! 'linear_salt' = -depressT * sss
    ! 'mushy' conforms with ktherm=2
    call NUOPC_CompAttributeGet(gcomp, name="tfreeze_option", value=tfrz_option_driver, &
         isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (.not. isPresent) then
       tfrz_option_driver = 'linear_salt'
    end if

    ! Advertise fields
    call ocn_advertise_fields(gcomp, importState, exportState, flds_scalar_name, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMLogMemInfo("Leaving "//trim(subname))

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    !-----------------------------------------------------------------------
    !  first initializaiton phase of pop2
    !  initialize the timers, communication routines, global reductions,
    !  domain decomposition, grid, and overflows
    !-----------------------------------------------------------------------
    use ESMF               , only: ESMF_VMGet
    use shr_const_mod      , only: shr_const_pi
!    use constants          , only: radius

    ! Initialize POP

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    !  local variables
    type(ESMF_VM)           :: vm
    type(ESMF_DistGrid)     :: distGrid
    type(ESMF_Mesh)         :: Emesh
    type(ESMF_Mesh)         :: tmpMesh
    type(ESMF_Grid)         :: Egrid 
    integer , allocatable   :: gindex_ocn(:)
    integer , allocatable   :: gindex_elim(:)
    integer , allocatable   :: gindex(:)
    !integer                 :: globalID
    character(CL)           :: cvalue
    integer                 :: num_elim_global
    !integer                 :: num_elim_local
    !integer                 :: num_elim
    !integer                 :: num_ocn
    !integer                 :: num_gidx
    !integer                 :: num_elim_gcells ! local number of eliminated gridcells
    !integer                 :: num_elim_blocks ! local number of eliminated blocks
    !integer                 :: num_total_blocks
    !integer                 :: my_elim_start
    !integer                 :: my_elim_end
    integer(i4)             :: lsize, lsizeL
    integer(i4)             :: shrlogunit      ! old values
    integer(i4)             :: iam
    character(len=32)       :: starttype
    integer                 :: n,i,j,l
    integer                 :: lbnum
    integer(i4)             :: errorCode       ! error code
    character(len=*), parameter  :: subname = "ocn_comp_nuopc:(InitializeRealize)"
    logical                 :: mpp_ocn
    integer(i4)             :: gsize
    integer(i4)             :: ier
    INTEGER                 :: inum            ! logical unit

    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS
    errorCode = 0 !POP_Success

    call ESMF_VMLogMemInfo("Entering "//trim(subname))

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
 
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    call ESMF_VMGet(vm, localPet=iam, PetCount=npes, rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return

!    call ESMF_VMGet(vm, pet=iam, peCount=nthreads, rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    if(nthreads==1) then
!       call NUOPC_CompAttributeGet(gcomp, "nthreads", value=cvalue, rc=rc)
!       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
!       read(cvalue,*) nthreads
!    endif

!$  call omp_set_num_threads(nThreads)

#if (defined _MEMTRACE)
    if (iam == 0) then
       lbnum=1
       call memmon_dump_fort('memmon.out','InitializeRealize:start::',lbnum)
    endif
#endif

    !-----------------------------------------------------------------------
    ! initialize the model run
    !-----------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !read(cvalue,*) runid
    read(cvalue,*) cn_exp 

    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) starttype

    if (trim(starttype) == trim('startup')) then
       runtype = "initial"
    else if (trim(starttype) == trim('continue') ) then
       runtype = "continue"
    else if (trim(starttype) == trim('branch')) then
       runtype = "continue"
    else
       write(numout,*) 'ocn_comp_nuopc: ERROR: unknown starttype'
       call shr_sys_abort('ocn_comp_nuopc: ERROR: unknown starttype')
    end if

    ! TODO: Set model_doi_url

!    call get_component_instance(gcomp, inst_suffix, inst_index, rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    inst_name = "OCN"

    !-----------------------------------------------------------------------
    !  first initializaiton phase of pop2
    !  initialize the timers, communication routines, global reductions,
    !  domain decomposition, grid, and overflows
    !-----------------------------------------------------------------------

    call t_startf ('nemo_init')

    call nemo_init(mpi_communicator_ocn)
    nproc = narea - 1

    ! check that all process are still there...
    ! If some process have an error, they will never enter in step and other
    ! processes will wait until the end of the cpu time!
    if ( lk_mpp ) call mpp_max('ocn_comp_nuopc', nstop )
    if (nstop /= 0) then
       if (lwp) then   ! error print
          write(numout,*) ' ===>>> : E R R O R'
          write(numout,*) nstop, ' error have been found after nemo_init from ocn_comp_nuopc'
       end if
       if (lk_mpp) call mppsync       ! sync PEs
       call shr_sys_abort('ocn_comp_nuopc: ocn_init_nuopc: '//SubName//': nstop>0 !')
    endif
    if (lwp) write(numout,cform_aaa)   ! Flag AAAAAAA

    !! Initialize the model time step
    !istp = nit000
    

    call t_stopf ('nemo_init')

    !---------------------------------------------------------------------------
    ! Determine the global index space needed for the distgrid
    !---------------------------------------------------------------------------

    ! Check if ocean domain decomposition only or not
    mpp_ocn = .false.
    if ( jpnij < jpni * jpnj ) then
       mpp_ocn = .true.
    endif

    ! Global domain size
    gsize = Ni0glo*Nj0glo   ! Number of grid points for the global domain (w/o halo)

    ! Local domain size
    lsize = Ni_0*Nj_0     ! Number of grid points for the local domain

    n = 0
    if ((lsize<=0).or.(lsize>gsize)) then
       n = 1
       write(message,FMT='(2(A,I))')     &
          'ocn_comp_nuopc: InitializeRealize: wrong local size: lsize=', lsize, &
          ' nproc=', nproc
    end if
    if (lk_mpp) call mpp_max('ocn_comp_nuopc',n)
    if (n>0) then
       if (lk_mpp) call mppsync       ! sync PEs
       call shr_sys_abort(message)
    end if

    n = lsize
    call mpp_sum('ocn_comp_nuopc',n)
    num_elim_global = gsize - n
    if (gsize /= n) then
       if (mpp_ocn ) then
          if (lwp) then
             write(numout,FMT='(A)') 'ocn_comp_nuopc: InitializeRealize: NEMO is using ocean only domains (mpp_init2)'
             write(numout,FMT='(2(A,I))') ' gsize=', gsize, ' mpp_sum(lsize)=', n
             write(numout,FMT='(2(A,I))') ' jpnij=', jpnij, ' < jpni * jpnj=', jpni * jpnj
          end if
       else
          write(numout,FMT='(A)') 'ocn_comp_nuopc: InitializeRealize: number of points in the global domain'
          write(numout,FMT='(2(A,I))') ' gsize=', gsize, ' mpp_sum(lsize)=', n
          if (lk_mpp) call mppsync       ! sync PEs
          call shr_sys_abort('ocn_comp_nuopc: InitializeRealize: mpp_sum(lsize) /= gsize !')
       end if
    else
       if (lwp) &
         write(numout,*) 'ocn_comp_nuopc: InitializeRealize: mpp_sum(lsize)=', n
    end if
    !
    allocate(gindex_ocn(lsize), stat=ier)
    ier = ABS(ier)
    if (lk_mpp) call mpp_max('ocn_comp_nuopc',ier)
    if (ier/=0) then
       if (lk_mpp) call mppsync       ! sync PEs
       call shr_sys_abort('ocn_comp_nuopc: InitializeRealize: failed allocation (gindex_ocn)!')
    end if
    gindex_ocn(:) = 0
    !
    n = 0
    do j=Njs0,Nje0
       l = (mjg0(j)-1)*Ni0glo
       do i=Nis0,Nie0
          n=n+1
          gindex_ocn(n) = l + mig0(i)
       enddo
    enddo

    n = 0
    if (any(gindex_ocn(:)<0)) then
       n = 1
    end if
    if (lk_mpp) call mpp_max('ocn_comp_nuopc',n)
    if (n>0) then
       if (lk_mpp) call mppsync       ! sync PEs
       call shr_sys_abort('ocn_comp_nuopc: InitializeRealize: in InitializeRealize gindex_ocn<0!')
     end if

    n = 0
    if (any(gindex_ocn(:)>gsize)) then
       n = 1
    end if
    if (lk_mpp) call mpp_max('ocn_comp_nuopc',n)
    if (n>0) then
       if (lk_mpp) call mppsync       ! sync PEs
       call shr_sys_abort('ocn_comp_nuopc: InitializeRealize: in InitializeRealize gindex_ocn>gsize!')
    end if

    ! Take into account suppressed land-only subdomains
    n = 0
    if (numsls>0) then
      !
      if (narea <= numsls) then    ! This process handles one of the suppressed land-only subdomains
         ! Number of grid points of the associated suppressed land-only subdomains
         lsizeL = (Nie0L-Nis0L+1)*(Nje0L-Njs0L+1)
      else
         lsizeL = 0
      end if
      lsizeL = MAX(0, lsizeL)
      ! compute total number of grid points belogning to suppressed land-only subdomains
      n = lsizeL
      if (lk_mpp) call mpp_sum('ocn_comp_nuopc',n)
      if (lwp) then
         write(numout,*) 'ocn_comp_nuopc: InitializeRealize: land points used in gindex_elim'
         write(numout,*) 'Expected value ', num_elim_global, ' computed value ', n
      end if

      n = 0
      if (narea <= numsls) then    ! This process handles one of the suppressed land-only subdomains
         ! allocate gindex_elim
         allocate(gindex_elim(lsizeL))
         gindex_elim(:) = 0
         ! add values to gindex_elim
         n = 0
         do j=Njs0L,Nje0L
            l = ((j + njmppL - 1 - nn_hls)-1)*Ni0glo
            do i=Nis0L,Nie0L
               n=n+1
               gindex_elim(n) = l + (i + nimppL - 1 - nn_hls)
            enddo
         enddo
         !
         n = 0
         if (any(gindex_elim(:)<0)) then
            n = 1
         end if
      end if

      if (lk_mpp) call mpp_max('ocn_comp_nuopc',n)
      if (n>0) then
         if (lk_mpp) call mppsync       ! sync PEs
         call shr_sys_abort('ocn_comp_nuopc: InitializeRealize: in InitializeRealize gindex_elim<0!')
      end if

      n = 0
      if (narea <= numsls) then    ! This process handles one of the suppressed land-only subdomains
         n = 0
         if (any(gindex_elim(:)>gsize)) then
            n = 1
         end if
      end if

      if (lk_mpp) call mpp_max('ocn_comp_nuopc',n)
      if (n>0) then
         if (lk_mpp) call mppsync       ! sync PEs
         call shr_sys_abort('ocn_comp_nuopc: InitializeRealize: in InitializeRealize gindex_elim>gsize!')
      end if
          !
      if (narea <= numsls) then    ! This process handles one of the suppressed land-only subdomains

          ! create a global index that includes both active and eliminated gridcells
          allocate(gindex(lsize + lsizeL))
          !
          gindex(:lsize) = gindex_ocn(:lsize)
          gindex(lsize+1:lsize+lsizeL) = gindex_elim(:lsizeL)
          !
          deallocate(gindex_elim)
          !
        else    ! This process does not handle suppressed land-only subdomains
          !
          allocate(gindex(lsize))
          !
          gindex(:) = gindex_ocn(:)

       end if

    else    ! No suppressed land-only subdomains
       !
       allocate(gindex(lsize))
       !
       gindex(:) = gindex_ocn(:)
       !
    end if
    !
    deallocate(gindex_ocn)
    !

    !---------------------------------------------------------------------------
    ! Create NEMO Mesh from SCRIP file and dist grid
    ! DP new order of commands for NEMO
    !---------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='mesh_ocn', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    Egrid = ESMF_GridCreate(filename=trim(cvalue),fileformat=ESMF_FILEFORMAT_SCRIP, addCornerStagger=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    tmpMesh = ESMF_MeshCreate(grid=Egrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridDestroy(Egrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
    Emesh = ESMF_MeshCreate(mesh=tmpMesh, elementDistGrid=DistGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MeshDestroy(tmpMesh, rc=rc) 
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate(gindex)

    if (nproc==0) then
       write(numout,*)'mesh file for nemo domain is ',trim(cvalue)
    end if

    !-----------------------------------------------------------------
    ! Realize the actively coupled fields
    !-----------------------------------------------------------------

    call ocn_realize_fields(gcomp, mesh=Emesh, flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMLogMemInfo("Leaving "//trim(subname))

  end subroutine InitializeRealize

  !===============================================================================

  subroutine DataInitialize(gcomp, rc)

    !-----------------------------------------------------------------------
    !  second initializaiton phase of pop2
    !  - initialize passive tracer modules -- after call init_forcing_coupled
    !  - initialize vertical mixing variables
    !  - initialize niw driven mixing
    !  - initialize geothermal heat flux
    !  - initialize horizontal mixing variables
    !  - initialize overflow regional values
    !  - initialize advection variables
    !  - initialize shortwave absorption
    !  - initialize estuary parameterization
    !  - partial coupling forcing initialization
    !  - initialize overflows output diagnostics filename
    !  - initialize output; subroutine init_output calls
    !  - initialize global budget diagnostics
    !  - initialize step timers
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)          :: clock
    type(ESMF_State)          :: importState
    type(ESMF_State)          :: exportState
    type(ESMF_TimeInterval)   :: timeStep        ! Model timestep
    type(ESMF_Time)           :: starttime
    character(CL)             :: cvalue
    integer(i4)               :: ocn_cpl_dt
!    integer(int_kind)         :: pop_cpl_dt
    integer(i4)               :: start_ymd
    integer(i4)               :: start_tod
    integer(i4)               :: start_year
    integer(i4)               :: start_day
    integer(i4)               :: start_month
    integer(i4)               :: start_hour
    integer(i4)               :: errorCode       ! error code
    integer(i4)               :: shrlogunit      ! old values
    integer(i4)               :: shrloglev       ! old values
!    integer                   :: ocnid
    character(len=*), parameter  :: subname = "ocn_comp_nuopc:(DataInitialize)"
!    character(len=32)  :: flux_epbal         ! infodata flux_epbal
    integer(i4)               :: nstp, n
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMLogMemInfo("Entering "//trim(subname))

    ! Initialize the model time step
    istp = nit000
    ndt05 = NINT( 0.5 * rn_Dt )

    !--------------------------------
    ! Reset shr logging to my log file
    !--------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (numout)

    !-----------------------------------------------------------------------
    ! register non-standard incoming fields
    !-----------------------------------------------------------------------

    ! query the Component for its importState, exportState and clock
    call ESMF_GridCompGet(gcomp, importState=importState, exportState=exportState, clock=clock, rc=rc)

    call t_startf ('nemo_cesm_init')

    call sbc_cpl_cesm_init
 
    call init_qflxice(tfrz_option_driver)

    call t_stopf ('nemo_cesm_init')

    !if (lwp) then
    !  write(numout,FMT='(A,I)') ' ocn_comp_nuopc:DataInitialize: NEMO ID ', OCNID
    !  write(numout,FMT='(A,I)') ' ocn_comp_nuopc:DataInitialize: NEMO MPICOM ', mpicom_o
    !end if

    call shr_file_getLogUnit (shrlogunit)
    !call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (numout)

    if (runtype == 'initial') then

       ! check for consistency of nemo ln_rstart and runtype from seq_infodata
       if (ln_rstart) then
          if (lwp) then
             write(numout,FMT='(A,I)') 'ocn_comp_nuopc: DataInitialize: '//&
                 'CESM runtype=initial BUT NEMO ln_rstart=.T. !'
          end if
       endif

       call ESMF_ClockGet( clock, startTime=startTime, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_TimeGet( startTime, yy=start_year, mm=start_month, dd=start_day, s=start_tod, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !call shr_cal_date2ymd(start_ymd,start_year,start_month,start_day)
       call shr_cal_ymd2date(start_year,start_month,start_day,start_ymd)

       if (nyear /= start_year) then
          if (lwp) then
             write(numout,fmt='(A,I)') 'DataInitialize: nyear      ', nyear
             write(numout,fmt='(A,I)') 'DataInitialize: start_year ', start_year
          endif
          call shr_sys_abort('ocn_comp_nuopc: DataInitialize: nyear does not match start_year')
       end if
       if (nmonth /= start_month) then
          if (lwp) then
             write(numout,fmt='(A,I)') 'DataInitialize: nmonth      ', nmonth
             write(numout,fmt='(A,I)') 'DataInitialize: start_month ', start_month
          endif
          call shr_sys_abort('ocn_comp_nuopc: DataInitialize: nmonth does not match start_month')
       end if
       if (nday /= start_day) then
          if (lwp) then
             write(numout,fmt='(A,I)') 'DataInitialize: nday      ', nday
             write(numout,fmt='(A,I)') 'DataInitialize: start_day ', start_day
          endif
       end if
#ifndef _HIRES 
       n = MOD(nsec_year-ndt05, NINT(rday))
       if (n /= start_tod) then
          if (lwp) then
             write(numout,fmt='(A,I)') 'DataInitialize: nsec      ', n
             write(numout,fmt='(A,I)') 'DataInitialize: start_tod ', start_tod
          end if
       end if
#endif
    else
       ! check for consistency of nemo ln_rstart and runtype from seq_infodata
       if (.not. ln_rstart) then
          if (lk_mpp) call mppsync       ! sync PEs
          call shr_sys_abort('ocn_comp_nuopc: DataInitialize: '//&
            'CESM runtype='//TRIM(runtype)//' BUT NEMO ln_rstart=.F. !')
       endif

    end if

    !-----------------------------------------------------------------
    ! Initialize MCT gsmaps and domains
    !-----------------------------------------------------------------

!    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    read(cvalue,*) ocnid  ! convert from string to integer
!
!    call nemo_mct_init(ocnid, mpi_communicator_ocn)
!    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    !-----------------------------------------------------------------------
    ! Initialize flags and shortwave absorption profile
    ! Note that these cpl_write_xxx flags have no freqency options
    ! set; therefore, they will retain a default value of .false.
    ! unless they are explicitly set .true.  at the appropriate times
    !-----------------------------------------------------------------------
!
!    call init_time_flag('cpl_write_restart',cpl_write_restart, owner = 'DataInitialize')
!    call init_time_flag('cpl_write_history',cpl_write_history, owner = 'DataInitialize')
!    call init_time_flag('cpl_write_tavg'   ,cpl_write_tavg,    owner = 'DataInitialize')
!    call init_time_flag('cpl_diag_global'  ,cpl_diag_global,   owner = 'DataInitialize')
!    call init_time_flag('cpl_diag_transp'  ,cpl_diag_transp,   owner = 'DataInitialize')
!
!    lsmft_avail = .true.
!    tlast_coupled = c0
!
    !-----------------------------------------------------------------------
    ! initialize necessary coupling info
    !-----------------------------------------------------------------------

    call ESMF_ClockGet(clock, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeIntervalGet( timeStep, s=ocn_cpl_dt, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (MOD(ocn_cpl_dt, INT(rn_Dt)) /= 0) then
       write(numout,*)' ocn_cpl_dt= ',ocn_cpl_dt, &
                      '      rn_dt= ',INT(rn_Dt)
       if (lk_mpp) call mppsync       ! sync PEs
       call shr_sys_abort('ocn_comp_nuopc: DataInitialize: ocn_cpl_dt must be an exact multiple of rn_dt')
    end if
    ! # of model times step in one coupling time step
    nn_ncpl = ocn_cpl_dt/INT(rn_Dt)
    if (nn_ncpl /= nn_fsbc) then
       write(numout,*)' nn_ncpl= ',nn_ncpl, '  nn_fsbc= ', nn_fsbc
       if (lk_mpp) call mppsync       ! sync PEs
       call shr_sys_abort('ocn_comp_nuopc: DataInitialize: nn_ncpl must be equal to nn_fsbc!')
    end if
    if (ldiag_cpl .AND. lwp) then
      write(numout,*) 'ocn_comp_nuopc: DataInitialize: coupling time step (sec)',  &
      ocn_cpl_dt
      write(numout,*) 'ocn_comp_nuopc: DataInitialize: coupling freq. (steps)  ',  &
      nn_ncpl
    end if

!    ! FIXME: used by POP when OCN_COUPLING=partial . Needed here ???
!    call NUOPC_CompAttributeGet(gcomp, name='flux_epbal', value=cvalue, rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    read(cvalue,*) flux_epbal
!!    call seq_infodata_GetData( infodata, flux_epbal=flux_epbal)
!    if (lwp) &
!      write(numout,FMT='(A)') 'ocn_comp_nuopc: DataInitialize: flux_epbal '//trim(flux_epbal)
!    IF (trim(flux_epbal) == 'ocn') THEN
!      if (lwp) &
!        write(numout,FMT='(A)') 'ocn_comp_nuopc: DataInitialize: send precip_fact to cpl'
!      ! From cesm1.2 precip_fact=1.0 ; Previous releases was precip_fact=1.0e6
!      ! !!!
!!      call NUOPC_CompAttributeSet(gcomp, name="precip_fact", value=1.0e6_wp, rc=rc)
!!      call seq_infodata_PutData( infodata, precip_fact=1.0e6_wp)
!    END IF

    !-----------------------------------------------------------------------
    ! send export state
    !-----------------------------------------------------------------------

    call ocn_sum_buffer(exportState, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ocn_export(exportState, flds_scalar_name, ldiag_cpl, errorCode, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (errorCode /= 0) then
       if (lk_mpp) call mppsync       ! sync PEs
       call shr_sys_abort('ocn_comp_nuopc: DataInitialize: in ocn_export')
    endif

!    call State_SetScalar(dble(nx_global), flds_scalar_index_nx, exportState, &
    call State_SetScalar(dble(Ni0glo), flds_scalar_index_nx, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

!    call State_SetScalar(dble(ny_global), flds_scalar_index_ny, exportState, &
    call State_SetScalar(dble(Nj0glo), flds_scalar_index_ny, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

!    if ( lsend_precip_fact ) then
!       call State_SetScalar(precip_fact, flds_scalar_index_precip_factor, exportState, &
!            flds_scalar_name, flds_scalar_num, rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    else
!       call State_SetScalar(1._r8, flds_scalar_index_precip_factor, exportState, &
!            flds_scalar_name, flds_scalar_num, rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    end if

!----------------------------------------------------------------------------
!
! If coupling frequency is higher than daily advance NEMO 1 coupling time step
! Needed to synchronize NEMO, which assumes that the starting time is 00:00.
! This applies to start_up and hybrid first simulation, while branch is like 
! a continuation run (see also setting of start_type in drv_in)
!----------------------------------------------------------------------------

    if ( runtype == 'initial' ) then
       ! Number of model time step per coupling time step
       ! NOTE: if daily coupling nstp=0, no need to synchronize NEMO
       nstp = MOD(ocn_cpl_dt, INT(rday))/INT(rn_Dt)

       if (nstp /= 0 .and. nstp /= nn_fsbc) then
          write(numout,*) ' Coupling frequency different from SBC/ICE frequency!'
          if (lk_mpp) call mppsync       ! sync PEs
          call shr_sys_abort('ocn_comp_nuopc: '//SubName//': nstp != nn_fsbc !')
       end if

       if (lwp .and. nstp/=0) then
          write(numout,*) ' Coupling frequency: ', nstp, ' model time steps'
          write(numout,*) ' Advance NEMO 1 coupling time step for coupling synchronization'
       end if

       do n=1, nstp
          lice_form_ts = lice_form(istp)
          lice_cpl_ts  = (MOD(istp-nit000+1, nn_ncpl)==0)
          call increment_tlast_ice(istp)
#if defined key_qco   ||   defined key_linssh
          call stp_MLF(istp)
#else
          call stp(istp)
#endif
          
          ! check for errors
          if (lk_mpp) call mpp_max('ocn_comp_nuopc', nstop )
          if (nstop /= 0) then
             if (lwp) then   ! error print
                write(numout,*) ' ===>>> : E R R O R'
                write(numout,*) nstop, ' error have been found in step from ocn_comp_nuopc'
             end if
             if (lk_mpp) call mppsync       ! sync PEs
             call shr_sys_abort('ocn_comp_nuopc: '//SubName//': nstop>0 !')
          endif

          call ocn_sum_buffer(exportState, rc)
          if (MOD(istp-nit000+1, nn_ncpl)==0) then
             call ocn_export(exportState, flds_scalar_name, ldiag_cpl, errorCode, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif

          istp = istp+1

       end do

    end if ! ln_rstart

#if (defined _MEMTRACE)
    if (iam  == 0) then
       lbnum=1
       call memmon_dump_fort('memmon.out','DataInitialize:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

    !-----------------------------------------------------------------------
    ! check whether all Fields in the exportState are "Updated"
    !-----------------------------------------------------------------------

    if (NUOPC_IsUpdated(exportState)) then
      call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="true", rc=rc)

      call ESMF_LogWrite("NEMO - Initialize-Data-Dependency SATISFIED!!!", ESMF_LOGMSG_INFO)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end if

    call ESMF_VMLogMemInfo("Leaving "//trim(subname))

  end subroutine DataInitialize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)

    !-----------------------------------------------------------------------
    ! Run POP for a coupling interval
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    !  local variables
    type(ESMF_State)             :: importState
    type(ESMF_State)             :: exportState
    type(ESMF_Clock)             :: clock
    type(ESMF_Alarm)             :: alarm
    type(ESMF_Time)              :: currTime
    type(ESMF_Time)              :: nextTime
    character(ESMF_MAXSTR)       :: cvalue
    integer(i4)                  :: errorCode  ! error flag
    integer                      :: ymd_sync   ! Sync clock current date (YYYYMMDD)
    integer                      :: tod_sync   ! Sync clcok current time of day (sec)
    integer                      :: lbnum
    integer                      :: yr_sync
    integer                      :: mon_sync
    integer                      :: day_sync
    integer                      :: shrlogunit ! old values
    integer                      :: shrloglev  ! old values
    logical                      :: first_time = .true.
    character(len=*), parameter  :: subname = "ocn_comp_nuopc: (ModelAdvance)"
    integer                      :: nitrst_old
    logical                      :: lsnd, lrcv
!DP+
    logical                      :: rstwr       ! .true. ==> write restart file before returning
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS
    !errorCode = POP_Success
    errorCode = 0

    call ESMF_VMLogMemInfo("Entering "//trim(subname))

    !-----------------------------------------------------------------------
    ! skip first coupling interval for an initial run
    !-----------------------------------------------------------------------

    ! NOTE: pop starts one coupling interval ahead of the rest of the system
    ! so to have it be in sync with the rest of the system, simply skip the first
    ! coupling interval for a initial run only
    if (first_time) then
       first_time = .false.
       if (runtype == 'initial') then
          if (nproc == 0) then
             write(numout,*)'Returning at first coupling interval'
          end if
          RETURN
       end if
    end if

!$  call omp_set_num_threads(nThreads)

nitrst_old = 0
nproc = narea - 1
#if (defined _MEMTRACE)
    !if(my_task == 0 ) then
    if(nproc == 0 ) then
       lbnum=1
       call memmon_dump_fort('memmon.out',subname//':start::',lbnum)
    endif
#endif

!!!DP it seems this check is not required with NEMO 
!    !--------------------------------------------------------------------
!    ! check that nemo internal clock is in sync with ESMF clock
!    !--------------------------------------------------------------------
!
!    ! NOTE: that in nuopc - the ESMF clock is updated at the end of the timestep
!    ! whereas in cpl7 it was updated at the beginning - so need to have the check
!    ! at the beginning of the time loop BEFORE pop updates its time step
!
!    ! nemo clock
!    ymd = iyear*10000 + imonth*100 + iday
!    tod = ihour*seconds_in_hour + iminute*seconds_in_minute + isecond
!
    ! model clock
    call NUOPC_ModelGet(gcomp, modelClock=clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet( currTime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr_sync, mon_sync, day_sync, ymd_sync)
!
!    ! check
!!$OMP MASTER
!    if ( (ymd /= ymd_sync) .or. (tod /= tod_sync) ) then
!       write(stdout,*)' pop2 ymd=',ymd     ,'  pop2 tod= ',tod
!       write(stdout,*)' sync ymd=',ymd_sync,'  sync tod= ',tod_sync
!       write(stdout,*)' Internal pop2 clock not in sync with Sync Clock'
!       call ESMF_LogWrite(subname//" Internal POP clock not in sync with ESMF model clock", ESMF_LOGMSG_INFO)
!       !call shr_sys_abort(subName// ":: Internal POP clock not in sync with ESMF model Clock")
!    end if
!!$OMP END MASTER
!    !-----------------------------------------------------------------------
!    !  start up the main timer
!    !-----------------------------------------------------------------------
!
!    call timer_start(timer_total)
!
    !-----------------------------------------------------------------------
    ! reset shr logging to my log file
    !----------------------------------------------------------------------------

 !   call shr_file_getLogUnit (shrlogunit)
 !   call shr_file_setLogUnit (stdout)
    call shr_file_getLogUnit (shrlogunit)
    !call shr_file_getLogLevel(shrloglev)
    call shr_sys_flush(shrlogunit)
    call shr_file_setLogUnit (numout)
    call shr_sys_flush(numout)

!    if (ldiag_cpl) then
!       call register_string ('info_debug_ge2')
!    endif
    if (lwp .AND. ldiag_cpl) write(numout,FMT='(A)') ' BEGINNING of ocn_init_run'

    !--------------------------------
    ! Query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! restart flag (rstwr) will assume only an eod restart for now
    !----------------------------------------------------------------------------

    ! Note this logic triggers off of the component clock rather than the internal pop time
    ! The component clock does not get advanced until the end of the loop - not at the beginning

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

!DP+
    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       rstwr = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       rstwr = .false.
    endif

!    if(rstwr) then
!       if (lwp) write(numout,*) 'nitrst=', nitrst, ' restore nitrst to ', nitrst_old
!       nitrst = nitrst_old  ! Restore original value
!
!!       call override_time_flag(cpl_write_restart, value=.true.)
!!       call ccsm_char_date_and_time ! set time_management module vars cyear, cmonth, ..
!!       write(message,'(6a)') 'driver requests restart file at eod  ', cyear,'/',cmonth,'/',cday
!!       call document ('ModelAdvance:', message)
!    endif

    !-----------------------------------------------------------------------
    ! advance the model in time over coupling interval
    ! write restart dumps and archiving
    !-----------------------------------------------------------------------

    ! Note that all ocean time flags are evaluated each timestep in time_manager
    ! tlast_coupled is set to zero at the end of ocn_export

    advance: do

       lsnd = .false.
       lrcv = .false.
       ! -----
       ! obtain import state data
       ! -----
       if (MOD(istp-nit000, nn_ncpl)==0) then 
!       if (check_time_flag(cpl_ts) .or. nsteps_run == 0) then

          call ocn_import(importState, flds_scalar_name, ldiag_cpl, errorCode, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

!          if (errorCode /= POP_Success) then
!             call POP_ErrorPrint(errorCode)
!             call exit_POP(sigAbort, 'ERROR in step')
!          endif
          if (errorCode /= 0) then
             if (lk_mpp) call mppsync       ! sync PEs
             call shr_sys_abort('ocn_comp_nuopc: '//SubName//': in ocn_import_nuopc')
          endif
          lrcv = .true.

!          ! update orbital parameters
!
!          call pop_orbital_update(clock, stdout, mastertask, orb_eccen, orb_obliqr, orb_lambm0, orb_mvelpp, rc)
!          if (ChkErr(rc,__LINE__,u_FILE_u)) return
!
!          call pop_set_coupled_forcing
       end if

!DP+
       ! open restart file if requested from driver
!      call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
!      if (ChkErr(rc,__LINE__,u_FILE_u)) return  
!
!      if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
!         if (ChkErr(rc,__LINE__,u_FILE_u)) return
!         call ESMF_AlarmRingerOff( alarm, rc=rc )
!         if (ChkErr(rc,__LINE__,u_FILE_u)) return

         if (rstwr .and. nn_ncpl>1 .and. (MOD(istp-nit000+2, nn_ncpl)==0)) then
!          if (nn_ncpl>1 .and. (MOD(istp-nit000+2, nn_ncpl)==0)) then
            if (lwp) write(numout,*) 'driver requests restart file at kt= ', istp+1
            if (istp+1 /= nitrst) THEN
               if (lwp) write(numout,*) 'nitrst=', nitrst, ' set nitrst to ', istp+1
               nitrst_old = nitrst  ! Save original value
               nitrst = istp+1
            endif
            if (lwp) write(numout,*) &
              'open restart file at kt= ', istp, ' ndastp=', ndastp
            if (lwp) call shr_sys_flush(numout)
          end if
!       end if

       ! Compute sea ice formation
       lice_form_ts = lice_form(istp)
!       lice_form_ts = .false.
       if (lice_form_ts) then
          if (lwp .AND. istp-nit000+1<48) then
            write(numout,*)  &
            SubName, ': time step (istp) ', istp, ' lice_form_ts=', lice_form_ts
          endif
       end if

       ! Compute sea ice formation/melting flux at coupling time step
       lice_cpl_ts = (MOD(istp-nit000+1, nn_ncpl)==0)
!       lice_cpl_ts = .false.
       if (lice_cpl_ts) then
          if (lwp .AND. istp-nit000+1<48) then
            write(numout,*)  &
            SubName, ': time step (istp) ', istp, ' lice_cpl_ts=', lice_cpl_ts
          endif
       end if
       call shr_sys_flush(numout)

       call increment_tlast_ice(istp)

       ! advance NEMO 1 model time step
#if defined key_qco   ||   defined key_linssh
       call stp_MLF(istp)
#else
       call stp(istp)
#endif
 
       if (lrcv .and. ldiag_cpl .and. lwp)     &
          write(numout,*)  &
          SubName, ': ocn_import called at time step (istp) ', istp, &
          ' ndastp=', ndastp

       if (ldiag_cpl .and. lwp)     &
          write(numout,*)  &
          SubName, ': stp() called at time step (istp) ', istp, ' ndastp=', ndastp

       ! check for errors
       if (lk_mpp) call mpp_max('ocn_comp_nuopc', nstop )
       if (nstop /= 0) then
          if (lwp) then   ! error print
             write(numout,*) ' ===>>> : E R R O R'
             write(numout,*) nstop, ' error have been found'
          end if
          if (lk_mpp) call mppsync       ! sync PEs
          call shr_sys_abort('ocn_comp_nuopc: '//SubName//': nstop>0 !')
       endif

       ! return export state to driver
       call ocn_sum_buffer(exportState, rc)
       if (MOD(istp-nit000+1, nn_ncpl)==0) then
          call ocn_export(exportState, flds_scalar_name, ldiag_cpl, errorCode, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return 
!          call ocn_export(o2x_o, errorCode)
          if (errorCode /= 0) then
             if (lk_mpp) call mppsync       ! sync PEs
             call shr_sys_abort('ocn_comp_nuopc: '//SubName//': in ocn_export')
          endif
          if (ldiag_cpl .and. lwp)     &
             write(numout,*)  &
             SubName, ': ocn_export called at time step (istp) ', istp, &
             ' ndastp=', ndastp
          istp = istp+1
          exit advance
       end if

       istp = istp+1

    enddo advance

    if(rstwr) then
       if (lwp) write(numout,*) 'nitrst=', nitrst, ' restore nitrst to ', nitrst_old
       nitrst = nitrst_old  ! Restore original value

!       call override_time_flag(cpl_write_restart, value=.true.)
!       call ccsm_char_date_and_time ! set time_management module vars cyear, cmonth, ..
!       write(message,'(6a)') 'driver requests restart file at eod  ', cyear,'/',cmonth,'/',cday
!       call document ('ModelAdvance:', message)
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)
    !call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(shrlogunit)

!    call timer_stop(timer_total)

#if (defined _MEMTRACE)
    !if(my_task == 0) then
    if(nproc == 0) then
       lbnum=1
       call memmon_dump_fort('memmon.out',subname//':end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

    call ESMF_VMLogMemInfo("Leaving  "//trim(subname))

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option ! Restart option units
    integer                  :: restart_n      ! Number until restart interval
    integer                  :: restart_ymd    ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    character(len=256)       :: stop_option    ! Stop option units
    integer                  :: stop_n         ! Number until stop interval
    integer                  :: stop_ymd       ! Stop date (YYYYMMDD)
    type(ESMF_ALARM)         :: stop_alarm
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname='ocn_comp_nuopc:(ModelSetRunClock)'
    !--------------------------------

    call ESMF_VMLogMemInfo("Entering  "//trim(subname))
    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart and stop alarms
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO)

       !----------------
       ! Restart alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------
       ! Stop alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="stop_option", value=stop_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="stop_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_n

       call NUOPC_CompAttributeGet(gcomp, name="stop_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_ymd

       call alarmInit(mclock, stop_alarm, stop_option, &
            opt_n   = stop_n,           &
            opt_ymd = stop_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_stop', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(stop_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    call ESMF_VMLogMemInfo("Leaving  "//trim(subname))

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelCheckImport(model, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)  :: clock
    type(ESMF_Time)   :: currTime
    integer(i4)       :: yy  ! current date (YYYYMMDD)
    integer(i4)       :: mon ! current month
    integer(i4)       :: day ! current day
    integer(i4)       :: tod ! current time of day (sec)
    character(len=*),parameter :: subname='ocn_comp_nuopc:(ModelCheckImport)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS
    nproc = narea - 1

    call ESMF_VMLogMemInfo("Entering  "//trim(subname))

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(currTime, yy=yy, mm=mon, dd=day, s=tod, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !if (mastertask) then
    !if (nproc == 0) then
    !   write(numout,*)' CheckImport nemo year = ',yy
    !   write(numout,*)' CheckImport nemo mon  = ',mon
    !   write(numout,*)' CheckImport nemo day  = ',day
    !   write(numout,*)' CheckImport nemo tod  = ',tod
    !end if

    call ESMF_VMLogMemInfo("Leaving  "//trim(subname))

  end subroutine ModelCheckImport

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)

    !--------------------------------
    ! NEMO finalization that shuts down NEMO gracefully (we hope).
    ! Exits the message environment and checks for successful execution.
    ! --------------------------------

!    use output          , only : final_output
!    use passive_tracers , only : passive_tracers_timer_print_all

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
!    integer (POP_i4) :: errorCode         ! error code
    character(len=*),parameter :: subname='ocn_comp_nuopc:(ModelFinalize)'
    !--------------------------------

    call ESMF_VMLogMemInfo("Entering  "//trim(subname))

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !                            !------------------------!
    !                            !==  finalize the run  ==!
    !                            !------------------------!
    if (lwp) write(numout,cform_aaa)   ! Flag AAAAAAA
    !
    if ( nstop /= 0 .and. lwp ) then   ! error print
       write(numout,*) ' ===>>> : E R R O R'
       write(numout,*) nstop, 'error have been found'
    endif
    !
    !CALL iom_context_finalize( cxios_context )   ! Finalize xios files
    !IF( ln_crs ) CALL iom_context_finalize( trim(cxios_context)//"_crs" ) !
    !
    IF ( ln_timing ) CALL timing_finalize
    !
    call nemo_closefile
#if defined key_xios
    CALL xios_finalize  ! end mpp communications with xios
#endif
    !
    ! Free allocated space
    call sbc_cpl_cesm_finalize

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    call ESMF_VMLogMemInfo("Leaving  "//trim(subname))

  end subroutine ModelFinalize

  !===============================================================================

 function lice_form(kt)

   integer, intent(in) :: kt

   logical :: lice_form

   !integer :: n

   !lice_form = .false.
   lice_form = .true.

   !do n=1,min(nn_nits,nn_ncpl)
   !   lice_form = lice_form .OR. (MOD(kt-nit000+n, nn_ncpl)==0)
   !   if (lice_form) exit
   !end do

 end function lice_form

  !===============================================================================

  subroutine ocn_advertise_fields(gcomp, importState, exportState, flds_scalar_name, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    type(ESMF_State)               :: importState
    type(ESMF_State)               :: exportState
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(out) :: rc

    ! local variables
    integer       :: n
    character(CS) :: stdname
    character(CS) :: cvalue
    character(CS) :: cname
    integer       :: ice_ncat
    logical       :: flds_i2o_per_cat  ! .true. => select per ocn thickness category
    logical       :: flds_co2a
    logical       :: flds_co2b
    logical       :: flds_co2c
    integer       :: ndep_nflds
    character(len=*), parameter :: subname='(ocn_advertise_fields)'
    !-------------------------------------------------------------------------------

    nproc = narea - 1
    rc = ESMF_SUCCESS

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !-----------------
    ! optional input from cice columns due to ice thickness categories
    !-----------------

    !-----------------
    ! advertise import fields
    !-----------------

    call fldlist_add(fldsToOcn_num, fldsToOcn, trim(flds_scalar_name))

    ! from ice
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Si_ifrac')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_melth')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_meltw')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_salt')

    ! from river
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_rofl')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_rofi')

    ! from mediator
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'So_duu10n')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_tauy')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_taux')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_lat')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_sen')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_lwup')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_evap')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_swnet')

    ! from atmosphere
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sa_pslv')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_lwdn')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_snow')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_rain')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_bcph'  , ungridded_lbound=1, ungridded_ubound=3)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_dstdry', ungridded_lbound=1, ungridded_ubound=4)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_dstwet', ungridded_lbound=1, ungridded_ubound=4)

    ! from atm co2 fields
    call NUOPC_CompAttributeGet(gcomp, name='flds_co2a', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2a
    !if (my_task == master_task) then
    if (nproc == 0) then
       write(numout,'(a)') trim(subname)//'flds_co2a = '// trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2b', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2b
    !if (my_task == master_task) then
    if (nproc == 0) then
       write(numout,'(a)') trim(subname)//'flds_co2b = '// trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2c', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2c
    !if (my_task == master_task) then
    if (nproc == 0) then
       write(numout,'(a)') trim(subname)//'flds_co2c = '// trim(cvalue)
    end if

    if (flds_co2a .or. flds_co2c) then
       call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sa_co2diag')
       call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sa_co2prog')
    end if

    do n = 1,fldsToOcn_num
       call NUOPC_Advertise(importState, standardName=fldsToOcn(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    !-----------------
    ! advertise export fields
    !-----------------

    ! Determine if ocn is sending temperature and salinity data to glc
    call NUOPC_CompAttributeGet(gcomp, name="ocn2glc_coupling", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ocn2glc_coupling
    !if (my_task == master_task) then
    if (nproc == 0) then
       write(numout,'(a,L1)') trim(subname) // 'ocn2glc coupling is ',ocn2glc_coupling
    end if

    ! Determine number of ocean levels and ocean level indices
    if (ocn2glc_coupling) then
       call NUOPC_CompAttributeGet(gcomp, name="ocn2glc_levels", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       num_ocn2glc_levels = shr_string_listGetNum(cvalue)
       allocate(ocn2glc_levels(num_ocn2glc_levels))
       do n = 1,num_ocn2glc_levels
          call shr_string_listGetName(cvalue, n, cname, rc)
          read(cname,*) ocn2glc_levels(n)
       end do
       !if (my_task == master_task) then
       if (nproc == 0) then
          write(numout,'(a,i0)') trim(subname)//' number of ocean levels sent to glc = ',num_ocn2glc_levels
          write(numout,*)' ',trim(subname)//' ocean level indices are ',ocn2glc_levels
       end if
    end if

    call fldlist_add(fldsFrOcn_num, fldsFrOcn, trim(flds_scalar_name))
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_omask')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_t')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_u')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_v')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_s')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_dhdx')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_dhdy')
    !call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_bldepth')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'Fioo_q')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'Faoo_fco2_ocn')
    if (ocn2glc_coupling) then
       call fldlist_add(fldsFrOcn_num , fldsFrOcn, 'So_t_depth', &
            ungridded_lbound=1, ungridded_ubound=num_ocn2glc_levels)
       call fldlist_add(fldsFrOcn_num , fldsFrOcn, 'So_s_depth', &
            ungridded_lbound=1, ungridded_ubound=num_ocn2glc_levels)
    end if

    do n = 1,fldsFrOcn_num
       call NUOPC_Advertise(exportState, standardName=fldsFrOcn(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ocn_advertise_fields

  !==============================================================================
  subroutine ocn_realize_fields(gcomp, mesh, flds_scalar_name, flds_scalar_num, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    type(ESMF_Mesh)  , intent(in)  :: mesh
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(in)  :: flds_scalar_num
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_State)      :: importState
    type(ESMF_State)      :: exportState
    type(ESMF_Field)      :: lfield
    integer               :: spatialDim
    integer               :: numOwnedElements
    integer               :: i,j,n !iblock,n
!    type(block)           :: this_block         ! block information for current block
    real(wp), allocatable :: mesh_areas(:)
    real(wp), allocatable :: model_areas(:)
    real(wp), pointer     :: dataptr(:)
    !integer               :: num_ocn
    real(wp)              :: max_mod2med_areacor
    real(wp)              :: max_med2mod_areacor
    real(wp)              :: min_mod2med_areacor
    real(wp)              :: min_med2mod_areacor
    real(wp)              :: max_mod2med_areacor_glob
    real(wp)              :: max_med2mod_areacor_glob
    real(wp)              :: min_mod2med_areacor_glob
    real(wp)              :: min_med2mod_areacor_glob
    real(wp), pointer     :: ownedElemCoords(:)
    real(wp), pointer     :: latModel(:), latMesh(:)
    real(wp), pointer     :: lonModel(:), lonMesh(:)
    real(wp)              :: diff_lon
    real(wp)              :: diff_lat
    logical               :: labort_lon, labort_lat
    type(ESMF_StateItem_Flag) :: itemflag
    character(len=*), parameter :: subname='(ocn_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    nproc = narea - 1

    ! Get import and export states
    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Realize import and export states
    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrOcn, &
         numflds=fldsFrOcn_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':NEMO_Export',&
         mesh=mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToOcn, &
         numflds=fldsToOcn_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':NEMO_Import',&
         mesh=mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine mesh lats and lons
    call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numownedelements))
    allocate(lonMesh(numOwnedElements))
    allocate(latMesh(numOwnedElements))
    allocate(lonModel(numOwnedElements))
    allocate(latModel(numOwnedElements))
    call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1,numOwnedElements
       lonMesh(n) = ownedElemCoords(2*n-1)
       latMesh(n) = ownedElemCoords(2*n)
    end do
    deallocate(ownedElemCoords)
    lonModel(:) = c0
    latModel(:) = c0

    ! Compare mesh lats/lons to model generated lats/lons
    labort_lon = .false.
    labort_lat = .false.
    n = 0
    do j=Njs0,Nje0  
      do i=Nis0,Nie0 
         n = n+1
         if (tmask(i,j,1)==1) then
             lonModel(n) = glamt(i,j)
             if (lonModel(n)<c0) then
               lonModel(n) = lonModel(n) + 360.0_wp
             end if
             latModel(n) = gphit(i,j)
             diff_lon = abs(lonMesh(n) - lonModel(n))
             if ( (diff_lon > 1.e2  .and. abs(diff_lon - 360.) > 1.e-1) .or.&
                  (diff_lon > 1.e-3 .and. diff_lon < c1) ) then
                write(numout,'(a,4(i6,3x),2(f21.13,3x),d21.5)') &
                     'ERROR: NEMO nproc, n, i, j, lonMesh, lonModel, diff_lon = ',&
                     nproc, n, i, j,lonMesh(n),lonModel(n), diff_lon
                labort_lon = .true.
             end if
             if (abs(latMesh(n) - latModel(n)) > 1.e-1) then
                write(numout,'(a,4(i6,3x),2(f21.13,3x),d21.5)') &
                     'ERROR: NEMO nproc, n, i, j, latMesh, latModel, diff_lat = ', &
                     nproc, n, i, j,latMesh(n),latModel(n), abs(latMesh(n)-latModel(n))
                labort_lat = .true.
             end if
          end if
       enddo
    enddo
    deallocate(lonMesh, latMesh, lonModel, latModel)

    n = 0
    if (labort_lon .or. labort_lat) then
       n = 1
       write(numout,'(a)') &
             'ERROR: NEMO nproc, numOwnedElements, jpi, jpj, Ni_0, Nj_0, Nis0, Nie0, Njs0, Nje0 = '
       write(numout,'(10(i8,1x))') &
              nproc, numOwnedElements, jpi, jpj, Ni_0, Nj_0, Nis0, Nie0, Njs0, Nje0  
       write(numout,'(a)') &
         '========================================================================================================'
       call flush(numout)
    endif
    if ( lk_mpp ) call mpp_max('ocn_comp_nuopc', n )
    if (n>0) then
       call mppsync()
       call shr_sys_abort(subname//': coordinates missmatch!')
    end if

    ! Determine mesh areas used in regridding
    lfield = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8 , meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridGetArea(lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(mesh_areas(numOwnedElements))
    mesh_areas(:) = dataptr(:)
    call ESMF_FieldDestroy(lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Determine flux correction factors (module variables)
    allocate(model_areas(numOwnedElements))
    allocate(mod2med_areacor(numOwnedElements))
    allocate(med2mod_areacor(numOwnedElements))
    mod2med_areacor(:) = 1._wp
    med2mod_areacor(:) = 1._wp
    n = 0
    do j=Njs0,Nje0
       do i=Nis0,Nie0
             n = n+1
             model_areas(n) = e1e2t(i,j)/(ra*ra)
             mod2med_areacor(n) = model_areas(n) / mesh_areas(n)
             med2mod_areacor(n) = mesh_areas(n) / model_areas(n)
       end do
    end do
    area = glob_sum('ocn_comp_nuopc', e1e2t(:,:))
    if (nproc==0) write(numout,'(a,1es28.19)') 'Global ocean area (m**2)', area

    min_mod2med_areacor = minval(mod2med_areacor)
    max_mod2med_areacor = maxval(mod2med_areacor)
    min_med2mod_areacor = minval(med2mod_areacor)
    max_med2mod_areacor = maxval(med2mod_areacor)
    call shr_mpi_max(max_mod2med_areacor, max_mod2med_areacor_glob, mpi_communicator_ocn)
    call shr_mpi_min(min_mod2med_areacor, min_mod2med_areacor_glob, mpi_communicator_ocn)
    call shr_mpi_max(max_med2mod_areacor, max_med2mod_areacor_glob, mpi_communicator_ocn)
    call shr_mpi_min(min_med2mod_areacor, min_med2mod_areacor_glob, mpi_communicator_ocn)

    if (nproc==0) then
       write(numout,'(2A,2g23.15,A )') trim(subname),' :  min_mod2med_areacor, max_mod2med_areacor ',&
            min_mod2med_areacor_glob, max_mod2med_areacor_glob, 'NEMO'
       write(numout,'(2A,2g23.15,A )') trim(subname),' :  min_med2mod_areacor, max_med2mod_areacor ',&
            min_med2mod_areacor_glob, max_med2mod_areacor_glob, 'NEMO'
    end if

    deallocate(model_areas)
    deallocate(mesh_areas)

!    call ESMF_StateGet(importState, 'Faxa_ndep', itemFlag, rc=rc)
!    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then
!       call ESMF_LogWrite(subname//' Faxa_ndep is in import state', ESMF_LOGMSG_INFO)
!    end if

  end subroutine ocn_realize_fields

  !==============================================================================
  subroutine ocn_import( importState, flds_scalar_name, ldiag_cpl, errorCode, rc )

    !-----------------------------------------------------------------------
    ! swnet  -- net short-wave heat flux                 (W/m2   )
    ! lwup   -- longwave radiation (up)                  (W/m2   )
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_State)   , intent(in)  :: importState
    character(len=*)   , intent(in)  :: flds_scalar_name
    logical            , intent(inout)  :: ldiag_cpl
    integer            , intent(out) :: errorCode
    integer            , intent(out) :: rc

    ! local variables
    character (lc)       :: label,  message
    real (wp)            :: work1(jpi,jpj)
    integer (i4)         :: i,j,k,n,ncol,nfld,nf  !iblock,nfld,nf
    real (wp)            :: m2percm2, gsum
    real (wp), pointer   :: dataptr1d(:)
    real (wp), pointer   :: dataptr2d(:,:)
    ! from mediator (virtual ocn)
    real (wp), pointer   :: foxx_swnet(:)
    real (wp), pointer   :: foxx_swnet_afracr(:)
    real (wp), pointer   :: foxx_taux(:)
    real (wp), pointer   :: foxx_tauy(:)
    real (wp), pointer   :: foxx_lwup(:)
    real (wp), pointer   :: foxx_sen(:)
    real (wp), pointer   :: foxx_lat(:)
    real (wp), pointer   :: foxx_evap(:)
    real (wp), pointer   :: foxx_rofl(:)
    real (wp), pointer   :: foxx_rofi(:)
    real (wp), pointer   :: so_duu10n(:)
    ! from atm
    real (wp), pointer   :: sa_pslv(:)
    real (wp), pointer   :: faxa_rain(:)
    real (wp), pointer   :: faxa_snow(:)
    real (wp), pointer   :: faxa_lwdn(:)
    real (wp), pointer   :: faxa_dstwet(:,:)
    real (wp), pointer   :: faxa_dstdry(:,:)
    real (wp), pointer   :: faxa_bcph(:,:)
    ! from ice
    real (wp), pointer   :: sf_afrac(:)
    real (wp), pointer   :: sf_afracr(:)
    real (wp), pointer   :: si_ifrac(:)
    real (wp), pointer   :: si_ifrac_n(:,:)
    real (wp), pointer   :: fioi_flxdst(:)
    real (wp), pointer   :: fioi_bcpho(:)
    real (wp), pointer   :: fioi_bcphi(:)
    real (wp), pointer   :: fioi_meltw(:)
    real (wp), pointer   :: fioi_melth(:)
    real (wp), pointer   :: fioi_salt(:)
    real (wp), pointer   :: fioi_swpen_ifrac_n(:,:)
    ! from wave
    !
    integer (i4)         :: fieldCount
    character (lc)       :: fldname
    type(ESMF_StateItem_Flag) :: itemflag
#ifdef _HIRES
    real (wp)            :: qsw_eps = -1.e-3_wp
#else
    real (wp)            :: qsw_eps = 0._wp
#endif
    character (lc), allocatable :: fieldNameList(:)
    character(len=*), parameter :: subname='(ocn_import_export:ocn_import)'
    real (wp)  :: sgn  
    real (wp)  :: rtmp1, rtmp2
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !-----------------------------------------------------------------------
    !  zero out padded cells
    !-----------------------------------------------------------------------

    nproc = narea - 1
    errorCode = 0
 
    ! lrecv should be reset to false in sbc_cpl_cesm_rcv
    if (lrecv) then
      errorCode = 1
      return
    end if

    work1 = c0
    !work2 = c0

    ! NOTE: if there are code changes associated with changing the names or
    !       the number of fluxes received from the coupler, then subroutine
    !       update_ghost_cells_coupler_fluxes will need to be modified also

    ! NOTE : RCALCT acts as a KMT mask, its 1 if KMT>=1 and 0 otherwise

    !-----------------------------------------------------------------------
    ! from mediator (virtual ocean)
    !-----------------------------------------------------------------------

    !  unpack and distribute wind stress, then convert to correct units
    !  and rotate components to local coordinates

    ! zonal wind and meridional wind stress  (W/m2)
    call state_getfldptr(importState, 'Foxx_taux', foxx_taux, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Foxx_tauy', foxx_tauy, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    taux_x2o(:,:) = c0
    tauy_x2o(:,:) = c0

    n = 0
    do j=Njs0,Nje0
       do i=Nis0,Nie0
             n = n + 1
             taux_x2o(i,j)  = foxx_taux(n) * med2mod_areacor(n)
             tauy_x2o(i,j)  = foxx_tauy(n) * med2mod_areacor(n)
       end do
    end do
!    ! rotate true zonal/meridional wind stress into local coordinates,
!    ! convert to dyne/cm**2, and shift SMFT to U grid
!    ! halo updates are performed in subroutine rotate_wind_stress,
!    ! following the rotation
!    call rotate_wind_stress(work1, work2)

    ! evaporation flux (kg/m2/s)
    call state_getfldptr(importState, 'Foxx_evap', foxx_evap, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! sensible heat flux (W/m2)
    call state_getfldptr(importState, 'Foxx_sen', foxx_sen, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! latent heat flux (W/m2)
    call state_getfldptr(importState, 'Foxx_lat', foxx_lat, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! long wave up flux  (W/m2)
    call state_getfldptr(importState, 'Foxx_lwup', foxx_lwup, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! shortwave net flux  (W/m2)
    call state_getfldptr(importState, 'Foxx_swnet', foxx_swnet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! 10m wind speed squared (m^2/s^2)
    call state_getfldptr(importState, 'So_duu10n', so_duu10n, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    evap_x2o(:,:)   = c0
    sen_x2o(:,:)    = c0
    lat_x2o(:,:)    = c0
    lwup_x2o(:,:)   = c0
    swnet_x2o(:,:)  = c0
    duu10n_x2o(:,:) = c0

    n = 0
    do j=Njs0,Nje0
       do i=Nis0,Nie0
             n = n + 1
             evap_x2o(i,j) = foxx_evap(n) * med2mod_areacor(n)
             sen_x2o(i,j) = foxx_sen(n) * med2mod_areacor(n)
             lat_x2o(i,j) = foxx_lat(n) * med2mod_areacor(n)
             lwup_x2o(i,j) = foxx_lwup(n) * med2mod_areacor(n)
             swnet_x2o(i,j) = foxx_swnet(n) * med2mod_areacor(n) 
             duu10n_x2o(i,j) = so_duu10n(n) 
       end do
    end do
    if (ANY(swnet_x2o < qsw_eps)) then
       do j=Njs0,Nje0
          do i=Nis0,Nie0
                write(6,*)'ERROR: j,i,swnet_x2o = ',&
                     j,i,swnet_x2o(i,j)
          enddo
       enddo
       call shr_sys_abort('(set_surface_forcing) ERROR: swnet_x2o < qsw_eps in set_surface_forcing')
    endif

    !-----------------------------------------------------------------------
    ! from atmosphere
    !-----------------------------------------------------------------------

    ! sea-level pressure (Pa)
    call state_getfldptr(importState, 'Sa_pslv', sa_pslv, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! snow
    call state_getfldptr(importState, 'Faxa_snow', faxa_snow, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! rain
    call state_getfldptr(importState, 'Faxa_rain', faxa_rain, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! longwave radiation (down) (W/m2)
    call state_getfldptr(importState, 'Faxa_lwdn', faxa_lwdn, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    pslv_x2o(:,:) = c0
    snow_x2o(:,:) = c0
    rain_x2o(:,:) = c0
    lwdn_x2o(:,:) = c0

    n = 0
    do j=Njs0,Nje0
       do i=Nis0,Nie0
             n = n+1
             pslv_x2o(i,j) = sa_pslv(n) 
             snow_x2o(i,j) = faxa_snow(n) * med2mod_areacor(n)
             rain_x2o(i,j) = faxa_rain(n) * med2mod_areacor(n) ! rain + snow
             lwdn_x2o(i,j) = faxa_lwdn(n) * med2mod_areacor(n)
       end do
    end do

    !-----------------------------------------------------------------------
    ! from sea-ice
    !-----------------------------------------------------------------------

    ! ice fraction
    call state_getfldptr(importState, 'Si_ifrac', si_ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! heat flux from sea ice snow & ice melt (W/m2)
    call state_getfldptr(importState, 'Fioi_melth', fioi_melth, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! salt from sea ice (kg(salt)/m2/s)
    call state_getfldptr(importState, 'Fioi_salt', fioi_salt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! snow melt flux from sea ice (kg/m2/s)
    call state_getfldptr(importState, 'Fioi_meltw', fioi_meltw, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ifrac_x2o(:,:) = c0
    melth_x2o(:,:) = c0
    meltw_x2o(:,:) = c0
    salt_x2o(:,:)  = c0

    n = 0
    do j=Njs0,Nje0
       do i=Nis0,Nie0
             n = n + 1
             ifrac_x2o(i,j) = si_ifrac(n) 
             melth_x2o(i,j) = fioi_melth(n) * med2mod_areacor(n)
             meltw_x2o(i,j) = fioi_meltw(n) * med2mod_areacor(n)
             salt_x2o(i,j) = fioi_salt(n) * med2mod_areacor(n)
       end do
    end do

    !-----------------------------------------------------------------------
    ! from wave
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! from river
    !-----------------------------------------------------------------------

    ! liquid runoff flux (kg/m2/s)
    call state_getfldptr(importState, 'Foxx_rofl', foxx_rofl, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! ice runoff flux (kg/m2/s)
    call state_getfldptr(importState, 'Foxx_rofi', foxx_rofi, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    roff_x2o(:,:) = c0
    ioff_x2o(:,:) = c0

    n = 0
    do j=Njs0,Nje0
       do i=Nis0,Nie0
             n = n+1
             roff_x2o(i,j) = foxx_rofl(n) * med2mod_areacor(n)
             ioff_x2o(i,j) = foxx_rofi(n) * med2mod_areacor(n)
       end do
    end do
    !roff_x2o = MAX(roff_x2o, 0.0_wp)
    !ioff_x2o = MAX(ioff_x2o, 0.0_wp)
    if (ANY(roff_x2o(:,:)*tmask(:,:,1) < c0)) then
       do j=Njs0,Nje0
          do i=Nis0,Nie0
             if (tmask(i,j,1)==1 .and. roff_x2o(i,j) < c0) then
                write(numout,*)'ERROR: j,i,roff_x2o = ',&
                     j,i,roff_x2o(i,j)
             end if
          enddo
       enddo
       call shr_sys_abort('(set_surface_forcing) ERROR: roff_x2o is negative')
    endif

    if (ANY(ioff_x2o(:,:)*tmask(:,:,1) < c0)) then
       do j=Njs0,Nje0
          do i=Nis0,Nie0
             if (tmask(i,j,1)==1 .and. ioff_x2o(i,j) < c0) then
                write(numout,*)'ERROR: j,i,ioff_x2o = ',&
                     j,i,ioff_x2o(i,j)
             end if
          enddo
       enddo
       call shr_sys_abort('(set_surface_forcing) ERROR: roff_x2o is negative')
    endif

    !-----------------------------------------------------------------------
    ! CO2 from atm
    !-----------------------------------------------------------------------

#if defined key_cpl_carbon_cycle
    call ESMF_StateGet(importState, 'Sa_co2prog', itemFlag, rc=rc)
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then
       call state_getimport(importState, 'Sa_co2prog', work1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       co2_x2o(:,:) = c0
       if (work1 > 0) then
          n = 0
          do j=Njs0,Nje0
            do i=Nis0,Nie0
              n = n + 1
              co2_x2o(i,j) = work1(i,j) 
            enddo
          enddo
       endif
    endif

    call ESMF_StateGet(importState, 'Sa_co2diag', itemFlag, rc=rc)
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then
       call state_getimport(importState, 'Sa_co2diag', work1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       co2_x2o(:,:) = c0
       if (work1 > 0) then
         n = 0
         do j=Njs0,Nje0
            do i=Nis0,Nie0
               n = n + 1
               co2_x2o(i,j) = work1(i,j) 
            enddo
         enddo
       endif
    endif
#endif

    ! coupling time step flag
    lrecv = .TRUE.

    !-----------------------------------------------------------------------
    !  diagnostics
    !-----------------------------------------------------------------------

    if (ldiag_cpl) then
      if (lwp) write(numout,*) 'nemo_recv_from_coupler'

       ! Determine all field names in import state
       call ESMF_StateGet(importState, itemCount=fieldCount, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(fieldNameList(fieldCount))
       call ESMF_StateGet(importState, itemNameList=fieldNameList, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! loop over all fields in import state
       ! from atm - black carbon deposition fluxes (3)
       ! (1) => Faxa_bcphidry, (2) => Faxa_bcphodry, (3) => Faxa_bcphiwet
       ! from atm - organic carbon deposition fluxes (3)
       ! (1) => Faxa_dstwet1, (2) => Faxa_dstwet2, (3) => Faxa_dstwet3, (4) => Faxa_dstwet4
       ! from atm - dry dust deposition frluxes (4 sizes)
       ! (1) => Faxa_dstdry1, (2) => Faxa_dstdry2, (3) => Faxa_dstdry3, (4) => Faxa_dstdry4

!       m2percm2  = mpercm*mpercm
       do nfld = 1, fieldCount
           if (fieldNameList(nfld) == trim(flds_scalar_name)) then
             CYCLE
           else
             call state_getfldptr(importState, trim(fieldNameList(nfld)), dataPtr1d, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             n = 0
             do j = Njs0,Nje0 !1, jpj !Njs0,Nje0
                do i = Nis0,Nie0
                    n = n + 1
                    work1(i,j) = dataPtr1d(n)
                end do
             end do
             sgn = 1._wp
             if (trim(fieldNameList(nfld))=='Foxx_taux' .or. trim(fieldNameList(nfld))=='Foxx_tauy') then
                sgn = -1._wp
             end if

             gsum = glob_sum('ocn_comp_nuopc', WORK1(:,:)*e1e2t(:,:))

             rtmp1=minval(WORK1, mask=(tmask(:,:,1)>0._wp))
             call mpp_min('ocn_comp_nuopc',rtmp1)
             rtmp2=maxval(WORK1, mask=(tmask(:,:,1)>0._wp))
             call mpp_max('ocn_comp_nuopc',rtmp2)

             if (nproc==0) then
                 write(numout,1100)'ocn','send', trim(fieldNameList(nfld)), gsum, gsum/area, &
                   rtmp1, rtmp2
             endif
           end if
       end do
       if (nproc==0) call shr_sys_flush(numout) 
    end if

1100 format ('comm_diag ', a3, 1x, a4, 1x, a16, 1x, 4es28.19:)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ocn_import

  !==============================================================================
  subroutine ocn_export(exportState, flds_scalar_name, ldiag_cpl, errorCode, rc)

    !-----------------------------------------------------------------------
    ! Create export state
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_State)                 :: exportState
    character(len=*)   , intent(in)  :: flds_scalar_name
    logical            , intent(inout)  :: ldiag_cpl
    integer (i4)       , intent(out) :: errorCode  ! pop error code
    integer            , intent(out) :: rc         ! returned error code

    ! local variables
    integer (i4)         :: n,i,j,nfld !,lev
    character (len=lc)   :: label
    real (wp)            :: work1(jpi,jpj)
    real (wp)            :: work2(jpi,jpj)
    real (wp)            :: work3(jpi,jpj)
    real (wp)            :: work4(jpi,jpj)
    real (wp)            :: worka(jpi,jpj)
    real (wp)            :: sgn
    real (wp)            :: gsum
    real (wp)            :: rtmp1, rtmp2
    real (wp), pointer   :: dataptr1(:)
    real (wp), pointer   :: dataptr2(:)
    real (wp), pointer   :: dataptr2d(:,:)
    integer (i4)         :: fieldCount
    character (len=lc), allocatable :: fieldNameList(:)
    character(len=*), parameter :: subname='(ocn_import_export:ocn_export)'
    logical              :: l_export
    logical, save        :: lfirst = .true.
    integer, dimension(jpi,jpj) :: iktop, ikbot          ! sea surface gradient
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS
    errorCode = 0 !POP_Success

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    nproc = narea - 1
    l_export = .true.
!    if (present(lexport)) l_export = lexport

    !-----------------------------------------------------------------------
    ! ocean mask
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_omask', dataPtr1, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    dataptr1(:) = c0
    n = 0
    do j=Njs0,Nje0
       do i=Nis0,Nie0
          n = n+1
          dataptr1(n) = tmask_i(i,j)
       enddo
    enddo

    !-----------------------------------------------------------------------
    ! interpolate onto T-grid points and rotate on T grid
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_u', dataPtr1, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'So_v', dataPtr2, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    dataptr1(:) = shr_const_spval
    dataptr2(:) = shr_const_spval

    work1(:,:) = c0
    work2(:,:) = c0
    work3(:,:) = c0
    work4(:,:) = c0

    ! Mask & average
    sbuff_sum_u(:,:) = sbuff_sum_u(:,:)*umask(:,:,1)/tlast_coupled
    sbuff_sum_v(:,:) = sbuff_sum_v(:,:)*vmask(:,:,1)/tlast_coupled
    ! Apply LBC
    call lbc_lnk('ocn_comp_nuopc', sbuff_sum_u(:,:), 'U', -1._wp, sbuff_sum_v(:,:), 'V', -1._wp)
    ! (U,V) -> T
    do j = 2, jpj
       do i = 2, jpi
          rtmp1 = umask(i,j,1)+umask(i-1,j,1)
          if (rtmp1>c0) then
             work1(i,j) = ( sbuff_sum_u(i,j)*umask(i,j,1) + &
                sbuff_sum_u(i-1,j)*umask(i-1,j,1) ) / rtmp1
          end if
          rtmp2 = vmask(i,j,1)+vmask(i,j-1,1)
          if (rtmp2>c0) then
             work2(i,j) = ( sbuff_sum_v(i,j)*vmask(i,j,1) + &
                sbuff_sum_v(i,j-1)*vmask(i,j-1,1) ) / rtmp2
          end if
       end do
    end do
    work1(:,:) = work1(:,:)*tmask(:,:,1)
    work2(:,:) = work2(:,:)*tmask(:,:,1)
    ! Apply LBC
    call lbc_lnk('ocn_comp_nuopc', work1, 'T', -1._wp, work2, 'T', -1._wp)
 
    if (.not. lfirst) then
       call iom_put("So_u_o2x", work1)
       call iom_put("So_v_o2x", work2)
    endif
 
    call rot_rep( work1, work2, 'T', 'ij->e', work3 )
    call rot_rep( work1, work2, 'T', 'ij->n', work4 )
 
    work1(:,:) = work3(:,:)*tmask(:,:,1)
    work2(:,:) = work4(:,:)*tmask(:,:,1)
 
    if (l_export) then
    n = 0
    do j=Njs0,Nje0
       do i=Nis0,Nie0
          n = n + 1
          dataptr1(n) = work1(i,j)
          dataptr2(n) = work2(i,j)
       enddo
    enddo
    end if


    !-----------------------------------------------------------------------
    ! convert and pack surface temperature
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_t', dataptr1, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    dataptr1(:) = shr_const_spval
 
    work1(:,:) = (sbuff_sum_t(:,:)/tlast_coupled + rt0)*tmask(:,:,1)

    if (l_export) then
    n = 0
    do j=Njs0,Nje0
       do i=Nis0,Nie0
          n = n + 1
          dataptr1(n) = work1(i,j)
       enddo
    enddo
    end if
! TO BE CORRECTED
!    if (ocn2glc_coupling) then
!       call state_getfldptr(exportState, 'So_t_depth', dataptr2d, rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!       dataptr2d(:,:) = c0
!       n = 0
!       call dom_zgr(iktop, ikbot)
!       do j = 1, jpj !1, jpj !Njs0,Nje0 !this_block%jb,this_block%je
!          do i = 1, jpi !1, jpi !Nis0,Nie0 !this_block%ib,this_block%ie
!             n = n + 1
!             do lev = 1,num_ocn2glc_levels
!                 if (ikbot(i,j) >= ocn2glc_levels(lev)) then
!                      dataptr2d(lev,n) = (sbuff_sum_t_depth(i,j,lev)/tlast_coupled + rt0)*tmask(:,:,1) 
!                 end if
!             end do
!          end do
!       end do
!    endif
    if (.not. lfirst) then
      call iom_put('So_t_o2x', WORK1)
    endif

    !-----------------------------------------------------------------------
    ! convert and pack salinity
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_s', dataptr1, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    dataptr1(:) = shr_const_spval
    work1(:,:) = sbuff_sum_s(:,:)*tmask(:,:,1)/tlast_coupled

    if (l_export) then
    n = 0
    do j=Njs0,Nje0
       do i=Nis0,Nie0
          n = n + 1
          dataptr1(n) = work1(i,j)
       enddo
    enddo
    end if
! TO BE CORRECTED 
!    if (ocn2glc_coupling) then
!       call state_getfldptr(exportState, 'So_s_depth', dataptr2d, rc)
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!       dataptr2d(:,:) = c0
!       n = 0
!       call dom_zgr(iktop, ikbot)
!       do j = 1, jpj !1, jpj !Njs0,Nje0 !this_block%jb,this_block%je
!          do i = 1, jpi !1, jpi !Nis0,Nie0 !this_block%ib,this_block%ie
!             n = n + 1
!             do lev = 1,num_ocn2glc_levels
!                 if (ikbot(i,j) >= ocn2glc_levels(lev)) then
!                      dataptr2d(lev,n) = sbuff_sum_s_depth(i,j,lev)*tmask(:,:,1)/tlast_coupled
!                 end if
!             end do
!          end do
!       end do
!    end if

    if (.not. lfirst) then
       call iom_put('So_s_o2x', work1)
    endif

    !-----------------------------------------------------------------------
    ! interpolate onto T-grid points and rotate on T grid
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_dhdx', dataPtr1, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'So_dhdy', dataPtr2, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    dataptr1(:) = shr_const_spval
    dataptr2(:) = shr_const_spval

    work1(:,:) = c0
    work2(:,:) = c0
    work3(:,:) = c0
    work4(:,:) = c0

    ! Apply LBC
    sbuff_sum_dhdx(:,:) = sbuff_sum_dhdx(:,:)*umask(:,:,1)/tlast_coupled
    sbuff_sum_dhdy(:,:) = sbuff_sum_dhdy(:,:)*vmask(:,:,1)/tlast_coupled
    call lbc_lnk('ocn_comp_nuopc', sbuff_sum_dhdx(:,:), 'U', -1._wp, sbuff_sum_dhdy(:,:), 'V', -1._wp)
    ! (U,V) -> T
    do j = 2, jpj
       do i = 2, jpi
          rtmp1 = umask(i,j,1)+umask(i-1,j,1)
          if (rtmp1>c0) then
             work1(i,j) = ( sbuff_sum_dhdx(i,j)*umask(i,j,1) + &
                sbuff_sum_dhdx(i-1,j)*umask(i-1,j,1) ) / rtmp1
          end if
          rtmp2 = vmask(i,j,1)+vmask(i,j-1,1)
          if (rtmp2>c0) then
             work2(i,j) = ( sbuff_sum_dhdy(i,j)*vmask(i,j,1) + &
                sbuff_sum_dhdy(i,j-1)*vmask(i,j-1,1) ) / rtmp2
          end if
       end do
    end do
    work1(:,:) = work1(:,:)*tmask(:,:,1)
    work2(:,:) = work2(:,:)*tmask(:,:,1)
    ! Apply LBC
    call lbc_lnk('ocn_comp_nuopc', work1, 'T', -1._wp, work2, 'T', -1._wp)
 
    if (.not. lfirst) then
       call iom_put("So_dhdx_o2x", work1)
       call iom_put("So_dhdy_o2x", work2)
    endif
 
    call rot_rep( work1, work2, 'T', 'ij->e', work3 )
    call rot_rep( work1, work2, 'T', 'ij->n', work4 )
 
    work1(:,:) = work3(:,:)*tmask(:,:,1)
    work2(:,:) = work4(:,:)*tmask(:,:,1)
 
    if (l_export) then
    n = 0
    do j=Njs0,Nje0
       do i=Nis0,Nie0
          n = n + 1
          dataptr1(n) = work1(i,j)
          dataptr2(n) = work2(i,j)
       enddo
    enddo
    end if

    !-----------------------------------------------------------------------
    ! pack heat flux due to freezing/melting (W/m^2)
    ! QFLUX computation and units conversion occurs in ice.F
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'Fioo_q', dataPtr1, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    dataptr1(:) = shr_const_spval
    if (l_export) then
    n = 0
    do j=Njs0,Nje0
       do i=Nis0,Nie0
          n = n + 1
          dataptr1(n) = QFLUX(i,j) * mod2med_areacor(n) 
       enddo
    enddo
    end if

    if (.not. lfirst) then
       call iom_put('So_qflux_o2x', MAX(c0, QFLUX))
    endif

    tlast_ice = c0
    AQICE     = c0
    QFLUX     = c0

    !-----------------------------------------------------------------------
    ! pack co2 flux, if requested (kg CO2/m^2/s)
    ! units conversion occurs where co2 flux is computed
    !-----------------------------------------------------------------------

    if ( State_FldChk(exportState, 'Faoo_fco2_ocn') .and. l_export) then
       call state_getfldptr(exportState, 'Faoo_fco2_ocn', dataPtr1, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       dataptr1(:) = shr_const_spval

       n = 0
       do j=Njs0,Nje0
         do i=Nis0,Nie0
            n = n + 1
            dataptr1(n) = (sbuff_sum_co2(i,j)*tmask(i,j,1)/tlast_coupled) * mod2med_areacor(n) 
         enddo
      enddo
    endif

    !-----------------------------------------------------------------------
    ! diagnostics
    !-----------------------------------------------------------------------

    if (ldiag_cpl .and. l_export) then
       if (lwp) write(numout,*)'nemo_send_to_mediator'

       ! Determine all field names in export state
       call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(fieldNameList(fieldCount))
       call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! loop over all fields in export state
       do nfld = 1, fieldCount
          if (trim(fieldNameList(nfld)) /= flds_scalar_name) then
             call state_getfldptr(exportState, trim(fieldNameList(nfld)), dataptr1, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             n = 0
             workA(:,:) = c0
             do j=Njs0,Nje0
                do i=Nis0,Nie0
                  n = n + 1
                  workA(i,j) = dataptr1(n)
                enddo
             enddo

             sgn = 1._wp
             if (trim(fieldNameList(nfld))=='So_u'    .or. trim(fieldNameList(nfld))=='So_v'     .or. &
                 trim(fieldNameList(nfld))=='So_dhdx' .or. trim(fieldNameList(nfld))=='So_dhdy') then
                sgn = -1._wp
             else if (trim(fieldNameList(nfld))=='Fioo_q') then
                workA(:,:)=MAX(c0, workA(:,:))
             end if
     
             gsum = glob_sum('ocn_comp_nuopc', workA(:,:)*e1e2t(:,:))
     
             rtmp1=minval(workA, mask=(tmask(:,:,1)>0._wp))
             call mpp_min('ocn_comp_nuopc',rtmp1)
             rtmp2=maxval(workA, mask=(tmask(:,:,1)>0._wp))
             call mpp_max('ocn_comp_nuopc',rtmp2)

             if (nproc==0) then
                 write(numout,1100)'ocn','send', trim(fieldNameList(nfld)), gsum, gsum/area, &
                   rtmp1, rtmp2
             endif
             
          end if
       end do
       if (nproc==0) call shr_sys_flush(numout)
    end if
1100 format ('comm_diag ', a3, 1x, a4, 1x, a16, 1x, 4es28.19:)

    tlast_coupled = c0
    lfirst = .false.

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ocn_export

  !===============================================================================
  subroutine fldlist_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    integer             , intent(inout) :: num
    type(fld_list_type) , intent(inout) :: fldlist(:)
    character(len=*)    , intent(in)    :: stdname
    integer, optional   , intent(in)    :: ungridded_lbound
    integer, optional   , intent(in)    :: ungridded_ubound

    ! local variables
    character(len=*), parameter :: subname='(fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call shr_sys_abort(trim(subname)//": ERROR num > fldsMax "//trim(stdname))
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

  end subroutine fldlist_add

  !===============================================================================
  subroutine fldlist_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    use NUOPC, only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(ocn_import_export:fldlist_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO)

             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

          else
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO)
             ! Create the field
             if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                                         ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                                         ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                                         gridToFieldMap=(/2/), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
          end if ! if not scalar field

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(fldlist_realize:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc) ! num of scalar values
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fldlist_realize

  !==============================================================================
  subroutine state_getimport(state, fldname, output, areacor, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(in)    :: state
    character(len=*)    , intent(in)    :: fldname
    real (wp)           , intent(inout) :: output(:,:) !output(:,:,:)
    real(wp) , optional , intent(in)    :: areacor(:)
    integer             , intent(out)   :: rc

    ! local variables
!    type(block)       :: this_block         ! block information for current block
    integer           :: i, j, n  !iblock, n   ! incides
    real(wp), pointer :: dataPtr1d(:)
    character(len=*), parameter :: subname='(import_export:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call state_getfldptr(state, trim(fldname), dataptr1d, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    n = 0
    if (present(areacor)) then
       do j=Njs0,Nje0
          do i=Nis0,Nie0
                n = n + 1
                output(i,j) = dataPtr1d(n) * areacor(n)
          end do
       end do
    else
       do j=Njs0,Nje0
          do i=Nis0,Nie0
                n = n + 1
                output(i,j) = dataPtr1d(n)
          end do
       end do
    end if
!    end do

  end subroutine state_getimport

  !===============================================================================
  subroutine State_GetFldPtr_1d(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get 1d pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)  , intent(in)     :: State
    character(len=*)  , intent(in)     :: fldname
    real(wp), pointer , intent(inout)  :: fldptr(:)
    integer, optional , intent(out)    :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    character(len=*),parameter :: subname='(ocn_import_export:State_GetFldPtr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine State_GetFldPtr_1d

  !===============================================================================
  subroutine State_GetFldPtr_2d(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get 2d pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)  , intent(in)     :: State
    character(len=*)  , intent(in)     :: fldname
    real(wp), pointer , intent(inout)  :: fldptr(:,:)
    integer, optional , intent(out)    :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    character(len=*),parameter :: subname='(ocn_import_export:State_GetFldPtr_2d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine State_GetFldPtr_2d

  !===============================================================================
  logical function State_FldChk(State, fldname)
    ! ----------------------------------------------
    ! Determine if field is in state
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State) , intent(in)  :: State
    character(len=*) , intent(in)  :: fldname

    ! local variables
    type(ESMF_StateItem_Flag) :: itemType
    ! ----------------------------------------------

    call ESMF_StateGet(State, trim(fldname), itemType)
    State_FldChk = (itemType /= ESMF_STATEITEM_NOTFOUND)

  end function State_FldChk

  !===============================================================================

  subroutine ocn_sum_buffer(exportState, rc)

    ! ----------------------------------------------
    ! Accumulates sums for averaging fields to be sent to mediator
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)     :: exportState
    integer, intent(out) :: rc

    ! local variables
!    type (block)             :: this_block ! local block info
    real (wp), dimension(jpi,jpj) :: work ! local work arrays
    real (wp), dimension(jpi,jpj) :: ssgu, ssgv          ! sea surface gradient
    integer,   dimension(jpi,jpj) :: iktop, ikbot          ! sea surface gradient
    real (wp)                :: delt                                      ! time interval since last step
    integer (i4)             :: sflux_co2_nf_ind = 0                      ! named field index of fco2
    integer                  :: ji,jj                                 ! indices
    logical, save            :: first = .true.                            ! only true for first call
    integer                  :: n, lev

    !-----------------------------------------------------------------------
    ! zero buffer if this is the first time after a coupling interval
    !-----------------------------------------------------------------------

    rc = 0

    if (tlast_coupled == c0) then
       if (.not. allocated(sbuff_sum_u)) then
          allocate(sbuff_sum_u (jpi,jpj))
       end if
       sbuff_sum_u    (:,:) = c0
       if (.not. allocated(sbuff_sum_v)) then
          allocate(sbuff_sum_v (jpi,jpj))
       end if
       sbuff_sum_v    (:,:) = c0
       if (.not. allocated(sbuff_sum_t)) then
          allocate(sbuff_sum_t (jpi,jpj))
       end if
       sbuff_sum_t    (:,:) = c0
       if (.not. allocated(sbuff_sum_s)) then
          allocate(sbuff_sum_s (jpi,jpj))
       end if
       sbuff_sum_s    (:,:) = c0
       if (.not. allocated(sbuff_sum_dhdx)) then
          allocate(sbuff_sum_dhdx (jpi,jpj))
       end if
       sbuff_sum_dhdx (:,:) = c0
       if (.not. allocated(sbuff_sum_dhdy)) then
          allocate(sbuff_sum_dhdy (jpi,jpj))
       end if
       sbuff_sum_dhdy (:,:) = c0
!       sbuff_sum_bld  (:,:,:) = c0
       if (.not. allocated(sbuff_sum_co2)) then
          allocate(sbuff_sum_co2 (jpi,jpj))
       end if
       sbuff_sum_co2  (:,:) = c0
       !if (.not. allocated(sbuff_sum_t_depth)) then
       !   allocate(sbuff_sum_t_depth (jpi,jpj,num_ocn2glc_levels))
       !end if
       !sbuff_sum_t_depth(:,:,:) = c0
       !if (.not. allocated(sbuff_sum_s_depth)) then
       !   allocate(sbuff_sum_s_depth (jpi,jpj,num_ocn2glc_levels))
       !end if
       !sbuff_sum_s_depth(:,:,:) = c0
    end if

    !-----------------------------------------------------------------------
    ! update time since last coupling
    !-----------------------------------------------------------------------

    delt = rn_Dt
    tlast_coupled = tlast_coupled + delt

    !-----------------------------------------------------------------------
    ! allow for fco2 field to not be registered on first call
    !    because init_forcing is called before init_passive_tracers
    ! use weight from previous timestep because flux used here is that
    !    computed during the previous timestep
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  accumulate sums of U,V,T,S and GRADP
    !  accumulate sum of co2 flux, if requested
    !     implicitly use zero flux if fco2 field not registered yet
    !  ice formation flux is handled separately in ice routine
    !-----------------------------------------------------------------------

    ssgu(:,:) = c0
    ssgv(:,:) = c0
    DO jj = 1, jpj-1             ! Sea surface gradient (now)
       DO ji = 1, jpi-1
            ssgu(ji,jj) = ( ssh(ji+1,jj,Nnn) - ssh(ji,jj,Nnn) ) / e1u(ji,jj)
            ssgv(ji,jj) = ( ssh(ji,jj+1,Nnn) - ssh(ji,jj,Nnn) ) / e2v(ji,jj)
       END DO
    END DO 
    
    sbuff_sum_u    (:,:) = sbuff_sum_u   (:,:) + delt * uu(:,:,1,Nnn) !UVEL(:,:,1,Nnn) !curtime,iblock)
    sbuff_sum_v    (:,:) = sbuff_sum_v   (:,:) + delt * vv(:,:,1,Nnn) !VVEL(:,:,1,Nnn) !curtime,iblock)
    WORK(:,:) = c0
    if ( ln_useCT ) then
        WORK(:,:) = eos_pt_from_ct(ts(:,:,1,jp_tem,Nnn),ts(:,:,1,jp_sal,Nnn))
    else
        WORK(:,:) = ts(:,:,1,jp_tem,Nnn)
    endif
    sbuff_sum_t    (:,:) = sbuff_sum_t   (:,:) + delt * WORK(:,:)
    sbuff_sum_s    (:,:) = sbuff_sum_s   (:,:) + delt * ts(:,:,1,jp_sal,Nnn) !TRACER(:,:,1,2,curtime,iblock)
    sbuff_sum_dhdx (:,:) = sbuff_sum_dhdx(:,:) + delt * ssgu(:,:)
    sbuff_sum_dhdy (:,:) = sbuff_sum_dhdy(:,:) + delt * ssgv(:,:)
!TO BE CORRECTED 
!    if (ocn2glc_coupling) then
!        call dom_zgr(iktop, ikbot)
!          do jj = 1, jpj !this_block%jb,this_block%je
!             do ji = 1, jpi !this_block%ib,this_block%ie
!                do n = 1,num_ocn2glc_levels
!                   lev = ocn2glc_levels(n)
!                   if (ikbot(ji,jj) >= lev) then
!                      sbuff_sum_t_depth(ji,jj,n) = sbuff_sum_t_depth(ji,jj,n) + &
!                           delt * ts(:,:,lev,jp_tem,Nnn)
!                      sbuff_sum_s_depth(ji,jj,n) = sbuff_sum_s_depth(ji,jj,n) + &
!                           delt * ts(:,:,lev,jp_sal,Nnn)
!                   end if
!                end do
!             end do
!          end do
!    end if

    first = .false.

  end subroutine ocn_sum_buffer

end module ocn_comp_nuopc

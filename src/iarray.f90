MODULE iarray
!
!  Contains all array var names.
!
    implicit none

    real, allocatable :: xface(:), yface(:), zface(:)
    real, allocatable :: rhokap(:,:,:)
    real, allocatable :: albedoar(:,:,:)
    real, allocatable :: jmean(:,:,:), jmeanGLOBAL(:,:,:), jme(:,:,:)
    real, allocatable :: refrac(:,:,:)

    real, allocatable :: si_step(:,:,:), si_step_GLOBAL(:,:,:)
    !real, allocatable :: si_step_flu(:,:,:), si_step_fluGLOBAL(:,:,:)
    real, allocatable ::  ua_ppix(:,:,:), ua_brain(:,:,:)

    real, allocatable:: percent_left(:), so_tot(:), o23_tot(:), o21_tot(:)

    real,allocatable:: white_matter(:,:,:), grey_matter(:,:,:), tumour(:,:,:), tumour_killed(:,:,:), tumour_resec(:,:,:)
    real,allocatable:: glial_matter(:,:,:), balloon(:,:,:),balloon_killed(:,:,:)

    real,allocatable::rhokap_gbm(:,:,:), albedoar_gbm(:,:,:), hggar(:,:,:)

    real,allocatable:: S_0(:,:,:), O2_3(:,:,:), O2_1(:,:,:)

    real,allocatable:: S0_slice(:,:,:), o23_slice(:,:,:), o21_slice(:,:,:), tumkill_slice(:,:,:), temp_slice(:,:,:)

    real,allocatable:: blob(:,:,:), blob_killed(:,:,:)

    real,allocatable:: temp(:,:,:), E_abs(:,:,:), temp_change(:,:,:)






end MODULE iarray

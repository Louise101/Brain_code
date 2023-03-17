MODULE iarray
!
!  Contains all array var names.
!
    implicit none

    real, allocatable :: xface(:), yface(:), zface(:)
    real, allocatable :: rhokap(:,:,:)
    real, allocatable :: albedoar(:,:,:)
    real, allocatable :: jmean(:,:,:), jmeanGLOBAL(:,:,:), jme(:,:,:), jm_av(:)
    real, allocatable :: jmean_flu(:,:,:), jmean_fluGLOBAL(:,:,:)
    real, allocatable :: refrac(:,:,:)

    real, allocatable :: si_step(:,:,:), si_step_GLOBAL(:,:,:)
    !real, allocatable :: si_step_flu(:,:,:), si_step_fluGLOBAL(:,:,:)
    real, allocatable :: nop_tot(:,:,:), time_tot(:,:,:)
    real, allocatable :: con_ppix(:,:,:), ua_ppix(:,:,:),con_test(:),con_test_GLOBAL(:)
    real, allocatable ::  p(:,:,:)
    integer, allocatable :: np_cell(:,:,:)

    real, allocatable :: flu_image(:,:)
    real, allocatable :: esc_flur(:), esc_flur_GLOBAL(:) ! array to bin escaped fluorescence
    real, allocatable ::alb_absor(:), pl_absor(:), pl_absor_GLOBAL(:) !arrays to compare the absorption vs depth for the albedo absorption and the path length absoption calculators
    real, allocatable:: phot_tot(:) !array to count total photons leaving the grid
    real,allocatable:: flu_phot_rel(:)
    real, allocatable:: flu_lay(:)
    real,allocatable:: obs_flur(:), obs_flur_GLOBAL(:)
    real,allocatable:: av_ppix_con_630(:,:) !array of average ppix concentration over each layer for each time step - used in Jaques validation of 420nm where 630nm does the photobleaching
    real,allocatable:: del_Q(:,:,:), del_pdd(:,:,:),PDD(:,:,:), PDD_GLOBAL(:,:,:)

    real,allocatable:: white_matter(:,:,:), grey_matter(:,:,:), tumour(:,:,:)
end MODULE iarray

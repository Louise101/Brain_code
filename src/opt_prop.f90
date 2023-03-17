MODULE opt_prop

    implicit none

    real    :: mua, mus, g2, hgg, kappa, albedo, wave, n1, n2
    real    :: e420, e630
    real    :: whitematter_mua, whitematter_mus, whitematter_g2, &
     whitematter_hgg, whitematter_kappa, whitematter_albedo
    real    :: greymatter_mua, greymatter_mus, greymatter_g2, greymatter_hgg, &
     greymatter_kappa, greymatter_albedo
    real*8:: csf_mua, csf_mus, csf_g2, csf_hgg, csf_kappa, csf_albedo
    real    :: glial_mua, glial_mus, glial_g2, glial_hgg, glial_kappa, glial_albedo
      real    :: gbm_mua, gbm_mus, gbm_g2, gbm_hgg, gbm_kappa, gbm_albedo
      real    :: mua_intra, mus_intra, g2_intra, hgg_intra, kappa_intra, albedo_intra
      real    :: mua_water, mus_water, g2_water, hgg_water, kappa_water, albedo_water

end MODULE opt_prop

MODULE ch_opt

implicit none

CONTAINS

  subroutine init_opt_630
    !brain optical properties at 630 nm

    use opt_prop
    use iarray
    use constants, only : nxg, nyg, nzg

    implicit none

    integer :: i,j,k

    !633nm

     hgg=0.9
    !white matter
    whitematter_hgg =0.85!0.84d0
    !g2  = hgg**2.
    whitematter_mua = 0.63!0.8d0
    whitematter_mus = 686.d0!409.d0

    whitematter_kappa  = whitematter_mus + whitematter_mua
    whitematter_albedo = whitematter_mus / whitematter_kappa

    !grey matter
    greymatter_hgg =0.85!0.89d0
    !greymatter_g2  = greymatter_hgg**2.
    greymatter_mua =0.99!0.2d0
    greymatter_mus =202.d0! 90.d0

    greymatter_kappa  = greymatter_mus + greymatter_mua
    greymatter_albedo = greymatter_mus / greymatter_kappa



  !CSF !660nm - dont have 633nm
    csf_hgg = 0.9d0
    csf_g2  = csf_hgg**2.
    csf_mua = 0.04d0
    csf_mus = 0.35d0

    csf_kappa  = csf_mus + csf_mua
    csf_albedo = csf_mus / csf_kappa

    !glial matter !using whtematter properties
    glial_hgg =0.85! 0.84d0
    glial_g2  = glial_hgg**2.
    glial_mua = 0.63!0.8d0
    glial_mus =686.d0! 409.d0

    glial_kappa  = glial_mus + glial_mua
    glial_albedo = glial_mus / glial_kappa

    !GBM cells
    gbm_hgg = 0.85!0.9d0
    !gbm_g2  = glial_hgg**2.
    gbm_mua = 1.3!2.1d0 !0.2
    gbm_mus = 218!34.9d0/(1.d0 - gbm_hgg)!160

    gbm_kappa  = gbm_mus + gbm_mua
    gbm_albedo = gbm_mus / gbm_kappa

    !salt water 630nm
    hgg_water = 0.9
    g2_water  = hgg_water**2.
    mua_water =0.003 !exploring the seafloor - feb 2017 - Garcia et al.
    mus_water =0.003

    kappa_water  = mus_water + mua_water
    albedo_water = mus_water / kappa_water


    !intralipid fluid
    hgg_intra = 0.875d0 !dupont parallel mcrt paper 2019
    g2_intra  = hgg_intra**2.
    mua_intra =0.001
    mus_intra =10.

    kappa_intra  = mus_intra + mua_intra
    albedo_intra = mus_intra / kappa_intra

   e630=0.0265
    ! calculate new ppix absorbtion coefficent and set new rhokap and albedo for each voxel
    do i= 1, nxg
      do j=1,nyg
        do k=1,nzg

          ua_ppix(i,j,k)=e630 * S_0(i,j,k)*tumour_resec(i,j,k) !ensures PpIX only present in tumour cells


          rhokap_gbm(i,j,k)=gbm_kappa + ua_ppix(i,j,k)
          albedoar_gbm(i,j,k)= gbm_mus / rhokap_gbm(i,j,k)

         end do
       end do
     end do


  end subroutine init_opt_630

  subroutine opt_jacques
! 630 nm
  use opt_prop
  use iarray, only :S_0, rhokap, albedoar, ua_ppix
  use constants, only : nxg, nyg, nzg

  implicit none

  integer:: i,j,k

  hgg = 0.9 !caclulated from paper using eq in fig 3 and tumour properties
  g2  = hgg**2.


  mua = 0.23
  mus = 21.d0/(1.d0 - hgg)

  kappa  = mus + mua
  albedo = mus/ kappa

  do i= 1, nxg
    do j=1,nyg
      do k=1,nzg

        ua_ppix(i,j,k)=e630 * S_0(i,j,k)

        rhokap(i,j,k)=mus + mua + ua_ppix(i,j,k)
        albedoar(i,j,k)= mus / rhokap(i,j,k)

       end do
     end do
   end do



end subroutine opt_jacques

subroutine opt_wang_m3
! 630 nm
use opt_prop
use iarray, only :S_0, rhokap, albedoar, ua_ppix
use constants, only : nxg, nyg, nzg

implicit none

integer:: i,j,k

hgg = 0.9 !caclulated from paper using eq in fig 3 and tumour properties
g2  = hgg**2.


mua = 1.51
mus = 12.29/(1.d0 - hgg)

kappa  = mus + mua
albedo = mus/ kappa

do i= 1, nxg
  do j=1,nyg
    do k=1,nzg

      ua_ppix(i,j,k)=e630 * S_0(i,j,k)

      rhokap(i,j,k)=mus + mua + ua_ppix(i,j,k)
      albedoar(i,j,k)= mus / rhokap(i,j,k)

     end do
   end do
 end do



end subroutine opt_wang_m3




   subroutine init_opt3
!
!  subroutine to set tissue optical properties 705nm
!

   use opt_prop

   implicit none

   hgg = 0.82445
   g2  = hgg**2.
   mua = 0.23d0
   mus = 17.d0/(1.d0 - hgg)

   kappa  = mus + mua
   albedo = mus / kappa


   end subroutine init_opt3


   subroutine sample(array, size_of, cdf, wave, iseed)
!
!  samples a random value from an array based upon its cdf
!
      implicit none

      integer, intent(IN)    :: iseed, size_of
      real,    intent(IN)    :: array(size_of, 2), cdf(size_of)
      real,    intent(OUT)   :: wave

      real :: ran2, value
      integer :: nlow

      value = ran2(iseed)

      call search_1D(size(cdf), cdf, nlow, value)
      call lin_inter_1D(array, cdf, value, size(cdf), nlow, wave)

   end subroutine sample

   subroutine lin_inter_1D(array, cdf, value, length, nlow, y)
!
!  linear interpolates between values for an array and its cdf
!
      implicit none

      real,    intent(OUT)  :: y
      integer, intent(IN)   :: length
      real,    intent(IN)   :: value,array(length,2),cdf(length-1)
      integer, intent(IN)   :: nlow

      y = array(nlow+1,1) + (array(nlow+2,1) - array(nlow+1,1)) * (value - cdf(nlow))/(cdf(nlow+1) - cdf(nlow))

   end subroutine lin_inter_1D

   subroutine lin_inter_2D(array,value,length,nlow,y)
!
!  linear interpolation for an array
!
      implicit none

      real,    intent(OUT)  :: y
      integer, intent(IN)   :: length
      real,    intent(IN)   :: value,array(length,2)
      integer, intent(IN)   :: nlow

      y = array(nlow,2) + (array(nlow+1,2) - array(nlow,2)) * (value - array(nlow,1))/(array(nlow+1,1) - array(nlow,1))

   end subroutine lin_inter_2D

   subroutine search_1D(length,array,nlow,value)
!
!  search by bisection for 1D array
!
      implicit none

      integer              :: nup,length,middle
      integer, intent(OUT) :: nlow
      real,    intent(in)  :: array(length),value

      nup = length
      nlow = 1
      middle = int((nup+nlow)/2.)

      do while((nup - nlow).gt.1)
         middle = int((nup + nlow)/2.)
         if(value.gt.array(middle))then
            nlow = middle
         else
            nup = middle
         end if
      end do
   end subroutine search_1D

   subroutine search_2D(length,array,nlow,value)
!
!  search by bisection for 2D array
!
      implicit none

      integer              :: nup,length,middle
      integer, intent(OUT) :: nlow
      real,    intent(in)  :: array(length,2),value

      nup = length
      nlow = 1
      middle = int((nup+nlow)/2.)

      do while((nup - nlow).gt.1)
         middle = int((nup + nlow)/2.)
         if(value.gt.array(middle,1))then
            nlow = middle
         else
            nup = middle
         end if
      end do
   end subroutine search_2D

   subroutine mk_cdf(array,cdf,length)
!
!  subroutine that creates cdf for an array of values.
!
      implicit none

      integer, intent(IN)    :: length
      real,    intent(IN)    :: array(length,2)
      real,    intent(INOUT) :: cdf(length)
      real                   :: summ
      integer                :: i,j

      do j=1,length-1
         summ=0.
         do i=1,j
            summ=summ+0.5*(array(i+1,2)+array(i,2))*(array(i+1,1)-array(i,1))
         end do
         cdf(j)=summ
      end do
      cdf=cdf/cdf(length-1)

   end subroutine mk_cdf
end module ch_opt

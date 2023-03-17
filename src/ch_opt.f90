MODULE ch_opt

implicit none

CONTAINS
  subroutine init_opt

    use opt_prop

    implicit none

    !633nm

    !white matter
    whitematter_hgg = 0.96d0
    whitematter_g2  = whitematter_hgg**2.
    whitematter_mua = 1.58d0
    whitematter_mus = 51

    whitematter_kappa  = whitematter_mus + whitematter_mua
    whitematter_albedo = whitematter_mus / whitematter_kappa

    !grey matter
    greymatter_hgg = 0.88d0
    greymatter_g2  = greymatter_hgg**2.
    greymatter_mua = 2.63d0
    greymatter_mus = 60.2

    greymatter_kappa  = greymatter_mus + greymatter_mua
    greymatter_albedo = greymatter_mus / greymatter_kappa

    !CSF !660nm - dont have 633nm
    csf_hgg = 0.9d0
    csf_g2  = csf_hgg**2.
    csf_mua = 0.04d0
    csf_mus = 0.35

    csf_kappa  = csf_mus + csf_mua
    csf_albedo = csf_mus / csf_kappa

    !glial matter !using whtematter properties
    glial_hgg = 0.96d0
    glial_g2  = glial_hgg**2.
    glial_mua = 1.58d0
    glial_mus = 51

    glial_kappa  = glial_mus + glial_mua
    glial_albedo = glial_mus / glial_kappa

    !GBM cells
    gbm_hgg = 0.875d0
    gbm_g2  = glial_hgg**2.
    gbm_mua = 0.2d0
    gbm_mus = 160

    gbm_kappa  = gbm_mus + gbm_mua
    gbm_albedo = gbm_mus / gbm_kappa







  end subroutine init_opt

   subroutine init_opt1
!
!  subroutine to set tissue optical properties 630nm
!
   use opt_prop
   use iarray, only :con_ppix, rhokap, albedoar, ua_ppix
   use constants, only : nxg, nyg, nzg

   implicit none

   integer:: i,j,k

   hgg = 0.8027
   g2  = hgg**2.
   mua = 0.23d0
   mus = 21.d0/(1.d0 - hgg)

   e630=0.0265

   kappa  = mus + mua
   albedo = mus / kappa

   ! calculate new ppix absorbtion coefficent and set new rhokap and albedo for each voxel
   do i= 1, nxg
     do j=1,nyg
       do k=1,nzg

         ua_ppix(i,j,k)=e630 * con_ppix(i,j,k)

         rhokap(i,j,k)=mus + mua + ua_ppix(i,j,k)
         albedoar(i,j,k)= mus / rhokap(i,j,k)

        end do
      end do
    end do
   end subroutine init_opt1

   subroutine init_opt2
!
!  subroutine to set tissue optical properties 420nm
!
   use opt_prop
   use iarray, only :con_ppix, rhokap, albedoar, ua_ppix
   use constants, only : nxg, nyg, nzg

   implicit none

   integer :: i,j,k

   hgg = 0.9
   g2  = hgg**2.
   mua = 1.8d0
   mus = 82.d0/(1.d0 - hgg)

   e420=0.105

   !kappa  = mus + mua
   !albedo = mus / kappa

   ! calculate new ppix absorbtion coefficent and set new rhokap and albedo for each voxel
   do i= 1, nxg
     do j=1,nyg
       do k=1,nzg

         ua_ppix(i,j,k)=e420 * con_ppix(i,j,k)
         rhokap(i,j,k)=mus + mua + ua_ppix(i,j,k)
         albedoar(i,j,k)= mus / rhokap(i,j,k)

        end do
      end do
    end do

   end subroutine init_opt2

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

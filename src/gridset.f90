MODULE gridset_mod

    implicit none

    private
    public :: gridset

    contains
        subroutine gridset(id,ts)

            use constants
            use iarray
            use opt_prop
            use ch_opt

            implicit none

            integer, intent(IN) :: id, ts

            integer             :: i, j, k
            real                :: x, y, z, taueq1, taupole1, taueq2, taupole2,rad,a, b, c,rc,len,rh, rad_small, rad_blob


            if(id == 0)then
                print*, ' '
                print *, 'Setting up density grid....'
            end if

            ! setup grid faces
            do i = 1, nxg + 1
                xface(i) = (i - 1) * 2. * xmax/nxg
                x = (i-1)*2.*xmax/nxg
            end do

            do i = 1, nyg + 1
                yface(i) = (i - 1) * 2. * ymax/nyg
                y=(j-1)*2.*ymax/nyg
            end do

            do i = 1, nzg + 2
                zface(i) = (i - 1) * 2. * zmax/nzg
                z=(k-1)*2.*zmax/nzg
            end do

            call init_opt_630
            refrac(:,:,:) = n2

            !set up optical properties grid
            do i = 1, nxg
                x = xface(i) - xmax + xmax/nxg
                do j = 1, nyg
                    y = yface(j) - ymax + ymax/nyg
                    do k = 1, nzg
                        z = zface(k) - zmax + zmax/nzg

                        !set tumour sticking out of brain due to balloon algorithm to zeros
                        if(white_matter(i,j,k) .eq. 0. .and. grey_matter(i,j,k) .eq. 0. &
                        .and. glial_matter(i,j,k) .eq. 0.)then
                          !tumour(i,j,k)=0.

                          if(tumour(i,j,k) .ne. 1.)then
                          tumour_killed(i,j,k)=0.
                          tumour_resec(i,j,k)=0.
                        endif
                      endif

                  !    rhokap(i,j,k) = white_matter(i,j,k)*abs(tumour_resec(i,j,k)-&
                  !    (white_matter(i,j,k)+grey_matter(i,j,k)+glial_matter(i,j,k)))*whitematter_kappa &
                  !    + grey_matter(i,j,k)*abs(tumour_resec(i,j,k)-(white_matter(i,j,k)+&
                  !     grey_matter(i,j,k)+glial_matter(i,j,k)))*greymatter_kappa &
                  !     +glial_matter(i,j,k)*abs(tumour_resec(i,j,k)-(white_matter(i,j,k)+&
                  !      grey_matter(i,j,k)+glial_matter(i,j,k)))*glial_kappa &
                  !     + tumour_resec(i,j,k)*rhokap_gbm(i,j,k)

                  if(white_matter(i,j,k) .ne. 0. .or. grey_matter(i,j,k) .ne. 0. &
                  .or. glial_matter(i,j,k) .ne. 0 .or. tumour_resec(i,j,k) .ne. 0.)then

                  rhokap(i,j,k)=((1-tumour_resec(i,j,k))*(white_matter(i,j,k)/(white_matter(i,j,k)&
                  +grey_matter(i,j,k)+glial_matter(i,j,k))))*whitematter_kappa &
                  +((1-tumour_resec(i,j,k))*(grey_matter(i,j,k)/(white_matter(i,j,k)&
                  +grey_matter(i,j,k)+glial_matter(i,j,k))))*greymatter_kappa &
                  + ((1-tumour_resec(i,j,k))*(glial_matter(i,j,k)/(white_matter(i,j,k)&
                  +grey_matter(i,j,k)+glial_matter(i,j,k))))*glial_kappa &
                  + tumour_resec(i,j,k)*rhokap_gbm(i,j,k)

                  albedoar(i,j,k)=((1-tumour_resec(i,j,k))*(white_matter(i,j,k)/(white_matter(i,j,k)&
                  +grey_matter(i,j,k)+glial_matter(i,j,k))))*whitematter_albedo &
                  +((1-tumour_resec(i,j,k))*(grey_matter(i,j,k)/(white_matter(i,j,k)&
                  +grey_matter(i,j,k)+glial_matter(i,j,k))))*greymatter_albedo &
                  + ((1-tumour_resec(i,j,k))*(glial_matter(i,j,k)/(white_matter(i,j,k)&
                  +grey_matter(i,j,k)+glial_matter(i,j,k))))*glial_albedo &
                  + tumour_resec(i,j,k)* albedoar_gbm(i,j,k)

                  hggar(i,j,k)=((1-tumour_resec(i,j,k))*(white_matter(i,j,k)/(white_matter(i,j,k)&
                  +grey_matter(i,j,k)+glial_matter(i,j,k))))*whitematter_hgg &
                  +((1-tumour_resec(i,j,k))*(grey_matter(i,j,k)/(white_matter(i,j,k)&
                  +grey_matter(i,j,k)+glial_matter(i,j,k))))*greymatter_hgg &
                  + ((1-tumour_resec(i,j,k))*(glial_matter(i,j,k)/(white_matter(i,j,k) &
                  +grey_matter(i,j,k)+glial_matter(i,j,k))))*glial_hgg &
                  + tumour_resec(i,j,k)* gbm_hgg


                endif

!print*, rhokap(i,j,k), albedoar(i,j,k), hggar(i,j,k)
                    ! albedoar(i,j,k) = white_matter(i,j,k)*abs(tumour_resec(i,j,k)-&
                     !(white_matter(i,j,k)+grey_matter(i,j,k)+glial_matter(i,j,k)))*whitematter_albedo &
                     !+ grey_matter(i,j,k)*abs(tumour_resec(i,j,k)-(white_matter(i,j,k)+&
                    !  grey_matter(i,j,k)+glial_matter(i,j,k)))*greymatter_albedo &
                    !  +glial_matter(i,j,k)*abs(tumour_resec(i,j,k)-(white_matter(i,j,k)+&
                    !   grey_matter(i,j,k)+glial_matter(i,j,k)))*glial_albedo &
                    !  + tumour_resec(i,j,k)*albedoar_gbm(i,j,k)


                  !  hggar(i,j,k)=((1-tumour(i,j,k))*(white_matter(i,j,k)/(white_matter(i,j,k) &
                  !  +grey_matter(i,j,k)+glial_matter(i,j,k))))*whitematter_hgg &
                  !  +((1-tumour(i,j,k))*(grey_matter(i,j,k)/(white_matter(i,j,k)&
                  !  +grey_matter(i,j,k)+glial_matter(i,j,k))))*greymatter_hgg &
                  !  + ((1-tumour(i,j,k))*(glial_matter(i,j,k)/(white_matter(i,j,k)&
                  !  +grey_matter(i,j,k)+glial_matter(i,j,k))))*glial_hgg &
                  !  + tumour(i,j,k)* gbm_hgg(i,j,k)

                      !hggar(i,j,k) =white_matter(i,j,k)*abs(tumour_resec(i,j,k)-&
                      !(white_matter(i,j,k)+grey_matter(i,j,k)+glial_matter(i,j,k)))*whitematter_hgg &
                      !+ grey_matter(i,j,k)*abs(tumour_resec(i,j,k)-(white_matter(i,j,k)+&
                      ! grey_matter(i,j,k)+glial_matter(i,j,k)))*greymatter_hgg &
                       !+glial_matter(i,j,k)*abs(tumour_resec(i,j,k)-(white_matter(i,j,k)+&
                      !  grey_matter(i,j,k)+glial_matter(i,j,k)))*glial_hgg &
                      ! + tumour_resec(i,j,k)*gbm_hgg

                      !sets refrac to n1 in materials that are not brain
                        if(rhokap(i,j,k) .eq. 10.001 .or. rhokap(i,j,k).eq.0.)then
                          refrac(i,j,k)=n1
                        endif
! x pos of light = 1/4 into grid - 236/4 = 59
!y pos of light = 0.166 into grid - 206*0.166 = 34.3 - 35
! z pos of light= 0.606 - 150*0.606 = 91

!                    end do
!                end do
!            end do

ua_brain(i,j,k) = (white_matter(i,j,k)*whitematter_mua)+(grey_matter(i,j,k)&
*greymatter_mua)+(glial_matter(i,j,k)*glial_mua)&
                   +(tumour_resec(i,j,k)*gbm_mua)

            !************ add cavity into brain ****************************************

            rad= 1.5 !sphere cavity radius for 3cm radius cavity
            rc=0.05 ! glass tube radius
            len=7. ! glass tube length
            rh=0.3
            rad_small=1.6



               !x^2 + y^2 + z^2 = r^2 is equation of a sphere

               !x^2/a^2 + y^2/b^2 + z^2/c^2 =1 is equation of an ellipse.
               !values taken from table 5.3 in Lille thesis for comparison
               !long in the y direction

               a=4.0/2.
               b=7.0/2.
               c=4.0/2.



          !   do i=1,nxg
          !    x = (i-1)*2.*xmax/nxg
          !     do j=1, nyg
          !     y=(j-1)*2.*ymax/nyg
          !       do k=1, nzg
          !      z=(k-1)*2.*zmax/nzg
            !insert space to clear brain at top of balloon
                !if(sqrt((x-(2.*xmax-3.9))**2+(y-(2.*ymax-1.5))**2+(z-(2.*zmax-6.37))**2) .le. sqrt(rad_small**2)) then
                !if(sqrt((x-(2.*xmax))**2+(y-(2.*ymax-1.5))**2+(z-(2.*zmax))**2) .le. sqrt(rad_small**2)) then
            !    if(sqrt((x-xmax+0.7)**2+(y-ymax+2.7)**2+(z-zmax)**2) .le. sqrt(rad_small**2)) then

            !       rhokap(i,j,k)=0.
            !       albedoar(i,j,k)=0.
            !       tumour_resec(i,j,k)=0.
            !       tumour_killed(i,j,k)=0.
            !       white_matter(i,j,k)=0.
            !       grey_matter(i,j,k)=0.
            !       ua_brain(i,j,k)=0.

            !       refrac(i,j,k)=n1
            !           hggar(i,j,k)=0.
                     !print*, 'done'
            !     end if

                !set resected cavity space to csf properties
                   if (tumour(i,j,k) .ge. 0.95)then
                     rhokap(i,j,k)=kappa_water
                     albedoar(i,j,k)=albedo_water
                     refrac(i,j,k)=1.33
                     tumour_resec(i,j,k)=0.
                     tumour_killed(i,j,k)=0.
                     white_matter(i,j,k)=0.
                     grey_matter(i,j,k)=0.
                     glial_matter(i,j,k)=0.
                     ua_brain(i,j,k)=0.
                   endif



                !if(sqrt((x-(2.*xmax-4.06))**2+(y-(2.*ymax-2.0))**2+(z-(2.*zmax-6.37))**2) .le. sqrt(rad**2)) then !sphere
                !if(((x-(2.*xmax-3.82))**2/a**2) + ((y-(2.*ymax-2.0))**2/b**2) +((z-(2.*zmax-6.37))**2/c**2) .le. 1.) then
                !!!if(((x-(2.*xmax-4.06))**2/a**2) + ((y-(2.*ymax-1.5))**2/b**2) +((z-(2.*zmax-6.37))**2/c**2) .le. 1.) then
                !if(((x-xmax)**2/a**2) + ((y-ymax)**2/b**2) +((z-zmax)**2/c**2) .le. 1.) then
                   ! values from van gemert paper - assuming lille intralipid is a 0.2%

                if(balloon(i,j,k) .eq. 100.)then
                   rhokap(i,j,k)=kappa_intra
                   albedoar(i,j,k)=albedo_intra
                   tumour_resec(i,j,k)=0.
                   tumour_killed(i,j,k)=0.
                   ua_brain(i,j,k)=0.

                    refrac(i,j,k)=1.33 !water refractive index

                    hggar(i,j,k)=hgg_intra
                     !print*, 'done'
                     !hgg=0.733 of intralipid at 635nm from p106 lille thesis
                 end if

                 rad_blob=0.2

                 !add in a blob
                 if(sqrt((x-(xmax+2.1))**2+(y-(ymax-0.7))**2+(z-(zmax))**2) .le. sqrt(rad_blob**2)) then !sphere
                   blob(i,j,k)=1.
                   rhokap(i,j,k)=gbm_kappa
                   albedoar(i,j,k)=gbm_albedo
                   hggar(i,j,k)=gbm_hgg
                   refrac(i,j,k)=n2
                   ua_brain(i,j,k)=gbm_mua

                   if(ts .eq.1)then
                     blob_killed(i,j,k) =1.
                     tumour_killed(i,j,k)=0.
                     tumour_resec(i,j,k)=0.
                   endif
                 endif


                 !if(((x-(xmax+1))**2/a**2) + ((y-(ymax+1.2))**2/b**2) +((z-(zmax))**2/c**2) .le. 1.) then
                  ! if(grey_matter(i,j,k).ne. 0. .and. tumour_resec(i,j,k) .ne. 0. .and. white_matter(i,j,k) .ne. 0. )then
              !     if(rhokap(i,j,k) .eq. 0. )then
                  !    rhokap(i,j,k)=kappa_intra
                  !     albedoar(i,j,k)=albedo_intra
                  !    white_matter(i,j,k)=0.
                  !    grey_matter(i,j,k)=0.
                  !    tumour_resec(i,j,k)=0.
                  !    tumour_killed(i,j,k)=0.

              !        refrac(i,j,k)=1.33 !water refractive index
!
!                      hggar(i,j,k)=hgg_intra
                  ! endif
                ! endif

                 !set all spaces to csf properties
                 if(rhokap(i,j,k) .eq. 0.)then
                   rhokap(i,j,k) = kappa_water!csf_kappa
                   albedoar(i,j,k)= albedo_water!csf_albedo
                   hggar(i,j,k)=hgg_water!csf_hgg
                   refrac(i,j,k)=1.33
                 endif


                 !if(sqrt((x-(2.*xmax-3.82))**2+(z-(2.*zmax-6.37))**2) .le. sqrt(rc**2) .and. y &
                 !if(sqrt((x-(2.*xmax-4.06))**2+(z-(2.*zmax-6.37))**2) .le. sqrt(rc**2) .and. y &
                 !.ge. (2.*ymax-1.5)+0.1-(len/2.) .and. y .le. (2.*ymax-1.5)+0.1+(len/2.)) then !sphere
                  ! rhokap(i,j,k)=50!1E-4 + 0.3
                   !albedoar(i,j,k)=0.3/(1E-4 + 0.3)
                  ! hggar(i,j,k)=0.9
                  ! tumour_resec(i,j,k)=0.
                   !tumour_killed(i,j,k)=0.

                    !refrac(i,j,k)=1.51 !borosilicat glass - from Lille parallel paper
                ! endif

                 !if(sqrt((x-(2.*xmax-3.82))**2+(z-(2.*zmax-6.37))**2) .le. sqrt(rh**2) .and. y &
                 !if(sqrt((x-(2.*xmax-4.06))**2+(z-(2.*zmax-6.37))**2) .le. sqrt(rh**2) .and. y &
                 !.ge. (2.*ymax-1.3)+0.1-(len/2.) .and. y .le. (2.*ymax-1.3)+0.1+(len/2.)) then !sphere
                  ! rhokap(i,j,k)=0.
                   !albedoar(i,j,k)=0.
                   !tumour_resec(i,j,k)=0.
                   !tumour_killed(i,j,k)=0.

                    !refrac(i,j,k)=n1 !air to make glass tube hollow
                ! endif

                    if(ts .eq. 1)then
                    S_0(i,j,k)=3.d0*(tumour_resec(i,j,k) + blob(i,j,k))
                    end if
                          !S_0(i,jj,k)=7.d0*(balloon(i,jj,k)/100.)





                   end do
                 end do
               end do

            !****************** Calculate equatorial and polar optical depths ****
            taueq1   = 0.
            taupole1 = 0.
            taueq2   = 0.
            taupole2 = 0.

            do i = 1, nxg
                taueq1 = taueq1 + rhokap(i,nyg/2,nzg/2)
            end do

            do i = 1, nzg
                taupole1 = taupole1 + rhokap(nxg/2,nyg/2,i)
            end do

            taueq1 = taueq1 * 2. * xmax/nxg
            taupole1 = taupole1 * 2. * zmax/nzg
            if(id == 0)then
                print'(A,F9.5,A,F9.5)',' taueq1 = ',taueq1,'  taupole1 = ',taupole1
            end if

            if(id == 0)then
                inquire(iolength=i)refrac(:,:,:nzg)
                open(newunit=j,file='refrac.dat',access='stream',form='unformatted',status='replace')
                write(j)refrac(:,:,:)
                close(j)
            end if
            ! call mpi_finalize()
            ! stop
        end subroutine gridset
end MODULE gridset_mod

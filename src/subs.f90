 MODULE subs

implicit none

!private
!public :: mcpolar_flu

    contains

        subroutine directory
        !  subroutine defines vars to hold paths to various folders
        !
        !
            use constants, only : cwd, homedir, fileplace, resdir

            implicit none

            !get current working directory

            call get_environment_variable('PWD', cwd)

            ! get 'home' dir from cwd
            homedir = trim(cwd(1:len(trim(cwd))-3))
            ! get data dir
            fileplace = trim(homedir)//'data/jmean/'
            ! get res dir
            resdir = trim(homedir)//'res/'

        end subroutine directory


        subroutine zarray

            use iarray

            !sets all arrays to zero
            implicit none


            jmean = 0.
            jmean_flu=0.
            xface = 0.
            yface = 0.
            zface = 0.
            rhokap = 0.
            albedoar=0.
            jmeanGLOBAL = 0.
            jmean_fluGLOBAL = 0.
            refrac = 0.

            si_step=0.
            si_step_GLOBAL=0.
            nop_tot=0.
            time_tot=0.
            con_ppix=0.
            ua_ppix=0.
            con_test=0.
            con_test_GLOBAL=0.

            esc_flur=0.
            esc_flur_GLOBAL=0.
            alb_absor=0.
            pl_absor=0.
            pl_absor_GLOBAL=0.
            flu_phot_rel=0.
            obs_flur=0.
            obs_flur_GLOBAL=0.
            jme=0.
            av_ppix_con_630=0.
            del_Q=0.
            PDD=0.
            PDD_GLOBAL=0.
            del_pdd=0.

            white_matter=0.
            grey_matter=0.
            tumour=0.



        end subroutine zarray


        subroutine alloc_array
        !  subroutine allocates allocatable arrays
        !
        !
            use iarray
            use constants,only : nxg,nyg,nzg,im_xvox,im_yvox, end_time, start_time

            implicit none

            allocate(xface(nxg+1), yface(nyg + 1), zface(nzg + 2))
            allocate(rhokap(nxg, nyg, nzg))
            allocate(albedoar(nxg, nyg, nzg+1))
            allocate(jmean(nxg, nyg, nzg+1), jmeanGLOBAL(nxg, nyg, nzg+1), jme(nxg,nyg,nzg+1), jm_av(nzg))
            allocate(jmean_flu(nxg, nyg, nzg+1), jmean_fluGLOBAL(nxg, nyg, nzg+1))
            allocate(refrac(nxg, nyg, nzg+1))

            allocate(si_step(nxg, nyg, nzg+1), si_step_GLOBAL(nxg, nyg, nzg+1))
            allocate(nop_tot(nxg,nyg,nzg+1), time_tot(nxg,nyg,nzg+1))
            allocate(con_ppix(nxg,nyg,nzg+1),ua_ppix(nxg,nyg,nzg+1), con_test((end_time-start_time)+1))
            allocate(con_test_GLOBAL((end_time-start_time)+1))
            allocate(p(nxg, nyg, nzg+1))
            allocate(np_cell(nxg, nyg, nzg+1))

            allocate(flu_image(im_xvox,im_yvox))

            allocate(esc_flur(nzg), esc_flur_GLOBAL(nzg))
            allocate(alb_absor(nzg), pl_absor(nzg), pl_absor_GLOBAL(nzg))
            allocate(flu_phot_rel(nzg))
            allocate(obs_flur((end_time-start_time)+1),obs_flur_GLOBAL((end_time-start_time)+1))
            allocate(av_ppix_con_630((end_time-start_time)+1,nzg))
            allocate(del_Q(nxg,nyg,nzg), PDD(nxg,nyg,nzg),PDD_GLOBAL(nxg,nyg,nzg),del_pdd(nxg,nyg,nzg))

            allocate(white_matter(nxg,nyg,nzg), grey_matter(nxg,nyg,nzg), tumour(nxg,nyg,nzg))



        end subroutine alloc_array

        subroutine time_step_calc(total_time,time_step,nphotons,xmax,ymax,zmax,numproc)

          use iarray, only : si_step_GLOBAL, nop_tot, jmeanGLOBAL, time_tot,ua_ppix
          use constants, only : nxg, nyg, nzg
          use opt_prop

          implicit none

          real, intent(OUT) :: time_step
          real, intent(INOUT) :: total_time
          real, intent(IN) :: xmax,ymax,zmax
          integer, intent(IN):: nphotons, numproc


          real :: q, nop, B, G, irad,p,ful_p
          integer :: i,j,k
          real :: h,c,w,Co,A,V
          real(16):: CC

      !  jmeanGLOBAL =jmeanGLOBAL * ((2.*xmax)**2./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg))) !convert to fluence

          w=630!*10**(-9)
          !h=6.626*10**(-34)
          !c=3*10**(8)
          !CC=1*10^-9/(h*c)
          CC=5.3*10**(15.)
          B=105
          e420=0.105
          e630=0.0265

          Co=0.2 !initial concentration of PpIX
          irad=0.1 ! irradiance of light
          A=(2.*xmax)*(2.*ymax) !irradiance area
          V=(2.*xmax)*(2.*ymax)*(2.*zmax) !voxel volume

          !I=1 , A=1 , N=500000, V=1, ua= updated absoption coefficent of ppix
          !number of photons = q*wavelength/h*c (h=placks constant, c=speed of light)
          ! total_time= - (ln(1-(Pdd/G)*B)/jmean
          ! G= (bECoB)
          !b=wavelength/(h*c)

          G= (w*CC)*e420*Co*B

         !ful_p=0.

          !do i= 1, nxg
          !  do j=1,nyg
          !    do k=1,nzg+1
                 !nop=0. !initialise number of photons absorbed in voxel for this step to 0.

          !       p(i,j,k)=ua_ppix(i,j,k)*si_step_GLOBAL(i,j,k)
          !       ful_p=ful_p + p(i,j,k)

                 !q=((irad*A)/(nphotons*numproc*V))*ua_ppix(i,j,k)*si_step_GLOBAL(i,j,k) !calculate energy absorbed in cell i,j,k by ppix
                 !nop= q*CC !convert energy absorbed into a number of photons
                 !nop_tot(i,j,k)=nop_tot(i,j,k) + nop !add number of photons absorbed in current step to total over all steps - aka the PDD

                 !time_tot(i,j,k)= -(LOG(1-(nop_tot(i,j,k)/G))*B)/jmeanGLOBAL(i,j,k) !calculates total time for each voxel after current time step
          !     end do
          !   end do
           !end do

           !do i= 1, nxg
            ! do j=1,nyg
            !   do k=1,nzg+1
                  !nop=0. !initialise number of photons absorbed in voxel for this step to 0.

            !      np_cell(i,j,k)=nphotons*numproc(p(i,j,k)/ful_p)

                  !q=((irad*A)/(nphotons*numproc*V))*ua_ppix(i,j,k)*si_step_GLOBAL(i,j,k) !calculate energy absorbed in cell i,j,k by ppix
                  !nop= q*CC !convert energy absorbed into a number of photons
                  !nop_tot(i,j,k)=nop_tot(i,j,k) + nop !add number of photons absorbed in current step to total over all steps - aka the PDD

                  !time_tot(i,j,k)= -(LOG(1-(nop_tot(i,j,k)/G))*B)/jmeanGLOBAL(i,j,k) !calculates total time for each voxel after current time step
            !    end do
            !  end do
            !end do


        !   do i= 1, nxg
        !     do j=1,nyg
        !       do k=1,nzg+1

        !         time_tot(i,j,k)=time_tot(i,j,k)-total_time !subtract overall total time from each calculated value to get individual time steps

!                end do
!              end do
!            end do



          !  do i= 1, nxg
          !    do j=1,nyg
          !      do k=1,nzg+1

                  !set zeros to the maximum value to remove them
          !        if (time_tot(i,j,k) .le. 0.) then
          !          time_tot(i,j,k)=10000000.!maxval(time_tot)
          !        end if

          !       end do
          !     end do
          !   end do



           !time_step=minval(time_tot) !choose the smallest time step from above as chosen time step

           !if(total_time .le. 180.)then
           time_step=0.2 !set time step to 1 second for testing
         !else
          ! time_step=60.
         !endif


           total_time= total_time + time_step
           print*, 'Total time', total_time

        !   print*, 'min time', minval(time_tot)
        !   print*, 'max time', maxval(time_tot)
        !   print*, 'timestep', time_step


      !still need to add check for total burn through of cell

    end subroutine time_step_calc!(total_time,time_step)

    subroutine ppix_con(total_time,nphotons,numproc,ts)

!calculates conentrarion of ppix at the start of each time step

      use iarray, only : jmean, con_ppix,jme, si_step, jm_av,con_test,av_ppix_con_630
      use constants, only : nxg, nyg, nzg,xmax,ymax,zmax

      implicit none

        real, intent(IN) :: total_time
        integer, intent(IN):: nphotons, numproc,ts

        real :: B,co,layav
        integer :: i,j,k

   !jmeanGLOBAL =jmeanGLOBAL * ((2.*xmax)**2./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg))) !convert to fluence


      co=0.2 !initial concentration of ppix
      B=105. !photobleaching constant
      jm_av=0.
print*, total_time,ts
      do k = 1, nzg
          do j = 1, nyg !fluence over whole layer...
              do i = 1, nxg
               jme(i,j,k)=jmean(i,j,k)*(0.1*(2.*xmax)**2./(nphotons*(ts)*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))

                con_ppix(i,j,k)=co*exp(-(jme(i,j,k)*(total_time))/B)

                if(i.eq.50 .and. j.eq.50 .and. k.eq.80.)then
                  con_test(ts)=con_ppix(i,j,k)
                endif

              end do
          end do
      end do

      !create array of average concentration at depth for each time step to use in 420nm simulation - where 630nm does photobleaching for jaques validation
      do k = 1, nzg
        layav=0.
          do j = 1, nyg
              do i = 1, nxg
                layav=layav+con_ppix(i,j,k)
              end do
          end do
          layav=layav/(nxg*nyg)
          av_ppix_con_630(ts,k)=layav
      end do

    ! print*, 'ppix_con',jmean(100,100,100),jme(100,100,100)
    !  jm_av = jm_av / (nxg*nyg)

    !  do i= 1, nxg
    !    do j=1,nyg
    !      do k=1,nzg

    !        con_ppix(i,j,k)=co*exp(-(jme(i,j,k)*total_time)/B)

    !       end do
    !     end do
    !   end do

       !print*, 'ppix', con_ppix(100,100,100)


    end subroutine ppix_con

    subroutine ppix_con_420(ts)

    !sets concentration from 630nm simulation

      use iarray, only : con_ppix, av_ppix_con_630
      use constants, only : nxg, nyg, nzg

      implicit none
        integer, intent(IN):: ts

        integer :: i,j,k

        do k = 1, nzg
            do j = 1, nyg
                do i = 1, nxg
                  con_ppix(i,j,k)=av_ppix_con_630(ts,k)
                end do
            end do
        end do

        end subroutine ppix_con_420



!updates photodynamic dose over the whole grid
    !subroutine pdd_calc(total_time)

    !  use iarray, only : jmeanGLOBAL, pdd
    !  use constants, only : nxg, nyg, nzg
    !  use opt_prop

    !  implicit none

    !    real, intent(IN) :: total_time

    !    integer :: w, i,j,k
    !    real :: CC, B, G, Co

    !    w=630

    !    CC=5.3*10**(15.)
    !    B=105.!photobleaching constant
    !    e420=0.105
    !    e630=0.0265
    !    Co=0.2 !initial concentration of PpIX

        ! G= (bECoB)
        !b=wavelength/(h*c)
        !h=6.626*10**(-34)
        !c=3*10**(8)  !h=6.626*10**(-34)
          !c=3*10**(8)
          !CC=1*10^-9/(h*c)
        !CC=1*10^-9/(h*c)

    !    G= (w*CC)*e630*Co*B

      !  do i= 1, nxg
      !    do j=1,nyg
      !      do k=1,nzg+1

      !        pdd(i,j,k)=G*(1-exp(-(jmeanGLOBAL(i,j,k)*total_time)/B))

      !       end do
      !     end do
      !   end do
!print*, 'pdd',pdd(200,100,100), pdd(150,100,100), pdd(100,100,100)

   !end subroutine pdd_calc

! call this at the end of the time step but before mpi stuff so that multiple cores can be used
   !subroutine flu_705
    ! use ch_opt

     !implicit none

     !can basically put a simplified and modified version of mcpolar here to run the fluorescne simulation
     !set up optical properties and grid for 705nm
     !remove end zero values from absoption position array obtained in main simulation (each core has its own one)
     !count to get total number of fluoresence photons
     !loop through remaining positions
     ! set initial photon weight to be ?0.1? and release it isotripically from positions
     !multiply weight by albedo at each interaction
     !set minimum weight as 1*10**(-3)
     !then either terminate it or do russian roulette (see mcrt lecture notes lecture 7)
     !kill if absorbed
     !if photon escapes out the top of the grid - count it and strore its inital postion (z depth)

     !call init_opt3 !changes optical properties to 705nm

  !end subroutine flu_705
  !**************************************************************************

  subroutine photon_cell_number(nphotons, numproc, iseed,total_time)

    use iarray
    use constants
    use opt_prop

    implicit none

      integer, intent(IN) :: nphotons, numproc, iseed
      real, intent(IN):: total_time

      real :: tt, dd, ran2
      integer :: i,j,k




            ! ful_p=0.
             np_cell=0.

              do i= 1, nxg
                do j=1,nyg
                  do k=1,nzg
                     !print*, ua_ppix(i,j,k),si_step(i,j,k)
                     !p(i,j,k)=ua_ppix(i,j,k)*si_step(i,j,k) !calculates power absorbed by ppix in each cell
                     !p(i,j,k)=(2.303*con_ppix(i,j,k)*refrac(i,j,k)*e630)*si_step(i,j,k)
                     p(i,j,k)=ua_ppix(i,j,k)*si_step(i,j,k)

                    if(total_time.eq.0.)then
                     ful_p=ful_p + p(i,j,k)
                   endif

                   end do
                 end do
               end do

               if(total_time.eq.1.)then
                print*, 'full_power',ful_p
               endif

               do i= 1, nxg
                 do j=1,nyg
                   do k=1,nzg


                      tt=nphotons*(p(i,j,k)/ful_p) ! calculates number of flourescence photons from each cell
                    !  print*, p(i,j,k)/ful_p, p(i,j,k), ful_p, ua_ppix(i,j,k)

                      dd=tt-int(tt)

                      if(ran2(iseed) .lt. dd .and. dd .gt. 0.) then  !helps for regions of low emmissivity
                      np_cell(i,j,k)=int(tt)+1
                      else
                      np_cell(i,j,k)=int(tt)
                      endif

                     pl_absor(k)= pl_absor(k)+np_cell(i,j,k)



                      !print*, i,j,k, np_cell(i,j,k)


                    end do
                  end do
                end do
                print*,'total photon', sum(np_cell)

  end subroutine photon_cell_number

 subroutine pdd_calc(total_time,nphotons,numproc)

   use iarray, only: p, del_Q,PDD,del_pdd
   use constants

   implicit none

   real, intent(IN):: total_time
   integer, intent(IN):: nphotons,numproc

   real:: b,c,del_t
   integer::i,j,k
   real::wl,h

   wl=630E-9
   c=3E8
   h=6.626E-34
   b=wl/(h*c)

   !if(total_time .le. 180.)then
   del_t=1. !set time step to 1 second for testing
  !else
   !del_t=60.
  !endif

   do i= 1, nxg
     do j=1,nyg
       do k=1,nzg

        del_Q(i,j,k)=p(i,j,k)*(0.1*(2.*xmax)**2./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))


        del_pdd(i,j,k)= del_Q(i,j,k)*b*del_t

        PDD(i,j,k)=PDD(i,j,k)+del_pdd(i,j,k)

        end do
      end do
    end do
!print*,'del_Q' ,del_Q(100,100,100), p(100,100,100), PDD(100,100,100)





 end subroutine pdd_calc



  subroutine mcpolar_flu(delta, flnscatt, id, iseed, nphotons,numproc, im_side,ts,total_time)

    use constants, only : nxg, nyg, nzg, xmax, ymax, zmax, TWOPI, PI, FOURPI
    use photon_vars
    use iarray, only : xface, yface, zface, albedoar, np_cell, flu_image, alb_absor, esc_flur, flu_phot_rel
    use ch_opt
    use stokes_mod
    use inttau2

    implicit none

    real, intent(IN) :: delta, im_side,total_time
    integer, intent(IN) :: id, nphotons, numproc,ts
    integer, intent(INOUT):: iseed
    double precision, intent(INOUT) :: flnscatt

    integer :: ic,jc,kc,j,ran, xcell, ycell, zcell, NP, u, z_orig, alab, rest
    logical :: tflag
    real :: ran2, tau_max, init_wgt, wgt1, nzp_tot

    !******************set observer plane direction***********************

  !  obs_cost= 1.d0 ! directly above top of grid
  !  obs_sint=sqrt(1. - obs_cost**2)
  !  obs_phi=0.
  !  obs_cosp=cos(obs_phi)
  !  obs_sinp=sin(obs_phi)


  !  obs_nxp = obs_sint * obs_cosp
  !  obs_nyp = obs_sint * obs_sinp
  !  obs_nzp = obs_cost
  iseed=-95648324+id+1
  iseed=-abs(iseed)

    call init_opt3 ! sets optical properties to 705nm
    call gridset_flu(id)
    call photon_cell_number(nphotons, numproc, iseed,total_time)

alab=0.
rest=0.
nzp_tot=0.

do ic= 1, nxg
  do jc=1,nyg
    do kc=1,nzg

      z_orig=kc


!print*, 'flu run',ic,jc,kc

      NP=np_cell(ic,jc,kc) !gets number of photons to emit from cell


    do j = 1, NP

      flu_phot_rel(z_orig)=flu_phot_rel(z_orig)+1

       tflag=.FALSE.

       if(mod(j,10000) == 0)then
          print *, j,' scattered photons completed on core: ',id
       end if

    !***** Release photon from random point in current grid cell ic,jc,kc and calculate distance to cell wall in initial direction*******************************
       call diffuse(xcell,ycell,zcell,ic,jc,kc,iseed, tflag, nzp_tot)
       !print*, xcell, ycell, zcell

       !***************** Find optical distance to top of grid (viewing plane)************************
       !call find_taumax(xcell,ycell,zcell,tau_max,iseed,tflag,obs_nxp,obs_nyp,obs_nzp)
       !print*, tau_max

       !*******************Set inital photon weight and bin into image***********************************
      ! wgt1= exp(-tau_max)/(2.*TWOPI)
       !print*, 'init_wgt',init_wgt

       !call bin_photons(im_side, init_wgt)
       !print*, 'BIN'

    !*********find maximum optical distance to edge of grid in chosen direction of travel and set photon weight*********

    !call find_taumax(xcell,ycell,zcell,tau_max,iseed,tflag,nxp,nyp,nzp)
    !print*, 'taumax'

    !if(tau_max .lt. 1.e-3)then

    !exit
    !end if

    !wgt1=(1.-exp(-tau_max))
    !wgt1=(exp(-tau_max))/(FOURPI*((tau_max/1.38)**2 ))!1.38=n
    !print*, wgt1

    !****** Find initial scattering location in direction of travel and move there

          call tauint_flu(xcell,ycell,zcell,tflag,iseed,delta,im_side,z_orig,wgt1,rest,ts)




    !******** Photon scatters in grid until it exits (tflag=TRUE)

          do while(tflag.eqv..FALSE.)

             ran = ran2(iseed)

             if(ran2(iseed) < albedoar(xcell,ycell,zcell))then!interacts with tissue
                !call peeloff(nxp,nyp,nzp,obs_nxp,obs_nyp,obs_nzp,TWOPI,xcell,ycell,zcell,iseed, wgt1, im_side,tflag,z_orig, tau_max)

                   call stokes(iseed)
                   flnscatt = flnscatt + 1

                else

                   tflag=.true.
                   alab=alab+1
                   exit
             end if


    !************ Find next scattering location

             call tauint_flu(xcell,ycell,zcell,tflag,iseed,delta,im_side,z_orig,wgt1,rest,ts)





          end do
          !print*, 'after while'


    end do      ! end loop over nph photons

    end do
    !  print*, 'here'
   end do
   !print*, 'here2'
  end do
  !print*, 'here3'

!  print*, 'alab',alab
!  print*, 'esc tot', sum(esc_flur)
!  print*, 'rest',rest
!  print*, 'nzp',nzp_tot

!print*, flu_image
!open(newunit=u,file='flu_image.dat',access='stream',status='REPLACE',form='unformatted')
!write(u) flu_image
!close(u)

   end subroutine mcpolar_flu

subroutine gridset_flu(id)

use constants, only : nxg, nyg, nzg, xmax, ymax, zmax
use iarray, only    : rhokap,xface,yface,zface, rhokap, refrac, albedoar, grey_matter, white_matter, albedoar
use opt_prop
use ch_opt

implicit none

integer, intent(IN) :: id

integer             :: i, j, k
real                :: x, y, z, taueq1, taupole1, taueq2, taupole2

if(id == 0)then
   print*, ' '
   print *, 'Setting up density grid....'
end if

! setup grid faces
!do i = 1, nxg + 1
!   xface(i) = (i - 1) * 2. * xmax/nxg
!end do

!do i = 1, nyg + 1
!   yface(i) = (i - 1) * 2. * ymax/nyg
!end do

!do i = 1, nzg + 2
!   zface(i) = (i - 1) * 2. * zmax/nzg
!end do

call init_opt3 ! sets grid to 705nm properties - already done at top of mcpolar_flu
refrac(:,:,:) = n1

!set up optical properties grid
do i = 1, nxg
  ! x = xface(i) - xmax + xmax/nxg
   do j = 1, nyg
  !     y = yface(j) - ymax + ymax/nyg
       do k = 1, nzg
  !         z = zface(k) - zmax + zmax/nzg
           rhokap(i,j,k) = (white_matter(i,j,k)*whitematter_kappa) + (grey_matter(i,j,k)*greymatter_kappa)

          albedoar(i,j,k) = (white_matter(i,j,k)*whitematter_albedo) &
          + (grey_matter(i,j,k)*greymatter_albedo)
           refrac(i,j,k) = n2
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

end subroutine gridset_flu


!emits photon from a random position in a grid cell - run instead of sourcephoton in mcpolar
subroutine diffuse(xcell,ycell,zcell,ic,jc,kc,iseed, tflag, nzp_tot)

use constants, only : nxg, nyg, nzg, xmax, ymax, zmax, TWOPI
use photon_vars
use iarray, only : xface, yface, zface
use inttau2

implicit none

integer, intent(OUT)   :: xcell, ycell, zcell
integer, intent(INOUT) :: iseed
integer, intent(IN) :: ic,jc,kc
logical, intent(INOUT):: tflag
real, intent(INOUT)::nzp_tot
real                :: ran2, xs,ys,zs, tau_max, thta, ran


!***** emit photon isotropically from grid cell ************
!print*, 'DIFFUSE'
xp=xface(ic)+ran2(iseed)*(xface(ic+1)-xface(ic))-xmax
yp=yface(jc)+ran2(iseed)*(yface(jc+1)-yface(jc))-ymax
zp=zface(kc)+ran2(iseed)*(zface(kc+1)-zface(kc))-zmax

!call random_number(ran)
thta=acos((2*ran2(iseed))-1.)

!call random_number(ran)
phi=TWOPI*ran2(iseed)


if(sin(thta) .le. 0.)then
  sint=0.
else
  sint=sqrt(sin(thta))
endif


!sint=(1.-cost*cost)
!if(sint.le.0.)then
!sint=0.
!else
!sint=sqrt(sint)
!endif

!phi=twopi*ran2(iseed)
!cosp=cos(phi)
!sinp=sin(phi)

!***** Set photon direction cosines for direction of travel *********
!nxp=sint*cosp
!nyp=sint*sinp
!nzp=cost

nxp=sint*cos(phi)
nyp=sint*sin(phi)
nzp=cos(thta)

!print*, nzp

nzp_tot=nzp_tot+nzp

!print*, 'start', xp,yp,zp,nxp,nyp,nzp

!***** Set Stokes fluxes ********************************************
!fi=1.
!fq=0.
!fu=0.
!fv=0.

!*************** Linear Grid *************************
xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
zcell=int(nzg*(zp+zmax)/(2.*zmax))+1

!*****************************************************
!call find_taumax(xcell,ycell,zcell, tau_max, iseed, tflag)

!print*,'tau_max=', tau_max
!return
end subroutine diffuse

subroutine peeloff(nxp,nyp,nzp,obs_nxp,obs_nyp,obs_nzp,TWOPI,xcell,ycell,zcell,iseed, wgt1, im_side, tflag, z_orig, tau_max)

  use ch_opt
  use opt_prop
  use iarray, only: albedoar, esc_flur
  use inttau2

  implicit none

  real, intent(IN) ::nxp,nyp,nzp,obs_nxp,obs_nyp,obs_nzp, TWOPI,wgt1, im_side, tau_max
  integer, intent(IN)::xcell,ycell,zcell,z_orig
  integer, intent(INOUT):: iseed
  logical, intent(INOUT):: tflag

  real :: calpha, hgfac, tau_max2,phot, theta, hg


 calpha=nxp*obs_nxp+nyp*obs_nyp+nzp*obs_nzp

 hgfac=(1.-g2)/(1.+g2-2.*hgg*calpha)**1.5/(2*TWOPI)

  !find maximum optical depth of observer direction using find tau_max again
  !print*,xcell,ycell,zcell,tau_max2,iseed,tflag,obs_nxp,obs_nyp,obs_nzp
  call find_taumax(xcell,ycell,zcell,tau_max2,iseed,tflag,obs_nxp,obs_nyp,obs_nzp)

  ! calculate phot (photon weight) and bin into image
  theta = abs(acos(nzp)-acos(obs_nzp))
  hg=(1-g2)/((1+g2-2*hgg*cos(theta))**(3/2))
  !phot=wgt1*hgfac*albedoar(xcell,ycell,zcell)*exp(-tau_max2)
  phot=albedoar(xcell, ycell, zcell) * (1-exp(-tau_max))*exp(-tau_max2)*hg

  esc_flur(z_orig)=esc_flur(z_orig) + phot

  call bin_photons(im_side,phot)

end subroutine peeloff






end MODULE subs

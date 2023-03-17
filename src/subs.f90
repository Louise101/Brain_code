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
            xface = 0.
            yface = 0.
            zface = 0.
            rhokap = 0.
            rhokap_gbm=0.
            albedoar=0.
            albedoar_gbm=0.
            jmeanGLOBAL = 0.
            refrac = 0.
            hggar=0.

            si_step=0.
            si_step_GLOBAL=0.
            ua_ppix=0.
            ua_brain=0.


            jme=0.

            white_matter=0.
            grey_matter=0.
            glial_matter=0.
            balloon=0.
            balloon_killed=0.
            tumour=0.
            tumour_killed=0.
            tumour_resec=0.
            percent_left=0.
            so_tot=0.
            o23_tot=0.
            o21_tot=0.

            S_0=0.
            O2_3=0.
            O2_1=0.

            S0_slice=0.
            o23_slice=0.
            o21_slice=0.
            tumkill_slice=0.
            temp_slice=0.

            blob=0.
            blob_killed=0.

            temp=0.
            temp_change=0.
            E_abs=0.










        end subroutine zarray


        subroutine alloc_array
        !  subroutine allocates allocatable arrays
        !
        !
            use iarray
            use constants,only : nxg,nyg,nzg,im_xvox,im_yvox, end_time, start_time

            implicit none

            allocate(xface(nxg+1), yface(nyg + 1), zface(nzg + 2))
            allocate(rhokap(nxg, nyg, nzg+1),rhokap_gbm(nxg, nyg, nzg+1))
            allocate(albedoar(nxg, nyg, nzg+1),albedoar_gbm(nxg, nyg, nzg+1))
            allocate(jmean(nxg, nyg, nzg+1), jmeanGLOBAL(nxg, nyg, nzg+1), jme(nxg,nyg,nzg+1))
            allocate(refrac(nxg, nyg, nzg+1), hggar(nxg, nyg, nzg+1) )

            allocate(si_step(nxg, nyg, nzg+1), si_step_GLOBAL(nxg, nyg, nzg+1))
            allocate(ua_ppix(nxg,nyg,nzg+1), ua_brain(nxg,nyg,nzg+1))

            allocate(white_matter(nxg,nyg,nzg), grey_matter(nxg,nyg,nzg), tumour(nxg,nyg,nzg), tumour_killed(nxg,nyg,nzg))
            allocate(tumour_resec(nxg,nyg,nzg), percent_left((end_time-start_time)+1), glial_matter(nxg,nyg,nzg))
            allocate(balloon(nxg,nyg,nzg),balloon_killed(nxg,nyg,nzg))
            allocate(so_tot((end_time-start_time)+1), o23_tot((end_time-start_time)+1), o21_tot((end_time-start_time)+1))

            allocate(S_0(nxg,nyg,nzg+1), O2_3(nxg,nyg,nzg+1),O2_1(nxg,nyg,nzg+1))

            allocate(S0_slice(nxg,nyg,(end_time-start_time)+1), o23_slice(nxg,nyg,(end_time-start_time)+1))
            allocate(o21_slice(nxg,nyg,(end_time-start_time)+1), tumkill_slice(nxg,nyg,(end_time-start_time)+1))
            allocate(temp_slice(nxg,nyg,(end_time-start_time)+1))

            allocate(blob(nxg,nyg,nzg+1), blob_killed(nxg,nyg,nzg+1))

            allocate(temp(nxg,nyg,nzg+1), E_abs(nxg,nyg,nzg+1), temp_change(nxg,nyg, nzg+1))

        end subroutine alloc_array

        subroutine time_step_calc(total_time,time_step)

          use iarray, only : si_step_GLOBAL, jmeanGLOBAL,ua_ppix
          use constants, only : nxg, nyg, nzg
          use opt_prop

          implicit none

          real, intent(OUT) :: time_step
          real, intent(INOUT) :: total_time

           time_step=1. !set time step to 1 second for testing


           total_time= total_time + time_step
           print*, 'Total time', total_time


    end subroutine time_step_calc!(total_time,time_step)

    subroutine ppix_con_wang(total_time,nphotons,numproc,ts)

!calculates conentrarion of ppix at the start of each time step

      use iarray
      use constants, only : nxg, nyg, nzg,xmax,ymax,zmax
      use opt_prop

      implicit none

        real, intent(IN) :: total_time
        integer, intent(IN):: nphotons, numproc,ts

        real :: B,co,layav, Be, gam, ep, n, del, kk
        real:: epsil, sig, beta, delta
        integer :: i,j,k

   !jmeanGLOBAL =jmeanGLOBAL * ((2.*xmax)**2./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg))) !convert to fluence


      !co=0.2 !initial concentration of ppix
      !B=105. !photobleaching constant
      !jm_av=0.

      Be=11.9
      gam=23E-6
      ep=21.d0
      n=90E-6
      del=128.d0
      kk=23.3

      !photofrin properties - change for ala!
      epsil= 3.7E-3
      sig=9E-5
      beta=11.9
      delta=33.d0



      do k = 1, nzg
          do j = 1, nyg
              do i = 1, nxg

                if(k .eq. 77)then
                  S0_slice(i,j,ts)=S_0(i,j,k)
                endif
               jme(i,j,k)=si_step_GLOBAL(i,j,k)*(2000.d0/(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))
               !jme(i,j,k)=si_step_GLOBAL(i,j,k)*(30.d0/(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))

           if (rhokap(i,j,k) .ne. 0 .and. rhokap(i,j,k) .ne. kappa_water .and. rhokap(i,j,k) .ne. kappa_intra)then
             !if (rhokap(i,j,k) .ne. 0 )then
               if(tumour_resec(i,j,k) + blob(i,j,k) .eq. 0.)then
            !   if(balloon(i,j,k) .eq. 0.)then
                 S_0(i,j,k)=0.1 !need to look at literature for this value
              !  S_0(i,j,k)= S_0(i,j,k)+ ((-epsil*sig*jme(i,j,k)*(S_0(i,j,k)+delta)*O2_3(i,j,k)&
              !                           /(O2_3(i,j,k)+beta))*S_0(i,j,k))
               else
                 S_0(i,j,k)= S_0(i,j,k)+ ((-epsil*sig*jme(i,j,k)*(S_0(i,j,k)+delta)*O2_3(i,j,k)&
                                          /(O2_3(i,j,k)+beta))*S_0(i,j,k))
              endif
            endif


                if(S_0(i,j,k).lt. 0.)then
                  S_0(i,j,k)=0.
                endif



              end do
          end do
      end do



    end subroutine ppix_con_wang

    subroutine trip_ox_wang(total_time,nphotons,numproc,ts) !changes triplet oxygen concentration with time

!calculates conentrarion of ppix at the start of each time step

      use iarray
      use constants, only : nxg, nyg, nzg,xmax,ymax,zmax
      use opt_prop

      implicit none

        real, intent(IN) :: total_time
        integer, intent(IN):: nphotons, numproc,ts

        real :: B,co,layav, Be, gam, ep, n, del, kk, gam_mm, k50, O23v
        real:: epsil, sig, beta, delta, gg, o2_3_t0, max_met, td, gtd
        integer :: i,j,k

   !jmeanGLOBAL =jmeanGLOBAL * ((2.*xmax)**2./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg))) !convert to fluence


      !co=0.2 !initial concentration of ppix
      !B=105. !photobleaching constant
      !jm_av=0.

      Be=11.9
      gam=23E-6
      ep=21.d0
      n=90E-6
      del=128.d0
      kk=23.3
      gam_mm=1.5
      k50=0.5
      O23v=80.d0

      epsil= 3.7E-3
      sig=9E-5
      beta=11.9
      delta=33.d0
      gg=21.6!10!37.4 !21.7! 0.022 ! see lab book notes for ref !subtracting metabolic rate of oxygen consumption see notes
      o2_3_t0=38!40!83!39.!83.d0


      td=(total_time-750)/632.1

      gtd=gg*((0.99*td**4+1.09*td**3+0.05*td**2+0.18*td+0.32)&
               /(td**4+1.16*td**3+0.18*td**2+0.24*td +0.31))


      do k = 1, nzg
          do j = 1, nyg
              do i = 1, nxg

                if(k .eq. 77)then
                  o23_slice(i,j,ts)=O2_3(i,j,k)
                endif
               jme(i,j,k)=si_step_GLOBAL(i,j,k)*(2000.d0/(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))
              ! jme(i,j,k)=si_step_GLOBAL(i,j,k)*(30.d0/(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))

               if (rhokap(i,j,k) .eq. 0 .or. rhokap(i,j,k) .eq. kappa_water .or. rhokap(i,j,k) .eq. kappa_intra)then
                 O2_3(i,j,k)=0.
               else
               O2_3(i,j,k)=O2_3(i,j,k) +gtd*(1-(O2_3(i,j,k)/o2_3_t0)) &
                                      - (epsil*jme(i,j,k)*S_0(i,j,k)/(O2_3(i,j,k)+beta))*O2_3(i,j,k)

                endif

              if(O2_3(i,j,k) .lt. 0.)then
                O2_3(i,j,k)=0.
              end if

            !  if(tumour_resec(i,j,k) .ne. 0)then
            !    print*, 'jme', jme(i,j,k)
            !  endif



              end do
          end do
      end do


end subroutine trip_ox_wang

subroutine sing_ox_wang(total_time,nphotons,numproc,ts) !changes triplet oxygen concentration with time

!calculates conentrarion of ppix at the start of each time step

  use iarray
  use constants, only : nxg, nyg, nzg,xmax,ymax,zmax
  use opt_prop

  implicit none

    real, intent(IN) :: total_time
    integer, intent(IN):: nphotons, numproc,ts

    real :: B,co,layav, Be, gam, ep, n, del, kk, gam_mm, k50, O23v, K0,a, D50
    real:: epsil, sig, beta, delta, gg, o2_3_t0
    integer :: i,j,k

!jmeanGLOBAL =jmeanGLOBAL * ((2.*xmax)**2./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg))) !convert to fluence


  !co=0.2 !initial concentration of ppix
  !B=105. !photobleaching constant
  !jm_av=0.

  Be=11.9
  gam=23E-6
  ep=21.d0
  n=90E-6
  del=128.d0
  kk=23.3
  gam_mm=1.5
  k50=0.5
  O23v=80.d0
  K0=0.037
  a=0.04
  D50=250.d0

  epsil= 3.7E-3
  sig=9E-5
  beta=11.9
  delta=33.d0



  do k = 1, nzg
      do j = 1, nyg
          do i = 1, nxg

            if(k .eq. 77)then
              o21_slice(i,j,ts)=O2_1(i,j,k)
            endif
           jme(i,j,k)=si_step_GLOBAL(i,j,k)*(2000.d0/(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))
           !jme(i,j,k)=si_step_GLOBAL(i,j,k)*(30.d0/(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))

           if (rhokap(i,j,k) .eq. 0 .or. rhokap(i,j,k) .eq. kappa_water .or. rhokap(i,j,k) .eq. kappa_intra)then
             O2_1(i,j,k)=0.
           else
           O2_1(i,j,k)=O2_1(i,j,k) + ((epsil*jme(i,j,k)*S_0(i,j,k)*O2_3(i,j,k))/(O2_3(i,j,k)+beta))
          endif


           if(O2_1(i,j,k) .lt. 0.)then
             O2_1(i,j,k)=0.
           end if



          end do
      end do
  end do


end subroutine sing_ox_wang

subroutine cell_kill(ts) !changes triplet oxygen concentration with time

!calculates conentrarion of ppix at the start of each time step

  use iarray
  use constants, only : nxg, nyg, nzg,xmax,ymax,zmax

  implicit none

  integer, intent(IN):: ts

    integer :: i,j,k



  do k = 1, nzg
      do j = 1, nyg
          do i = 1, nxg

            if(k .eq. 77)then
              tumkill_slice(i,j,ts)=tumour_killed(i,j,k) + blob_killed(i,j,k)
            endif
        if(O2_1(i,j,k) .ge. 560)then
        !if(O2_1(i,j,k) .ge. 1000.)then
          tumour_killed(i,j,k)=0.
          blob_killed(i,j,k)=0.
          !balloon_killed(i,j,k)=0.
        end if



          end do
      end do
  end do


end subroutine cell_kill

subroutine percent_killed_calc(total_time,ts)

!calculates percentage of tumour left at end of time step

  use iarray, only: tumour_killed, tumour_resec, percent_left, balloon, balloon_killed, blob, blob_killed
  use constants, only: nxg,nyg,nzg


  implicit none

  real, intent(IN) :: total_time
  integer, intent(IN):: ts

  integer :: i,j,k
  real :: per_kill, per_left

  per_left= (sum(tumour_killed+blob_killed)/sum(tumour_resec+blob))*100.
!  per_left= (sum(balloon_killed)/sum(balloon))*100.
  per_kill = 100 - per_left

  percent_left(ts)=per_left

end subroutine percent_killed_calc

subroutine summing_calc(total_time,ts)

!calculates percentage of tumour left at end of time step

  use iarray, only: S_0, O2_3, O2_1, so_tot, o23_tot, o21_tot

  implicit none

  real, intent(IN) :: total_time
  integer, intent(IN):: ts




so_tot(ts)=sum(S_0)
o23_tot(ts)=sum(O2_3)
o21_tot(ts)=sum(O2_1)

end subroutine summing_calc

subroutine temp_calc(total_time,ts, numproc, nphotons)

!calculates temperature change in each time step - no heat diffusion or cooling

  use iarray, only: temp, E_abs, ua_brain, temp_change, si_step_GLOBAL, temp_slice
  use constants, only : nxg, nyg, nzg,xmax,ymax,zmax

  implicit none

  real, intent(IN) :: total_time
  integer, intent(IN):: ts, numproc, nphotons

  real :: density, c, mass
  integer :: i,j,k

  density= 1.046 !g/cm3
  c = 3.630 ! average specific heat capacity J/g/degC
  mass=density * (2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg) !g

  if (ts .eq. 1)then
    temp=37.
  endif

  do k = 1, nzg
      do j = 1, nyg
          do i = 1, nxg

            if(k .eq. 77)then
              temp_slice(i,j,ts)=temp(i,j,k)
            endif

            E_abs(i,j,k)=si_step_GLOBAL(i,j,k)*ua_brain(i,j,k)&
            *(2.d0/(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg))) !j/cm3

            temp_change(i,j,k)=(E_abs(i,j,k))*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)/(mass * c)

            temp(i,j,k)= temp(i,j,k)+ temp_change(i,j,k)

            if(si_step_GLOBAL(i,j,k)*ua_brain(i,j,k) .gt. 0.)then
          !  print*, mass,si_step_GLOBAL(i,j,k), E_abs(i,j,k), temp_change(i,j,k), temp(i,j,k)
          !  print*, (2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg), xmax, nxg, nphotons, numproc

          endif

          end do
      end do
  end do









end subroutine temp_calc



end MODULE subs

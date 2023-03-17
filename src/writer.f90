module writer_mod

implicit none

    contains
        subroutine writer(xmax, ymax, zmax, nphotons, numproc,ts,total_time)

            use constants, only : nxg,nyg,nzg,fileplace, end_time, start_time
            use iarray!,   only : jmean, esc_flur , alb_absor, pl_absor, np_cell, flu_phot_rel, obs_flur,jme, con_test!), si_step_GLOBAL

            implicit none

            integer           :: nphotons, i, u, numproc, j,k
            real              :: xmax, ymax, zmax, jm(nzg)
            character(len=70) :: fn

            integer, intent(IN):: ts
            real, intent(IN):: total_time


!jmeanGLOBAL=jmeanGLOBAL*(33/(nphotons*(ts)*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))
jmeanGLOBAL=jmeanGLOBAL*(2000.d0/(nphotons*(ts)*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg))) !power =2W, 2000mW
!jmeanGLOBAL=jmeanGLOBAL*(30.d0/(nphotons*(ts)*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg))) !power =2W, 2000mW

!print*, 'writer', jmean(100,100,100), jme(100,100,100)

          !  print*, 'step=',si_step_GLOBAL(190,190,190)

            inquire(iolength=i)jmeanGLOBAL

            open(newunit=u,file=trim(fileplace)//'jmean_brain_code26.dat',access='stream'&
            ,status='REPLACE',form='unformatted')
            write(u) jmeanGLOBAL
            close(u)

            open(newunit=u,file=trim(fileplace)//'rhokap_brain_code26.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) rhokap
            close(u)

            open(newunit=u,file=trim(fileplace)//'albedo_brain_code26.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) albedoar
            close(u)

            open(newunit=u,file=trim(fileplace)//'hgg_brain_code26.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) hggar
            close(u)

          !  open(newunit=u,file=trim(fileplace)//'jmean/np_cell_full.dat',access='stream',status='REPLACE',form='unformatted')
          !  write(u) np_cell
          !  close(u)

          !  open(newunit=u,file=trim(fileplace)//'esc_flur.dat',access='stream',status='REPLACE',form='unformatted')
          !  write(u) esc_flur_GLOBAL
          !  close(u)

          !  open(newunit=u,file=trim(fileplace)//'pl_absor.dat',access='stream',status='REPLACE',form='unformatted')
          !  write(u) pl_absor_GLOBAL
          !  close(u)


            open(newunit=u,file=trim(fileplace)//'percent_left_code26.dat',access='stream',&
            status='REPLACE',form='unformatted')
            write(u) percent_left
            close(u)

            open(newunit=u,file=trim(fileplace)//'so_tot_code26.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) so_tot
            close(u)

            open(newunit=u,file=trim(fileplace)//'o21_tot_code26.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) o21_tot
            close(u)

            open(newunit=u,file=trim(fileplace)//'o23_tot_code26.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) o23_tot
            close(u)

            open(newunit=u,file=trim(fileplace)//'so_slice26.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) S0_slice
            close(u)

            open(newunit=u,file=trim(fileplace)//'o23_slice26.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) o23_slice
            close(u)

            open(newunit=u,file=trim(fileplace)//'o21_slice26.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) o21_slice
            close(u)

            open(newunit=u,file=trim(fileplace)//'tumkill_slice26.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) tumkill_slice
            close(u)

            open(newunit=u,file=trim(fileplace)//'temp_slice26.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) temp_slice
            close(u)


            !open(newunit=u,file=trim(fileplace)//'av_ppix_con_630.dat',access='stream',status='REPLACE',form='unformatted')
            !write(u) av_ppix_con_630
            !close(u)

          !  open(newunit=u,file=trim(fileplace)//'jmean/alb_absor.dat',access='stream',status='REPLACE',form='unformatted')
          !  write(u) alb_absor
          !  close(u)

          !  open(newunit=u,file=trim(fileplace)//'jmean/flu_phot_rel.dat',access='stream',status='REPLACE',form='unformatted')
          !  write(u) flu_phot_rel
          !  close(u)

          !if(ts .eq. 1)then
      !    if(ts .eq.1 .or. ts .eq. 30 .or. ts .eq. 60 .or. ts .eq. 120 .or.&
      !     ts .eq. 180 .or. ts .eq. 240 .or. ts .eq. 300 &!) then! &
      !    .or. ts .eq. 420 .or. ts .eq. 480 .or. ts .eq. 540 .or. ts .eq. 558 .or. ts .eq. 600 &!.or. ts .eq. 1200 .or. ts .eq. 1800 &
      !    .or. ts .eq. 720 .or. ts .eq. 900 .or. ts .eq. 1020 .or. ts .eq. 1200  )then
      !    .or. ts .eq. 2400 .or. ts .eq. 3000 .or. ts .eq. 3600)then
      !    write(fn,fmt='(i0,a)') ts, 'so.dat'
      !   open(newunit=u,file=trim(fileplace)//fn, status="replace", access='stream', form='unformatted')
      !    write(u) S_0
      !   close(u)

      !   write(fn,fmt='(i0,a)') ts, 'o23.dat'
      !  open(newunit=u,file=trim(fileplace)//fn, status="replace", access='stream', form='unformatted')
      !   write(u) O2_3
      !  close(u)

      !  write(fn,fmt='(i0,a)') ts, 'o21_code14.dat'
      ! open(newunit=u,file=trim(fileplace)//fn, status="replace", access='stream', form='unformatted')
    !    write(u) O2_1
    !   close(u)

    !   write(fn,fmt='(i0,a)') ts, 'o23_code14.dat'
    !  open(newunit=u,file=trim(fileplace)//fn, status="replace", access='stream', form='unformatted')
    !   write(u) O2_3
    !  close(u)

    !  write(fn,fmt='(i0,a)') ts, 's0_code14.dat'
     !open(newunit=u,file=trim(fileplace)//fn, status="replace", access='stream', form='unformatted')
      !write(u) S_0
     !close(u)

    !   write(fn,fmt='(i0,a)') ts, 'jmean.dat'
    !  open(newunit=u,file=trim(fileplace)//fn, status="replace", access='stream', form='unformatted')
    !   write(u) jmeanGLOBAL
    !  close(u)

      !write(fn,fmt='(i0,a)') ts, 'tum_kill_code15.dat'
     !open(newunit=u,file=trim(fileplace)//fn, status="replace", access='stream', form='unformatted')
    !  write(u) tumour_killed
      !  write(u) balloon_killed
     !close(u)




      !  endif

            jm = 0.d0
            do k = 1, nzg
                do j = 1, nyg
                    do i = 1, nxg
                        jm(k) = jm(k) + jmeanGLOBAL(i,j,k)
                    end do
                end do
            end do

            jm = jm / (nxg*nyg)

            open(newunit=u,file=trim(fileplace)//"validation-630.dat",status="replace")
            do i = nzg,1,-1
                write(u,*)real(nzg-i)*(2./nzg),jm(i)
            end do
            close(u)
        end subroutine writer
end module writer_mod

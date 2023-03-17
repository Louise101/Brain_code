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

!jmeanGLOBAL =jmeanGLOBAL * ((2.*xmax)**2./(nphotons*numproc*((end_time-start_time)+1)*(2.*xmax/nxg) *(2.*ymax/nyg)*(2.*zmax/nzg)))
jmeanGLOBAL=jmeanGLOBAL*(0.1*(2.*xmax)**2./(nphotons*(ts-1)*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))

!print*, 'writer', jmean(100,100,100), jme(100,100,100)

          !  print*, 'step=',si_step_GLOBAL(190,190,190)

            inquire(iolength=i)jmeanGLOBAL

            open(newunit=u,file=trim(fileplace)//'jmean-test.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) jmeanGLOBAL
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

            open(newunit=u,file=trim(fileplace)//'obs_flur.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) obs_flur_GLOBAL
            close(u)

            open(newunit=u,file=trim(fileplace)//'con_test.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) con_test
            close(u)

            open(newunit=u,file=trim(fileplace)//'whitematter.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) white_matter
            close(u)

            open(newunit=u,file=trim(fileplace)//'greymatter.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) grey_matter
            close(u)

            open(newunit=u,file=trim(fileplace)//'tumour.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) tumour
            close(u)

            open(newunit=u,file=trim(fileplace)//'rhokap.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) rhokap
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

          !if(ts .eq. 60 .or. ts .ge. 120 .and. ts .le. 128 .or. ts .eq. 138 .or. ts .eq. 148)then
          if(ts .eq. 300 .or. ts .eq. 600 )then
        !  if(ts .eq. 60 .or. ts .eq. 120 .or. ts .eq. 180 .or. ts .eq. 240 .or. ts .eq. 300 .or. ts .eq. 360 &
        !  .or. ts .eq. 420 .or. ts .eq. 480 .or. ts .eq. 540 .or. ts .eq. 600 .or. ts .eq. 1200 .or. ts .eq. 1800 &
        !  .or. ts .eq. 2400 .or. ts .eq. 3000 .or. ts .eq. 3600)then
          write(fn,fmt='(i0,a)') ts, 'pdd.dat'

        !   open it with a fixed unit number
          open(newunit=u,file=trim(fileplace)//fn, status="replace", access='stream', form='unformatted')

          ! write something
          write(u) PDD_GLOBAL

          ! close it
          close(u)

        endif

            jm = 0.d0
            do k = 1, nzg
                do j = 1, nyg
                    do i = 1, nxg
                        jm(k) = jm(k) + jmeanGLOBAL(i,j,k)
                    end do
                end do
            end do

            jm = jm / (nxg*nyg)

            open(newunit=u,file=trim(fileplace)//"validation-630_1.dat",status="replace")
            do i = nzg,1,-1
                write(u,*)real(nzg-i)*(2./nzg),jm(i)
            end do
            close(u)
        end subroutine writer
end module writer_mod

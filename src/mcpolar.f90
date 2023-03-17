program mcpolar

 use mpi

!shared data
use constants
use photon_vars
use iarray
use opt_prop

!subroutines
use subs
use gridset_mod
use sourceph_mod
use inttau2
use ch_opt
use stokes_mod
use writer_mod

implicit none

integer          :: nphotons, iseed, j, xcell, ycell, zcell
logical          :: tflag
double precision :: nscatt, flnscatt
real             :: ran, delta, start,finish,ran2, im_side

integer :: id, error, numproc, ts,u!, P
real    :: nscattGLOBAL, total_time, time_step


call cpu_time(start)

!set directory paths
call directory

!allocate and set arrays to 0
call alloc_array
call zarray

!init MPI
 call MPI_init(error)
 call MPI_Comm_size(MPI_COMM_WORLD, numproc, error)
 call MPI_Comm_rank(MPI_COMM_WORLD, id, error)
!id = 0
!numproc = 1


!**** Read in parameters from the file input.params
open(10,file=trim(resdir)//'input.params',status='old')
   read(10,*) nphotons
   read(10,*) xmax
   read(10,*) ymax
   read(10,*) zmax
   read(10,*) im_side
   read(10,*) n1
   read(10,*) n2
   close(10)

!****************** Read in brain tumour data grids*************************************
open(newunit=u,file='grey_matter_cut.txt',status='old')!status='old')
!read(u,end=102) grey_matter
read(u,*,end=102) grey_matter
102 close(u)

open(newunit=u,file='white_matter_cut.txt',status='old')
read(u,*,end=101) white_matter
101 close(u)

open(newunit=u,file='tumour_cut.txt',status='old')
read(u,*,end=103) tumour
103 close(u)

!print*, 'whitematter',maxval(white_matter)
!print*, 'greymatter',size(grey_matter)
!print*, 'tumour',maxval(tumour)
!read in ppix_con file for 420nm simulation
!open(100,file=trim(fileplace)//'av_ppix_con_630.dat',access='stream')
!read(100)av_ppix_con_630
!close(100)
!print*,'ppix con 630', av_ppix_con_630(120,100)

! set seed for rnd generator. id to change seed for each process
iseed=-95648323+id
iseed=-abs(iseed)  ! Random number seed must be negative for ran2

total_time=0.
ts=1.

call ppix_con(total_time,nphotons,numproc,ts) !set initial ppix concentration
!call ppix_con_420(ts)


!*********** NEW TIME STEP STARTS HERE ****************************
do ts=start_time,end_time ! start loop over time step
   !set seed for rnd generator. id to change seed for each process, ts to change seed for each time step
  iseed=-95648324+id+ ts
  !iseed=-abs(iseed)

!flu_image=0. !set the flouresscence image back to 0. for the new time step
call init_opt

if(id == 0)then
   print*, ''
   print*,'# of photons to run',nphotons*numproc
end if

!***** Set up density grid *******************************************
call gridset(id)

!***** Set small distance for use in optical depth integration routines
!***** for roundoff effects when crossing cell walls
delta = 1.e-8*(2.*zmax/nzg)
nscatt=0
flnscatt=0.

!set escaped flouroescnce to 0 for each new time step
!esc_flur=0.
!esc_flur_GLOBAL=0.
!pl_absor=0.
!pl_absor_GLOBAL=0.

! call cpu_time(start)

!loop over photons
print*,'Photons now running on core: ',id

!set path length arrays to 0 at the beginning of the new time step
si_step=0.
si_step_GLOBAL=0.
del_Q=0.
del_pdd=0.
!time_tot=0.


do j = 1, nphotons


   tflag=.FALSE.

   if(mod(j,10000) == 0)then
      print *, j,' scattered photons completed on core: ',id
   end if

!***** Release photon from point source *******************************
   call evenDis(xcell,ycell,zcell,iseed)
!****** Find scattering location

      call tauint1(xcell,ycell,zcell,tflag,iseed,delta)
!******** Photon scatters in grid until it exits (tflag=TRUE)

!here after else - set up new optical property grid with 705nm and run stokes and tauint1 until it absorbes or escapes the grid.
      do while(tflag.eqv..FALSE.)

         ran = ran2(iseed)

         if(ran < albedoar(xcell,ycell,zcell))then!interacts with tissue
               call stokes(iseed)
               nscatt = nscatt + 1
            else !is absorbed and releases a fluorescence photon at 705nm - add 1 to fluoes list and store starting position.
             !store position of absorption in array size N (can remove zeros from end later from where photons escape grid)
               tflag=.true.
               alb_absor(zcell)=alb_absor(zcell)+1
               exit
         end if

!************ Find next scattering location

         call tauint1(xcell,ycell,zcell,tflag,iseed,delta)

      end do
end do      ! end loop over nph photons

call mcpolar_flu(delta, flnscatt, id, iseed, nphotons, numproc, im_side,ts,total_time) !runs fluoresence simulation

jmeanGLOBAL = jmean
nscattGLOBAL = nscatt
! collate fluence from all processes
 !call mpi_reduce(jmean, jmeanGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
 !call mpi_reduce(nscatt,nscattGLOBAL,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
 !call mpi_reduce(si_step, si_step_GLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
 !call mpi_reduce(esc_flur, esc_flur_GLOBAL, (nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
 !call mpi_reduce(pl_absor, pl_absor_GLOBAL, (nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
 !call mpi_reduce(obs_flur, obs_flur_GLOBAL, ((end_time-start_time)+1),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)


 !!mpi reduce flourescne stuff before writing results to file!!
 !calculate time step length
 call time_step_calc(total_time,time_step,nphotons,xmax,ymax,zmax,numproc)

 !calculate new ppix concentration over whole grid
 call ppix_con(total_time,nphotons,numproc,ts)
 !call ppix_con_420(ts)

 !calculate photodynamic dose over the whole grid
 call pdd_calc(total_time,nphotons,numproc)

 !write out files
! call writer(xmax,ymax,zmax,nphotons, numproc,ts)






!jmeanGLOBAL = jmean
!nscattGLOBAL = nscatt

!jmean=jmean*(0.1*(2.*xmax)**2./(nphotons*(ts)*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))
! collate fluence from all processes
 !call mpi_reduce(jmean, jmeanGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
 call mpi_reduce(nscatt,nscattGLOBAL,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
! call mpi_reduce(si_step, si_step_GLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
 call mpi_reduce(esc_flur, esc_flur_GLOBAL, (nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
 call mpi_reduce(pl_absor, pl_absor_GLOBAL, (nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
 call mpi_reduce(obs_flur, obs_flur_GLOBAL, ((end_time-start_time)+1),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
 !call mpi_reduce(con_test, con_test_GLOBAL, ((end_time-start_time)+1),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)

  call mpi_reduce(PDD, PDD_GLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)

 !open(newunit=u,file=trim(fileplace)//'obs_flur.dat',access='stream',status='REPLACE',form='unformatted')
 !write(u) obs_flur!_GLOBAL
 !close(u)

 if(id == 0)then
    print*,'Average # of scatters per photon:',nscattGLOBAL/(nphotons*numproc)
    call writer(xmax,ymax,zmax,nphotons, numproc,ts,total_time)
    print*,'write done'
 endif

end do ! end loop over time step

!jmeanGLOBAL = jmean
!nscattGLOBAL = nscatt
! collate fluence from all processes
 call mpi_reduce(jmean, jmeanGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
 !call mpi_reduce(nscatt,nscattGLOBAL,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
! call mpi_reduce(si_step, si_step_GLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
 !call mpi_reduce(esc_flur, esc_flur_GLOBAL, (nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
 !call mpi_reduce(pl_absor, pl_absor_GLOBAL, (nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
 !call mpi_reduce(obs_flur, obs_flur_GLOBAL, ((end_time-start_time)+1),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
 !call mpi_reduce(con_test, con_test_GLOBAL, ((end_time-start_time)+1),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)

 if(id == 0)then
    print*,'Average # of scatters per photon:',nscattGLOBAL/(nphotons*numproc)
    call writer(xmax,ymax,zmax,nphotons, numproc,ts,total_time)
    print*,'write done'
 endif

!if(id == 0)then
!   print*,'Average # of scatters per photon:',nscattGLOBAL/(nphotons*numproc)
!   call writer(xmax,ymax,zmax,nphotons, numproc,ts,total_time)
!   print*,'write done'
!endif

call cpu_time(finish)
if(finish-start.ge.60. .and. id==0)then
   print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
else
   if(id==0)print*, 'time taken ~',floor(finish-start/60.),'s'
end if

 call MPI_Finalize(error)
end program mcpolar

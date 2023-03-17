MODULE constants
!
! Module containing constants:
!         PI,TWOPI, the number of grid elements in each direction n%g,
!         the various bin # params,
!         and the vars for the filepaths.


    implicit none
    save


    integer, parameter :: nxg=236, nyg=206, nzg=150, im_xvox=200, im_yvox=200, start_time=1, end_time=1 !time is actually just time step number, 153 steps should produce full jaques plot

    !integer, parameter :: nxg=320, nyg=320, nzg=176, im_xvox=200, im_yvox=200, start_time=1, end_time=1 !time is actually just time step number, 153 steps should produce full jaques plot
    real,    parameter :: PI=4.*atan(1.), TWOPI=2.*PI, OFFSET=1.e-2*(2.*.5/nxg), FOURPI=4.*PI
    real               :: xmax, ymax, zmax, v(3), costim, sintim, cospim, sinpim,ful_p
    character(len=255) :: cwd, homedir, fileplace, resdir

end MODULE constants

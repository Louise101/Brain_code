module sourceph_mod

    implicit none

    contains
      subroutine fibre_lille(xcell, ycell, zcell, iseed)
  ! Emits photons evenly across the top of the grid

      use constants, only : nxg, nyg, nzg, xmax, ymax, zmax, TWOPI, PI
      use photon_vars

      implicit none

      integer, intent(OUT)   :: xcell, ycell, zcell
      integer, intent(INOUT) :: iseed

      real :: ran2, theta, dd,d, ranran


          dd=ran2(iseed)*3.9
          d=dd/2.

          ranran=ran2(iseed)

        !  if(ranran .le. 0.5)then
        !  xp=(xmax) -( 4.06)!0-(2.*xmax-4.06)!0.
        !  yp=(ymax) - (2.4) +d!(0.1-d)!-(2.*ymax-2.0)!0.1 -d
        !  zp=zmax-6.37!0-(2.*zmax-6.37)!0.
        !else
        !  xp=(xmax) -( 4.06)!0-(2.*xmax-4.06)!0.
        !  yp=(ymax) - (2.4) -d!(0.1+d)!-(2.*ymax-2.0)!0.1 +d
        !  zp=zmax-6.37!0-(2.*zmax-6.37)!0.
        !endif

        if(ranran .le. 0.5)then !see moving oxygen lab notes for new dimension calcs
        xp=0-0.2!(xmax) !-( 1.4)!0-(2.*xmax-4.06)!0.
        yp=0 - 0.5 +d!(ymax)! - (0.88) +d!(0.1-d)!-(2.*ymax-2.0)!0.1 -d
        zp=0!zmax!-2.22!0-(2.*zmax-6.37)!0.
      else
        xp=0-0.2!(xmax)! -( 1.4)!0-(2.*xmax-4.06)!0.
        yp=0-0.5-d!(ymax)! - (0.88) -d!(0.1+d)!-(2.*ymax-2.0)!0.1 +d
        zp=0!zmax!-2.22!0-(2.*zmax-6.37)!0.
      endif



           cost=2.*ran2(iseed)-1.
           sint=(1.-cost*cost)
           if(sint .le. 0.)then
             sint=0.
           else
             sint=sqrt(sint)
           end if

           phi=TWOPI*ran2(iseed)
           cosp=cos(phi)
           sinp=sin(phi)



              nxp = sint * cosp
              nyp = sint * sinp
              nzp = cost


          !*************** Linear Grid *************************
          xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
          ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
          zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
          !*****************************************************



      end subroutine fibre_lille
        subroutine sourceph(xcell, ycell, zcell, iseed)
        ! get intial photon position


            use constants, only : nxg, nyg, nzg, xmax, ymax, zmax
            use photon_vars

            implicit none


            integer, intent(OUT)   :: xcell, ycell, zcell
            integer, intent(INOUT) :: iseed
            real                   :: ran2

            zp = zmax - epsilon(1.d0)
            xp = 2.*xmax*ran2(iseed)-xmax!ranu(-xmax, xmax, iseed)
            yp = 2.*xmax*ran2(iseed)-xmax!ranu(-ymax, ymax, iseed)

            phi = 0.
            cosp = 0.d0
            sinp = 0.d0
            cost = -1.d0
            sint =  0.d0

            nxp = sint * cosp
            nyp = sint * sinp
            nzp = cost

            !*************** Linear Grid *************************
            xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
            ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
            zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
        !*****************************************************
        end subroutine sourceph

        subroutine evenDis(xcell, ycell, zcell, iseed)
   	! Emits photons evenly across the top of the grid

   	use constants, only : nxg, nyg, nzg, xmax, ymax, zmax, TWOPI
        use photon_vars

        implicit none

        integer, intent(OUT)   :: xcell, ycell, zcell
        integer, intent(INOUT) :: iseed

        real :: ran2, theta

        	 !changed zp as photons were being inputted from the bottom of the grid instead of the top.
        	zp = zmax + (1.d-5 * (2.d0*zmax/1190))



        	if(ran2(iseed) .gt. 0.5)then
        	xp=-ran2(iseed)*xmax
        	else
        	xp=ran2(iseed)*xmax
        	end if


        	if(ran2(iseed) .gt. 0.5)then
        	yp=-ran2(iseed)*ymax
        	else
        	yp=ran2(iseed)*ymax
        	end if


   		   phi = TWOPI * ran2(iseed)
            cosp = cos(phi)
            sinp = sin(phi)
            cost = -1.d0 !direct irradiation
       !     cost=ran2(iseed) !diffuse irradiation
            sint = sqrt(1. - cost**2)

            nxp = sint * cosp
            nyp = sint * sinp
            nzp = cost

            !*************** Linear Grid *************************
            xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
            ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
            zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
            !*****************************************************



        end subroutine evenDis

        subroutine fibre(xcell, ycell, zcell, iseed)
! Emits photons evenly across the top of the grid

use constants, only : nxg, nyg, nzg, xmax, ymax, zmax, TWOPI, PI
use photon_vars

implicit none

integer, intent(OUT)   :: xcell, ycell, zcell
integer, intent(INOUT) :: iseed

real :: ran2, theta, dd,d,rot_theta, ranran



    dd=ran2(iseed)*1. !1cm in length
    d=dd/2.

    ranran=ran2(iseed)

    if(ran2(iseed) .le. 0.5)then
    xp=0.
    yp=0.
    zp=0.5 -d
  else
    xp=0.
    yp=0.
    zp=0.5 +d
  endif

  !xp=0
  !yp=0
  !zp=0

     cost=2.*ran2(iseed)-1.
     sint=(1.-cost*cost)
     if(sint .le. 0.)then
       sint=0.
     else
       sint=sqrt(sint)
     end if

     phi=TWOPI*ran2(iseed)
     cosp=cos(phi)
     sinp=sin(phi)


        nxp = sint * cosp
        nyp = sint * sinp
        nzp = cost




    !*************** Linear Grid *************************
    xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
    ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
    zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
    !*****************************************************



end subroutine fibre


        real function ranu(a, b, iseed)

            implicit none


            real, intent(IN)       :: a, b
            integer, intent(INOUT) :: iseed

            real :: ran2

            ranu = a + ran2(iseed) * (b - a)

        end function ranu
end module sourceph_mod

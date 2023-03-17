module sourceph_mod

    implicit none

    contains
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
        	zp = zmax - (1.d-5 * (2.d0*zmax/1190))



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


        real function ranu(a, b, iseed)

            implicit none


            real, intent(IN)       :: a, b
            integer, intent(INOUT) :: iseed

            real :: ran2

            ranu = a + ran2(iseed) * (b - a)

        end function ranu
end module sourceph_mod

MODULE gridset_mod

    implicit none

    private
    public :: gridset

    contains
        subroutine gridset(id)

            use constants, only : nxg, nyg, nzg, xmax, ymax, zmax
            use iarray, only    : rhokap,xface,yface,zface, rhokap, refrac, grey_matter, white_matter, albedoar
            use opt_prop,
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
            do i = 1, nxg + 1
                xface(i) = (i - 1) * 2. * xmax/nxg
            end do

            do i = 1, nyg + 1
                yface(i) = (i - 1) * 2. * ymax/nyg
            end do

            do i = 1, nzg + 2
                zface(i) = (i - 1) * 2. * zmax/nzg
            end do

            call init_opt
            refrac(:,:,:) = n1

            !set up optical properties grid
            do i = 1, nxg
                x = xface(i) - xmax + xmax/nxg
                do j = 1, nyg
                    y = yface(j) - ymax + ymax/nyg
                    do k = 1, nzg
                        z = zface(k) - zmax + zmax/nzg
                      !  rhokap(i,j,k) = grey_matter(i,j,k) !already set in ch_opt.f90
                        refrac(i,j,k) = n2

                        rhokap(i,j,k) = (white_matter(i,j,k)*whitematter_kappa) + (grey_matter(i,j,k)*greymatter_kappa)

                        albedoar(i,j,k) = (white_matter(i,j,k)*whitematter_albedo) &
                        + (grey_matter(i,j,k)*greymatter_albedo)

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

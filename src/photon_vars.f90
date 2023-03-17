MODULE photon_vars
!
! Module containing photon vars:
!           %p:current position of photon
!           n%p:current direction of photon
!           sin/cos(%)/phi:various angles related to photons flight
!           angles measured from ''north'' in thetas case.

    implicit none

    real :: xp,yp,zp,nxp,nyp,nzp,sint,cost,sinp,cosp,phi
    real ::obs_cost,obs_sint,obs_phi, obs_cosp,obs_sinp, obs_nxp, obs_nyp, obs_nzp
    real :: leave_theta, leave_phi

end MODULE photon_vars

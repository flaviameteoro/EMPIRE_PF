module hadcm3_config
  Integer, parameter :: a_nxn=96
  integer, parameter :: a_nyn=73
  integer, parameter :: a_levels=19
  integer, parameter :: a_ntimesteps=48

  Integer, parameter :: o_nxn=290
  integer, parameter :: o_nyn=144
  integer, parameter :: o_levels=20
  integer, parameter :: o_ntimesteps=24
end module hadcm3_config

module hadcm3_data
  real(kind=kind(1.0D0)), allocatable, dimension(:,:,:,:) :: u,v,qt,thetal
  real(kind=kind(1.0D0)), allocatable, dimension(:,:,:) :: pstar
  real(kind=kind(1.0D0)), allocatable, dimension(:,:,:,:) :: t_o,sal,b_u,b_v
  real(kind=kind(1.0D0)), allocatable, dimension(:,:,:) :: mld
end module hadcm3_data

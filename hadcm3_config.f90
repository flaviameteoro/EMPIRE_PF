module hadcm3_config
  Integer, parameter :: nxn=96
  integer, parameter :: nyn=73
  integer, parameter :: levels=19
  integer, parameter :: ntimesteps=48
end module hadcm3_config

module hadcm3_data
  real(kind=kind(1.0D0)), allocatable, dimension(:,:,:)   :: u,v,qt,thetal
  real(kind=kind(1.0D0)), allocatable, dimension(:,:)   :: pstar
end module hadcm3_data

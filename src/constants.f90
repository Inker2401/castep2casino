module constants
  !===============================================================================!
  !                   ____                _              _                        !
  !                  / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___                  !
  !                 | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __|                 !
  !                 | |__| (_) | | | \__ \ || (_| | | | | |_\__ \                 !
  !                  \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/                 !
  ! This module contains the various physical and mathematical constants          !
  ! used throughout this program.                                                 !
  ! ------------------------------------------------------------------------------!
  ! Modules used                                                                  !
  ! standard Fortran environment: iso_fortran_env                                 !
  ! ------------------------------------------------------------------------------!
  ! Public variables:                                                             !
  ! See below                                                                     !
  !===============================================================================!
  use iso_fortran_env, only : real64,int64

  implicit none

  ! Define kind parameter for double precision
  integer,parameter::dp=real64
  integer,parameter::int_dp=int64
  complex(kind=dp),public,parameter::cmplx_0=cmplx(0.0,kind=dp) ! complex 0 for double precision

  real(kind=dp),parameter :: pi = acos(-1.0_dp)   ! Ratio of circumference of a circle to its diameter
  real(kind=dp),parameter :: bohr_radius=0.529177210903_dp ! Bohr radius in Angstroms
  real(kind=dp),parameter :: eV=27.211386245_dp ! eV in 1 Hartree
end module constants

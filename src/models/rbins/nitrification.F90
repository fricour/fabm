#include "fabm_driver.h"

module rbins_nitrification
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_rbins_nitrification
      type (type_state_variable_id) :: id_no3, id_nh4
      real(rk) :: nitriRate, o2, kso2 ! o2 constant because we work at saturation?
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

  subroutine initialize(self, configunit)
    class (type_rbins_nitrification), intent(inout), target :: self
    integer,                          intent(in)            :: configunit

    call self%register_implemented_routines((/source_do/))

    call self%get_parameter(self%nitriRate, 'nitriRate', 'day-1', 'nitrification rate', default=0.1_rk)
    call self%get_parameter(self%o2, 'o2', 'mmol m-3', 'oxygen concentration in total system', default=300._rk)
    call self%get_parameter(self%kso2, 'kso2', 'mmol m-3', 'half-saturation constant of oxygen', default=3._rk)

    call self%register_state_variable(self%id_no3, 'no3', 'mmolN m-3', 'nitrate concentration', minimum=0.0_rk) ! where are the initial values?
    call self%register_state_variable(self%id_nh4, 'nh4', 'mmolN m-3', 'ammonium concentration', minimum=0.0_rk)
  end subroutine initialize

  subroutine do(self, _ARGUMENTS_DO_)
    class (type_rbins_nitrification), intent(in) :: self
    _DECLARE_ARGUMENTS_DO_

    real(rk) :: no3, nh4

    _LOOP_BEGIN_
        ! Obtain concentration of nutrients
        _GET_(self%id_no3, no3)
        _GET_(self%id_no4, no4)

        ! Compute nitrification
        nitrif = self%nitriRate * nh4 * (self%o2/(self%o2+self%kso2))

        ! Send rates of change to FABM.
        _ADD_SOURCE_(self%id_no3,nitrif)
        _ADD_SOURCE_(self%id_nh4,-nitrif)

        ! Send the value of diagnostic variables to FABM.
        _SET_DIAGNOSTIC_(self%id_mu,mu)

        _LOOP_END_
  end subroutine do

end module

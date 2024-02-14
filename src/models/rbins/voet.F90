#include "fabm_driver.h"

module rbins_voet
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_rbins_voet
      type (type_state_variable_id) :: id_p, id_d
      real(rk) :: k_f, k_d
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

  subroutine initialize(self, configunit)
    class (type_rbins_voet), intent(inout), target :: self
    integer,                  intent(in)            :: configunit

    call self%register_implemented_routines((/source_do/))

    call self%get_parameter(self%k_f, 'k_f', 'day-1', 'fixation rate', default=0.2_rk, scale_factor=1.0_rk/86400.0_rk)
    call self%get_parameter(self%k_d, 'k_d', 'day-1', 'dissolution rate', default=0.3_rk, scale_factor=1.0_rk/86400.0_rk)

    call self%register_state_variable(self%id_p, 'p', 'mmolC m-3', 'particulate concentration', initial_value=200.0_rk, minimum=0.0_rk) 
    call self%register_state_variable(self%id_d, 'd', 'mmolC m-3', 'dissolved concentration', minimum=0.0_rk)
    
  end subroutine initialize

  subroutine do(self, _ARGUMENTS_DO_)
    class (type_rbins_voet), intent(in) :: self
    _DECLARE_ARGUMENTS_DO_

    real(rk) :: p, d

    _LOOP_BEGIN_
        ! Obtain concentration of nutrients
        _GET_(self%id_p, p)
        _GET_(self%id_d, d)

        ! Send rates of change to FABM.
        _ADD_SOURCE_(self%id_p, -self%k_d*p + self%k_f*d)
        _ADD_SOURCE_(self%id_d, -self%k_f*d + self%k_d*p)

        _LOOP_END_
  end subroutine do

end module

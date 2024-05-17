#include "fabm_driver.h"

module rbins_base_model
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_rbins_base_model ! it is a coincidence that I chose the same base_model as the type_base_model but it's not at all related...
      type (type_state_variable_id) :: id_p, id_d 
      type (type_dependency_id) :: id_I_0
      real(rk) :: k_p, k_p_par_n
   contains      
      procedure :: initialize
      procedure :: do
   end type

contains

  subroutine initialize(self, configunit)
    class (type_rbins_base_model), intent(inout), target :: self
    integer,                  intent(in)            :: configunit

    real(rk), parameter :: secs_per_day = 86400.0_rk
    real(rk), parameter :: hours_per_day = 3600.0_rk 
    real(rk) :: w_p ! sinking speed

    ! parameters
    call self%get_parameter(self%k_p, 'k_p', 'h-1', 'precipitation rate', default=0.006_rk, scale_factor=1.0_rk/hours_per_day)
    call self%get_parameter(self%k_p_par_n, 'k_p_par_n', '[h-1/(mumol photons/m^2/s)]', 'linear PAR dependency of precipitation rate', default=0.00000176_rk, scale_factor=1.0_rk/hours_per_day)
    call self%get_parameter(w_p, 'w_p', 'm day-1', 'vertical velocity (<0 for sinking)', default=-1.0_rk, scale_factor=1.0_rk/secs_per_day)

    ! state variables 
    call self%register_state_variable(self%id_p, 'p', 'mmolC m-3', 'particulate concentration', initial_value=200.0_rk, minimum=0.0_rk, vertical_movement=w_p) 
    call self%register_state_variable(self%id_d, 'd', 'mmolC m-3', 'dissolved concentration', minimum=0.0_rk)
    
    ! dependencies
    call self%register_dependency(self%id_I_0, standard_variables%downwelling_photosynthetic_radiative_flux) ! units in W/m^2 
  end subroutine initialize

  subroutine do(self, _ARGUMENTS_DO_)
    class (type_rbins_base_model), intent(in) :: self
    _DECLARE_ARGUMENTS_DO_

    real(rk) :: d
    real(rk) :: I_0
    real(rk) :: conversion_factor = 4.6_rk ! conversion factor between W/m^2 and micromol photon/m^2/s  

    _LOOP_BEGIN_
        _GET_(self%id_d, d)
        _GET_(self%id_I_0, I_0)

        ! Send rates of change to FABM.
        !_ADD_SOURCE_(self%id_p, (self%k_p*d + self%k_p_par_n*I_0)*d)
        _ADD_SOURCE_(self%id_p, (self%k_p + self%k_p_par_n*I_0*conversion_factor)*d) 
        _ADD_SOURCE_(self%id_d, -self%k_p*d)

    _LOOP_END_
  end subroutine do

end module

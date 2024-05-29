#include "fabm_driver.h"

module rbins_base_model
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_rbins_base_model ! it is a coincidence that I chose the same base_model as the type_base_model but it's not at all related...
      type (type_state_variable_id) :: id_p, id_d  
      type (type_bottom_state_variable_id) :: id_p_benthos
      type (type_dependency_id) :: id_I_0
      real(rk) :: k_p, k_p_par, reminpart, burialpart, w_p, airconc, piston_velocity 
   contains      
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
      procedure :: do_surface
   end type

contains

  subroutine initialize(self, configunit)
    class (type_rbins_base_model), intent(inout), target :: self
    integer,                  intent(in)            :: configunit

    real(rk), parameter :: secs_per_day = 86400.0_rk
    real(rk), parameter :: secs_per_hour = 3600.0_rk 

    ! parameters (always put a default value)
    call self%get_parameter(self%k_p, 'k_p', 'h-1', 'precipitation rate', default=0.006_rk, scale_factor=1.0_rk/secs_per_hour)
    call self%get_parameter(self%k_p_par, 'k_p_par', '[h-1/(mumol photons/m^2/s)]', 'linear PAR dependency of precipitation rate', default=0.0000017_rk, scale_factor=1.0_rk/secs_per_hour)
    call self%get_parameter(self%w_p, 'w_p', 'm day-1', 'vertical velocity (<0 for sinking)', default=-1.0_rk, scale_factor=1.0_rk/secs_per_day)
    call self%get_parameter(self%reminpart, 'reminpart', 'day-1', default=0.1_rk, scale_factor=1.0_rk/secs_per_day)
    call self%get_parameter(self%burialpart, 'burialpart', 'day-1', default=0.1_rk, scale_factor=1.0_rk/secs_per_day)
    call self%get_parameter(self%airconc, 'airconc', 'mol m-3', 'air concentration of dissolved', default=0.02_rk)
    call self%get_parameter(self%piston_velocity, 'piston_velocity', 'm day-1', 'piston velocity', default=0.1_rk, scale_factor=1.0_rk/secs_per_day)

    ! state variables 
    call self%register_state_variable(self%id_p, 'p', 'molC m-3', 'particulate concentration', initial_value=200.0_rk, minimum=0.0_rk, vertical_movement=self%w_p) 
    call self%register_state_variable(self%id_d, 'd', 'molC m-3', 'dissolved concentration', minimum=0.0_rk)
    call self%register_state_variable(self%id_p_benthos, 'p_benthos', 'molC m-2', 'benthic particulate concentration', initial_value=10.0_rk, minimum=0.0_rk)
 
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
        _ADD_SOURCE_(self%id_p, (self%k_p + self%k_p_par*I_0/conversion_factor)*d) ! we divide by a conversion factor to get back to micromol photon/m^2/s
        _ADD_SOURCE_(self%id_d, -self%k_p*d)

    _LOOP_END_
  end subroutine do

  subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
     class (type_rbins_base_model), intent(in) :: self
     _DECLARE_ARGUMENTS_DO_BOTTOM_

     real(rk) :: p_benthos, p

     _BOTTOM_LOOP_BEGIN_
        _GET_HORIZONTAL_(self%id_p_benthos, p_benthos)
        _GET_(self%id_p, p)
           
        _ADD_BOTTOM_SOURCE_(self%id_p_benthos, -self%reminpart*p_benthos -self%burialpart*p_benthos - self%w_p*p) ! the last term impacts the content of the benthic state variable and the sign is "-" because self%w_p is negative 
        _ADD_BOTTOM_FLUX_(self%id_p, -self%w_p*p) 
        _ADD_BOTTOM_FLUX_(self%id_d, self%reminpart*p_benthos)
 
     _BOTTOM_LOOP_END_
  end subroutine do_bottom 

  subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
     class (type_rbins_base_model), intent(in) :: self
     _DECLARE_ARGUMENTS_DO_SURFACE_

     real(rk) :: d

     _SURFACE_LOOP_BEGIN_
       _GET_(self%id_d, d)

       _ADD_SURFACE_FLUX_(self%id_d, -self%piston_velocity*(d - self%airconc)) 
 
     _SURFACE_LOOP_END_ 

  end subroutine do_surface

end module

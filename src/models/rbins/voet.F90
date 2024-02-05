#include "fabm_driver.h"

module rbins_voet

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_rbins_voet
      ! Add variable identifiers and parameters here.
   contains
      procedure :: initialize
      ! Reference model procedures here.
   end type

contains

   subroutine initialize(self, configunit)
      class (type_rbins_voet), intent(inout), target :: self
      integer,                          intent(in)            :: configunit
 
      ! Register model parameters and variables here.
   end subroutine initialize

   ! Add model subroutines here.

end module

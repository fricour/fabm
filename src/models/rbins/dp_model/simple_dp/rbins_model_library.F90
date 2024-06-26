module rbins_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   ! Add use statements for new models here
   use rbins_voet

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: rbins_model_factory

contains

   subroutine create(self, name, model)

      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         ! Add case statements for new models here
         case ('voet'); allocate(type_rbins_voet::model)
      end select

   end subroutine create

end module

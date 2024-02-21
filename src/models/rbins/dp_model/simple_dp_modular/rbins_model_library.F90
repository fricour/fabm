module rbins_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   ! Add use statements for new models here
   !use rbins_voet
   use rbins_base_model
   use rbins_dissolution	

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
         !case ('voet'); allocate(type_rbins_voet::model)
         case ('base_model'); allocate(type_rbins_base_model::model)
         case ('dissolution'); allocate(type_rbins_dissolution::model)
      end select

   end subroutine create

end module

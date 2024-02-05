module rbins_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use rbins_voet
   ! Add use statements for new models here

   ...

contains

   subroutine create(self, name, model)

      ...

      select case (name)
         case ('voet'); allocate(type_rbins_voet::model)
         ! Add new case statements for new models here
      end select

   end subroutine create

end module

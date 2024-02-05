add_library(fabm_models_rbins OBJECT
            rbins_model_library.F90
            voet.F90
           )

add_dependencies(fabm_models_rbins fabm_base)

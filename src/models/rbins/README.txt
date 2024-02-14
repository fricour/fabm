When you want to create a new biogeochemical model, follow the instructions at https://github.com/fabm-model/fabm/wiki/Developing-a-new-biogeochemical-model UNTIL "Including the model during compilation".

When you want to test your empty model, you need to run FABM. For this, I chose to use the PYTHON host therefore you firstly need to BUILD (https://github.com/fabm-model/fabm/wiki/Building-and-installing) the configuration with the selected host (here: python) with
cmake <FABMDIR> -DFABM_HOST=<HOSTNAME> where <FABMDIR> is where you have the source code and then you build and install with cmake --build . --target install.

If you modify the FABM source code, you only need to repeat the last step (however, I do not know for now if I need to recompile this when I update the files in the new_model folder, we'll see...).

Lastly, it is possible that you need to reinstall pyfabm (when your host is python) so python -m pip install <FABMDIR> when you modify the FABM source again (again, we'll see if it needs to be done every time we add new files in the new_model folder).

--- PYTHON HOST ---

So basically, when you want to use fabm, use the right conda environment such as "conda activate fabm" then build the code with your new biogeochemical model python -m pip install <FABMDIR> then in a python terminal, run the wiki example https://github.com/fabm-model/fabm/wiki/python

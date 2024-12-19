# Template for CMT GiD homework

Use this template as starting point for your GiD Problemtype. 

Inside the annex folder, you can find an example of input files for kratos

Please add your implementation to the following files:

- simplecfd_default.spd: Where the data tree is defined.
- simplecfd.xml: Where the problemtype is defined.
- write.tcl: Place the scripts to generate the mdpa file.
- write_parameters.tcl: Place the scripts to generate the ProjectParameters.json
- write_materials.tcl: Place the scipts to generate the FluidMaterials.json

Check all this files and fill all the TODO sections with your implementation
ORCA Input Generator

This a Django based web application which prepares .inp files for the ORCA computational chemistry software based upon the user's selection of parameters. Ideal for beginners no additional software needs to be downloaded, and a simple input selection page allows for ease of use. In future I plan to add a help page with explanations of various parameters.

Current features:
-Supports molecular energy, optimization, vibrational frequency, NMR, and UV-vs input generation
-Accepts either SMILES or InChi strings as molecular input
-Choice from different groupings of methods and basis sets (easily expandable in future)
-Allows for specification of charge and spin multiplicity of the molecule
-Presents user with a depiction of the molecule to confirm the structure generated is correct
-Downloadable .inp file

Note that the current deployment on dycl1.pythonanywhere.com is behind the version here on github, as I only have the free PythonAnywhere plan.
# Trnsys_Type653_HX_effectiveness
 Contains a python script for the calculation of the HX effectiveness as input to Trnsys Type 653 and an example for the coupling with Type 3157.

This repository contains three files:
* *effectiveness.py* is a script to calculate the effectiveness as an input for the floor heat exchanger Type 653
* *pipeLengthCalculator.py* is a script to calculate the pipe length of a given heat exchanger design
* *example.tpf* is an example Trnsys project where the effectiveness is implemented via Type 3157

**Type 653**: Thornton, J.; Bradley, D.; McDowell, T.; Blair, N.; Duffy, M.; LaHam, N.; Naik, A. TESSLibs 17 Ground Coupling LibraryMathematical Reference. Thermal Energy System Specialists, LLC, 2012.
**Type 3157**: Bernier, N.; Marcotte, B.; Kummert, M. Calling Python from TRNSYS with CFFI. Polytechnique Montréal, 2022. https://doi.org/10.5281/zenodo.6523078.
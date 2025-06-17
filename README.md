# DCDMcode

## Demographic Carbon Distribution Model
**Repository for**: _Modeling Carbon Dynamics Using Free-air CO2 Enrichment Data from Duke Forest_ <br>
Corresponding Author: Daniel Poll (polldb@cofc.edu) <br>
Contributing Authors: Allan Strand, Seth Pritchard, Lucas Vaughn, Katilyn Beidler, and Benton Taylor <br>


## Overview

This repository contains the required data and MATLAB scripts used to generate the figures and tables presented in the paper _Modeling Carbon Dynamics Using Free-air CO2 Enrichment Data from Duke Forest_. <br>

The repository includes:
- [face_dataset_paper.csv](./face_dataset_paper.csv): CSV file containing empirical data sampled from DUKE Forest
- [model_data.csv](./model_data.csv): CSV file containing ecosystem estimates from LANDCARB and CENTURY models, (McCormack 2015)
- [figure1.m](./figure1.m): MATLAB script used to generate panels B-D of Figure 1
- [figure2.m](./figure2.m): MATLAB script used to generate Figure 2
- [figure3.m](./figure3.m): MATLAB script used to generate Figure 3
- [table1.m](./table1.m): MATLAB script used to generate parameter estimates given in Table 1
- figureS1.m: MATLAB script to generate Supplemental Figure 1
- figureS2.m: MATLAB script to generate Supplemental Figure 2


## Software Requirements

- MATLAB (tested on R2024b)
- The following MATLAB toolboxes:
    - Optimization Toolbox
    - Statistics and Machine Learning Toolbox
    - Financial Toolbox
    - Curve Fitting Toolbox
 
Some MATLAB scripts may not require all toolboxes. <br>


## How to Use

1. Clone the repository or download the relevant .csv and .m files

2. Open MATLAB, and verify the folder location where the downloaded files are stored

3. In MATLAB, run the desired script. For example, if generating figure 3, type
~~~
run('figure3.m')
~~~
<br>

## Citing This Work

To be updated. <br> <br>

## License

This code is licensed under the MIT License. See the [LICENSE](./LICENSE) file for details. <br>

## Contact Information
For questions, please contact: <br>

Daniel Poll (polldb@cofc.edu) <br>
College of Charleston <br> 
Department of Mathematics <br>

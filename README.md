# Multiple-Linear-Regression
A tool for multiple linear regression analysis. 

#### Background:
[Multiple linear regression](https://www.statology.org/multiple-linear-regression/)
(MLS) estimates linear relationships  
between a 'dependent' variable and several 'independent' variables.

For photomicography, instead of fitting [high order polynomials](https://github.com/blekenbleu/CorrCA)
or [cubic with radial symmetry](https://lensfun.github.io/calibration-tutorial/lens-tca.html),  
consider multiple linear regression, replacing its linear kernel with a bicubic:  
 &emsp; `CAxy = a + b*x + c*x*x + d*x*x*x + e*y + f*y*y + g*y*y*y + h*abs(y*x)`

While more complex than [typical radially symmetric Transverse Chromatic
Aberration](https://lensfun.github.io/calibration-tutorial/lens-tca.html) (TCA) models:  
 &emsp; `CAr = a * r^4 + b * r^3 + c * r^2 + v * r`  
.. this MLS anticipates photomicography issues:
- possible stage tilt, introducing some LoCA, presenting as asymmetric TCA
- imperfect, misaligned and uncentered optics

#### How does this tool work?:
The program reads a sample data file into numeric matrices (implemented as a C struct);  
it contains dependent and associated independent variable values.  
[Ordinary Least Squares](https://www.datacamp.com/tutorial/ols-regression) (OLS) generates coeffiecients  
relating dependent variables to remaining independent variables.  
Model statistics are then calculated for significance and accuracy of this relationship. 


#### Example output:
```
    ******************************************************************************
    *                                                                            *
    * Title: Multiple Linear Regression                                          *
    * Description:  This program performs multiple linear regression.            *
    * Author:       Oscar Zealley                                                *
    * Instructions: append path to data file on the command line, e.g.           *
    *               caBicubic.exe data/health_data.txt                           *
    *               File must be formatted like this:                            *
    *               - one row of comma-separated variable names (no spaces)      *
    *               - rows of corresponding comma-separated numbers, e.g.        *
    *                                                                            *
    *  Dependent_Variable,Independent_Variable_1,Indendent_Variable_2,...        *
    *  85,143,56.3,...                                                           *
    *  112.5,191,62.5,...                                                        *
    *  .,.,.,...                                                                 *
    *                                                                            *
    * NB: Max number of variables is 10.                                         *
    *                                                                            *
    * Press enter to see results for example data/health_data.txt                *
    *                                                                            *
    ******************************************************************************
```
...*after pressing* `enter`, [statistics are generated and displayed](data/health_data_stats.txt)

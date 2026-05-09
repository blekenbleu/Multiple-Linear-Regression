# Multiple-Linear-Regression
A tool for multiple linear regression analysis. 

#### Background:
[Multiple linear regression](https://www.statology.org/multiple-linear-regression/) (MLS) estimates linear relationships  
between a 'dependent' variable and several 'independent' variables.

For photomicography, instead of fitting [high order polynomials](https://github.com/blekenbleu/CorrCA) or [cubic with radial symmetry](https://lensfun.github.io/calibration-tutorial/lens-tca.html),  
consider multiple linear regression, replacing its linear kernel with a bicubic:  
 &emsp; `CAxy = a + b*x + c*x*x + d*x*x*x + e*y + f*y*y + g*y*y*y + h*y*x`

While more complex than [typical radially symmetric Transverse Chromatic Aberration](https://lensfun.github.io/calibration-tutorial/lens-tca.html) (TCA) models:  
 &emsp; `CAr = a * r^4 + b * r^3 + c * r^2 + v * r`  
.. MLS anticipates photomicography issues:
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
    * Description:  This program takes an inputted data file and performs        *
    *               multiple linear regression analysis on the data.             *
    * Author:       Oscar Zealley                                                *
    * Instructions: Put your data in a .txt file in the same directory as        *
    *               this program. Data must be formatted like this:              *
    *                                                                            *
    *  Dependent Variable, Independent Variable 1, Indendent Variable 2,...      *
    *                   4,                      5,                    8,...      *
    *                   7,                     12,                    5,...      *
    *                   .                       .                     .          *
    *                   .                       .                     .          *
    *                                                                            *
    * NB: Max number of variables is 10.                                         *
    *                                                                            *
    * Press enter to see a demonstration using sample data.                      *
    *                                                                            *
    ******************************************************************************
    
    
    Regression Equation:
    Weight = -127.82 +0.24 Age +3.09 Height
    
      Source |       SS       df       MS                  Number of obs =     237
    ---------+------------------------------               F(  2,   234) =  199.61
       Model |  56233.2543     2  28116.6271               Prob > F      =  0.0000
    Residual |  32960.7605   234  140.857951               R-squared     =  0.6305
    ---------+------------------------------               Adj R-squared =  0.6273
       Total |  89194.0148   236  377.940741               Root MSE      =  140.07
    
    ------------------------------------------------------------------------------
     Weight  |      Coef.   Std. Err.       t     P>|t|       [95% Conf. Interval]
    ---------+--------------------------------------------------------------------
     Age     |  0.2402749  0.05510303    4.3605   0.000       0.132273   0.3482769
     Height  |   3.090048   0.2573415    12.008   0.000       2.585659    3.594437
     Const   |  -127.8199      12.099   -10.565   0.000      -151.5339   -104.1059
    ------------------------------------------------------------------------------
```

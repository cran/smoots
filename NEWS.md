# smoots 1.1.0

- the README file has been adjusted.
- new functions were added for graphical testing of linearity assumptions and 
  for forecasting.
- a rescaling function was added to simplify the transfer of the obtained 
  derivative estimates on the time interval [0, 1] to the actual observation  
  time points.
- the package now uses compiled C++ code within selected functions and needs to
  import the Rcpp package for this purpose.
- due to the implemented C++ code the performance of the functions of version
  1.0.0 has been slightly improved for this version.
- the documentation has been adjusted with respect to the newly introduced 
  functions.
- typos and false or missing information in the documentation have been 
  corrected. 
- the documentation has generally been edited with a stronger focus on text 
  formatting so that actual mathematical formulae are visible in the package's 
  PDF-manual
- the S3 method plot.smoots() has been reworked and allows for more arguments 
  of the standard plot function and thus for more flexibility.
- improved feedback on incorrect input.
- some default values for arguments of functions existing since package version 
  1.0.0 have been adjusted, e.g. p = v + 1 is now the default in gsmooth() as 
  opposed to p = 1 in previous versions.

# smoots 1.0.1

- the Imports field in the DESCRIPTION file was adjusted.
- minor performance improvements were made.
- the example in the README file was adjusted.
- minor changes in the documentation were made.
- minor changes in the print method were made.

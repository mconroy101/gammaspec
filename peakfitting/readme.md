## Peakfit.py

### Instructions:
To run, type ```python peakfit.py <file_name> <fit_mode>```.
If either varable is not specified, they can be input with text prompts.

A plot will appear showing the spectrum to fit. Use the Matplotlib zoom and navigation tools to locate the peak of interest.
Then, **right click** to the left and right of the fitting region to define a region of interest. This is defined by black dashed lines.
Once a ROI has been selected, **right click** on the centroids of the peak(s) to be fitted. There is no hard limit of peaks to fit, but too many may cause the fit to fail.
When all peaks have been selected, press **ENTER** to fit.

The fitted spectrum will then be displayed and the fit parameters output in the terminal.
If you are happy with the fit and wish to save the fit parameters, press **SPACE**.
To save the figure, use the inbuild Matplotlib controls, or press **S**.
To close the fitted spectrum and return to the peak selection for another fit, press **ENTER**.

To close the program at any point press **ESCAPE**.


### Fit Modes:

0) Quadratic background with a step function, skewed Gaussian peaks

  $$F = Q+S+G_{sk}$$
2) Quadratic background with a step function, Gaussian peaks with a skewed Gaussian low energy tail
3) Quadratic background, Gaussian peaks


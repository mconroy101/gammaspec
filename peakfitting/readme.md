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

---

### Fit Modes:

0) Quadratic background with a step function, skewed Gaussian peaks. This is a very flexible fitting routine, but is less physically motivated than 1.
  
   $$F = Q+S+G_{sk}$$

   where

   $$Q=ax^2+bx+c$$
   
   $$S = \frac{N}{\sigma\sqrt{2\pi}} \frac{stp}{1+e^{\frac{x-x_0}{\sigma}}}$$
   
   $$G_{sk} = \frac{N}{\sigma\sqrt{2\pi}} e^{\frac{(x-x_0)^2}{2\sigma^2}} \left(1+erf\left(\frac{skw(x-x_0)}{\sigma\sqrt(2)}\right)\right)$$
  
1) Quadratic background with a step function, Gaussian peaks with a skewed Gaussian low energy tail. This is analagous to RadWare, where each component is physically motivated.

   $$F=Q+S+G+T$$

   where

   $$G=\frac{N}{\sigma\sqrt{2\pi}} e^{\frac{(x-x_0)^2}{2\sigma^2}}$$

   $$T= N e^{\frac{x-x_0}{beta}} erfc\left(\frac{x-x_0}{\sigma\sqrt{2}} + \frac{\sigma}{\beta\sqrt{2}}\right)$$
   
2) Quadratic background, Gaussian peaks. This is a simple routine, good for comparison to Buffit (UoB).

   $$F=Q+G$$


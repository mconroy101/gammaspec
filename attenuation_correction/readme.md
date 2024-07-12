# Foil Correction Factor Calculation
For extended sources, it is necessary to consider the difference in solid angle a detector subtends compared to a point source. Additionally, for thick, high Z sources, low energy gamma rays can be self attenuated in the source and not make it to the detector. A correction factor can be calculated which accounts for both of these effects via a Monte Carlo approach.

Here is a ROOT macro and a Python script which can both compute the correction factor for square foils of holmium and dsyprosium. Other foils can be modelled by using the relevant mass attenuation coefficientes from NIST.

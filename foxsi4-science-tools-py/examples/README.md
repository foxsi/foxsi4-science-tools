# `foxsi-4-science-tools-py` Examples

This code will provide simple examples of making FOXSI images, lightcurves, spectra, and spectrograms from FOXSI-4 Level 3 CdTe (HXR) data. Currently only the image example is included.

The code is quite basic because in all cases this amounts to making a histogram in whichever dimensions you want -- for example, an image is produced by making a 2D histogram in X and Y spatial dimensions.

In each case, be careful to limit to a chosen range in each dimension that you are NOT histogramming.  For example, if you're making an image, you should choose a time and energy range.  If you do not restrict the energy range you'll get all photons across all energies, which is often fine. However, if you do not restrict the time range, you might end up 
including a little bit of pointing drift at the very start and end of the observation.

# Pat's Comps - MEATS

My comps project was the MUSLCES Extension for Atmospheric Transmission Spectroscopy (MEATS). It is an extension of Kevin France's MUSCLES treasury survey which created x-ray - IR SEDs for M and K dwarfs to be used as inputs exoplanet atmospheric models. My project was essentially to create the same data products for hotter, more distant stars which have guaranteed JWST observations.

This directory contains my [comps powerpoint presentation](comps_presentation_v3.pptx), the [submitted version of my paper](MUSCLES_Extension_comps_edition.pdf) (since updated), and a fun animation showing the [absorption of EUV photons](absorption.mp4) due to interstellar H (never included in presentation, sadly).

[`assemble_spectrum.ipynb`](assemble_spectrum.ipynb): This jupyter notebook handles the all the main processes of reading in my data and splicing together all the various pieces. This is a reconstruction of a previous notebook that was lost in a fatal harddrive malfunction (backup your work, folks). Not all of the original functionality has been replicated but it's close!

[`spfuncs.py`](spfuncs.py): This contains some functions that I used often for convolving, splicing, and otherwise manipulating spectral data

[`UV_absorption`](UV_absorption.ipynb): This jupyter notebook is used to create the EUV absorption animation. It calculates the extinction, including the Gaunt factor, of EUV photons due (only) to neutral hydrogen.

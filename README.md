# GS_AT
This repository contains codes for in-line holography muti-height Gerchberg-Saxton phase retrieval with automatic affine transform preprocessing. <br>
Follow the 'main_GSwithAT.m' code-file for full guide how to use this algorithm.

# Processing path
Firstly, each in-line hologram is reconstructed (propagated to the object plane), then 'AutoAffineTransform.m' function is applied to extract and match features between n-th and last reconstruction. Basing on these features, affine transforms (AT) are estimated and applied to the input holograms to correct the xy shift and magnification mismatch between holograms (gif below) <br> <br>
![](https://github.com/MRogalski96/GS_AT/blob/main/github_images/vid.gif) <br> <br>
After this preprocessing step, Gerchberg-Saxton (GS) multi-height algorithm is performed to retrieve object phase with minimized twin image noise comparing to single-frame angular spectrum backpropagation (fig. below)  <br> <br>
![](https://github.com/MRogalski96/GS_AT/blob/main/github_images/ASGS_phase.png) <br> <br>
# Cite as
Affine transform-based twin-image suppression for in-line lensless digital holographic microscopy (M. J. Marzejon, M. Rogalski, M. Trusiak) in Rosen, J., Alford, S., Allan, B., et al, "Roadmap on computational methods in optical imaging and holography [invited]," Appl. Phys. B 130, 166 (2024). https://doi.org/10.1007/s00340-024-08280-3

# Created by
Mikołaj Rogalski, <br>
mikolaj.rogalski.dokt@pw.edu.pl <br>
Institute of Micromechanics and Photonics, <br>
Warsaw University of Technology, Poland <br>

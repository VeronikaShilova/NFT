# NFT
Direct and inverse nonlinear Fourier transform (NFT)

This program performs Nonlinear Fourier Transfrom with an optical signal or potential through solving direct and inverse Zaharov-Shabat problems.

Class Potential contains most commonly used signals in Photonics ("Potential.h"): Satsuma-Yajima, Rectangular, Multisoliton, Dark Soliton signals. These signals are used to perform verification of numerical methods for NFT.

Boffeta-Osborne method is used to find continouos spectrum a signal ("BOsborne.h"). 
To find discrete spectrum, new hybrid method of contour integration is used ("HybridContourIntegration.h").

For inverse INFT two different approaches are used: inner bordering algorithm (TIB) ("InnerBordering.h"), algorithm which is based on transition to differential equations ("DifferentialEquationsINFT.h").

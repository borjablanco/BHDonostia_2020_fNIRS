Title:
Task-based connectivity analysis for functional NIRS data

Project lead:

Borja Blanco, twitter.com/borja_blanco4

Irene Arrieta, twitter.com/irenearrieta3

Project collaborators:

César Caballero-Gaudes, twitter.com/CaballeroGaudes Mattermost: @CesarCaballeroGaudes

Eneko Uruñuela, twitter.com/eurunuela Mattermost: @eurunuela

Registered Brainhack Global 2020 Event:

Brainhack Donostia 2020, San Sebastián-Donostia

Project Description:

The aim of this project is to learn, comprehend and implement implement two types of task-based connectivity analyses for functional near infrared spectroscopy (fNIRS) data, namely generalized psycho-physiological interactions (gPPI) and Dynamic Causal Modelling (DCM). These approaches will be evaluated in fNIRS data collected in 4-month-old infants while they listened to forward and backward speech sentences during sleep. Coding will mostly be in MATLAB (similar to the main programs for fNIRS data analysis), although implementation in Python can be explored.

Data to use:

Link to project repository/sources:

Goals for Brainhack Global 2020:

The goals of the project will be split into two different parts and days.
Days 1-3:

Implement gPPI algorithms based on the current implementation using Parametric Empirical Bayes estimation available in SPM12, and adapt this formulation to deal with the fNIRS data structure. Moreover, alternative implementation of the deconvolution algorithm based on stability selection will be explored also based on current code available in Matlab. Milestone: Compute gPPI at the channel-level and global-level for several datasets.
Days 4 and 5:

Learn and understand the implementation of DCM for fNIRS available in SPM12 (see chapter 46 of SPM12 manual, and related articles).
Using the gPPI results, formulate and implement different models within the DCM framework.
Milestone: Perform DCM data analysis in several datasets and interpret the results. Comparison between gPPI and DCM.
Good first issues:

Recommended Readings:

gPPI:
Gitelman, D.R., Penny, W.D., Ashburner, J. and Friston, K.J., 2003. Modeling regional and psychophysiologic interactions in fMRI: the importance of hemodynamic deconvolution. Neuroimage, 19(1), pp.200-207. https://doi.org/10.1016/S1053-8119(03)00058-2
McLaren, D.G., Ries, M.L., Xu, G. and Johnson, S.C., 2012. A generalized form of context-dependent psychophysiological interactions (gPPI): a comparison to standard approaches. Neuroimage, 61(4), pp.1277-1286. https://doi.org/10.1016/j.neuroimage.2012.03.068
Hassanpour, M.S., Eggebrecht, A.T., Peelle, J.E. and Culver, J.P., 2017. Mapping effective connectivity within cortical networks with diffuse optical tomography. Neurophotonics, 4(4), p.041402. https://doi.org/10.1117/1.NPh.4.4.041402
Gerchen, M.F., Bernal‐Casas, D. and Kirsch, P., 2014. Analyzing task‐dependent brain network changes by whole‐brain psychophysiological interactions: A comparison to conventional analysis. Human brain mapping, 35(10), pp.5071-5082. https://doi.org/10.1002/hbm.22532

DCM:
Tak, S., Kempny, A., Friston, K.J., Leff, A.P. and Penny, W.D., 2015. Dynamic causal modelling for functional near-infrared spectroscopy. Neuroimage, 111, pp.338-349. https://doi.org/10.1016/j.neuroimage.2015.02.035
Bulgarelli, C., Blasi, A., Arridge, S., Powell, S., de Klerk, C.C., Southgate, V., Brigadoi, S., Penny, W., Tak, S. and Hamilton, A., 2018. Dynamic causal modelling on infant fNIRS data: A validation study on a simultaneously recorded fNIRS-fMRI dataset. NeuroImage, 175, pp.413-424. https://doi.org/10.1016/j.neuroimage.2018.04.022
Chapter 46 of SPM12 Manual. https://www.fil.ion.ucl.ac.uk/spm/doc/spm12_manual.pdf

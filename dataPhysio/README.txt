This repository contains characterizations of individual neuronal firing rates in response to stimuli with different binocular disparities. The data herein were recorded across a range of different studies and are largely previously published. Their sources and structure are summarized below.

Dataset: V1V2
This dataset comprises a single Matlab .mat file containing a set of disparity responses recorded from macaque V1 and V2 over the course of a decade in Bruce Cumming’s laboratory at the NIH. The data collection methods are most thoroughly described in:

Read, JCA and Cumming, BG (2003) Testing Quantitative Models of Binocular Disparity Selectivity in Primary Visual Cortex. Journal of Neurophysiology 90: 2795-2817

For each neuron, the data structure contains the following fields:

--visualarea: Either 1 (V1) or 2 (V2)

--dx: The stimulus disparity, in terms of displacement on the monitor in visual degrees. The values are not adjusted for the geometrical horopter. If you plan to make a horopter adjustment, keep in mind that for many of the eccentric RFs we had to move the fixation point to the opposite side of the screen so that the RF was not off screen or too close to the edge. The location of the fixation point is recorded in the “fp” field.

--fp: Location of fixation point on screen. These are mostly the screen center (0,0), but sometimes the fixation was moved to the edge of the screen. Some of these are NaNs, indicating that we don't have the recorded value–in these cases it is likely 0,0. 

--rf: Estimated X,Y coordinates of the RF on screen. These were computed as atan(x/d), where x is displacement on screen, d is distance to nodal point from the screen. Positive x values are leftward and positive y is upward. The RF estimates are mostly based on hand measurements online. You may notice there are a large number of cells with an identical RF location very close to the fovea.  These were all recorded with a Utah array, so the recorded "hand estimate" of the RF is the same for all of them.

--resp: Measured response in units of mean spike count, not divided by duration

--counts: Number of repeated trials

--pval: P-value for a one-way ANOVA. We have already removed all cells that did not pass a one-way ANOVA with p < 0.01.

--duration: Duration of recording. This was 0.5, 0.47 or 0.42 sec (different studies).  


--or: The RF orientation, which defines the axis along which disparities were displaced. If the RF orientation is 90, then this is horizontal disparity. When RF or is -90 the disparity is still horizontal, but the direction is reversed.  So now negative values are uncrossed. In experiments where horizontal disparity was applied, ignoring orientation, these are set to NaN.   

--exptype: If this is a 3 or a 4 it means that the disparity was applied orthogonal to the RF orientation. 



Dataset: MT
This dataset comprises a Matlab .mat file for each of 501 neurons recorded from macaque MT. The data were originally described in:

DeAngelis GC, Uka T. Coding of horizontal disparity and velocity by MT neurons in the alert macaque. J Neurophysiol. 2003 Feb;89(2):1094-111. 

Inside the .mat files, there is a data structure that contains a few variables:
--cellnum: the filename for that recording (which is also in the name of the .mat file)
--disparity: the stimulus disparity for each trial
--speed: the stimulus speed of motion for each trial
--firing_rate: the firing rate for each trial

For the disparity variable, -99 indicates a left-eye (monocular) control condition, +99 indicates a right-eye control condition, and 98 indicates a condition in which binocularly uncorrelated dots were presented. Also, trials that have a -9999 value for either the disparity or speed variables correspond to "null" control trials in which only the fixation point was presented.

There is also an .xlsx file that indicates the angular location of each receptive field in x and y. Positive x is right in the visual field, positive y is up in the visual field.



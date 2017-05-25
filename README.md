# Two Line Element Reader / Orbit Plotter #

Given multiple NORAD Two-Line-Element (TLE) files, this matlab code plots the orbits of the satellites around Earth. This is meant to be a simple orbit plotting / visulization tool. This was developed for orbit visualization during my PhD thesis undertaken in the GPS Research Lab in the Department of Aeronautics and Astronautics at Stanford University: 

[1]	T. G. R. Reid, "Orbital Diversity for Global Navigation Satellite Systems," Doctor of Philosophy, Aeronautics and Astronautics, Stanford University, Stanford, CA, 2017.

This thesis is available at the following link: https://purl.stanford.edu/dc409wn9227

## How to Use ##

'MAIN_TLE_READ_Multifile.m' reads in specified TLE files and plots their orbits around Earth. TLE orbit data be downloaded from the Celestrak website: http://www.celestrak.com/NORAD/elements/master.asp

'MAIN_TLE_READ_Multifile_ALL_TLE.m' does the same as 'MAIN_TLE_READ_Multifile.m' only is preset to read the TLE data for all operational orbits and plot them. 

## More Info ## 

For more info, please refer to the Stanford University GPS Research Lab in the Department of Aeronautics and Astronautics: https://gps.stanford.edu 
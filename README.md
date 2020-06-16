# spectro-lammps

This is a code written primarily by Matthew Mostrom, Dr. Preston Moore, and Andres Saez in Python for the purpose of using classical dynamics to characterize
spectroscopic phenomena of molecular systems. Input currently is made to
resemble LAMMPS input scripts for the purpose of being able to quickly and
easily check the validity of our results with the folks at Sandia for any system
we care to look at.

Running MD:
  Like LAMMPS, our code takes the name of the program to be run, followed by input
  file. A sample input line would look like this:

    python mdlammps.py morse.in

  Sample inputs are kept in the "Models" directory. Currently only NVE has been
  exhaustively tested to match LAMMPS output, and NVT is being worked on.

Correlator:
  One of the main functionalities of this code is to compare Instantaneous Normal
  Mode frequency histograms as a result of the system's Hessian matrix to the 
  observed frequencies revealed by the oscillation of some physical observable.
  Getting these properties involves taking the correlation function of the
  observable over the length of the simulation, and then taking its Fourier
  transform. In the current testing phase, we're only considering the motion of
  the atoms (and not the dipole moment, for example), and are thus correlating
  their velocities. In the future, a similar technique might be used to
  reproduce the IR or Raman spectra by computing dipole moments and
  polarizabilities to make a similar comparison.

  The correlator, vac.py, takes the following arguments:
  
    python vac.py [input_file] [output_file]

  The correlator uses Einstein Summation to arrive at its statistical average,
  and currently has an exponential interpolation to ensure that it is completely
  decorrelated by the end of the simulation window. The FFT used only uses a
  basis of cosine functions, as we are only interested in non-complex
  frequencies. Because of this, it is important to have it start at 1 and end at
  zero, otherwise the DCT signal is left with a residual ring that will make the
  spectrum difficult to interpret.

Distribution Moments:
  For whatever reason, you might want to analyze a set of data by looking at the
  first few moments of its distribution. I made a script that takes a minimum
  and maximum domain window, that takes a noise threshold, and calculates the
  first two moments. Easily adapted to take the third or fourth if needed.

  Syntax:
    
    moment.py [input]

  It outputs the result into moment.dat. I currently have the output data look
  like this:

    XXX [first moment] [second moment]

  so that it can be used in conjunction with run.sh, letting the user run the
  script in multiple directories and concatonate the data with whatever label
  they want by string-replacing "XXX".

Scanning Different Initial Conditions:
  "bash run.sh" is the command I use when I want to run the same simulation with
  one parameter being changed incrementally. It should be invoked in the same
  directory that you want your output data to be copied to, like moment.dat,
  replays, histograms, etc.

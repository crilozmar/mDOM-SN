# mDOM-SN
mDOM - SNe and solar neutrinos

Tools to simulate SNe and solar neutrinos interactions in ice. 

The world is simulated as a cylinder of ice, where the mDOM is deployed in the middle and the neutrinos are coming from one of its faces. This simplifies the weights calculation.

The different classes of primary generator allows to simulate different interactions (mdomPrimaryGeneratorActionX.cc/hh, this is chosen from input parameter --SNgun X):
- 0: Use gps file.
- 1: Simulate elastic scattering of electronic neutrinos from SNe.
- 2: Simulate inverse beta decay of electronic antineutrinos from SNe.
- 3: Simulate elastic scattering of electronic neutrinos from the Sun (#this is not updated!).

mdomPrimaryGeneratorAction 1,2,3 use the class mdomSNTools, where all common stuff are calculated.

A description of the weights, fluxes, cross section, energy distribution used here and many more can be found in https://inspirehep.net/literature/1870591 (to be changed with the latest review comments after acceptation). All used cross section are also referred in the corresponding functions where the distributions are built.

# mdomPrimaryGeneratorAction1/2 (with gfixmeanenergy == false)
Here the simulation of the SN neutrinos takes place. Neutrinos are obviously not directly simulated but only the resulting particle from its interaction (so a positron for inverse beta decay and an electron for elastic scattering). No simulation of neutron capture after inverse beta decay is included (NOTE: it was included some commits ago, I deleted it because I wasn't using it. It was part of a project to study the impact of modifying the ice with Gd to increase the energy of the emitted gamma after the neutron capture.).

Important! Many stuff in these files assume that the SN is in the z-direction, so that's where the neutrinos are coming from. So if it is not, one should change the corresponding here. Also the world should be a cylinder for the weights calculation

Both classes work the same way:
 * In the constructor, the information from the files is read. The files are located in the folder "SNandSolarStuff" and all of them have 4 columns: time in s, flux in s-1 m-2 (but units here do not matter), mean energy <E> in MeV, and mean squared energy <E^2> in MeV^2. The constructor also defines a couple of constant to be used later.
  * The member function "GeneratePrimaries" do everything else:
  
  ** First a particle is created (e- or e+) and a random position in the world is obtained. This position is checked to be contained in the world volume and not inside a module.
  
  ** Then a time of the spectrum is sampled from the fluxes read before, using the intensity of the flux. This and every other sampling in the functions is done with the InverseCumulAlgorithm member function in mdomSNTools.
  
  ** The corresponding <E> and <E^2> at the chosen time t are selected from the read data. For that, a linear interpolation is used.
  
  ** Then the energy distribution is created using <E> and <E^2>. A single energy value is then sampled from this energy distribution. This is the nu/nubar energy.
  
  ** The nu_energy/nu_bar_energy is used to build the angular distribution of produced electrons/positrons. Then a costheta is chosen from this distribution. This is then assuming SN neutrinos coming from the z-direction. A random phi is given.
  
  ** Using costheta and nu_energy/nubar_energy, the energy of the resulting electron/positron is calculated.
  
  ** Later the cross section and weight are obtained. The weight is given by sigma*NTargets*(2*gHeight), where sigma is the total cross section; NTargets is the number of targets for the chosen interaction in ice and gHeight is the height (from the center) of the cylinder.
  
  ** Then the particle is generated and the chosen quantities are given to the analysis manager.

# mdomPrimaryGeneratorAction1/2 (with gfixmeanenergy == true)
  * The energy and alpha parameters of the distribution are given by the user using input parameters, then some parts of the code change but the working principle remains the same.

# How to use this thing
  
  * Compile it (geant4.10.02.p03, GNU 9.3.0). If you use a different compiler you might need to change small stuff here and there.
  * Then run it. "./mdom --help" would list the input parameters that are available (Some input parameters are old/not used anymore. #TODO: Clean it)
  * In order to run SN simulations, you should write something like the following:
  
  ** ./mdom --QE --SNgun 2 --depthpos 88 --nmdoms 1 --SN 4 -w 20 -n 10000 -o youroutputfilepath.txt
  
  *** --QE: Quantum effiency ON.
  
  *** --depthpos X: select a position from the tabulated depths (up to 110). Corresponding depth and absorption/scattering lenghts can be seen directly from the vector in mdomDetectorConstruction.cc around lines 100 (#TODO: Write this into separated files)
  
  *** --nmdoms X: how many mdoms do you want to simulate. Centered by default.
  
  *** -n X: number of particles to simulate
  
  *** -o filename: output filename. 
  
  *** -- SNgun X: Chooses the interaction, or the mdomPrimaryGeneratorAction that will be use (1 = enes, 2 = ibd)
  
  *** -- SN X: Chooses the fluxes/SNe to be used:
  
    **** 0: type II 27 progenitor solar mass SN (sometimes refered as heavy SN in the code)
  
    **** 1: type II 9.6 progenitor solar mass SN (sometimes refered as light SN in the code)
  
    **** 2: type I SN based on DDT scenario
  
    **** 3: type I SN based on GCD scenario
  
    **** 4: type II SN (long tailed)
  
  ** With the flag -v you can open the visualizer. Better do this without tne -n X, and only after the visualizer is opened run/beamOn 1 (only 1 or you'll see nothing just photons)
  
Type I fluxes are a bit weird and I prefer to do not use them (#TODO: To be studied and replaced if necessary).

# How does the output look like (when using --SNgun 1 or 2)

The outputs are written by the mdomAnalysisManager class. Lets say the output file that was given to geant4 is -o data.txt. 2 files will be created: data.txt and data_info.txt.
  * data_info.txt: Contains the general information of the simulated event
  
  ** Time of Flux [s] | Mean energy of nubar | nubar energy | costheta of e+ from z dir | e+ energy | event weigh | Vertex Position (X, Y, Z) [m] | Primary direction (Px,Py,Pz
  
  * data.txt: Contains the hit information. Notice that this file must be read by lines, since each line can have different column size depending on the detected photons.
  
  ** Total hits | PMTs hit |...for each PMT hit...| Module number | PMT number | Hits in that PMT |...for each hit in that PMT... hit time |
  
  * Only detected events are saved. Each line in data.txt correspond to the same line in data_info.txt

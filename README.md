# mDOM-SN
mDOM - SNe and solar neutrinos

Tools to simulate SNe and solar neutrinos interactions in ice. 

The world is simulated as a cylinder of ice, where the mDOM is deployed in the middle and the neutrinos are coming from one of its faces.

The different classes of primary generator allows to simulate different kind of particles:
- 0: Use gps file.
- 1: Simulate elastic scattering of electronic neutrinos from SNe.
- 2: Simulate inverse beta decay of electronic antineutrinos from SNe.
- 3: Simulate elastic scattering of electronic neutrinos from the Sun.

Two type II SNe can be chosen (lighter and heavier ones) based on the ls220 EoS. For the type I SNe, 2 different models (DDT and GCD scenarios) can be chosen.

# LED-IceColumn branch

Used with --SNGun 0. 

Active LEDs are chosen like --LED 1001101010 (more than one is possible at the same time, although maybe not recomended).

IceColumn properties and position can be tunned with input parameters, although for some of them the visualization might not work.

One can set no ice properties for the bulk ice with --bulk_ice 0 and therefore there would be scattering and absorption only for hole ice.

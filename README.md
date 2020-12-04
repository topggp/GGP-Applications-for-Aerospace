# GGP-Applications-for-Aerospace

Based off from the recently published article in Ref[1], GGP is a MATLAB based code that has successfully integrated major geometric projection based topology optimizer, i.e., Geometric Projection (GP), MMC (Moving Morphable Components) and MNA (Moving Node Approach).

## How To Use

Running the code is simple enough! Navigate to the directory in MATLAB file browser where the files have been downloaded, and use the following command in the MATLAB console to run structural topology optimizer:

```bash
GGP_main(100,100,0.3,'MBB','GP')
```
The first two arguements dictate the resolution of the design space (100 x 100 elements in this case), the third argument provides the maximum volume fraction constraint, the 4th argument mentions the boundary conditions and then the method of solving comes last. Read below to know more about the last two arguments.

For heat based application, use the following command:

```bash
GGP_heat(100,100,0.3,'GP')
```
The code is seperated from the main file due to changes in the K-matrix for heat conduction problems.

### Boundary Conditions

The 4th argument for the structural topology optimization is for the boundary condition, i.e., the problem set-up. Of course custom BCs can be used, but for reference, there are several example BCs already used namely:

'MBB' : The classic Messerschmitt-Bolkow-Blohm beam.

'L_Shape' : The L-Shaped Cantilever beam.

'Short_Cantilever' : The simple short Cantilever beam.

'Compliant' : Boundary Conditions and objective functions tailored to provide a simple compliant structure with Linear Mechanics.

'RIB' : A simple application to provide design for a Eppler 420 wing rib, not accounted for multiple load cases seen practically.

'LW' : BC, Objective and Constraint function implemented to obtain a Morphing winglet using Linear mechanics. Uses Objective and Constraint (and their sensitivities) from objective_fcn.m and constraint_fcn.m. Also dedicates 1 component to provide optimum placement of actuator. A work in progress.

### Method
The final argument for both structural and thermal boundary conditions is the method of solving. This can be changed to 3 methods:

'GP' : Also known as Geometric Projection, first envisioned by Julian Norato (Ref[2]).

'MMC' : Moving Morphable Components, by Xu Guo and team (Ref[3])

'MNA' : Moving Node Approach, first conceptualised by J.T. Overvelde (Ref[4]), further developed in house at ISAE SUPAERO, Toulouse (France).

## Note:

1) While GP and MNA can be used largly unchanged, the parameters for MMC must be carefully changed depending upon the resolution. If not, there is a possibility of obtained large fraction of gray areas.
2) Lines 240-247 provides code for loading a file if previously present. If not, It'll create a new one. Can be replaced using a simple fopen() command.
3) Line 248 can be commented if use does not desire to store the variable post code execution.
4) If the user desires to read the output directly on the console, in Line 351, delete the section "f1," in
```bash
fprintf(f1,'It.:%5i ........);
```
5) The program automatically creates a folder with basic parameters of any run, and stores the density and component plot of selected iterations. The files and variables recorded also are stored here.
6) The user is adviced to run a smaller resolution if th intention is to mock execute the code, say 50 x 50. A full run with 100 x 100 elements takes well over 2 hours.

## References

[1. Coniglio, S., Morlier, J., Gogu, C., and Amargier, R., 2019.“Generalized geometry projection:  A unified approach forgeometric feature based topology optimization”.](https://link.springer.com/article/10.1007/s11831-019-09362-8)

[2. Norato, J., Bell, B., and Tortorelli, D. A., 2015.  “A geom-etry projection method for continuum-based topology op-timization with discrete elements”.](https://www.sciencedirect.com/science/article/pii/S0045782515001711?casa_token=Xr892VegDc0AAAAA:89vzo5j0SLHYUh81j6ct9CI6nLxcAElsgHH-j3wqz5d1toX4X8BYiRwC3ZdPUg8Lu_Wyf3BtltM)

[3. Zhang,   W.,   Yuan,   J.,   Zhang,   J.,   and  Guo,   X.,   2016.“A  new  topology  optimization  approach  based  on  mov-ing  morphable  components  (mmc)  and  the  ersatz  mate-rial model”.](https://idp.springer.com/authorize/casa?redirect_uri=https://link.springer.com/content/pdf/10.1007/s00158-015-1372-3.pdf&casa_token=iAKD3Y2P-30AAAAA:yMzRxgj07Jrk8lFPfZERQh7l05SX_PkJFCOmzNqBWRilfAOllY0mJ0dcDsOG7wX5qjq-66Ap8BkrqI2p)

[4. Overvelde,  J.  T.,  2012.    “The  moving  node  approach  intopology optimization”.](https://repository.tudelft.nl/islandora/object/uuid:86c056d8-f368-4239-893f-07ca3a22e112/datastream/OBJ1/download)

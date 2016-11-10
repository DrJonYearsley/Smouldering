# ABC_CA_3layer

Matlab scripts to use sequential Approximate Batesian Computation to infer parameter values for the 3-layer cellular auomata that simulates smouldering combustion of moist peat.

Run the file ABC_moistfire_SMC to run the sequential ABC parameter estimation 

Files are:
+ **ABC_moistfire_SMC.m** run the sequential ABC parameter estimation
+ **ABC_moistfire_simulation.m** run one simulation and save results
+ **ABC_moistfire_setup.m** initialize the parameters for a simulation
+ **ABC_distance.m** calculate the distance between a simulation and the experimental data


## Sequential Monte-Carlo Algorithm

1. Set tolerance criterion to accept a simulation
2. Simulate picking parameters from a prior until tolerance reached.
3. Repeat steps 1 & 2 until n simulations (chains) created
4. Reduce tolerance (advance one time step)
5. Use simulations to create a prior
6. Return to step 1

 This algorithm is described in [Toni et al (2009), Approximate Bayesian computation scheme for parameter inference and model
 selection in dynamical systems, Interface, Vol. 6, p187-202.](http://rsif.royalsocietypublishing.org/content/6/31/187)
 doi:10.1098/rsif.2008.0172

# CA of smouldering combustion that incorporates the heat equation

FIREOX2_Temp_implicit.m is the version of the FIREOX cellular automata that has 2 layers: the temperature and the state of combustion.

The cellular automaton has two layers:

1. A layer giving the temperature of each grid square
2. A layer giving the combustion state of each grid square

All non-local effects are governed by heat diffusion in the temperature layer.  This is slightly different from standard CA which define a neighbourhood around each grid square. The dynamics of the temperature layer are driven by the heat equation with
parameters suitable for peat. Boundary conditionas are that heat can be lost through the sides and from the top of each grid square. The sides and top are considered to be a heat bath of constant temperature. To try and be computationally efficient
yet stable integration scheme the code uses operator splitting with an alternating direction implicit method 

There are 5 combustion states: unburned peat -> pyrolysis -> char -> oxidation -> ash

The transition probabilities from one state to another are a function of the temperature of the grid square. The function is sigmoidal with two parameters, Thalf and A

The transition probabilities have the form

prob of transition = X / (1+X)

where 
+ X = exp(A * (T-Thalf) / Thalf)
+ T is the temperature of the grid square (Kelvin), 
+ Thalf is the temperature when the transition probability is 0.5
+ A is a parameter that determines how quickly the temperature causes the prob of transition to increase.

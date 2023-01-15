# Cahn-Hilliard Systems with dynamic boundary conditions

The Cahn-Hilliard equation describes phase separations in a alloy with two parts.

In this project we will solve Cahn-Hilliard Systems with differnt kinds of boundary values. 
We will look at systems with basic homogenous Neumann values (HN) as well as systems with dynamic boundary values. 
In case of dynamic boundary values this code implements systems with boundary values of Allen-Cahn Type (AC), Liu Wu (LW) and 
Goldstein, Miranville and Schimperna (GMS) boundary conditions.

We will apply the FEM to handle these systems and use the structure of each problem to efficently solve it with the newton method.
With these solutions we can create animations out of the solutions or study the properties like mass or energy of the system we are looking at.

Some of this code uses the special structure of the decoupled weak formulation to set up different spartial and time related discretizations.
This helps to handle the different rates of phase separation of the bulk and the boundary and get a better result.

Finally we want to look at convergence of the systems. We expect order 1 for HN1, AC1, LW1 and GMS1. We obtain these plots by solving a 
inhomogenous problem and compare the solution with the one the code computes.

To understand the code better, there is a operation manual in the repository.

Have fun!

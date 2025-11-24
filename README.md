# CoordFreeNP: Coordinate-free approach to NP calculations

This work is funded by TUBITAK 1001 Program, Grant Number 123R114.

Copy the file **NPDerivatives.sage** into your directory and load it by

**load("NPDerivatives.sage")**

Please see the file **Coordinate_Free_Example.ipynb** as an example.

* Defined parameters:

-- alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau
  
-- lambda, defined as "lambdaa". Displays as LaTeX $\lambda$.
  
-- Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Lambda
  
-- Psi0, Psi1, Psi2, Psi3, Psi4

* Commands:

- General:
  
-- "showall" function: The function to show all NP parameters.

-- "opconj" function: Special function for conjugation. (If conjugate is used, derivatives in the l and n directions also become complex).

- Derivatives:
  
-- Defined as Dl (D), Dn (Capital Delta), Dm (small delta), and Dmbar (small deltabar). 

-- Derivatives can be nested.

- Commutators, such that they equal zero (X input):
  
-- Dn_Dl_com(X): [Dn,Dl]

-- Dm_Dl_com(X): [Dm,Dl]

-- Dm_Dn_com(X): [Dm, Dn]

-- Dmbar_Dm_comNP(X): [Dmbar, Dm]

- NPeqs():
  
-- Calculates and writes the NP equations such that they equal zero. 

-- Outputs can be used in the form NP1, .., NP18.

- Bianchi():

-- Calculates and writes Bianchi identities such that they equal 0. 

-- Outputs can be used in the form BI1, ..., BI11.

- Transformations defined in the Carmeli and Kaye paper (z input):

-- CKAdo(z): A. Null Rotation about l (Do transformation and the parameters will change)

-- CKBdo(z): B. Boost in l-n Plane and Spatial Rotation in m-mbar Plane (Do transformation and the parameters will change)

-- CKCdo(z): C. Null Rotation about n (Do transformation and the parameters will change)

-- CKAsee(z): A. Null Rotation about l (See the results of the transformation and the parameters will not change)

-- CKBsee(z): B. Boost in l-n Plane and Spatial Rotation in m-mbar Plane (See the results of the transformation and the parameters will not change)

-- CKCsee(z): C. Null Rotation about n (See the results of the transformation and the parameters will not change)

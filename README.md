# CoordFreeNP: Coordinate-free approach to NP calculations

This work is funded by TUBITAK 1001 Program, Grant Number 123R114.

Please see the file Kinnersley_NUT_May_5_2024.ipynb as an example.

# * Defined parameters:

- alpha,beta,epsilon,gamma,kappa,mu,nu,pi,rho,sigma,tau
- lambda is defined as "lambdaa". Displays as LaTeX $\lambda$.
- Phi00,Phi01,Phi02,Phi10,Phi11,Phi12,Phi20,Phi21,Phi22,Lambda
- Psi0,Psi1,Psi2,Psi3,Psi4

# * Commands:

- Derivatives:
  
-- Dl (D), Dn (capital Delta), Dm (small delta) and Dmbar (small deltabar). 

--Derivatives can be nested.

- "opconj" function: A special function for conjugation. 
(If "conjugate" is used, l and n-directed derivatives are also complex).

- Commutators, equals zero (X input):

-- Dn_Dl_com(X): [Dn,Dl]

-- Dm_Dl_com(X): [Dm,Dl]

-- Dm_Dn_com(X): [Dm, Dn]

-- Dmbar_Dm_comNP(X): [Dmbar, Dm]

- NPeqs(): 

-- Calculate and write NP equations with =0. 

-- Outputs can be used as NP1, ..., NP18.

- Bianchi():

-- Calculate and write Bianchi equations with =0. 

-- Outputs can be used as BI1, ..., BI11.

- Transformations as defined in the Carmeli and Kaye paper (z input):

-- CKA(z): A. Null Rotation about l

-- CKB(z): B. Boost in l-n Plane and Spatial Rotation in m-mbar Plane

-- CKC(z): C. Null Rotation about n

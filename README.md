# CoordFreeNP: Coordinate-free approach to NP calculations

Please see the file Kinnersley_NUT_May_5_2024.ipynb as an example.

# * Tanimli olan parametreler:

# - alpha,beta,epsilon,gamma,kappa,mu,nu,pi,rho,sigma,tau
# - lambda, "lambdaa" olarak tanimli. Yazarken LaTeX $\lambda$ olarak gösteriyor.
# - Phi00,Phi01,Phi02,Phi10,Phi11,Phi12,Phi20,Phi21,Phi22,Lambda
# - Psi0,Psi1,Psi2,Psi3,Psi4

# * Kullanilabilecek komutlar:

# - Turevler:
# -- Dl (D), Dn (buyuk Delta), Dm (kucuk delta) ve Dmbar (kucuk deltabar) olarak tanimli. 
# -- Turevler ic ice kullanilabilir.

# "opconj" fonksiyonu: Eslenik almak icin ozel fonksiyon. 
# (conjugate kullanilirsa l ve n yonundeki turevler de kompleks oluyor).

# - Komutatorler, esittir sifir olacak sekilde (X girdi):
# -- Dn_Dl_com(X): [Dn,Dl]
# -- Dm_Dl_com(X): [Dm,Dl]
# -- Dm_Dn_com(X): [Dm, Dn]
# -- Dmbar_Dm_comNP(X): [Dmbar, Dm]

# - NPeqs(): 
# -- NP denklemlerini =0 olacak sekilde hesaplar ve yazar. 
# -- Ciktilar NP1, .., NP18 seklinde kullanilabilir.

# - Bianchi():
# -- Bianchi ozdesliklerini =0 olacak sekilde hesaplar ve yazar. 
# -- Ciktilar BI1, .., BI11 seklinde kullanilabilir.

# - Carmeli ve Kaye makalesinde tanimlanan donusumler (z girdi):
# -- CKA(z): A. Null Rotation about l
# -- CKB(z): B. Boost in l-n Plane and Spatial Rotation in m-mbar Plane
# -- CKC(z): C. Null Rotation about n


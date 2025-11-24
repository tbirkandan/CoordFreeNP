# * Defined parameters:

# - alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau
# - lambda, defined as "lambdaa". Displays as LaTeX $\lambda$.
# - Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Lambda
# - Psi0, Psi1, Psi2, Psi3, Psi4

# * Commands:

# - General: 
# -- "showall" function: The function to show all NP parameters.
# -- "opconj" function: Special function for conjugation. (If conjugate is used, derivatives in the l and n directions also become complex).

# - Derivatives:
# -- Defined as Dl (D), Dn (Capital Delta), Dm (small delta), and Dmbar (small deltabar). 
# -- Derivatives can be nested.

# - Commutators, such that they equal zero (X input):
# -- Dn_Dl_com(X): [Dn,Dl]
# -- Dm_Dl_com(X): [Dm,Dl]
# -- Dm_Dn_com(X): [Dm, Dn]
# -- Dmbar_Dm_comNP(X): [Dmbar, Dm]

# - NPeqs(): 
# -- Calculates and writes the NP equations such that they equal zero. 
# -- Outputs can be used in the form NP1, .., NP18.

# - Bianchi():
# -- Calculates and writes Bianchi identities such that they equal 0. 
# -- Outputs can be used in the form BI1, .., BI11.

# - Transformations defined in the Carmeli and Kaye paper (z input):
# -- CKAdo(z): A. Null Rotation about l (Do transformation and the parameters will change)
# -- CKBdo(z): B. Boost in l-n Plane and Spatial Rotation in m-mbar Plane (Do transformation and the parameters will change)
# -- CKCdo(z): C. Null Rotation about n (Do transformation and the parameters will change)
# -- CKAsee(z): A. Null Rotation about l (See the results of the transformation and the parameters will not change)
# -- CKBsee(z): B. Boost in l-n Plane and Spatial Rotation in m-mbar Plane (See the results of the transformation and the parameters will not change)
# -- CKCsee(z): C. Null Rotation about n (See the results of the transformation and the parameters will not change)

##################################################################
from sage.symbolic.operators import add_vararg, mul_vararg
##################################################################

##################################################################
def opconjfunc(self,X):
    if X in RR or X in CC:
        return conjugate(X)
    try:
        myop=X.operator()
    except:
        pass
    try:
        myoprnd=X.operands()[0]
    except:
        myoprnd=X

    if X.operator() == exp:
        return exp(opconj(log(X)))

    if X.operator() == log:
        return log(opconj(exp(X)))

    if X.operator() == operator.pow:
        try:
            a, b = X.operands()
            return (opconj(a))**b
        except:
            if X.operator() == mul_vararg:
                oplist = X.operands()
                return prod(opconj(op) for op in oplist)
            elif X in RR or X in CC:
                    return conjugate(X)

    if X.operator() == mul_vararg:
        oplist=X.operands()
        try:
            return prod(opconj(op) for op in oplist)
        except:
            if X.operator() == add_vararg:
                oplist = X.operands()
                return sum(opconj(op) for op in oplist)

    if myoprnd.is_symbol() and myop==conjugate:
        return conjugate(X)

    if X.operator() == add_vararg:
        oplist = X.operands()
        try:
            return sum(opconj(op) for op in oplist)
        except:
            if X.operator() == mul_vararg:
                oplist = X.operands()
                return prod(opconj(op) for op in oplist)

    
    if myop==conjugate:
        return(conjugate(myoprnd))
    if myop==Dl:
        return Dl(conjugate(myoprnd))
    if myop==Dn:
        return Dn(conjugate(myoprnd))
    if myop==Dm:
        return Dmbar(conjugate(myoprnd))
    if myop==Dmbar:
        return Dm(conjugate(myoprnd))
    if X.is_symbol():
        return conjugate(X)
    
opconj = function("opconj",eval_func=opconjfunc)
##################################################################

##################################################################
# D derivative (Command: Dl)
##################################################################
def efDl(self, x):
    # x is a number:
    if x in RR or (I*x) in RR:
        return 0
    
    # Exponential exp(a)
    elif x.operator()==exp:
        return Dl((log(x)).simplify())*x
    
    # Logarithm log(a)
    elif x.operator()==ln:
        return (1/exp(x))*Dl(exp(x))
        
    # Single symbolic input
    elif x.is_symbol():
        if x.is_integer() or (I*x).is_integer() or x in RR or (I*x) in RR:
            return 0
        else:
            pass
    elif (-1*(x)).is_symbol():
        if x.is_integer() or (I*x).is_integer() or x in RR or (I*x) in RR:
            return 0
        else:
            return -Dl(-x)
        
    # Power in input with plus (a1^c1)
    elif x.operator() == operator.pow:
        try:
            a, b = x.operands()
            if a.is_integer() or a in RR or (I*a) in RR:
                return a*(b^(a-1))*Dl(b)
            else:
                return b*(a^(b-1))*Dl(a)
        except:
            a, b, c = x.operands()
            if a.is_integer() or (I*a).is_integer() or a in RR or (I*a) in RR:
                return a*(b^(a-1))*Dl(b)
            else:
                return b*(a^(b-1))*Dl(a)
    
    # Power in input with minus (-a1^c1)
    elif (-1*x).operator() == operator.pow:
        try:
            a, b = (-1*x).operands()
            if a.is_integer() or (I*a).is_integer() or a in RR or (I*a) in RR:
                return -a*(b^(a-1))*Dl(b)
            else:
                return -b*(a^(b-1))*Dl(a)
        except:
            a, b, c = (-1*x).operands()
            if a.is_integer() or (I*a).is_integer() or a in RR or (I*a) in RR:
                return -a*(b^(a-1))*Dl(b)
            else:
                return -b*(a^(b-1))*Dl(a)
            
    # Summation in input
    elif x.operator() == add_vararg:
        oplist = x.operands()
        return sum(Dl(op) for op in oplist)
    
    # Multiplication in input
    elif x.operator() == mul_vararg:
        ops=x.operands()
        result=0
        for i in range(len(ops)):
            prod=1
            for j in range(len(ops)):
                if i!=j:
                    prod=prod*ops[j]
            result=result+Dl(ops[i])*prod
        return result

##################################################################
# Delta derivative (Command: Dn)
##################################################################
def efDn(self, x):
    # x is a number:
    if x in RR or (I*x) in RR:
        return 0
    
    # Exponential exp(a)
    elif x.operator()==exp:
        return Dn((log(x)).simplify())*x
    
    # Logarithm log(a)
    elif x.operator()==ln:
        return (1/exp(x))*Dn(exp(x))
        
    # Single symbolic input
    elif x.is_symbol():
        if x.is_integer() or (I*x).is_integer() or x in RR or (I*x) in RR:
            return 0
        else:
            pass
    elif (-1*(x)).is_symbol():
        if x.is_integer() or (I*x).is_integer() or x in RR or (I*x) in RR:
            return 0
        else:
            return -Dn(-x)
        
    # Power in input with plus (a1^c1)
    elif x.operator() == operator.pow:
        try:
            a, b = x.operands()
            if a.is_integer() or a in RR or (I*a) in RR:
                return a*(b^(a-1))*Dn(b)
            else:
                return b*(a^(b-1))*Dn(a)
        except:
            a, b, c = x.operands()
            if a.is_integer() or (I*a).is_integer() or a in RR or (I*a) in RR:
                return a*(b^(a-1))*Dn(b)
            else:
                return b*(a^(b-1))*Dn(a)
    
    # Power in input with minus (-a1^c1)
    elif (-1*x).operator() == operator.pow:
        try:
            a, b = (-1*x).operands()
            if a.is_integer() or (I*a).is_integer() or a in RR or (I*a) in RR:
                return -a*(b^(a-1))*Dn(b)
            else:
                return -b*(a^(b-1))*Dn(a)
        except:
            a, b, c = (-1*x).operands()
            if a.is_integer() or (I*a).is_integer() or a in RR or (I*a) in RR:
                return -a*(b^(a-1))*Dn(b)
            else:
                return -b*(a^(b-1))*Dn(a)
            
    # Summation in input
    elif x.operator() == add_vararg:
        oplist = x.operands()
        return sum(Dn(op) for op in oplist)
    
    # Multiplication in input
    elif x.operator() == mul_vararg:
        ops=x.operands()
        result=0
        for i in range(len(ops)):
            proDn=1
            for j in range(len(ops)):
                if i!=j:
                    proDn=proDn*ops[j]
            result=result+Dn(ops[i])*proDn
        return result

##################################################################
# delta derivative (Command: Dm)
##################################################################
def efDm(self, x):
    # x is a number:
    if x in RR or (I*x) in RR:
        return 0
    
    # Exponential exp(a)
    elif x.operator()==exp:
        return Dm((log(x)).simplify())*x
    
    # Logarithm log(a)
    elif x.operator()==ln:
        return (1/exp(x))*Dm(exp(x))
        
    # Single symbolic input
    elif x.is_symbol():
        if x.is_integer() or (I*x).is_integer() or x in RR or (I*x) in RR:
            return 0
        else:
            pass
    elif (-1*(x)).is_symbol():
        if x.is_integer() or (I*x).is_integer() or x in RR or (I*x) in RR:
            return 0
        else:
            return -Dm(-x)
        
    # Power in input with plus (a1^c1)
    elif x.operator() == operator.pow:
        try:
            a, b = x.operands()
            if a.is_integer() or a in RR or (I*a) in RR:
                return a*(b^(a-1))*Dm(b)
            else:
                return b*(a^(b-1))*Dm(a)
        except:
            a, b, c = x.operands()
            if a.is_integer() or (I*a).is_integer() or a in RR or (I*a) in RR:
                return a*(b^(a-1))*Dm(b)
            else:
                return b*(a^(b-1))*Dm(a)
    
    # Power in input with minus (-a1^c1)
    elif (-1*x).operator() == operator.pow:
        try:
            a, b = (-1*x).operands()
            if a.is_integer() or (I*a).is_integer() or a in RR or (I*a) in RR:
                return -a*(b^(a-1))*Dm(b)
            else:
                return -b*(a^(b-1))*Dm(a)
        except:
            a, b, c = (-1*x).operands()
            if a.is_integer() or (I*a).is_integer() or a in RR or (I*a) in RR:
                return -a*(b^(a-1))*Dm(b)
            else:
                return -b*(a^(b-1))*Dm(a)
            
    # Summation in input
    elif x.operator() == add_vararg:
        oplist = x.operands()
        return sum(Dm(op) for op in oplist)
    
    # Multiplication in input
    elif x.operator() == mul_vararg:
        ops=x.operands()
        result=0
        for i in range(len(ops)):
            proDm=1
            for j in range(len(ops)):
                if i!=j:
                    proDm=proDm*ops[j]
            result=result+Dm(ops[i])*proDm
        return result

##################################################################
# deltabar derivative (Command: Dmbar)
##################################################################
def efDmbar(self, x):
    # x is a number:
    if x in RR or (I*x) in RR:
        return 0
    
    # Exponential exp(a)
    elif x.operator()==exp:
        return Dmbar((log(x)).simplify())*x
    
    # Logarithm log(a)
    elif x.operator()==ln:
        return (1/exp(x))*Dmbar(exp(x))
        
    # Single symbolic input
    elif x.is_symbol():
        if x.is_integer() or (I*x).is_integer() or x in RR or (I*x) in RR:
            return 0
        else:
            pass
    elif (-1*(x)).is_symbol():
        if x.is_integer() or (I*x).is_integer() or x in RR or (I*x) in RR:
            return 0
        else:
            return -Dmbar(-x)
        
    # Power in input with plus (a1^c1)
    elif x.operator() == operator.pow:
        try:
            a, b = x.operands()
            if a.is_integer() or a in RR or (I*a) in RR:
                return a*(b^(a-1))*Dmbar(b)
            else:
                return b*(a^(b-1))*Dmbar(a)
        except:
            a, b, c = x.operands()
            if a.is_integer() or (I*a).is_integer() or a in RR or (I*a) in RR:
                return a*(b^(a-1))*Dmbar(b)
            else:
                return b*(a^(b-1))*Dmbar(a)
    
    # Power in input with minus (-a1^c1)
    elif (-1*x).operator() == operator.pow:
        try:
            a, b = (-1*x).operands()
            if a.is_integer() or (I*a).is_integer() or a in RR or (I*a) in RR:
                return -a*(b^(a-1))*Dmbar(b)
            else:
                return -b*(a^(b-1))*Dmbar(a)
        except:
            a, b, c = (-1*x).operands()
            if a.is_integer() or (I*a).is_integer() or a in RR or (I*a) in RR:
                return -a*(b^(a-1))*Dmbar(b)
            else:
                return -b*(a^(b-1))*Dmbar(a)
            
    # Summation in input
    elif x.operator() == add_vararg:
        oplist = x.operands()
        return sum(Dmbar(op) for op in oplist)
    
    # Multiplication in input
    elif x.operator() == mul_vararg:
        ops=x.operands()
        result=0
        for i in range(len(ops)):
            proDmbar=1
            for j in range(len(ops)):
                if i!=j:
                    proDmbar=proDmbar*ops[j]
            result=result+Dmbar(ops[i])*proDmbar
        return result


##################################################################
# Derivative operators:
#######################
Dl = function("Dl",latex_name=r"D",eval_func=efDl)
Dn = function("Dn",latex_name=r"\Delta",eval_func=efDn)
Dm = function("Dm",latex_name=r"\delta",eval_func=efDm)
Dmbar = function("Dmbar",latex_name=r"\bar{\delta}",eval_func=efDmbar)
##################################################################

##################################################################
# Commutators 
# Page 77, Eq.(7.6)
var('alpha,beta,epsilon,gamma,kappa,mu,nu,pi,rho,sigma,tau',domain='complex')
var('lambdaa',latex_name=r"\lambda",domain='complex')
var('Phi01,Phi02,Phi10,Phi12,Phi20,Phi21,Lambda',domain='complex')
var('Phi00,Phi11,Phi22',domain='real')
var('Psi0,Psi1,Psi2,Psi3,Psi4',domain='complex')

# Commutations as =0.
def Dn_Dl_com(X):
    comresult=(gamma+conjugate(gamma))*Dl(X)+(epsilon+conjugate(epsilon))*Dn(X)-(tau+conjugate(pi))*Dmbar(X)-(conjugate(tau)+pi)*Dm(X)-Dn(Dl(X))+Dl(Dn(X))
    return comresult

def Dm_Dl_com(X):
    comresult=(conjugate(alpha)+beta-conjugate(pi))*Dl(X)+kappa*Dn(X)-sigma*Dmbar(X)-(conjugate(rho)+epsilon-conjugate(epsilon))*Dm(X)-Dm(Dl(X))+Dl(Dm(X))
    return comresult

def Dm_Dn_com(X):
    comresult=-conjugate(nu)*Dl(X)+(tau-conjugate(alpha)-beta)*Dn(X)+conjugate(lambdaa)*Dmbar(X)+(mu-gamma+conjugate(gamma))*Dm(X)-Dm(Dn(X))+Dn(Dm(X))
    return comresult

def Dmbar_Dm_com(X):
    comresult=(conjugate(mu)-mu)*Dl(X)+(conjugate(rho)-rho)*Dn(X)-(conjugate(alpha)-beta)*Dmbar(X)-(conjugate(beta)-alpha)*Dm(X)-Dmbar(Dm(X))+Dm(Dmbar(X))
    return comresult  
    
def Dmbar_Dl_com(X):
    comresult=(alpha+conjugate(beta)-pi)*Dl(X)+conjugate(kappa)*Dn(X)-conjugate(sigma)*Dm(X)-(rho+conjugate(epsilon)-epsilon)*Dmbar(X)-Dmbar(Dl(X))+Dl(Dmbar(X))
    return comresult

def Dmbar_Dn_com(X):
    comresult=-nu*Dl(X)+(conjugate(tau)-alpha-conjugate(beta))*Dn(X)+lambdaa*Dm(X)+(conjugate(mu)-conjugate(gamma)+gamma)*Dmbar(X)-Dmbar(Dn(X))+Dn(Dmbar(X))
    return comresult
    
##################################################################
##################################################################
def NPeqs():
    global NP1,NP2,NP3,NP4,NP5,NP6,NP7,NP8,NP9
    global NP10,NP11,NP12,NP13,NP14,NP15,NP16,NP17,NP18
    NP1=-Dl(rho)+Dmbar(kappa)+(rho^2+sigma*conjugate(sigma))+(epsilon+conjugate(epsilon))*rho-conjugate(kappa)*tau-kappa*(3*alpha+conjugate(beta)-pi)+Phi00
    NP2=-Dl(sigma)+Dm(kappa)+sigma*(3*epsilon-conjugate(epsilon)+rho+conjugate(rho))+kappa*(conjugate(pi)-tau-3*beta-conjugate(alpha))+Psi0
    NP3=-Dl(tau)+Dn(kappa)+(tau+conjugate(pi))*rho+(conjugate(tau)+pi)*sigma+(epsilon-conjugate(epsilon))*tau-(3*gamma+conjugate(gamma))*kappa+Psi1+Phi01
    NP4=-Dl(alpha)+Dmbar(epsilon)+(rho+conjugate(epsilon)-2*epsilon)*alpha+beta*conjugate(sigma)-conjugate(beta)*epsilon-kappa*lambdaa-conjugate(kappa)*gamma+(epsilon+rho)*pi+Phi10
    NP5=-Dl(beta) + Dm(epsilon) +(alpha + pi)*sigma + (conjugate(rho) - conjugate(epsilon))*beta - (mu + gamma)*kappa - (conjugate(alpha) - conjugate(pi))*epsilon + Psi1
    NP6=-Dl(gamma) + Dn(epsilon) + (tau + conjugate(pi))*alpha + (conjugate(tau)+pi)*beta - (epsilon + conjugate(epsilon))*gamma - (gamma + conjugate(gamma))*epsilon + tau*pi - nu*kappa + Psi2 + Phi11 - Lambda
    NP7=-Dl(lambdaa) + Dmbar(pi) + (rho*lambdaa + conjugate(sigma)*mu) + pi^2 + (alpha - conjugate(beta))*pi - nu*conjugate(kappa) - (3*epsilon - conjugate(epsilon))*lambdaa + Phi20
    NP8=-Dl(mu) + Dm(pi) + (conjugate(rho)*mu + sigma*lambdaa) +pi*conjugate(pi) - (epsilon + conjugate(epsilon))*mu - (conjugate(alpha) - beta)*pi - nu*kappa + Psi2 + 2*Lambda
    NP9=-Dl(nu) + Dn(pi) +(pi + conjugate(tau))*mu + (conjugate(pi) + tau)*lambdaa + (gamma - conjugate(gamma))*pi - (3*epsilon + conjugate(epsilon))*nu + Psi3 + Phi21
    NP10=-Dn(lambdaa) + Dmbar(nu) - (mu + conjugate(mu))*lambdaa - (3*gamma - conjugate(gamma))*lambdaa + (3*alpha + conjugate(beta) + pi - conjugate(tau))*nu - Psi4
    NP11=-Dm(rho) + Dmbar(sigma) + (conjugate(alpha) + beta)*rho - (3*alpha - conjugate(beta))*sigma + (rho - conjugate(rho))*tau + (mu - conjugate(mu))*kappa - Psi1 + Phi01
    NP12=-Dm(alpha) + Dmbar(beta) + (mu*rho - lambdaa*sigma) + alpha*conjugate(alpha) + beta*conjugate(beta) - 2*alpha*beta + (rho-conjugate(rho))*gamma + (mu-conjugate(mu))*epsilon - Psi2 + Phi11 + Lambda
    NP13=-Dm(lambdaa) + Dmbar(mu) + (rho - conjugate(rho))*nu + (mu - conjugate(mu))*pi + (alpha + conjugate(beta))*mu + (conjugate(alpha) - 3*beta)*lambdaa - Psi3 + Phi21
    NP14=-Dm(nu) + Dn(mu) +(mu^2 + lambdaa*conjugate(lambdaa)) + (gamma + conjugate(gamma))*mu - conjugate(nu)*pi + (tau - 3*beta - conjugate(alpha))*nu + Phi22
    NP15=-Dm(gamma) + Dn(beta) + (tau - conjugate(alpha) - beta)*gamma + mu*tau - sigma* nu - epsilon*conjugate(nu) - (gamma - conjugate(gamma) - mu)*beta + alpha*conjugate(lambdaa) + Phi12
    NP16=-Dm(tau) + Dn(sigma) + (mu*sigma + conjugate(lambdaa)*rho) + (tau + beta - conjugate(alpha))*tau - (3*gamma - conjugate(gamma))*sigma - kappa*conjugate(nu)  + Phi02
    NP17=-Dn(rho) + Dmbar(tau) - (rho*conjugate(mu) + sigma*lambdaa) + (conjugate(beta) - alpha - conjugate(tau))*tau + (gamma + conjugate(gamma))*rho + nu*kappa - Psi2 - 2*Lambda
    NP18=-Dn(alpha) + Dmbar(gamma) + (rho + epsilon)*nu - (tau + beta)*lambdaa + (conjugate(gamma)- conjugate(mu))*alpha + (conjugate(beta)- conjugate(tau))*gamma - Psi3
    show("NP1=0=",NP1)
    show("NP2=0=",NP2)
    show("NP3=0=",NP3)
    show("NP4=0=",NP4)
    show("NP5=0=",NP5)
    show("NP6=0=",NP6)
    show("NP7=0=",NP7)
    show("NP8=0=",NP8)
    show("NP9=0=",NP9)
    show("NP10=0=",NP10)
    show("NP11=0=",NP11)
    show("NP12=0=",NP12)
    show("NP13=0=",NP13)
    show("NP14=0=",NP14)
    show("NP15=0=",NP15)
    show("NP16=0=",NP16)
    show("NP17=0=",NP17)
    show("NP18=0=",NP18)
    print("NP equations are ready to use as NP1(=0), NP2(=0), ..., NP18(=0).")
##################################################################
##################################################################
def Bianchi():
    global BI1, BI2, BI3,BI4,BI5,BI6,BI7,BI8,BI9,BI10,BI11
    BI1=-(Dmbar(Psi0)-Dl(Psi1)+Dl(Phi01)-Dm(Phi00))+(4*alpha-pi)*Psi0-2*(2*rho+epsilon)*Psi1+3*kappa*Psi2+(conjugate(pi)-2*conjugate(alpha)-2*beta)*Phi00+2*(epsilon+conjugate(rho))*Phi01+2*sigma*Phi10-2*kappa*Phi11-conjugate(kappa)*Phi02
    BI2=-Dn(Psi0) + Dm(Psi1) - Dl(Phi02) + Dm(Phi01) + (4*gamma-mu)*Psi0 - 2*(2*tau + beta)*Psi1 + 3*sigma*Psi2 + (2*epsilon - 2*conjugate(epsilon) + conjugate(rho))*Phi02 + 2*(conjugate(pi) - beta)*Phi01 + 2*sigma*Phi11 - 2*kappa*Phi12 - conjugate(lambdaa)*Phi00
    BI3=-Dmbar(Psi3) + Dl(Psi4) - Dmbar(Phi21) + Dn(Phi20) + (4*epsilon - rho)*Psi4 - 2*(2*pi + alpha)*Psi3 + 3*lambdaa*Psi2 + (2*gamma - 2*conjugate(gamma) + conjugate(mu))*Phi20 + 2*(conjugate(tau) - alpha)*Phi21 + 2*lambdaa*Phi11 - 2*nu*Phi10 - conjugate(sigma)*Phi22
    BI4=-Dn(Psi3) + Dm(Psi4) - Dmbar(Phi22) + Dn(Phi21) + (4*beta - tau)*Psi4 - 2*(2*mu + gamma)*Psi3 + 3*nu*Psi2 + (conjugate(tau) - 2*conjugate(beta) - 2*alpha)*Phi22 + 2*(gamma+conjugate(mu))*Phi21 + 2*lambdaa*Phi12 - 2*nu*Phi11 - conjugate(nu)*Phi20
    BI5=-Dl(Psi2) + Dmbar(Psi1) - Dn(Phi00) + Dmbar(Phi01) - 2*Dl(Lambda) - lambdaa*Psi0 + 2*(pi-alpha)*Psi1+ 3*rho*Psi2 - 2*kappa*Psi3 + (2*gamma + 2*conjugate(gamma) - conjugate(mu))*Phi00 - 2*(conjugate(tau) + alpha)*Phi01 - 2*tau*Phi10 + 2*rho*Phi11 + conjugate(sigma)*Phi02
    BI6=-Dn(Psi2) + Dm(Psi3) - Dl(Phi22) + Dm(Phi21) - 2* Dn(Lambda) + sigma*Psi4 + 2*(beta-tau)*Psi3 - 3*mu*Psi2 + 2*nu*Psi1 + (conjugate(rho) - 2*epsilon - 2*conjugate(epsilon))*Phi22 + 2*(conjugate(pi) + beta)*Phi21 + 2*pi*Phi12 - 2*mu*Phi11 - conjugate(lambdaa)*Phi20
    BI7=-Dl(Psi3) + Dmbar(Psi2) + Dl(Phi21) -Dm(Phi20) + 2*Dmbar(Lambda) - kappa*Psi4 + 2*(rho-epsilon)*Psi3 + 3*pi*Psi2 - 2*lambdaa*Psi1 + (2*conjugate(alpha)-2*beta - conjugate(pi))*Phi20 - 2*(conjugate(rho)-epsilon)*Phi21 - 2*pi*Phi11 + 2*mu*Phi10 + conjugate(kappa)*Phi22
    BI8=-Dn(Psi1) + Dm(Psi2) + Dn(Phi01) - Dmbar(Phi02) + 2*Dm(Lambda) + nu*Psi0 + 2*(gamma-mu)*Psi1 - 3*tau*Psi2 + 2*sigma*Psi3 + (conjugate(tau) - 2*conjugate(beta) + 2*alpha)*Phi02 + 2*(conjugate(mu) - gamma)*Phi01 + 2*tau*Phi11-2*rho*Phi12-conjugate(nu)*Phi00
    BI9=-Dl(Phi11) + Dm(Phi10) + Dmbar(Phi01) - Dn(Phi00) - 3*Dl(Lambda) + (2*gamma - mu + 2*conjugate(gamma) - conjugate(mu))*Phi00 + (pi - 2*alpha - 2*conjugate(tau))*Phi01 + (conjugate(pi) - 2*conjugate(alpha) - 2*tau)*Phi10 + 2*(rho+conjugate(rho))*Phi11 + conjugate(sigma)*Phi02 + sigma*Phi20 - conjugate(kappa)*Phi12 - kappa*Phi21
    BI10=-Dl(Phi12) + Dm(Phi11) + Dmbar(Phi02) - Dn(Phi01) - 3*Dm(Lambda) + (-2*alpha + 2*conjugate(beta) + pi - conjugate(tau))*Phi02 + (conjugate(rho) + 2*rho - 2*conjugate(epsilon))*Phi12 + 2*(conjugate(pi) - tau)*Phi11 + (2*gamma - 2*conjugate(mu) - mu)*Phi01 + conjugate(nu)*Phi00 - conjugate(lambdaa)*Phi10 + sigma*Phi21 - kappa*Phi22
    BI11=-Dl(Phi22) + Dm(Phi21) + Dmbar(Phi12) - Dn(Phi11) - 3*Dn(Lambda) + (rho + conjugate(rho) - 2*epsilon - 2*conjugate(epsilon))*Phi22 + (2*conjugate(beta) + 2*pi - conjugate(tau))*Phi12 + (2*beta + 2*conjugate(pi) - tau)*Phi21 - 2*(mu+conjugate(mu))*Phi11 + nu*Phi01 + conjugate(nu)*Phi10 - conjugate(lambdaa)*Phi20 - lambdaa*Phi02
    show("BI1=0=",BI1)
    show("BI2=0=",BI2)
    show("BI3=0=",BI3)
    show("BI4=0=",BI4)
    show("BI5=0=",BI5)
    show("BI6=0=",BI6)
    show("BI7=0=",BI7)
    show("BI8=0=",BI8)
    show("BI9=0=",BI9)
    show("BI10=0=",BI10)
    show("BI11=0=",BI11)
    print("Bianchi identities are ready to use as BI1(=0), BI2(=0), ..., BI11(=0).")
##################################################################
def showall():
    global alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau, lambdaa, Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Lambda, Psi0, Psi1, Psi2, Psi3, Psi4
    show("alpha=",alpha)
    show("beta=",beta)
    show("epsilon=",epsilon)
    show("gamma=",gamma)
    show("kappa=",kappa)
    show("mu=",mu)
    show("nu=",nu)
    show("pi=",pi)
    show("rho=",rho)
    show("sigma=",sigma)
    show("tau=",tau)
    show("lambda=",lambdaa)
    show("Phi00=",Phi00)
    show("Phi01=",Phi01)
    show("Phi02=",Phi02)
    show("Phi10=",Phi10)
    show("Phi11=",Phi11)
    show("Phi12=",Phi12)
    show("Phi20=",Phi20)
    show("Phi21=",Phi21)
    show("Phi22=",Phi22)
    show("Lambda=",Lambda)
    show("Psi0=",Psi0)
    show("Psi1=",Psi1)
    show("Psi2=",Psi2)
    show("Psi3=",Psi3)
    show("Psi4=",Psi4)

##################################################################
#Carmeli and Kaye Article Transformations
##################################################################
##################################
def simplifyPhis(alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau, lambdaa, Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Lambda, Psi0, Psi1, Psi2, Psi3, Psi4):

    s01, s02, s10, s12, s20, s21 = SR.var('Phi01, Phi02, Phi10, Phi12, Phi20, Phi21')

    subs_dict = {
        conjugate(s01):s10,
        conjugate(s10):s01,
        conjugate(s02):s20,
        conjugate(s12):s21,
        conjugate(s20):s02,
        conjugate(s21):s12        
    }

    input_list = [alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau, lambdaa, Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Lambda, Psi0, Psi1, Psi2, Psi3, Psi4]
    output_list = []

    for expr in input_list:
        temp = expr.subs(subs_dict)
        final_val = temp.canonicalize_radical()
        
        output_list.append(final_val)

    return output_list
##################################
def CKAdo(z):
    global alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau, lambdaa, Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Lambda, Psi0, Psi1, Psi2, Psi3, Psi4

    var('newalpha, newbeta, newepsilon, newgamma, newkappa, newmu, newnu, newpi, newrho, newsigma, newtau, newlambdaa, newPhi01, newPhi02, newPhi10, newPhi12, newPhi20, newPhi21, newPsi0, newPsi1, newPsi2, newPsi3, newPsi4',domain='complex')

    var('newPhi00, newPhi11, newPhi22',domain='real')

    print("Null Rotation about l")
    #A. Null Rotation about l

    newrho=rho+z*kappa
    newalpha=alpha+z*(rho+epsilon)+z^2*kappa
    newlambdaa=lambdaa+z*(pi+2*alpha)+z^2*(rho+2*epsilon)+z^3*kappa+Dmbar(z)+z*Dl(z)
    newkappa=kappa
    newepsilon=epsilon+z*kappa
    newpi=pi+2*z*epsilon+z^2*kappa+Dl(z)
    newsigma=sigma+conjugate(z)*kappa
    newbeta=beta+z*sigma+conjugate(z)*epsilon+z*conjugate(z)*kappa
    newmu=mu+2*z*beta+conjugate(z)*pi+z^2*sigma+2*z*conjugate(z)*epsilon+z^2*conjugate(z)*kappa+Dm(z)+conjugate(z)*Dl(z)
    newtau=tau+z*sigma+conjugate(z)*rho+z*conjugate(z)*kappa
    newgamma=gamma+z*(tau+beta)+conjugate(z)*alpha+z^2*sigma+z*conjugate(z)*(rho+epsilon)+z^2*conjugate(z)*kappa
    newnu=nu+z*(mu+2*gamma)+conjugate(z)*lambdaa+z^2*(2*beta+tau)+z*conjugate(z)*(2*alpha+pi)+z^2*conjugate(z)*(2*epsilon+rho)+z^3*sigma+z^3*conjugate(z)*kappa+Dn(z)+z*Dm(z)+conjugate(z)*Dmbar(z)+z*conjugate(z)*Dl(z)

    newPsi0=Psi0
    newPsi1=z*Psi0+Psi1
    newPsi2=z^2*Psi0+2*z*Psi1+Psi2
    newPsi3=z^3*Psi0+3*z^2*Psi1+3*z*Psi2+Psi3
    newPsi4=z^4*Psi0+4*z^3*Psi1+6*z^2*Psi2+4*z*Psi3+Psi4

    newPhi00=(Phi00)
    newPhi01=(conjugate(z)*Phi00+Phi01)
    newPhi02=(conjugate(z)^2)*Phi00+2*conjugate(z)*Phi01+Phi02
    newPhi11=(z*conjugate(z)*Phi00+z*Phi01+conjugate(z)*Phi10+Phi11)
    newPhi12=(z*conjugate(z)^2*Phi00+2*conjugate(z)*z*Phi01+z*Phi02+2*conjugate(z)*Phi11+conjugate(z)^2*Phi10+Phi12)
    newPhi22=(z^2*conjugate(z)^2*Phi00+2*conjugate(z)*z^2*Phi01+z^2*Phi02+2*conjugate(z)^2*z*Phi10+4*conjugate(z)*z*Phi11+2*z*Phi12+conjugate(z)^2*Phi20+2*conjugate(z)*Phi21+Phi22)

    alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau, lambdaa, Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Psi0, Psi1, Psi2, Psi3, Psi4=newalpha, newbeta, newepsilon, newgamma, newkappa, newmu, newnu, newpi, newrho, newsigma, newtau, newlambdaa, newPhi00, newPhi01, newPhi02, newPhi10, newPhi11, newPhi12, newPhi20, newPhi21, newPhi22, newPsi0, newPsi1, newPsi2, newPsi3, newPsi4

    Phi10=opconj(Phi01)
    Phi20=opconj(Phi02)
    Phi21=opconj(Phi12)

    (alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau, lambdaa, Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Lambda, Psi0, Psi1, Psi2, Psi3, Psi4)=simplifyPhis(alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau, lambdaa, Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Lambda, Psi0, Psi1, Psi2, Psi3, Psi4)
    

#####################################
def CKAsee(z):
    print("Null Rotation about l")
    #A. Null Rotation about l

    var('l,l1,n,n1,m,m1')
    var('l1',latex_name=r"l'")
    var('n1',latex_name=r"n'")
    var('m1',latex_name=r"m'")

    #Null Tetrad Transforms

    #Here comp show that component
    #Here l1,n1,m1 equal to l',n',m' respectively
    show("Null Tetrad Transforms")
    show("comp1=0=",l-l1)
    show("comp2=0=",conjugate(z)*l+m-m1)
    show("comp3=0=",z*conjugate(z)*l+z*m+conjugate(z)*conjugate(m)+n-n1)

    show("########################################################")
    #Spin Coefficients Transform

    var('lambdaa1',latex_name=r"\lambda'")
    var('rho1,alpha1,kappa1,epsilon1,pi1,sigma1,beta1,mu1,tau1,nu1,gamma1')
    var('alpha1',latex_name=r"\alpha'")
    var('beta1',latex_name=r"\beta'")
    var('kappa1',latex_name=r"\kappa'")
    var('epsilon1',latex_name=r"\epsilon'")
    var('pi1',latex_name=r"\pi'")
    var('sigma1',latex_name=r"\sigma'")
    var('mu1',latex_name=r"\mu'")
    var('tau1',latex_name=r"\tau'")
    var('nu1',latex_name=r"\nu'")
    var('gamma1',latex_name=r"\gamma'")
    var('rho1',latex_name=r"\rho'")

    show("Spin Coefficients Transform")

    show("rhotr=0=",rho+z*kappa-rho1)
    show("alphatr=0=",alpha+z*(rho+epsilon)+z^2*kappa-alpha1)
    show("lambdatr=0=",lambdaa+z*(pi+2*alpha)+z^2*(rho+2*epsilon)+z^3*kappa+Dmbar(z)+z*Dl(z)-lambdaa1)
    show("kappatr=0=",kappa-kappa1)
    show("epsilontr=0=",epsilon+z*kappa-epsilon1)
    show("pitr=0=",pi+2*z*epsilon+z^2*kappa+Dl(z)-pi1)
    show("sigmatr=0=",sigma+conjugate(z)*kappa-sigma1)
    show("betatr=0=",beta+z*sigma+conjugate(z)*epsilon+z*conjugate(z)*kappa-beta1)
    show("mutr=0=",mu+2*z*beta+conjugate(z)*pi+z^2*sigma+2*z*conjugate(z)*epsilon+z^2*conjugate(z)*kappa+Dm(z)+conjugate(z)*Dl(z)-mu1)
    show("tautr=0=",tau+z*sigma+conjugate(z)*rho+z*conjugate(z)*kappa-tau1)
    show("gammatr=0=",gamma+z*(tau+beta)+conjugate(z)*alpha+z^2*sigma+z*conjugate(z)*(rho+epsilon)+z^2*conjugate(z)*kappa-gamma1)
    show("nutr=0=",nu+z*(mu+2*gamma)+conjugate(z)*lambdaa+z^2*(2*beta+tau)+z*conjugate(z)*(2*alpha+pi)+z^2*conjugate(z)*(2*epsilon+rho)+z^3*sigma+z^3*conjugate(z)*kappa+Dn(z)+z*Dm(z)+conjugate(z)*Dmbar(z)+z*conjugate(z)*Dl(z)-nu1)

    show("##################################################################")
    #Weyl Tensor Transform
    var('Psi00,Psi11,Psi22,Psi3,Psi44')
    var('Psi00',latex_name=r"\Psi_0'")
    var('Psi11',latex_name=r"\Psi_1'")
    var('Psi22',latex_name=r"\Psi_2'")
    var('Psi33',latex_name=r"\Psi_3'")
    var('Psi44',latex_name=r"\Psi_4'")

    show("Weyl Tensor Transform")

    show("Psi0tr=0=",Psi0-Psi00)
    show("Psi1tr=0=",z*Psi0+Psi1-Psi11)
    show("Psi2tr=0=",z^2*Psi0+2*z*Psi1+Psi2-Psi22)
    show("Psi3tr=0=",z^3*Psi0+3*z^2*Psi1+3*z*Psi2+Psi3-Psi33)
    show("Psi4tr=0=",z^4*Psi0+4*z^3*Psi1+6*z^2*Psi2+4*z*Psi3+Psi4-Psi44)

    show("####################################################################")

    #Ricci Tensor Transform
    var('Phi00x,Phi01x,Phi02x,Phi11x,Phi12x,Phi22x')
    var('Phi00x',latex_name=r"\Phi_{00}'")
    var('Phi01x',latex_name=r"\Phi_{01}'")
    var('Phi02x',latex_name=r"\Phi_{02}'")
    var('Phi11x',latex_name=r"\Phi_{11}'")
    var('Phi12x',latex_name=r"\Phi_{12}'")
    var('Phi22x',latex_name=r"\Phi_{22}'")
    show("Ricci Tensor Transform")

    show("Phi00tr=0=",Phi00-Phi00x)
    show("Phi01tr=0=",conjugate(z)*Phi00+Phi01-Phi01x)
    show("Phi02tr=0=",conjugate(z)^2*Phi00+2*conjugate(z)*Phi01+Phi02-Phi02x)
    show("Phi11tr=0=",z*conjugate(z)*Phi00+z*Phi01+conjugate(z)*Phi10+Phi11-Phi11x)
    show("Phi12tr=0=",z*conjugate(z)^2*Phi00+2*conjugate(z)*z*Phi01+z*Phi02+2*conjugate(z)*Phi11+conjugate(z)^2*Phi10+Phi12-Phi12x)
    show("Phi22tr=0=",z^2*conjugate(z)^2*Phi00+2*conjugate(z)*z^2*Phi01+z^2*Phi02+2*conjugate(z)^2*z*Phi10+4*conjugate(z)*z*Phi11+2*z*Phi12+conjugate(z)^2*Phi20+2*conjugate(z)*Phi21+Phi22-Phi22x)
##################################################################
def CKBdo(z):
    global alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau, lambdaa, Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Lambda, Psi0, Psi1, Psi2, Psi3, Psi4

    var('newalpha, newbeta, newepsilon, newgamma, newkappa, newmu, newnu, newpi, newrho, newsigma, newtau, newlambdaa, newPhi01, newPhi02, newPhi10, newPhi12, newPhi20, newPhi21, newPsi0, newPsi1, newPsi2, newPsi3, newPsi4',domain='complex')

    var('newPhi00, newPhi11, newPhi22',domain='real')

    print("Boost in l-n Plane and Spatial Rotation in m-mbar Plane")
    #B. Boost in l-n Plane and Spatial Rotation in m-mbar Plane

    newrho=z*conjugate(z)*rho
    newalpha=z^(-1)*conjugate(z)*(alpha-z*Dmbar(z^(-1)))
    newlambdaa=z^(-3)*conjugate(z)*lambdaa
    newkappa=z^(3)*conjugate(z)*kappa
    newepsilon=z*conjugate(z)*(epsilon-z*Dl(z^(-1)))
    newpi=z^(-1)*conjugate(z)*pi
    newsigma=z^(3)*conjugate(z)^(-1)*sigma
    newbeta=z*conjugate(z)^(-1)*(beta-z*Dm(z^(-1)))
    newmu=z^(-1)*conjugate(z)^(-1)*mu
    newtau=z*conjugate(z)^(-1)*tau
    newgamma=z^(-1)*conjugate(z)^(-1)*(gamma-z*Dn(z^(-1)))
    newnu=z^(-3)*conjugate(z)^(-1)*nu

    newPsi0=z^(4)*Psi0
    newPsi1=z^(2)*Psi1
    newPsi2=Psi2
    newPsi3=z^(-2)*Psi3
    newPsi4=z^(-4)*Psi4

    newPhi00=(z^(2)*conjugate(z)^2*Phi00)
    newPhi01=(z^(2)*Phi01)
    newPhi02=(z^(2)*conjugate(z)^(-2)*Phi02)
    newPhi11=(Phi11)
    newPhi12=(conjugate(z)^(-2)*Phi12)
    newPhi22=(z^(-2)*conjugate(z)^(-2)*Phi22)

    alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau, lambdaa, Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Psi0, Psi1, Psi2, Psi3, Psi4=newalpha, newbeta, newepsilon, newgamma, newkappa, newmu, newnu, newpi, newrho, newsigma, newtau, newlambdaa, newPhi00, newPhi01, newPhi02, newPhi10, newPhi11, newPhi12, newPhi20, newPhi21, newPhi22, newPsi0, newPsi1, newPsi2, newPsi3, newPsi4

    Phi10=opconj(Phi01)
    Phi20=opconj(Phi02)
    Phi21=opconj(Phi12)

    (alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau, lambdaa, Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Lambda, Psi0, Psi1, Psi2, Psi3, Psi4)=simplifyPhis(alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau, lambdaa, Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Lambda, Psi0, Psi1, Psi2, Psi3, Psi4)

#####################################
def CKBsee(z):
    print("Boost in l-n Plane and Spatial Rotation in m-mbar Plane")
    #B. Boost in l-n Plane and Spatial Rotation in m-mbar Plane

    var('l,l1,n,n1,m,m1')
    var('l1',latex_name=r"l'")
    var('n1',latex_name=r"n'")
    var('m1',latex_name=r"m'")

    #Null Tetrad Transforms

    #Here comp show that component
    #Here l1,n1,m1 equal to l',n',m' respectively
    show("Null Tetrad Transforms")
    show("comp1=0=",z*conjugate(z)*l-l1)
    show("comp2=0=",z/conjugate(z)*m-m1)
    show("comp3=0=",1/z*1/conjugate(z)*n-n1)

    show("########################################################")
    #Spin Coefficients Transform

    var('lambdaa1',latex_name=r"\lambda'")
    var('rho1,alpha1,kappa1,epsilon1,pi1,sigma1,beta1,mu1,tau1,nu1,gamma1')
    var('alpha1',latex_name=r"\alpha'")
    var('beta1',latex_name=r"\beta'")
    var('kappa1',latex_name=r"\kappa'")
    var('epsilon1',latex_name=r"\epsilon'")
    var('pi1',latex_name=r"\pi'")
    var('sigma1',latex_name=r"\sigma'")
    var('mu1',latex_name=r"\mu'")
    var('tau1',latex_name=r"\tau'")
    var('nu1',latex_name=r"\nu'")
    var('gamma1',latex_name=r"\gamma'")
    var('rho1',latex_name=r"\rho'")

    show("Spin Coefficients Transform")

    show("rhotr=0=",z*conjugate(z)*rho-rho1)
    show("alphatr=0=",z^(-1)*conjugate(z)*(alpha-z*Dmbar(z^(-1)))-alpha1)
    show("lambdatr=0=",z^(-3)*conjugate(z)*lambdaa-lambdaa1)
    show("kappatr=0=",z^(3)*conjugate(z)*kappa-kappa1)
    show("epsilontr=0=",z*conjugate(z)*(epsilon-z*Dl(z^(-1)))-epsilon1)
    show("pitr=0=",z^(-1)*conjugate(z)*pi-pi1)
    show("sigmatr=0=",z^(3)*conjugate(z)^(-1)*sigma-sigma1)
    show("betatr=0=",z*conjugate(z)^(-1)*(beta-z*Dm(z^(-1)))-beta1)
    show("mutr=0=",z^(-1)*conjugate(z)^(-1)*mu-mu1)
    show("tautr=0=",z*conjugate(z)^(-1)*tau-tau1)
    show("gammatr=0=",z^(-1)*conjugate(z)^(-1)*(gamma-z*Dn(z^(-1)))-gamma1)
    show("nutr=0=",z^(-3)*conjugate(z)^(-1)*nu-nu1)

    show("##################################################################")
    #Weyl Tensor Transform
    var('Psi00,Psi11,Psi22,Psi3,Psi44')
    var('Psi00',latex_name=r"\Psi_0'")
    var('Psi11',latex_name=r"\Psi_1'")
    var('Psi22',latex_name=r"\Psi_2'")
    var('Psi33',latex_name=r"\Psi_3'")
    var('Psi44',latex_name=r"\Psi_4'")

    show("Weyl Tensor Transform")

    show("Psi0tr=0=",z^(4)*Psi0-Psi00)
    show("Psi1tr=0=",z^(2)*Psi1-Psi11)
    show("Psi2tr=0=",Psi2-Psi22)
    show("Psi3tr=0=",z^(-2)*Psi3-Psi33)
    show("Psi4tr=0=",z^(-4)*Psi4-Psi44)

    show("####################################################################")
    #Ricci Tensor Transform
    var('Phi00x,Phi01x,Phi02x,Phi11x,Phi12x,Phi22x')
    var('Phi00x',latex_name=r"\Phi_{00}'")
    var('Phi01x',latex_name=r"\Phi_{01}'")
    var('Phi02x',latex_name=r"\Phi_{02}'")
    var('Phi11x',latex_name=r"\Phi_{11}'")
    var('Phi12x',latex_name=r"\Phi_{12}'")
    var('Phi22x',latex_name=r"\Phi_{22}'")
    show("Ricci Tensor Transform")

    show("Phi00tr=0=",z^(2)*conjugate(z)^2*Phi00-Phi00x)
    show("Phi01tr=0=",z^(2)*Phi01-Phi01x)
    show("Phi02tr=0=",z^(2)*conjugate(z)^(-2)*Phi02-Phi02x)
    show("Phi11tr=0=",Phi11-Phi11x)
    show("Phi12tr=0=",conjugate(z)^(-2)*Phi12-Phi12x)
    show("Phi22tr=0=",z^(-2)*conjugate(z)^(-2)*Phi22-Phi22x)
##################################################################
def CKCdo(z):
    global alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau, lambdaa, Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Lambda, Psi0, Psi1, Psi2, Psi3, Psi4

    var('newalpha, newbeta, newepsilon, newgamma, newkappa, newmu, newnu, newpi, newrho, newsigma, newtau, newlambdaa, newPhi01, newPhi02, newPhi10, newPhi12, newPhi20, newPhi21, newPsi0, newPsi1, newPsi2, newPsi3, newPsi4',domain='complex')

    var('newPhi00, newPhi11, newPhi22',domain='real')
    
    print("Null Rotation about n")
    #C. Null Rotation about n

    newrho=rho + 2*z*alpha + conjugate(z)*tau + z^2*lambdaa + 2*z*conjugate(z)*gamma + z^2*conjugate(z)*nu - Dmbar(z) - conjugate(z)*Dn(z)
    newalpha=alpha + z*lambdaa + conjugate(z)*gamma + z*conjugate(z)*nu
    newlambdaa=lambdaa + conjugate(z)*nu
    newkappa=kappa + z*(rho + 2*epsilon) + conjugate(z)*sigma + z^2*(2*alpha + pi) + z*conjugate(z)*(2*beta + tau) + z^2*conjugate(z)*(2*gamma + mu) + z^3*lambdaa + z^3*conjugate(z)*nu - Dl(z)- z*Dmbar(z)- conjugate(z)*Dm(z)- z*conjugate(z)*Dn(z)
    newepsilon=epsilon + z*(pi + alpha) + conjugate(z)*beta + z^2*lambdaa + z*conjugate(z)*(mu + gamma) + z^2*conjugate(z)*nu
    newpi=pi + 2*lambdaa + conjugate(z)*mu + z*conjugate(z)*nu
    newsigma=sigma + z*(tau + 2*beta) + z^2*(mu + 2*gamma)+ z^3*nu - Dm(z) - z*Dn(z)
    newbeta=beta + z*(mu + gamma) + z^2*nu
    newmu=mu + z*nu
    newtau=tau + 2*gamma*z + z^2*nu - Dn(z)
    newgamma=gamma + z*nu
    newnu=nu

    newPsi0=Psi0 + 4*z*Psi1 + 6*z^2*Psi2 + 4*z^3*Psi3 + z^4*Psi4
    newPsi1=Psi1 + 3*z*Psi2 + 3*z^2*Psi3 + z^3*Psi4
    newPsi2=Psi2 + 2*z*Psi3 + z^2*Psi4
    newPsi3=Psi3 + z*Psi4
    newPsi4=Psi4

    newPhi00=(Phi00 + 2*conjugate(z)*Phi01 + 2*z*Phi10 + 4*z*conjugate(z)*Phi11 + conjugate(z^2)*Phi02 + z^2*Phi20 + 2*z^2*conjugate(z)*Phi21 + 2*conjugate(z^2)*Phi12 + conjugate(z^2)*z^2*Phi22)
    newPhi10=(Phi10 + 2*conjugate(z)*Phi11 + z*Phi20 + 2*z*conjugate(z)*Phi21 + conjugate(z^2)*Phi12 + conjugate(z^2)*z*Phi22)
    newPhi11=(Phi11 + conjugate(z)*Phi12 + z*Phi21 + conjugate(z)*z*Phi22)
    newPhi20=(Phi20 + 2*conjugate(z)*Phi21 + conjugate(z^2)*Phi22)
    newPhi21=(Phi21 + conjugate(z)*Phi22)
    newPhi22=(Phi22)

    alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau, lambdaa, Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Psi0, Psi1, Psi2, Psi3, Psi4=newalpha, newbeta, newepsilon, newgamma, newkappa, newmu, newnu, newpi, newrho, newsigma, newtau, newlambdaa, newPhi00, newPhi01, newPhi02, newPhi10, newPhi11, newPhi12, newPhi20, newPhi21, newPhi22, newPsi0, newPsi1, newPsi2, newPsi3, newPsi4

    Phi01=opconj(Phi10)
    Phi02=opconj(Phi20)
    Phi12=opconj(Phi21)

    (alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau, lambdaa, Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Lambda, Psi0, Psi1, Psi2, Psi3, Psi4)=simplifyPhis(alpha, beta, epsilon, gamma, kappa, mu, nu, pi, rho, sigma, tau, lambdaa, Phi00, Phi01, Phi02, Phi10, Phi11, Phi12, Phi20, Phi21, Phi22, Lambda, Psi0, Psi1, Psi2, Psi3, Psi4)

#####################################
def CKCsee(z):
    print("Null Rotation about n")
    #C. Null Rotation about n

    var('l,l1,n,n1,m,m1')
    var('l1',latex_name=r"l'")
    var('n1',latex_name=r"n'")
    var('m1',latex_name=r"m'")

    #Null Tetrad Transforms

    #Here comp show that component
    #Here l1,n1,m1 equal to l',n',m' respectively
    show("Null Tetrad Transforms")
    show("comp1=0=", l + conjugate(z)*m + z*conjugate(m) + z*conjugate(z)*n - l1)
    show("comp2=0=", m + z*n - m1)
    show("comp3=0=", n- n1)

    show("########################################################")
    #Spin Coefficients Transform

    var('lambdaa1',latex_name=r"\lambda'")
    var('rho1,alpha1,kappa1,epsilon1,pi1,sigma1,beta1,mu1,tau1,nu1,gamma1')
    var('alpha1',latex_name=r"\alpha'")
    var('beta1',latex_name=r"\beta'")
    var('kappa1',latex_name=r"\kappa'")
    var('epsilon1',latex_name=r"\epsilon'")
    var('pi1',latex_name=r"\pi'")
    var('sigma1',latex_name=r"\sigma'")
    var('mu1',latex_name=r"\mu'")
    var('tau1',latex_name=r"\tau'")
    var('nu1',latex_name=r"\nu'")
    var('gamma1',latex_name=r"\gamma'")
    var('rho1',latex_name=r"\rho'")

    show("Spin Coefficients Transform")

    show("rhotr=0=", rho + 2*z*alpha + conjugate(z)*tau + z^2*lambdaa + 2*z*conjugate(z)*gamma + z^2*conjugate(z)*nu - Dmbar(z) - conjugate(z)*Dn(z) - rho1)
    show("alphatr=0=",alpha + z*lambdaa + conjugate(z)*gamma + z*conjugate(z)*nu - alpha1)
    show("lambdatr=0=", lambdaa + conjugate(z)*nu - lambdaa1)
    show("kappatr=0=", kappa + z*(rho + 2*epsilon) + conjugate(z)*sigma + z^2*(2*alpha + pi) + z*conjugate(z)*(2*beta + tau) + z^2*conjugate(z)*(2*gamma + mu) + z^3*lambdaa + z^3*conjugate(z)*nu - Dl(z)- z*Dmbar(z)- conjugate(z)*Dm(z)- z*conjugate(z)*Dn(z) -kappa1 )
    show("epsilontr=0=", epsilon + z*(pi + alpha) + conjugate(z)*beta + z^2*lambdaa + z*conjugate(z)*(mu + gamma) + z^2*conjugate(z)*nu - epsilon1)
    show("pitr=0=", pi + 2*lambdaa + conjugate(z)*mu + z*conjugate(z)*nu - pi1)
    show("sigmatr=0=", sigma + z*(tau + 2*beta) + z^2*(mu + 2*gamma)+ z^3*nu - Dm(z) - z*Dn(z) - sigma1)
    show("betatr=0=", beta + z*(mu + gamma) + z^2*nu - beta1)
    show("mutr=0=", mu + z*nu - mu1)
    show("tautr=0=", tau + 2*gamma*z + z^2*nu - Dn(z) - tau1)
    show("gammatr=0=", gamma + z*nu - gamma1)
    show("nutr=0=", nu - nu1)


    show("##################################################################")
    #Weyl Tensor Transform
    var('Psi00,Psi11,Psi22,Psi3,Psi44')
    var('Psi00',latex_name=r"\Psi_0'")
    var('Psi11',latex_name=r"\Psi_1'")
    var('Psi22',latex_name=r"\Psi_2'")
    var('Psi33',latex_name=r"\Psi_3'")
    var('Psi44',latex_name=r"\Psi_4'")


    show("Weyl Tensor Transform")

    show("Psi0tr=0=", Psi0 + 4*z*Psi1 + 6*z^2*Psi2 + 4*z^3*Psi3 + z^4*Psi4 - Psi00)
    show("Psi1tr=0=", Psi1 + 3*z*Psi2 + 3*z^2*Psi3 + z^3*Psi4 - Psi11)
    show("Psi2tr=0=", Psi2 + 2*z*Psi3 + z^2*Psi4 - Psi22)
    show("Psi3tr=0=", Psi3 + z*Psi4 - Psi33)
    show("Psi4tr=0=", Psi4 -Psi44)


    show("####################################################################")

    #Ricci Tensor Transform
    var('Phi00x,Phi01x,Phi02x,Phi11x,Phi12x,Phi22x')
    var('Phi00x',latex_name=r"\Phi_{00}'")
    var('Phi10x',latex_name=r"\Phi_{10}'")
    var('Phi21x',latex_name=r"\Phi_{21}'")
    var('Phi11x',latex_name=r"\Phi_{11}'")
    var('Phi20x',latex_name=r"\Phi_{20}'")
    var('Phi22x',latex_name=r"\Phi_{22}'")
    show("Ricci Tensor Transform")

    show("Phi00tr=0=", Phi00 + 2*conjugate(z)*Phi01 + 2*z*Phi10 + 4*z*conjugate(z)*Phi11 + conjugate(z^2)*Phi02 + z^2*Phi20 + 2*z^2*conjugate(z)*Phi21 + 2*conjugate(z^2)*Phi12 + conjugate(z^2)*z^2*Phi22 - Phi00x)
    show("Phi10tr=0=", Phi10 + 2*conjugate(z)*Phi11 + z*Phi20 + 2*z*conjugate(z)*Phi21 + conjugate(z^2)*Phi12 + conjugate(z^2)*z*Phi22 - Phi10x )
    show("Phi11tr=0=", Phi11 + conjugate(z)*Phi12 + z*Phi21 + conjugate(z)*z*Phi22 - Phi11x)
    show("Phi20tr=0=", Phi20 + 2*conjugate(z)*Phi21 + conjugate(z^2)*Phi22 - Phi20x)
    show("Phi21tr=0=", Phi21 + conjugate(z)*Phi22 - Phi21x)
    show("Phi22tr=0=", Phi22 -Phi22x)
##################################################################

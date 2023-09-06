#!/usr/bin/python

import math

class Jaffe_Diode:
  """Exact solution of the cold diode problem.

  """

  def __init__(self, d, J, W, V1):
    """Precompute as much of the solution as possible.

    Then hopefully functions which return field variables will be fast.
    """

    PI   = 3.141592654
    c    = 299792458.0    # (m/s)
    mu0  = 4.0*PI*1.0E-7  # Permeability =      N/A^2 (Newton/Ampere^2)
    eps0 = 1.0/(c*c*mu0)  # Permittivity =      (A^2.s^4/kg.m^3) = C^2/(N.m^2) = (s^2.C^2)/m^3.kg)
    qe   = 1.60217662e-19 # elementary charge = (e = 1.60217662e-19) C
    me   = 9.10938356e-31 # electron mass     = (m_e = 9.10938356e-31) kG

    v0    = math.sqrt( 2.0*W/me ) 
    J0    = (16.0/9.0)*eps0*math.sqrt( 2.0*qe/me )*math.pow(W/qe,1.5)/(d*d)

    # Non-dimensionalize
    kappa = (16.0/9.0)*J/J0
    xi0   = 2.0*math.sqrt(kappa)
    eta0  = qe*V1/W
    u0    = math.sqrt(1.0 + eta0)

    alpha1       = 0.0
    alpha2       = 0.0
    delta        = 0.0
    etamin_exact = 0.0
    umin_exact   = 0.0

    case_type = 0
    if (V1 > 0.0):
      case_type = 1
      alpha1       = self.__compute_alpha1_case1(xi0,u0)
      alpha2       = self.__compute_alpha2(alpha1)
    elif (V1 == 0.0):
      case_type = 2
      alpha1       = self.__compute_alpha1_case2(xi0,u0)
      delta        = self.__compute_delta(alpha1)
      etamin_exact = -(1.0 - alpha1*alpha1)
      umin_exact   = math.sqrt(1.0 + etamin_exact)
    else:
      print("Case not supported!")
      exit()

#    print ("Jaffe_Diode: ", xi0, u0, alpha1, delta, etamin_exact, v0)
#    print ("J,J0,kappa: ", J, J0, kappa)
#    print (" ")

    self.__setInput(case_type, d, W, J, V1, kappa, v0, u0, xi0, alpha1, alpha2, delta, etamin_exact, umin_exact)

  def __setInput(self, case_type, d, W, J, V1, kappa, v0, u0, xi0, alpha1, alpha2, delta, etamin_exact, umin_exact):

    """Make all the user inputs object data."""

    self.case_type    = case_type    # (1-Branch 1: V1 > 0, 2-Branch 2: V1=0)
    self.d            = d     # (m)
    self.W            = W     # (eV)
    self.J            = J     # (Amps/m^2/s)
    self.V1           = V1    # (Volts)
    self.kappa        = kappa
    self.v0           = v0    # (m/s)
    self.u0           = u0
    self.xi0          = xi0
    self.alpha1       = alpha1
    self.alpha2       = alpha2
    self.delta        = delta
    self.etamin_exact = etamin_exact
    self.umin_exact   = umin_exact

    self.PI   = 3.141592654
    self.c    = 299792458.0                  # (m/s)
    self.mu0  = 4.0*self.PI*1.0E-7           # Permeability =      N/A^2 (Newton/Ampere^2)
    self.eps0 = 1.0/(self.c*self.c*self.mu0) # Permittivity =      (A^2.s^4/kg.m^3) = C^2/(N.m^2) = (s^2.C^2)/m^3.kg)
    self.qe   = 1.60217662e-19               # elementary charge = {e = 1.60217662e-19} C
    self.me   = 9.10938356e-31               # electron mass     = {m_e = 9.10938356e-31} kG

    return

# ------------------------------------------------------------------------------

  def compute_phi(self, x):
    xi  = 2.0*math.sqrt(self.kappa)*x/self.d
    u   = self.__compute_u(xi)
    eta = u*u - 1.0;
    phi = eta*self.W/self.qe
    return phi

  def compute_phimin(self):
    phimin_exact = self.etamin_exact*self.W/self.qe
    return phimin_exact

  def compute_eta(self, x):
    xi  = 2.0*math.sqrt(self.kappa)*x/self.d
    u   = self.__compute_u(xi)
    eta = u*u - 1.0;
    return eta

  def compute_velocity(self, x):
    xi  = 2.0*math.sqrt(self.kappa)*x/self.d
    u   = self.__compute_u(xi)
    eta = u*u - 1.0
    v   = self.v0*math.sqrt( 1.0 + eta )
    return v

  def compute_ne(self, x):
    xi  = 2.0*math.sqrt(self.kappa)*x/self.d
    u   = self.__compute_u(xi)
    eta = u*u - 1.0
    v   = self.v0*math.sqrt( 1.0 + eta )
    ne  = self.J/(v*self.qe)
    return ne

  def compute_E(self, x):
    if   (self.case_type == 1):
      E = self.__compute_E_case1(x) 
    elif (self.case_type == 2):
      E = self.__compute_E_case2(x) 
    return E

  def __compute_E_case1(self,x):
    xi  = 2.0*math.sqrt(self.kappa)*x/self.d
    u   = self.__compute_u(xi)
    eta = u*u - 1.0;
    sgn = +1.0
    detadxi = math.sqrt(math.sqrt(1.0+eta) + self.alpha1)
    E       = -2.0*sgn*math.sqrt(self.kappa)*(self.W/self.qe)*detadxi/self.d
    return E

  def __compute_E_case2(self,x):
    xi  = 2.0*math.sqrt(self.kappa)*x/self.d
    u   = self.__compute_u(xi)
    eta = u*u - 1.0;
    if xi <= self.delta:
      sgn = -1.0 
    else:
      sgn = +1.0 
    detadxi = math.sqrt(math.sqrt(1.0+eta) + self.alpha1)
    E       = -2.0*sgn*math.sqrt(self.kappa)*(self.W/self.qe)*detadxi/self.d
    return E

# ------------------------------------------------------------------------------

  def __compute_u(self,xi):

    if   (self.case_type == 1):
      u   = self.__compute_u_case1(xi)
    elif (self.case_type == 2):
      u   = self.__compute_u_case2(xi)
    return u
 
  def __compute_u_case1(self, xi):

    # Solves eqn. (9) for u
    it    = 0
    MAXIT = 100
    TOL   = 1.0E-10
    error = 1.0E+10

    ul   = 1.0
    uh   = 10.0

    sgna = -1.0
    sgnx = -1.0

    um = 0.5*(ul + uh)
    fl = self.__fu_case1(ul,xi,sgna,sgnx)
    fh = self.__fu_case1(uh,xi,sgna,sgnx)
    fm = self.__fu_case1(um,xi,sgna,sgnx)

    while (error > TOL and it < MAXIT):

      f0 = fm

      if ( fm*fl > 0.0 ):
       ul = um
       fl = fm
      else:
       uh = um
       fh = fm

      um = 0.5*(ul + uh)
      fm = self.__fu_case1(um,xi,sgna,sgnx)

      error = abs(fm - f0)
      it += 1

    if (error > TOL or it == MAXIT):
      print ("compute_u_case1 failed: it = ", it, " error = ", error ," xi = ", xi, " um = ", um)

    return um

  def __compute_u_case2(self, xi):

    # Solves eqn. (14-15) for u
    it    = 0
    MAXIT = 100
    TOL   = 1.0E-10
    error = 1.0E+10

    if xi < self.delta:
      sgnd = -1.0
      sgnx = +1.0
      ul   = 1.0 
      uh   = self.umin_exact
    elif xi >= self.delta:
      sgnd = +1.0
      sgnx = -1.0
      ul   = self.umin_exact
      uh   = 1.0

    um = 0.5*(ul + uh)
    fl = self.__fu_case2(ul,xi,sgnd,sgnx)
    fh = self.__fu_case2(uh,xi,sgnd,sgnx)
    fm = self.__fu_case2(um,xi,sgnd,sgnx)

    while (error > TOL and it < MAXIT):

      f0 = fm

      if ( fm*fl > 0.0 ):
       ul = um
       fl = fm
      else:
       uh = um
       fh = fm

      um = 0.5*(ul + uh)
      fm = self.__fu_case2(um,xi,sgnd,sgnx)

      error = abs(fm - f0)
      it += 1

    if (error > TOL or it == MAXIT):
      print ("compute_u_case2 failed: it = ", it, " error = ", error ," xi = ", xi, " um = ", um)

    return um

  def __fu_case1(self, u, xi, sgna, sgnx):

    # Solves eqn. (9) for u
    f = (4.0/3.0)*(u - 2.0*self.alpha1)*math.sqrt(u + self.alpha1) + sgna*self.alpha2 + sgnx*xi
    return f

  def __fu_case2(self, u, xi, sgnd, sgnx):

    # Solves eqn. (14-15) for u
    f = (4.0/3.0)*(u - 2.0*self.alpha1)*math.sqrt(u + self.alpha1) + sgnd*self.delta + sgnx*xi
    return f

  def __compute_alpha1_case2(self, xi0, u0):

    # Solves eqn. (19,13) 
    it    = 0
    MAXIT = 100
    TOL   = 1.0E-12
    error = 1.0E+12

    al = -0.999
    ah = -0.601
    am = 0.5*(al + ah)

    fl = self.__fa_case2(al,xi0,u0)
    fh = self.__fa_case2(ah,xi0,u0)
    fm = self.__fa_case2(am,xi0,u0)

    while (error > TOL and it < MAXIT):

      f0 = fm

      if ( fm*fl > 0.0 ):
       al = am
       fl = fm
      else:  
       ah = am
       fh = fm

      am = 0.5*(al + ah)
      fm = self.__fa_case2(am,xi0,u0)

      error = abs(fm - f0)
      it += 1

    if (error > TOL or it == MAXIT):
      print ("compute_alpha1_case2 failed: it = ", it, " error = ", error, " xi0 = ", self.xi0, " am = ", am)

    return am

  def __compute_alpha1_case1(self, xi0, u0):

    # Solves eqn. (19,13) 
    it    = 0
    MAXIT = 100
    TOL   = 1.0E-12
    error = 1.0E+12

    al = -0.999
    ah =  5.0
    am = 0.5*(al + ah)

    fl = self.__fa_case1(al,xi0,u0)
    fh = self.__fa_case1(ah,xi0,u0)
    fm = self.__fa_case1(am,xi0,u0)

    while (error > TOL and it < MAXIT):

      f0 = fm

      if ( fm*fl > 0.0 ):
       al = am
       fl = fm
      else:
       ah = am
       fh = fm

      am = 0.5*(al + ah)
      fm = self.__fa_case1(am,xi0,u0)

      error = abs(fm - f0)
      it += 1

    if (error > TOL or it == MAXIT):
      print ("compute_alpha1_case1 failed: it = ", it, " error = ", error, " xi0 = ", self.xi0, " am = ", am)

    return am

  def __fa_case2(self, alpha1, xi0, u0):
    # Solves eqn. (19,13) 
    sgnd  = +1.0
    sgnx  = -1.0
    delta = self.__compute_delta(alpha1)
    f     = (4.0/3.0)*(u0 - 2.0*alpha1)*math.sqrt(u0 + alpha1) + sgnd*delta + sgnx*xi0
    return f

  def __fa_case1(self, alpha1, xi0, u0):
    # Solves eqn. (11,12) 
    sgna   = -1.0
    sgnx   = -1.0
    alpha2 = self.__compute_alpha2(alpha1)
    f     = (4.0/3.0)*(u0 - 2.0*alpha1)*math.sqrt(u0 + alpha1) + sgna*alpha2 + sgnx*xi0
    return f

  def __compute_delta(self, alpha1):
    delta = (4.0/3.0)*(1.0 - 2.0*alpha1)*math.sqrt(1.0 + alpha1)
    return delta 

  def __compute_alpha2(self, alpha1):
    alpha2 = (4.0/3.0)*(1.0 - 2.0*alpha1)*math.sqrt(1.0 + alpha1)
    return alpha2

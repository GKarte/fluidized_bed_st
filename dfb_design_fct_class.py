# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 15:23:42 2023

@author: Gregor Karte
"""

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import CoolProp.CoolProp as CP # he1 = CP.PropsSI("H", "P", p, "T", data["Te1"][i]+273.15, "Water")
from CoolProp.HumidAirProp import HAPropsSI
from dataclasses import dataclass
from scipy.optimize import fsolve, minimize_scalar
# from CoolProp.HumidAirProp import HAPropsSI

# constants
g = 9.806 # gravitational acceleration
sigma_boltz = 5.67*10**(-8) # Stefan-Boltzmann Konstante
R_gas = 8.314 # Gaskonstante
# units
units = {'d_p':'m', 'Phi':'-','d_sv':'m', 'rho_p':'kg/m^3', 'subst':'-', 'T':'°C', 'p':'Pa', 'rho_f':'kg/m^3', 'mu':'Pa*s',
         'D':"m", 'A':"m^2", 'h':"m", 'A_q':"m^2", 'U':"m/s", 'Vn_dot':"Nm^3/h", 'V_dot':"m^3/h", 'm_dot':"kg/h",
         'U_to_Umf':"-", 'U_to_Uc':"-", 'U_to_Use':"-", 'Ar':"-", 'Re_p':"-", 'Fr_D':"-", 'Fr_p':"-", 'DR':"-",
         'U_star':"-", 'd_p_star':"-", 'eps_mf':"-", 'U_mf':"m/s", 'U_c_av':"m/s", 'U_t':"m/s", 'U_se':"m/s",
         'Re_mf':"-", 'Re_c_av':"-", 'Re_se':"-"        
        }
df_units = pd.DataFrame(units.values(),units.keys(),columns=["unit"])
# print(df_units)

"""classes"""


@dataclass
class BedMaterial:
    d_p: float
    rho: float
    Phi: float
    
    def __post_init__(self):
        self.d_sv = self.d_p*self.Phi
        
    def __add__(self,other):
        d_p = self.d_p+other.d_p
        v = 1/self.rho+1/other.rho
        Phi = self.Phi+other.Phi
        return self.__class__(d_p,1/v,Phi)
    def __mul__(self,K):
        d_p = self.d_p*K
        v = (1/self.rho)*K
        Phi = self.Phi*K
        return self.__class__(d_p,1/v,Phi)
        

    
@dataclass
class FluidState:
    subst: str
    T: float
    p: float
    rho: float = None
    mu: float = None
    nu: float = None
    # M: float = None
    # lam: float = None
    # Pr: float = None
    # cp: float = None
    
    def __post_init__(self):
        if self.rho==None:
            self.rho = CP.PropsSI("D", "P", self.p, "T", self.T+273.15, self.subst)
            self.mu = CP.PropsSI("V", "P", self.p, "T", self.T+273.15, self.subst)
            self.nu = self.mu/self.rho
            # self.M = CP.PropsSI("M", "P", self.p, "T", self.T+273.15, self.subst)
            # self.lam = CP.PropsSI("L", "P", self.p, "T", self.T+273.15, self.subst)
            # self.Pr = CP.PropsSI("PRANDTL", "P", self.p, "T", self.T+273.15, self.subst)
            # self.cp = CP.PropsSI("C", "P", self.p, "T", self.T+273.15, self.subst)
        
    def __add__(self,other):
        rho = self.rho+other.rho
        mu = self.mu+other.mu
        nu = self.nu+other.nu
        # M = self.M+other.M
        return self.__class__(self.subst+"+"+other.subst,self.T,self.p,rho,mu,nu)
    def __mul__(self,K):
        rho = self.rho*K
        mu = self.mu*K
        nu = self.nu*K
        # M = self.M*K
        return self.__class__(self.subst,self.T,self.p,rho,mu,nu)
    
    
def gas_mix(fluids, mol_fracts):
    M_mix = sum([fl.M*mol_fracts[ind] for ind, fl in enumerate(fluids)])
    mass_fracts = np.array([fl.M*mol_fracts[ind] for ind, fl in enumerate(fluids)])/M_mix
    rho_mix = sum([mass_fracts[ind]/fl.rho for ind, fl in enumerate(fluids)]) ** (-1)
    
    # viscosity
    return mass_fracts, M_mix, rho_mix
        
        
    
        
class FluidizedBed():
    
      def __init__(self, bed, fluid):
          """
          bed --> BedMaterial dataclass
          fluid --> FluidState material dataclass
          """
          self.bed = bed
          self.fluid = fluid
          props = FluidBedProperties(bed, fluid)
          self.Ar, self.Re_mf, self.U_mf, self.Re_se, self.U_se, self.Re_c_av, self.U_c_av, self.DR, self.d_p_star, self.U_t = props
          undefined = 15
          self.U_to_Use, self.U_to_Uc, self.m_dot, self.Vn_dot, self.U_to_Umf, self.Re_p, self.U_star, self.A, self.Fr_D, self.Fr_p, self.U, self.D, self.h, self.A_q, self.V_dot = [None]*undefined
          self.est_eps_mf_ergun()
          
      def est_eps_mf_ergun(self):
          x0 = [0.4]
          param = (self.Re_mf, self.Ar, self.bed.Phi)
          self.eps_mf = fsolve(sys_eps_mf_ergun, x0, args=param)[0]
          
      def set_geometry(self, D=0, h=0, A=None):
          if A == None:
              self.D =          D # reactor diameter
              self.A =          (D/2)**2 * np.pi # cross sectional area
          else:
              self.A = A
              self.D = np.sqrt( self.A*4/np.pi )
          self.h =          h # reactor height
          self.A_q =        self.D*np.pi*self.h # lateral surface for heat exchange with surrounding
          self.Fr_D =       None
          self.Fr_p =       None
          self.U =          None
          self.U_star =     None
          self.Re_p =       None
          self.U_to_Umf =   None
          self.U_to_Uc =    None
          self.U_to_Use =   None
          
      def set_U(self,U):
          if self.D == None:
              print("Diameter not defined. Use \"set_geometry()\" first.")
          else:
              self.U = U
              self.Re_p =       U*self.bed.d_sv / self.fluid.nu # Re
              self.U_star =     self.Re_p/self.d_p_star
              self.Fr_D =       self.U**2/(self.D*g) # Froude
              self.Fr_p =       self.U**2/(self.bed.d_sv*g) # Froude particle
              self.V_dot =      U*self.A*3600
              self.Vn_dot =     V_to_Vn(self.V_dot, self.fluid.T, self.fluid.p)
              self.m_dot =      self.V_dot*self.fluid.rho
              self.U_to_Umf =   self.U/self.U_mf
              self.U_to_Uc =   self.U/self.U_c_av
              self.U_to_Use =   self.U/self.U_se
          
      def set_V_dot(self,V_dot):
          if self.D == None:
              print("Diameter not defined. Use \"set_geometry()\" first.")
          else:
              self.V_dot =      V_dot
              self.Vn_dot =     V_to_Vn(self.V_dot, self.fluid.T, self.fluid.p)
              self.m_dot =      self.V_dot*self.fluid.rho
              self.U =          V_dot/self.A/3600
              self.Fr_D =       self.U**2/(self.D*g) # Froude
              self.Fr_p =       self.U**2/(self.bed.d_sv*g) # Froude
              self.Re_p =       self.U*self.bed.d_sv / self.fluid.nu # Re
              self.U_star =     self.Re_p/self.d_p_star
              self.U_to_Umf =   self.U/self.U_mf
              self.U_to_Uc =   self.U/self.U_c_av
              self.U_to_Use =   self.U/self.U_se
              
      def set_Vn_dot(self,Vn_dot):
          if self.D == None:
              print("Diameter not defined. Use \"set_geometry()\" first.")
          else:
              self.Vn_dot =     Vn_dot
              self.V_dot =      Vn_to_V(self.Vn_dot, self.fluid.T, self.fluid.p)
              self.m_dot =      self.V_dot*self.fluid.rho
              self.U =          self.V_dot/self.A/3600
              self.Fr_D =       self.U**2/(self.D*g) # Froude
              self.Fr_p =       self.U**2/(self.bed.d_sv*g) # Froude
              self.Re_p =       self.U*self.bed.d_sv / self.fluid.nu # Re
              self.U_star =     self.Re_p/self.d_p_star
              self.U_to_Umf =   self.U/self.U_mf
              self.U_to_Uc =   self.U/self.U_c_av
              self.U_to_Use =   self.U/self.U_se
          
      def diameter_boundaries(self, V_dot, regime="BFB"):
          if regime == "CFB":
              Umin = self.U_se
              Umax = self.U_se*10**4
          elif regime == "BFB":
              Umin = self.U_mf
              Umax = self.U_c_av
          elif regime == "TFB":
              Umin = self.U_c_av
              Umax = self.U_se
          A_min =          V_dot/Umax/3600
          D_min =          np.sqrt(A_min/np.pi)*2
          A_max =          V_dot/Umin/3600
          D_max =          np.sqrt(A_max/np.pi)*2
          text = f"D_min = {D_min*1000} mm \nD_max = {D_max*1000} mm"
          print(text)
          return D_min, D_max
          
      def operational_boundaries(self, regime="BFB"):
          if self.D == None:
              print("Diameter not defined. Use \"set_geometry()\" first.")
          else:
              if regime == "CFB":
                  Umin = self.U_se
                  Umax = self.U_se*10**4
              elif regime == "BFB":
                  Umin = self.U_mf
                  Umax = self.U_c_av
              elif regime == "TFB":
                  Umin = self.U_c_av
                  Umax = self.U_se
              V_dot_min =      self.A*Umin*3600
              V_dot_max =      self.A*Umax*3600
              text = f"V_dot_min = {V_dot_min} m^3/h \nV_dot_max = {V_dot_max} m^3/h"
              print(text)
          return V_dot_min, V_dot_max
      
      def dsv_bounderies(self):
          d_ar = np.linspace(self.bed.d_sv*0.01, self.bed.d_sv*50, 2000)
          d_min = None
          d_max = None
          for d in d_ar:
              bed_i = BedMaterial(d, self.bed.rho, 1)
              U_t_i = TerminalVelocity(bed_i, self.fluid)
              U_mf_i = FluidBedProperties(bed_i, self.fluid)[2]
              if d_min == None and U_t_i > self.U:
                  d_min = d
              elif d_max == None and U_mf_i > self.U:
                  d_max = d
              elif d_max != None and d_min != None:
                  break
          return d_min, d_max
              
                        
          
      def calc_operating_points(self,V_dot: np.ndarray,plot=False,ax=None):
          if self.D == None:
              print("Diameter not defined. Use \"set_geometry()\" first.")
          else:
              U =          V_dot/self.A/3600 # U
              Re_p =       U*self.bed.d_sv / self.fluid.nu # Re
              U_star =     Re_p/self.d_p_star
              d_p_star =   np.array([self.d_p_star]*len(V_dot))
              if plot:
                  ax.scatter(d_p_star, U_star ,label="operating points",marker="x",color="chocolate")
              return U, d_p_star, U_star
          
      def dp_bed(self, h_bed, eps_bed=None): # Kaiser 2003
          if eps_bed == None:
              eps_bed = self.eps_bed()
          dp_bed = (1-eps_bed)*self.bed.rho*h_bed*g
          m_bed = dp_bed*self.A / g
          return dp_bed, m_bed, eps_bed
      
      def eps_bed(self): # Kaiser 2003
          delta_b_bed = delta_b(self.U, self.U_mf, self.bed.d_p)
          eps_bed = delta_b_bed + (1-delta_b_bed)*self.eps_mf
          return eps_bed

      def dp_bottom(self, h_bed, dp_bed=None):
          if dp_bed==None:
              dp_bed=self.dp_bed(h_bed)[0]
          dp_bottom = dp_bed * ( 0.01 + 0.2*( 1 - np.exp(-self.D/(2*h_bed)) ) )
          return dp_bottom
      
      def Vn_p_corr(self, p_corr):
          return V_to_Vn(self.Vn_dot, 0, p_corr)
      
      def Ut_theory(self):
          Ut_th = TerminalVelocity(self.bed, self.fluid)
          return Ut_th
          
      def downscale_glicksmanSimpl(self, fluid_cm, ratio):
          """
          CFM scaling based on Fr_D, U/U_mf, density ratio=const. (simplified Glicksman rule):
              rho_p_cm, d_p_cm, U_cm, D_cm = f(fluid_cm, geom. ratio)
              point in Grace != const.
              rho_p usually not changeable
          """
          h_cm = self.h*ratio
          D_cm = self.D*ratio
          rho_p_cm = self.DR*fluid_cm.rho
          param = (fluid_cm.nu, fluid_cm.rho, self.U_to_Umf, self.Fr_D, D_cm, rho_p_cm)
          x0 =    [self.Ar, self.Re_mf, self.U_mf, self.U, self.bed.d_sv]        
          Ar_cm, Re_mf_cm, U_mf_cm, U_cm, d_sv_cm = fsolve(sys_glicksmanSimpl, x0, args=param)          
          # define CFM object
          bed_cm = BedMaterial(d_p=d_sv_cm/self.bed.Phi, rho=rho_p_cm, Phi=self.bed.Phi)
          obj = self.__class__(bed_cm,fluid_cm)
          obj.set_geometry(D_cm, h_cm)
          obj.set_U(U_cm)         
          # comparison 
          df = self.createTable()
          df.rename(columns = {df.columns[1]:'hot'}, inplace = True)
          df_cold = obj.createTable()
          df["cold"] = df_cold[df_cold.columns[1]]
          df=df.drop(["h", "A_q"])
          # print(df)
          df["ratio"] = df["hot"]/df["cold"]          
          return obj, df
      
      def downscale_glicksmanFull(self, fluid_cm):
          """
          CFM scaling based on Fr_p, Ar, density ratio, D/d_sv=const. (full Glicksman rule):
              rho_p_cm, d_p_cm, U_cm, D_cm = f(fluid_cm, rho_p_cm, geom. ratio)
              point in Grace != const.
              rho_p usually not changeable
          """
          rho_p_cm = self.DR*fluid_cm.rho
          param = (fluid_cm.nu, fluid_cm.rho, self.Ar, self.Fr_p, self.bed.d_sv, self.D, rho_p_cm)
          x0 =    [self.U, self.bed.d_sv, self.D]        
          U_cm, d_sv_cm, D_cm = fsolve(sys_glicksmanFull, x0, args=param)
          ratio = D_cm/self.D
          h_cm = self.h*ratio
          # define CFM object
          bed_cm = BedMaterial(d_p=d_sv_cm/self.bed.Phi, rho=rho_p_cm, Phi=self.bed.Phi)
          obj = self.__class__(bed_cm, fluid_cm)
          obj.set_geometry(D_cm, h_cm)
          obj.set_U(U_cm)         
          # comparison 
          df = self.createTable()
          df.rename(columns = {df.columns[1]:'hot'}, inplace = True)
          df_cold = obj.createTable()
          df["cold"] = df_cold[df_cold.columns[1]]
          df=df.drop(["h", "A_q"])
          df["ratio"] = df["hot"]/df["cold"]          
          return obj, df
      
      def downscale_marx(self, fluid_cm, rho_p_cm, Phi_cm, ratio): # marx diss
          """
          CFM scaling based on Ar=const. and Re_p=const. (point in Grace = const.):
              d_p_cm, U_cm, D_cm = f(fluid_cm, rho_p_cm, geom. ratio)
              Fr_D, density ratio != const. (violation of Glicksman)
              similar to Marx Dissertation
          """
          h_cm = self.h*ratio
          D_cm = self.D*ratio        
          d_sv_cm = (self.Ar * (g*(rho_p_cm-fluid_cm.rho) / (fluid_cm.nu**2*fluid_cm.rho))**(-1))**(1/3)
          U_cm = self.Re_p*fluid_cm.nu / d_sv_cm
          bed_cm = BedMaterial(d_p=d_sv_cm/Phi_cm, rho=rho_p_cm, Phi=Phi_cm)
          # define CFM object
          obj = self.__class__(bed_cm,fluid_cm)
          obj.set_geometry(D_cm, h_cm)
          obj.set_U(U_cm)
          # comparison 
          df = self.createTable()
          df.rename(columns = {df.columns[1]:'hot'}, inplace = True)
          df_cold = obj.createTable()
          df["cold"] = df_cold[df_cold.columns[1]]
          df=df.drop(["h", "A_q"])
          df["ratio"] = df["hot"]/df["cold"]            
          return obj, df
      
      def downscale_detU(self, fluid_cm, bed_cm, ratio, param="Fr_p"): # proell paper
          """
          CFM scaling based on Fr_p=const.
              d_p_cm, U_cm, D_cm = f(fluid_cm, rho_p_cm, geom. ratio)
              Fr_D, density ratio != const. (violation of Glicksman)
              similar to Marx Dissertation
          """
          h_cm = self.h*ratio
          D_cm = self.D*ratio
          obj = self.__class__(bed_cm,fluid_cm)
          obj.set_geometry(D_cm, h_cm)
          if param == "Fr_p":
              U_cm = np.sqrt(self.Fr_p*g*bed_cm.d_sv)
          elif param == "Re_p":
              U_cm = self.Re_p*fluid_cm.nu/bed_cm.d_sv
          elif param == "Umf":
              U_cm = obj.U_mf*self.U_to_Umf
          elif param == "min":
              param = (self.Re_p, self.Fr_D, self.Fr_p, self.U_to_Umf, bed_cm.d_sv, fluid_cm.nu, obj.U_mf, D_cm)
              res = minimize_scalar(obj_detU_min, args=param)
              print(res)
              U_cm = res.x
          elif param == "Use":
              U_cm = obj.U_se*self.U/self.U_se
          elif param == "Uc":
              U_cm = obj.U_c_av*self.U/self.U_c_av
          elif param == "Ut":
              U_cm = obj.U_t*self.U/self.U_t
              
          # set U of CFM object
          obj.set_U(U_cm)
          # comparison 
          df = self.createTable()
          df.rename(columns = {df.columns[1]:'hot'}, inplace = True)
          df_cold = obj.createTable()
          df["cold"] = df_cold[df_cold.columns[1]]
          df=df.drop(["h", "A_q"])
          df["ratio"] = df["hot"]/df["cold"]            
          return obj, df
      
      def downscale_set(self, fluid_cm, bed_cm, ratio, Vn_dot):
          """
          CFM scaling based on Fr_p=const.
              d_p_cm, U_cm, D_cm = f(fluid_cm, rho_p_cm, geom. ratio)
              Fr_D, density ratio != const. (violation of Glicksman)
              similar to Marx Dissertation
          """
          h_cm = self.h*ratio
          D_cm = self.D*ratio
          obj = self.__class__(bed_cm,fluid_cm)
          obj.set_geometry(D_cm, h_cm)
          obj.set_Vn_dot(Vn_dot)
          # comparison 
          df = self.createTable()
          df.rename(columns = {df.columns[1]:'hot'}, inplace = True)
          df_cold = obj.createTable()
          df["cold"] = df_cold[df_cold.columns[1]]
          df=df.drop(["h", "A_q"])
          df["ratio"] = df["hot"]/df["cold"]            
          return obj, df
      
           
      
      def createTable(self):
          attr = vars(self)
          attr.update(vars(self.bed))
          attr.update(vars(self.fluid))
          attr["rho_p"] = self.bed.rho
          attr["rho_f"] = self.fluid.rho
          df = pd.DataFrame(columns=["unit","FB"], index=units.keys())
          df["unit"] = df_units 
          for key in units.keys():
              df["FB"][key] = attr[key]
          df=df.drop(["subst"])
          return df
      
      def createTable_op(self, full=False):
          ind = ["U", "U_to_Umf", "Vn_dot", "V_dot"]
          df = pd.DataFrame(columns=["unit","value"], index=ind)
          for i in ind:
              df["unit"][i] = units[i]
          df["value"]["U"] = self.U
          df["value"]["U_to_Umf"] = self.U_to_Umf
          df["value"]["Vn_dot"] = self.Vn_dot
          df["value"]["V_dot"] = self.V_dot
          return df
              
      def grace(self, ax, label="label", marker='*', color="red"):
          ax.scatter(self.d_p_star, self.U_star,label=label, marker=marker, color=color)
      def grace_ccc(self, ax, k=0.771):
          d = np.log10(self.U_star)-k*np.log10(self.d_p_star)
          dp_ar = np.logspace(np.log10(3), np.log10(20), num=20, base=10)
          U_ar = 10 ** (k*np.log10(dp_ar)+d)
          ax.plot(dp_ar, U_ar, "--", label="ccc-iso", color="black")

          
          
      
"""functions"""

def Archimedes(d,rho_p,rho_f,nu_f):
    Ar =         d**3*g*(rho_p-rho_f) / (nu_f**2*rho_f)
    return Ar

def FluidBedProperties(bed, fluid):
    Ar =         Archimedes(bed.d_sv,bed.rho,fluid.rho,fluid.nu)
    Re_mf =      np.sqrt(27.2**2 + 0.0408*Ar)-27.2
    U_mf =       Re_mf*fluid.nu / bed.d_sv
    Re_se =      1.53*np.sqrt(Ar) # Re min. fluid.
    U_se =       Re_se*fluid.nu / bed.d_sv # U min. fluid.
    Re_c_av =    Ar**(19/30) / (0.85+0.85*Ar**(1/5))
    U_c_av =     Re_c_av*fluid.nu / bed.d_sv # U min. fluid.
    DR =         (bed.rho-fluid.rho)/fluid.rho # density ratio
    d_p_star =   Ar**(1/3) # d_p_star
    Re_t =       (18/d_p_star**2 + (2.335-1.744*0.95)/d_p_star**0.5) ** (-1) * d_p_star
    U_t =        Re_t*fluid.nu / bed.d_sv # U t
    return Ar, Re_mf, U_mf, Re_se, U_se, Re_c_av, U_c_av, DR, d_p_star, U_t

def TerminalVelocity(bed, fluid):
    U_lam = (bed.rho-fluid.rho)*bed.d_sv**2*g / (18*fluid.nu*fluid.rho)
    Re_lam = bed.d_sv*U_lam/fluid.nu
    U_turb = np.sqrt(4/3 * (bed.rho-fluid.rho)/fluid.rho * bed.d_sv*g/0.43)
    Re_turb = bed.d_sv*U_turb/fluid.nu
    # print(Re_lam, U_lam, Re_turb, U_turb)

    if Re_lam <= 0.2:
        U_t = U_lam
    elif Re_turb >= 1000:
        U_t = U_turb
    else:
        param_trans = (bed.d_sv, fluid.nu, bed.rho, fluid.rho)
        # x0 = [3, 10, 2]
        x0 = [(Re_lam+Re_turb)/2, 8, 2]
        Re_trans, C_w_trans, U_trans = fsolve(sys_U_t_trans, x0, args=param_trans)
        U_t = U_trans
        # print(Re_trans,C_w_trans, U_trans)
    return U_t

def sys_U_t_trans(x, *parameter):
    Re, C_w, U_t = x
    d, nu, rho_p, rho_f = parameter
    return [U_t     - np.sqrt(4/3 * (rho_p-rho_f)/rho_f * d*g/C_w),
            Re      - d*U_t/nu,
            C_w     - (24/Re + 4/np.sqrt(Re) + 0.4) 
            ]


def sys_glicksmanSimpl(x, *parameter):
    Ar_cm, Re_mf_cm, U_mf_cm, U_cm, d_sv_cm = x
    nu_f_cm, rho_f_cm, U_to_Umf, Fr_D, D_cm, rho_p_cm = parameter
    return [Ar_cm -      d_sv_cm**3*g*(rho_p_cm-rho_f_cm) / (nu_f_cm**2*rho_f_cm),
            Re_mf_cm -      (np.sqrt(27.2**2 + 0.0408*Ar_cm)-27.2),
            U_mf_cm -       Re_mf_cm*nu_f_cm / d_sv_cm,
            U_to_Umf -      U_cm/U_mf_cm,
            Fr_D -       U_cm**2/(D_cm*g)]

def sys_glicksmanFull(x, *parameter):
    U_cm, d_sv_cm, D_cm = x
    nu_f_cm, rho_f_cm, Ar, Fr_p, d_sv, D, rho_p_cm = parameter    
    return [Ar -      d_sv_cm**3*g*(rho_p_cm-rho_f_cm) / (nu_f_cm**2*rho_f_cm),
            D_cm/d_sv_cm - D/d_sv,
            Fr_p -       U_cm**2/(d_sv_cm*g)]

def sys_scale_marx(x, *parameter):
    Re_mf_cm, U_mf_cm, U_cm, d_sv_cm = x
    nu_f_cm, rho_f_cm, U_to_Umf, Ar_cm, Re_p, D_cm, rho_p_cm = parameter
    
    return [Ar_cm -         d_sv_cm**3*g*(rho_p_cm-rho_f_cm) / (nu_f_cm**2*rho_f_cm), 
            Re_mf_cm -      (np.sqrt(27.2**2 + 0.0408*Ar_cm)-27.2),
            U_mf_cm -       Re_mf_cm*nu_f_cm / d_sv_cm,
            Re_p -          U_cm*d_sv_cm/nu_f_cm]

def obj_detU_min(U, *parameter):
    Re, Fr_D, Fr, X, d_sv_cm, nu_f_cm, U_mf_cm, D_cm = parameter
    Re_cm = U*d_sv_cm/nu_f_cm
    Fr_cm = U**2/(g*d_sv_cm)
    Fr_D_cm = U**2/(g*D_cm)
    X_cm = U/U_mf_cm
    dRe = abs(Re_cm-Re)/Re
    dFr = abs(Fr_cm-Fr)/Fr
    dFr_D = abs(Fr_D_cm-Fr_D)/Fr_D
    dX = abs(X_cm-X)/X
    res = np.sqrt(dRe + dFr_D + dFr + dX)
    return res

def sys_eps_mf_ergun(x, *parameter): # Kaiser 2003
    eps_mf = x[0]
    Re_mf, Ar, Phi = parameter
    K1 = 150*Re_mf / Ar
    K2 = (150*Re_mf+1.75*Re_mf**2) / Ar
    # K1 = 150*Re_mf / Ar / Phi**2
    # K2 = (150*Re_mf+1.75*Re_mf**2*Phi) / Ar / Phi**2
    return [eps_mf**3 + eps_mf*K1 - K2]

def createGrace(title="no"):
    pct_name = 'grace_schmid.png'
    grace_png = plt.imread(pct_name)
    fig, ax = plt.subplots(1, 1, figsize=(6,7.2), dpi=100)
    fig.tight_layout(pad=0.5)
    xmin, xmax, ymin, ymax = 0.5, 100, 0.005, 20
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('$d_p^* = Ar^{1/3}$')
    ax.set_ylabel('$U^* = Re/Ar^{1/3}$')
    ax.set_zorder(2)
    ax.set_facecolor('none')
    ax_tw_x = ax.twinx()
    ax2 = ax_tw_x.twiny()
    ax_tw_x.axis('off')
    ax2.axis('off')
    ax2.imshow(grace_png, extent=[xmin, xmax, ymin, ymax], aspect='auto')
    if title != "no":
        ax.set_title(title)
    return fig, ax

def Vn_to_V(Vn, T, p):
    V = Vn * 1.01325*10**5/p * (T+273.15)/273.15
    return V

def V_to_Vn(V, T, p):
    Vn = V * p/(1.01325*10**5) * 273.15/(T+273.15)
    return Vn

def delta_b(U, Umf, d_p):
    delta_b = ( 1 + 1.3*(0.15+U-Umf)**0.33*(U-Umf)**(-0.8) / (0.26+0.7*np.exp(-3.3*d_p)) ) ** (-1)
    return delta_b

def pg_prop_T(T): # Bickel DA, extrapolation
    rho = 0.21*(766+273)/(T+273)
    mu = 3.3*10**(-5) + (T-766)/(965-766)*(3.8-3.3)*10**(-5)
    nu = mu/rho
    return rho, mu, nu


def alpha_BFB(FB, T_surf, cp_bed, emiss_surf=0.9, media=None):
    """
    Berechnung des Waermeuebergangskoeffizienten in einer Wirbelschicht. gemäss VDI-Waermeatlas
    Funktioniert nur fuer folgende Medien: H2, H2O, Luft, CO2
    """
    fluid = FB.fluid
    bed = FB.bed
    d_sv = bed.d_sv
    T_bed = fluid.T
    eps = FB.eps_bed()
    # radiation
    T_m = (T_surf+T_bed)/2
    alpha_r = 4*emiss_surf*sigma_boltz* (T_m+273.15)**3
    
    # gas convection
    Nu_g = 0.009 * fluid.Pr**0.33 * FB.Ar**0.5
    alpha_g = Nu_g*fluid.lam/d_sv
    
    # particle convection
    C_A_dict = {"H2":10000, "Water":3.62, "H2O":3.62, "Air":2.8, "air":2.8, "CO2":2.8}
    if media == None:
        media = FB.fluid.subst
    C_A = C_A_dict[media]
    C_K = 2.6
    gamma = ( 10**(0.6-(1000/(T_bed+273.15)+1)/C_A) + 1 )**(-1)
    Lambda = np.sqrt(2*np.pi*R_gas*(T_bed+273.15)/fluid.M) * fluid.lam/(fluid.p*(2*fluid.cp-R_gas/fluid.M))
    l = 2*Lambda*(2/gamma-1)
    Nu_WP_max = 4 * ( (1+2*l/d_sv)*np.log(1+d_sv/(2*l)) -1 )
    Z = 1/6 * bed.rho*cp_bed/fluid.lam * np.sqrt ( g*d_sv**3*(eps-FB.eps_mf) / ( 5*(1-eps)*(1-FB.eps_mf) ) )
    N = Nu_WP_max / (C_K*Z)
    Nu_p = (1-eps)*Z*(1-np.exp(-N))
    alpha_p = Nu_p*fluid.lam/d_sv
    alpha = alpha_p + alpha_g + alpha_r
    return alpha, alpha_r, alpha_g, alpha_p, Lambda, gamma, l, Z, N, Nu_WP_max
    # return alpha
    
# def sys_nozzles(x, *param):
#     d_c, U_c = x
#     V_dot, dp_bottom, dp_tube, d_i, N, xi = param
#     U = V_dot/N / (d_i**2*np.pi/4)
#     return [U*(d_i**2*np.pi/4) - U_c*(d_c**2*np.pi/4),
#             dp_c - (dp_bottom-dp_tube),
            
# def Re(U,L,fluid):
#     return U*L/fluid.nu
# def dp_gas(xi, fluid, U):
#     return xi*fluid.rho*U**2/2          

# def calc_nozzles(fluid, V_dot, dp_bottom, L_tube, d_tube, d_i, N, xi_c):
#     U_tube = V_dot / (d_tube**2*np.pi/4)
#     Re_tube = Re(U_tube,d_tube,fluid)
#     xi_tube = (1.8*np.log10(Re_tube)-1.5)**(-2) * L_tube/d_tube
#     dp_tube = dp_gas(xi_tube, fluid, U_tube)
#     U_i = V_dot/N / (d_i**2*np.pi/4)
#     xi_i = 1.25
#     dp_i = dp_gas(xi_i, fluid, U_i)
#     dp_c = dp_bottom-dp_tube-dp_i
#     U_c = np.sqrt(2*dp_c/(xi_c*fluid.rho))
#     A_c = V_dot/N / U_c
#     d_c = np.sqrt(A_c*4/np.pi)
#     return d_c, Re_tube, dp_i, dp_c, dp_tube, U_i, U_c
            
            
            
# def dp_blende(fluid, V_dot, d_i, D):
#     phi = (d_i/D)**2
#     alpha = 0.6+0.4*phi
#     # xi = 1.8
#     xi = (1/alpha - phi)**2
#     if xi < 1:
#         xi=1
#     Ui = V_dot / (d_i**2*np.pi/4)
#     U = V_dot / (D**2*np.pi/4)
#     dp = dp_gas(xi,fluid,Ui)
#     print(Ui, U, xi)
#     return dp
    

"""old"""


if __name__ == '__main__':
    bed1 = BedMaterial(100*10**(-6), 8900, 1)
    fluid = FluidState("Water", 965, 1.01325*10**5)
    air20 = FluidState("Air", 20, 1.01325*10**5)
    OBJ = FluidizedBed(bed1, air20)
    OBJ.set_geometry(0.036)
    OBJ.set_U(15*OBJ.U_mf)
    print(OBJ.createTable())
    
    GF_o = FluidizedBed(bed1, fluid)
    GF_o.set_geometry(0.08), GF_o.set_Vn_dot(24)
    GF_o.downscale_glicksmanSimpl(air20, 0.5)
    p = 1.01325*10**5
    T=965
    h2o = FluidState("Water", T, p)
    n2 = FluidState("N2", T, p)
    co2 = FluidState("CO2", T, p)
    h2 = FluidState("H2", T, p)
    ch4 = FluidState("CH4", T, p)
    d_co = CP.PropsSI("D", "P", p, "T", T+273.15, "CO")
    co = FluidState("CH4", T, p)
    co_rho = d_co
    # fl = (h2o*0.3 + co2*0.07 + h2*0.4 + ch4*0.08 + co*0.15)
    # dd = h2o*0.7+co2*0.3
    # fluid = FluidState("Water", 965, 1.01325*10**5)
    # test = gas_mix([h2o,co2,h2,ch4,co,n2],[0.3572,0.2205*(1-0.3572),0.5037*(1-0.3572),0.0863*(1-0.3572),0.177*(1-0.3572),0.0111*(1-0.3572)])
    # print(test)
    # rho_test, mu_test, nu_test = pg_prop_T(965)
    # print(rho_test, mu_test, nu_test)
    
    # test eps: Hoeltl Diss. S. 53
    bed_eps = BedMaterial(300*10**(-6), 2700, 0.9)
    fluid_eps = FluidState("Air", 700, 1.01325*10**5)
    OBJ_eps = FluidizedBed(bed_eps, fluid_eps)
    OBJ_eps.set_geometry(0.036)
    U_eps = np.linspace(OBJ_eps.U_mf*1.05, OBJ_eps.U_se, 50)
    UtoUmf_eps = U_eps/OBJ_eps.U_mf
    eps_eps = []
    for U in U_eps:
        OBJ_eps.set_U(U)
        eps_eps.append(OBJ_eps.eps_bed())
    feststoffanteil_eps = np.ones(len(eps_eps)) - np.array(eps_eps)
    fig, ax = plt.subplots(1, 1, figsize=(5,5), dpi=200)
    # ax.set_xlim((0,10))
    # ax.set_ylim((0,0.04))
    ax.plot(UtoUmf_eps, feststoffanteil_eps)
    ax.set_ylabel("Feststoffanteil (1-eps)")
    ax.set_xlabel("U/U_mf")
    ax.grid()
    
    # test Vn_p_corr()
    Vn_corr = OBJ_eps.Vn_p_corr(1.1*10**5)
    print(f"Pressure correction of Vn: {OBJ_eps.Vn_dot}, {Vn_corr}")
    
    
    # heat transfer
    # bed = BedMaterial(300*10**(-6), 2700, 1)
    # fluid = FluidState("Water", 730, 1.01325*10**5)
    # FB = FluidizedBed(bed, fluid)
    # FB.set_geometry(0.1), FB.set_U(4*FB.U_mf)
    # res = alpha_BFB(FB, T_surf=400, cp_bed=800, emiss_surf=0.6)
    # print(res)
    
    # dp_bottom
    
    # dp_bot = OBJ.dp_bottom(0.15)
    # dp_rat = dp_bot/OBJ.dp_bed(0.3)[0]
    # dp_bot = OBJ.dp_bed(0.15)[0]*0.2
    # # res_dp_bot = calc_nozzles(OBJ.fluid, OBJ.V_dot/3600, dp_bot, 2, 0.01, 0.008, 4, 1.5)
    # # print(dp_bot, dp_rat)
    # # print(res_dp_bot)
    # N_ar = [1,2,3,4,5]
    # dp_ar = []
    # for N in N_ar:
    #     V_doti = OBJ.V_dot / 3600 / N
    #     print(N)
    #     dp_ar.append((N,dp_blende(OBJ.fluid, V_doti, 0.0025, 0.009)))
    # print(dp_ar)
    # print(dp_bot)
    
    
    

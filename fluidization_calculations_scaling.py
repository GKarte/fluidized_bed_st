# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 10:00:10 2024

@author: Gregor
"""

from matplotlib import pyplot as plt
import CoolProp.CoolProp as CP # he1 = CP.PropsSI("H", "P", p, "T", data["Te1"][i]+273.15, "Water")
import numpy as np
import pandas as pd
from scipy import integrate
from dfb_design_fct_class import BedMaterial, FluidState, FluidizedBed, createGrace, pg_prop_T, alpha_BFB
import additional_functions as gf
pd.options.display.precision = 3


# geometry
# gasifier
Ac_to_A = 0.3
D_ccc = 0.11
D_c = D_ccc*np.sqrt(Ac_to_A)
L_f = 0.3
W_f = 0.3
L_u = 0.15
W_u = 0.1
# segregator
L_seg=0.230
W_seg=0.07
D_seg_o = 0.028
# converter
L_conv=0.15
W_conv=0.08
# riser
D_r = 0.08
# connections
D_con = 0.048
# loop seals
L_ls = 0.12
W_ls = 0.05

# Querschnittsflächen für rechteckige Wirbelschichten
A_seg = W_seg*L_seg
A_conv = W_conv*L_conv
A_f = L_f*W_f
A_u = L_u*W_u
A_ls = L_ls*W_ls

# operating conditions
p = 1.01325 * 10**5

# media
steam650 = FluidState("Water", T=650 , p=p)
steam700 = FluidState("Water", T=700 , p=p)
steam750 = FluidState("Water", T=750 , p=p)
steam800 = FluidState("Water", T=800 , p=p)
steam900 = FluidState("H2O", T=900 , p=p)
air20 = FluidState("Air", T=20 , p=p)
air22 = FluidState("Air", T=22 , p=p)
air800 = FluidState("Air", T=800 , p=p)
air900 = FluidState("Air", T=900 , p=p)
pg650 = FluidState("Water", T=650 , p=p)
pg650.subst = "pg"
pg650.rho, pg650.mu, pg650.nu = pg_prop_T(650)
pg700 = FluidState("Water", T=700 , p=p)
pg700.subst = "pg"
pg700.rho, pg700.mu, pg700.nu = pg_prop_T(700)
pg800 = FluidState("Water", T=800 , p=p)
pg800.subst = "pg"
pg800.rho, pg800.mu, pg800.nu = pg_prop_T(800)
pg900 = FluidState("Water", T=900 , p=p)
pg900.subst = "pg"
pg900.rho, pg900.mu, pg900.nu = pg_prop_T(900)

# bed material
rho_CaO = 1500 # Benedikt Diss
rho_CaCO3 = 2650 # Benedikt Diss
dp = 300
rho_calcined = rho_CaO
rho_carbonated = 1/(0.2*1/rho_CaCO3 + 0.8*1/rho_CaO)
bed_carbonated = BedMaterial(d_p=dp*10**(-6), rho=rho_carbonated, Phi=1)
bed_calcined = BedMaterial(d_p=dp*10**(-6), rho=rho_calcined, Phi=1)
bronze = BedMaterial(d_p=72*10**(-6), rho=8720, Phi=1)
# bed = CaO

# fluidization criteria
# U_r = 7.
# UtoUt_seg_o = 1.
UtoUmf_seg_u = 2
UtoUmf_ls = 5.
# UtoUmf_gf_in = 7.
# UtoUmf_conv = 7.

Vn_conv_in = 5
Vn_conv_out = 8.6

Vn_gr_in = 6.8
Vn_gr_out = 15.4



"""fluidized beds"""
# segregator
SEG_u = FluidizedBed(bed, steam700) # unten
SEG_u.set_geometry(A=A_seg), SEG_u.set_U(SEG_u.U_mf*UtoUmf_seg_u)

SEG_o = FluidizedBed(bed, steam700) # oben
SEG_o.set_geometry(D=D_seg_o), SEG_o.set_U(UtoUt_seg_o*SEG_o.U_t)

SEG_m = FluidizedBed(bed, steam700) # mitte
SEG_m.set_geometry(A=A_seg), SEG_m.set_Vn_dot(SEG_o.Vn_dot)

# converter
CONV_in = FluidizedBed(bed, steam900) # unten vor gsf
CONV_in.set_geometry(A=A_conv), CONV_in.set_Vn_dot(Vn_conv_in)

CONV_out = FluidizedBed(bed, pg900) # nach gsf
CONV_out.set_geometry(A=A_conv), CONV_out.set_Vn_dot( Vn_conv_out )

# gasifier
GF_BFB_in = FluidizedBed(bed, steam700) # unten vor gsf
GF_BFB_in.set_geometry(A=A_u), GF_BFB_in.set_Vn_dot( Vn_gr_in ) 

GF_BFB_out = FluidizedBed(bed, pg700) # unten nach gsf
GF_BFB_out.set_geometry(A=A_u), GF_BFB_out.set_Vn_dot( Vn_gr_out )


Vn_ccc = Vn_gr_out+Vn_conv_out+SEG_o.Vn_dot


GF_ccc = FluidizedBed(bed, pg800)
GF_ccc.set_geometry(D=D_ccc), GF_ccc.set_Vn_dot(Vn_ccc)
GF_c = FluidizedBed(bed, pg800)
GF_c.set_geometry(D=D_c), GF_c.set_Vn_dot(Vn_ccc)

# riser
R = FluidizedBed(bed, steam900) # oben, diluted
R.set_geometry(D=D_r), R.set_U(U_r)

# loop seals
LLS = FluidizedBed(bed, steam700)
LLS.set_geometry(A=A_ls), LLS.set_U(LLS.U_mf*UtoUmf_ls)


CFM, df_CFM = CONV_in.downscale_detU(air22, bronze, 2, param="Umf")
CFM_GF_in, df_CFM_GF_in = GF_BFB_in.downscale_detU(air22, bronze, 2, param="Umf")
CFM_R, df_CFM_R = R.downscale_detU(air22, bronze, 2, param="Use")

"""diagrams"""
# Grace
# fig,ax = createGrace()
# GF_BFB_in.grace(ax,"GR BFB inlet","o","blue")
# GF_BFB_out.grace(ax,"GR BFB after gasif.","s","blue")
# GF_ccc.grace(ax,"CCC","^","blue")
# GF_c.grace(ax,"CCC constr.","x","blue")
# SEG_u.grace(ax,"SEG lower zone","o","green")
# SEG_m.grace(ax,"SEG upper zone","s","green")
# SEG_o.grace(ax,"SEG transport zone","^","green")
# CONV_in.grace(ax,"CONV inlet","o","orange")
# CONV_out.grace(ax,"CONV after gasif.","s","orange")
# R.grace(ax,"Riser","o","red")
# LLS.grace(ax,"Loop seals","o","magenta")
# # ax.axvline(CFM.d_p_star, label="CFM", linestyle="--", color='purple')
# # CFM.grace(ax,"CFM: CONV inlet","v",'purple')
# # CFM_R.grace(ax,"CFM: Riser","*",'purple')
# ax.legend(loc="upper right")

fig,ax = createGrace()
GF_BFB_in.grace(ax,"GR BFB inlet","o","blue")
GF_BFB_out.grace(ax,"GR BFB after gasif.","s","blue")
GF_ccc.grace(ax,"CCC","^","blue")
GF_c.grace(ax,"CCC constr.","x","blue")
# SEG_u.grace(ax,"SG lower zone","o","green")
# SEG_m.grace(ax,"SG upper zone","s","green")
# SEG_o.grace(ax,"SG outlet","^","green")
CONV_in.grace(ax,"CV inlet","o","orange")
CONV_out.grace(ax,"CV outlet","s","orange")
R.grace(ax,"Riser","o","red")
LLS.grace(ax,"Loop seals","o","magenta")
ax.axvline(CFM.d_p_star, label="CFM", linestyle="--", color='purple')
CFM_GF_in.grace(ax,"CFM: GR BFB inlet","v",'purple')
CFM_R.grace(ax,"CFM: Riser","*",'purple')
ax.legend(loc="upper right")




# heat transfer
# alpha = alpha_BFB(LLS, 750, cp_bed=800, emiss_surf=0.6, media=None)[0]
# print(alpha)




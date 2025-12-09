# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 18:51:55 2023

@author: Gregor Karte
"""

from io import BytesIO
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas as pd
# import thermo

def fig2img(fig):
    """
    Konvertiert matplotlib figure in ein png.

    Parameters
    ----------
    fig : matplotlib.figure
        Matplotlib Plot.

    Returns
    -------
    img : png
        png des plots.

    """
    buf = BytesIO()
    fig.savefig(buf, format='png')
    img = buf
    img.seek(0)
    # img = Image.open(buf)
    return img

def str_date_time():
    # using now() to get current time
    current_time = datetime.datetime.now()
    string = f'{current_time.day}{current_time.month}{current_time.year}_{current_time.hour}{current_time.minute}'
    return string


def read_ipse_txt(txt_file):
    df = pd.read_csv(txt_file, sep="\t\t", header=None, index_col=0, engine='python')
    return df





# def calc_gas_mixture_props(subst_dict, T, p):
    
#     # visc_cond_methods = dict([('nitrogen', "VDI_PPDS"),
#     #                     ('oxygen', "VDI_PPDS"),
#     #                     ('methane', "VDI_PPDS"),
#     #                     ('carbon dioxide', "VDI_PPDS"),
#     #                     ('carbon monoxide', "VDI_PPDS"),
#     #                     ('water', "REFPROP_FIT"),
#     #                     ('ethene', "VDI_PPDS"),
#     #                     ('ethane', "VDI_PPDS"),
#     #                     ('hydrogen', "VDI_PPDS")                        
#     #                     ])
    
#     visc_cond_methods = dict([('nitrogen', "COOLPROP"),
#                         ('methane', "VDI_PPDS"),
#                         ('carbon dioxide', "COOLPROP"),
#                         ('carbon monoxide', "VDI_PPDS"),
#                         ('water', "COOLPROP"),
#                         ('ethene', "VDI_PPDS"),
#                         ('ethane', "VDI_PPDS"),
#                         ('hydrogen', "DIPPR_PERRY_8E"),
#                         ('oxygen', "COOLPROP")                        
#                         ])
    
    
#     comp = list(subst_dict.keys())
#     zs = np.array(list(subst_dict.values()))
#     zs = zs/sum(zs)
#     m = thermo.Mixture(IDs=comp, zs=zs, T=T, P=p)

#     visc_gases = []
#     cond_gases = []
#     for ind in range(len(m.names)):
#         visc_gas = thermo.viscosity.ViscosityGas(CASRN=m.CASs[ind],method=visc_cond_methods[m.names[ind]])
#         visc_gases.append(visc_gas)
#         cond_gas = thermo.thermal_conductivity.ThermalConductivityGas(CASRN=m.CASs[ind],method=visc_cond_methods[m.names[ind]])
#         cond_gases.append(cond_gas)
    
#     m_visc_mixt_obj = thermo.viscosity.ViscosityGasMixture(MWs=m.MWs, molecular_diameters=m.molecular_diameters, Stockmayers=m.Stockmayers, CASs=m.CASs, ViscosityGases=visc_gases)
#     m_cond_mixt_obj = thermo.thermal_conductivity.ThermalConductivityGasMixture(MWs=m.MWs, Tbs=m.Tbs, CASs=m.CASs, ThermalConductivityGases=cond_gases, ViscosityGases=visc_gases)
    
#     mu = m_visc_mixt_obj.calculate(T=T, P=p, zs=m.zs, ws=m.ws, method="BROKAW") # BROKAW
#     lambd = m_cond_mixt_obj.calculate(T=T, P=p, zs=m.zs, ws=m.ws, method="LINDSAY_BROMLEY") # LINDSAY_BROMLEY
#     rho = m.rho
#     cp = m.Cp
#     Pr = cp*mu/lambd

#     nu = mu/rho
#     return mu, rho, nu, cp, lambd, Pr


# if __name__ == "__main__":
    
#     # print(read_ipse_txt(r"C:\Users\Gregor\Downloads\PG_SER_20MW.txt"))
    
#     comp_dict = dict([  ('N2', 0),
#                         ('O2', 0),
#                         ('CH4', 0),
#                         ('CO2', 1.),
#                         ('carbon monoxide', 0.0),
#                         ('H2O', 0),
#                         ('C2H4', 0.000),
#                         ('C2H6', 0.000),
#                         ('H2', 0.000)
#                         ])
    
#     print(calc_gas_mixture_props(comp_dict, T=700+273.15, p=1.01325E5))
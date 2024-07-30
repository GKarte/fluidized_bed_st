# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 13:47:14 2023

@author: Gregor Karte
"""

import streamlit as st
from dfb_design_fct_class import BedMaterial, FluidState, FluidizedBed, createGrace, pg_prop_T
import general_functions as gf
import pandas as pd
import numpy as np
from io import BytesIO
from matplotlib import pyplot as plt


def create_plot(figsize=(6, 5), dpi=200, x_range=(0,1), y_range=(0,1), x_label="x", y_label="y", second_ax=False, y2_range=None, y2_label="y2", title=None, grid=True, grid_fine=True):
    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)
    fig.tight_layout()
    xmin, xmax = x_range
    ymin, ymax = y_range
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if grid:
        ax.grid()
    if grid_fine:
        ax.grid(which='minor', color='lightgray', linestyle=':', linewidth=0.5)
        ax.minorticks_on()
    if title:
        ax.set_title(title)
    if second_ax:
        ax2 = ax.twinx()
        y2min, y2max = y2_range
        ax2.set_ylim(y2min, y2max)
        ax2.set_ylabel(y2_label)
        return fig, (ax, ax2)
    else:
        return fig, ax

# pd.options.display.precision = 2

def scaling(bed_cm, fluid_cm, criteria="Umf", ratio=2):
    fullGlicks = FB.downscale_glicksmanFull(fluid_cm)
    simpleGlicks = FB.downscale_glicksmanSimpl(fluid_cm, 1/ratio)
    used = FB.downscale_detU(fluid_cm, bed_cm, 1/ratio, criteria)
    
    df_cm = fullGlicks[1]
    df_cm.rename(columns = {'cold':'full Glicksman', 'ratio':'f. Gl. ratio'}, inplace = True)
    df_cm["simple Glicksman"] = simpleGlicks[1]["cold"]
    df_cm["s. Gl. ratio"] = simpleGlicks[1]["ratio"]
    df_cm["used"] = used[1]["cold"]
    df_cm["used ratio"] = used[1]["ratio"]
    
    fig_cm, ax_cm = createGrace()
    FB.grace(ax_cm, "Hot", "o", "red")
    fullGlicks[0].grace(ax_cm, "full Glicksman", ".", "blue")
    simpleGlicks[0].grace(ax_cm, "simple Glicksman", "x", "green")
    used[0].grace(ax_cm, "used", "*", "magenta")
    ax_cm.legend(loc="upper right")
    
    return df_cm, fig_cm


def createEpsDiagram():
    fig, ax = create_plot(figsize=(6, 5), dpi=200, x_range=(0,40), y_range=(0,1), x_label="$U/U_{mf}$", y_label="$\epsilon$")
    FB_dummy = initialize_FB(bed, fluid, geom_val, OP_val)
    eps_list = []
    U_list = list(np.linspace(FB_dummy.U_mf, FB_dummy.U_mf*50, 51))
    for U in U_list:
        FB_dummy.set_U(U)
        eps_list.append(FB_dummy.eps_bed())
    eps_list.insert(0, eps_list[0])
    U_list.insert(0, 0)
    U_ar = np.array(U_list)
    ax.plot(U_ar/FB_dummy.U_mf, eps_list, label="$\epsilon$")
    ax.legend()
    return fig, ax
        
    

def initialize_FB(bed, fluid, geom_val, OP_val):
    FB = FluidizedBed(bed, fluid)
    if st.session_state.geom == "D / mm":
        FB.set_geometry(D=geom_val/1000)
    elif st.session_state.geom == "A / mm^2":
        FB.set_geometry(A=geom_val/1000000)
        
    if st.session_state.OP == "U_to_Umf / -":
        FB.set_U(OP_val*FB.U_mf)
    elif st.session_state.OP == "U_to_Ut / -":
        FB.set_U(OP_val*FB.U_t)
    elif st.session_state.OP == "U_to_Use / -":
        FB.set_U(OP_val*FB.U_se)
    elif st.session_state.OP == "U_to_Uc / -":
        FB.set_U(OP_val*FB.U_c_av)
    elif st.session_state.OP == "U / (m/s)":
        FB.set_U(OP_val)
    elif st.session_state.OP == "V_dot / (m^3/h)":
        FB.set_V_dot(OP_val)
    elif st.session_state.OP == "Vn_dot / (Nm^3/h)":
        FB.set_Vn_dot(OP_val)
    return FB


st.title("Fluidized Bed Calculations")

col1, col2 = st.columns(2)

with st.sidebar:
    st.header("Settings")
    st.markdown("**Bed material:**")
    dp = st.number_input(label="particle diameter / 10^-6 m", value=300) * 10**(-6)
    rho_s = st.number_input(label="density / (kg/m^3)", value=3000)
    Phi = st.number_input(label="sphericity / -", value=0.9)
    
    
    st.markdown("**Media:**")
    subst = st.selectbox("Substance", ("Water", "Air", "CO2", "N2", "product gas (approx.)", "manual"))
    if subst == "manual":
        pf = st.number_input(label="pressure / bar", value=1.01325, format="%0.5f") * 10**(5)
        Tf = st.number_input(label="temperature / °C", value=650.)
        fluid = FluidState("Water", T=Tf , p=pf)
        fluid.subst = "manual"
        fluid.rho = st.number_input(label="density / (kg/m^3)", value=1., format="%0.5f")
        fluid.nu = st.number_input(label="kin. viscosity / (m^2/s)", value=0.000001, format="%0.7f")
        fluid.mu = fluid.nu*fluid.rho
    elif subst == "product gas (approx.)":
        pf = st.number_input(label="pressure / bar", value=1.01325, format="%0.5f") * 10**(5)
        Tf = st.number_input(label="temperature / °C", value=650.)
        fluid = FluidState("CO2", T=Tf , p=pf)
        fluid.subst = "pg"
        fluid.rho, fluid.mu, fluid.nu = pg_prop_T(Tf)
    else:
        pf = st.number_input(label="pressure / bar", value=1.01325, format="%0.5f") * 10**(5)
        Tf = st.number_input(label="temperature / °C", value=650.)
        fluid = FluidState(subst, T=Tf , p=pf)
    
    
    st.markdown("**Geometry:**")
    geom_switch = st.selectbox("Set parameter", ("D / mm", "A / mm^2"), key="geom")
    geom_val = st.number_input("Set value", value=70.)
    
    
    st.markdown("**Fludization:**")
    OP_switch = st.selectbox("Set parameter", ("U_to_Umf / -", "U_to_Ut / -", "U_to_Use / -", "U_to_Uc / -",
                                               "U / (m/s)", "V_dot / (m^3/h)", "Vn_dot / (Nm^3/h)"), key="OP" )
    OP_val = st.number_input("Set value", value=8.0)  



bed = BedMaterial(d_p=dp, rho=rho_s, Phi=Phi)
FB = initialize_FB(bed, fluid, geom_val, OP_val)


tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(["Grace", "Summary", "Pressure drop", "Scaling", "Particle boundaries", "Diagrams"])

with tab1:
    if st.button("Create Grace diagram"):
        fig, ax = createGrace()
        FB.grace(ax, "Operating point", "o", "blue")
        st.pyplot(fig)
        # st.markdown("Diagram taken from: Schmid, J. C. (2014). Development of a novel dual fluidized bed gasification system for increased fuel flexibility [Dissertation, Technische Universität Wien]. reposiTUm. https://doi.org/10.34726/hss.2014.25397")
    with st.expander("References"):
        st.write('''
            Diagram taken from: Schmid, J. C. (2014). Development of a novel dual fluidized bed gasification system for increased fuel flexibility [Dissertation, Technische Universität Wien]. reposiTUm. https://doi.org/10.34726/hss.2014.25397
            ''')

with tab2:
    style_dict1 = {"FB":"{:.2e}"}
    df = FB.createTable()
    st.table(df.style.format(style_dict1))
    # download button
    output = BytesIO()
    writer = pd.ExcelWriter(output, engine = 'xlsxwriter')
    df.to_excel(writer, sheet_name="Fluidized bed", float_format="%.5f", startrow=0)
    writer.close()
    filename = f'FluidizedBedNumbers_{gf.str_date_time()}'
    st.download_button(
        label="Download as xlsx-file",
        data=output.getvalue(),
        file_name= f'{filename}.xlsx'
    )
    
with tab3:
    h_bed = st.number_input(label="bed height in fluidized state / m", value=0.2)
    eps_bed = st.number_input(label="porosity in fluidized state (estimated value based on set U/Umf) / -", value=FB.eps_bed())
    dp_bed = FB.dp_bed(h_bed, eps_bed)[0]
    dp_bot = FB.dp_bottom(h_bed, dp_bed=dp_bed)
    st.write(f"""
             Bed: $\Delta p$ = {dp_bed:.4} Pa\n
             Bottom (required): $\Delta p$ = {dp_bot:.4} Pa\n
             Pressure drop ratio (bottom/bed): {dp_bot/dp_bed*100:.4} %\n
             $\epsilon$ = {eps_bed:.4}
             """)
    if st.button("Create diagram"):
        d_rel_op = FB.D/(2*h_bed)
        dp_rel_op = 0.01 + 0.2*( 1 - np.exp(-d_rel_op) )
        d_rel = np.arange(0, 10, 0.2)
        dp_rel = 0.01 + 0.2*( 1 - np.exp(-d_rel) )
        fig_dp_bottom, ax_dp_bottom = create_plot(figsize=(6, 5), dpi=200, x_range=(0,10), y_range=(0,0.25), x_label="$D/(2H)$", y_label="$dp_{bottom}/dp_{bed}$")
        ax_dp_bottom.plot(d_rel,dp_rel)
        ax_dp_bottom.scatter([d_rel_op],[dp_rel_op], c="r", label="operating point")
        ax_dp_bottom.legend()
        st.pyplot(fig_dp_bottom)
             
with tab4:  
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**Cold bed material**")
        dp_cm = st.number_input(label="cold particle diameter / 10^-6 m", value=80) * 10**(-6)
        rho_s_cm = st.number_input(label="cold density / (kg/m^3)", value=8900)
        Phi_cm = st.number_input(label="cold sphericity / -", value=0.9)
    with col2:
        st.markdown("**Cold media**")
        subst_cm = st.selectbox("cold substance", ("Air", "CO2"))
        p_cm = st.number_input(label="cold pressure / bar", value=1.01325, format="%0.5f") * 10**(5)
        T_cm = st.number_input(label="cold temperature / °C", value=20.)
    
    st.markdown("**Scaling settings**")
    ratio = st.number_input(label="Geometric ratio", value=2.)
    criteria = st.selectbox("Scaling parameter", ("Umf", "Ut", "Use", "Uc","Re_p", "Fr_p"), key="crit")
    
    fluid_cm = FluidState(subst=subst_cm, T=T_cm , p=p_cm)
    bed_cm = BedMaterial(d_p=dp_cm, rho=rho_s_cm, Phi=Phi_cm)
    if st.button("Scale"):
        df_cm, fig_cm = scaling(bed_cm, fluid_cm, criteria, ratio)
        style_dict2 = {"hot":"{:.2e}", "full Glicksman":"{:.2e}", "simple Glicksman":"{:.2e}", "used":"{:.2e}",
                       "f. Gl. ratio":"{:.2f}", "s. Gl. ratio":"{:.2f}", "used ratio":"{:.2f}"}
        st.dataframe(df_cm.style.format(style_dict2))
        
        # download button
        output = BytesIO()
        writer = pd.ExcelWriter(output, engine = 'xlsxwriter')
        df_cm.to_excel(writer, sheet_name="Scaled FBs", float_format="%.5f", startrow=0)
        writer.close()
        filename = f'ScaledFBs_{gf.str_date_time()}'
        st.download_button(
            label="Download as xlsx-file",
            data=output.getvalue(),
            file_name= f'{filename}.xlsx'
        )
        
        st.pyplot(fig_cm)
        
with tab5:
    if st.button("Calculate"):
        try:
            dmin, dmax = FB.dsv_bounderies()
            st.write(f"d_sv_min = {dmin*10**6:.0f} * 10^-6 m (smaller particles get entrained)")
            st.write(f"d_sv_max = {dmax*10**6:.0f} * 10^-6 m (bigger particles aren't fluidized)")
        except:
            st.write("Calculation not successful")
            
            
with tab6:
    if st.button("Create diagrams"):
        st.write("Porosity of fluidized bed ( m^3_fluid / (m^3_fluid+m^3_particles) ) as function of U/Umf (Source: Kaiser 2003, DOI:10.1016/S0009-2509(03)00233-1):")
        fig_eps, ax_eps = createEpsDiagram()
        st.pyplot(fig_eps)
        

        
        
        

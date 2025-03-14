# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 13:47:14 2023

@author: Gregor Karte
"""

import streamlit as st
from dfb_design_fct_class import BedMaterial, FluidState, FluidizedBed, createGrace, pg_prop_T
import additional_functions as af
import pandas as pd
import numpy as np
from io import BytesIO
from matplotlib import pyplot as plt


# colors = ["red", "green", "blue", "magenta", "black", "orange", "cyan"]
colors = ['#377eb8', '#ff7f00', '#4daf4a',
          '#f781bf', '#a65628', '#984ea3',
          '#999999', '#e41a1c', '#dede00']
markers = ["x","o", "v", "^","*","+", "."]

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
    simpleGlicksLowRe = FB.downscale_glicksmanSimplLowRe(fluid_cm, 1/ratio, bed_cm.rho)
    used = FB.downscale_detU(fluid_cm, bed_cm, 1/ratio, criteria)
    
    df_cm = fullGlicks[1]
    df_cm.rename(columns = {'cold':'full Glicksman', 'ratio':'f. Gl. ratio'}, inplace = True)
    df_cm["simple Glicksman"] = simpleGlicks[1]["cold"]
    df_cm["s. Gl. ratio"] = simpleGlicks[1]["ratio"]
    df_cm["simple Glicksman (low Re)"] = simpleGlicksLowRe[1]["cold"]
    df_cm["s. Gl. ratio (low Re)"] = simpleGlicksLowRe[1]["ratio"]
    df_cm["used"] = used[1]["cold"]
    df_cm["used ratio"] = used[1]["ratio"]
    
    fig_cm, ax_cm = createGrace()
    FB.grace(ax_cm, "Hot", "o", "red")
    fullGlicks[0].grace(ax_cm, "full Glicksman", ".", "blue")
    simpleGlicks[0].grace(ax_cm, "simple Glicksman", "x", "green")
    simpleGlicksLowRe[0].grace(ax_cm, "simple Glicksman (low Re)", "^", "orange")
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

def export_dfb_design_settings(df_beds, df_fluids, df_zones):
    output = BytesIO()
    writer = pd.ExcelWriter(output, engine = 'xlsxwriter')
    df_zones.to_excel(writer, sheet_name="FB_zones", float_format="%.5f", startrow=0, index=False)
    df_beds.to_excel(writer, sheet_name="bed_materials", float_format="%.5f", startrow=0, index=False)
    df_fluids.to_excel(writer, sheet_name="fluids", float_format="%.5f", startrow=0, index=False)
    writer.close()
    return output.getvalue()

def export_dfb_design_results(df_beds, df_fluids, df_zones, df_results):
    output = BytesIO()
    writer = pd.ExcelWriter(output, engine = 'xlsxwriter')
    df_results.to_excel(writer, sheet_name="results", float_format="%.5f", startrow=0)
    df_zones.to_excel(writer, sheet_name="FB_zones", float_format="%.5f", startrow=0, index=False)
    df_beds.to_excel(writer, sheet_name="bed_materials", float_format="%.5f", startrow=0, index=False)
    df_fluids.to_excel(writer, sheet_name="fluids", float_format="%.5f", startrow=0, index=False)
    writer.close()
    return output.getvalue()

def calculate_dfb_design(df_beds, df_fluids, df_zones):
    fig, ax = createGrace()
    FB_zones = {}
    df_results = None
    reactors = {}
    for ind in df_zones.index:
        zone = df_zones.loc[ind]
        name = df_zones.loc[ind, "name"]
        bed_name = df_zones.loc[ind, "bed material"]
        fluid_name = df_zones.loc[ind, "fluid"]
        bed = df_beds.loc[df_beds['name'] == bed_name]
        fluid = df_fluids.loc[df_fluids['name'] == fluid_name]
        bed_zone = BedMaterial(d_p=float(bed["d_sv / 10^-6m"])/10**6, rho=float(bed["density / kg/m^3"]), Phi=1)
        if str(fluid_name) in known_fluids:
            fluid_zone = FluidState(subst=fluid_name, T=float(zone["T / °C"]) , p=float(zone["p / bar"])*10**5)
        else:
            fluid_zone = FluidState("Water", T=float(zone["T / °C"]) , p=float(zone["p / bar"])*10**5)
            fluid_zone.subst = fluid_name
            fluid_zone.rho = float(fluid["density / kg/m^3"])
            fluid_zone.nu = float(fluid["kin. viscosity / mm^2/s"])/10**6
            fluid_zone.mu = fluid_zone.nu*fluid_zone.rho
        if zone["shape"] == "rectangular":
            A_zone = zone["W / mm"] * zone["L / mm"] / 10**6
        else:
            A_zone = zone["D / mm"]**2*np.pi/4 / 10**6
            
        FB_zone = def_FB_zone(bed_zone, fluid_zone, A=A_zone)
        
        if zone["Fluid. parameter"] == "U_to_Umf / -":
            FB_zone.set_U(float(zone["Fluid. value"])*FB_zone.U_mf)
        elif zone["Fluid. parameter"] == "U_to_Ut / -":
            FB_zone.set_U(float(zone["Fluid. value"])*FB_zone.U_t)
        elif zone["Fluid. parameter"] == "U_to_Use / -":
            FB_zone.set_U(float(zone["Fluid. value"])*FB_zone.U_se)
        elif zone["Fluid. parameter"] == "U_to_Uc / -":
            FB_zone.set_U(float(zone["Fluid. value"])*FB_zone.U_c_av)
        elif zone["Fluid. parameter"] == "U / (m/s)":
            FB_zone.set_U(float(zone["Fluid. value"]))
        elif zone["Fluid. parameter"] == "V_dot / (m^3/h)":
            FB_zone.set_V_dot(float(zone["Fluid. value"]))
        elif zone["Fluid. parameter"] == "Vn_dot / (Nm^3/h)":
            FB_zone.set_Vn_dot(float(zone["Fluid. value"]))
        
        FB_zones[name] = FB_zone
        if df_results is not None:
            df_zone_i = FB_zone.createTable()
            df_results[name] = df_zone_i[df_zone_i.columns[1]]  
        else:
            df_results = FB_zone.createTable()
            df_results.rename(columns = {df_results.columns[1]:name}, inplace = True)
        
        reactor = name.split(" ")[0]
        if reactor not in reactors.keys():
            reactors[reactor] = 0
        else:
            reactors[reactor] += 1
        
        FB_zone.grace(ax, name, markers[reactors[reactor]], colors[list(reactors.keys()).index(reactor)])
    
    ax.legend()
    st.pyplot(fig)   
    st.dataframe(df_results.style.format(subset=list(df_results.columns)[1:], formatter="{:.2e}",thousands="", decimal="."), use_container_width=True)
    
    return df_results




def def_FB_zone(bed, fluid, A):
    FB_zone = FluidizedBed(bed, fluid)
    FB_zone.set_geometry(A=A)
    return FB_zone

def load_dfb_settings(df_beds, df_fluids, df_zones, ipse_items):
    edited_df_beds = placeholder_beds.data_editor(df_beds, num_rows="dynamic")
    if ipse_items is None:
        edited_df_fluids = placeholder_fluids.data_editor(df_fluids,
                                                          column_config={       
                                                              "IPSE stream": st.column_config.SelectboxColumn(required=False,
                                                                                                              options=[]),
                                                              "IPSE object": st.column_config.SelectboxColumn(required=False,
                                                                                                              options=[])
                                                              },
                                                              num_rows="dynamic")
    else:
        edited_df_fluids = placeholder_fluids.data_editor(df_fluids,
                                                          column_config={       
                                                              "IPSE stream": st.column_config.SelectboxColumn(required=False,
                                                                                                              options=ipse_items),
                                                              "IPSE object": st.column_config.SelectboxColumn(required=False,
                                                                                                              options=ipse_items)
                                                              },
                                                              num_rows="dynamic")
    edited_df_zones = placeholder_zones.data_editor(df_zones,
                               column_config={
                                   "bed material": st.column_config.SelectboxColumn(required=True,
                                    options=list(edited_df_beds["name"])),
                                   "fluid": st.column_config.SelectboxColumn(required=True,
                                    options=known_fluids+list(edited_df_fluids["name"])),
                                   "shape": st.column_config.SelectboxColumn(required=True,
                                    options=["circular","rectangular"]),
                                   "Fluid. parameter": st.column_config.SelectboxColumn(required=True,
                                    options=fluid_parameters)
                                   },
                               num_rows="dynamic")
    return edited_df_beds, edited_df_fluids, edited_df_zones
    
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
    subst = st.selectbox("Substance", ("Water", "Air", "CO2", "N2", "gas mixture", "product gas (approx.)", "manual"))
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
    elif subst == "gas mixture":
        pf = st.number_input(label="pressure / bar", value=1.01325, format="%0.5f") * 10**(5)
        Tf = st.number_input(label="temperature / °C", value=650.)
        
        H2O = st.number_input(label="y_H2O / -", value=0.45, format="%0.3f")
        N2 = st.number_input(label="y_N2_dry / -", value=0.0, format="%0.3f")
        H2 = st.number_input(label="y_H2_dry / -", value=0.42, format="%0.3f")
        CH4 = st.number_input(label="y_CH4_dry / -", value=0.12, format="%0.3f")
        CO2 = st.number_input(label="y_CO2_dry / -", value=0.21, format="%0.3f")
        CO = st.number_input(label="y_CO_dry / -", value=0.22, format="%0.3f")
        C2H4 = st.number_input(label="y_C2H4_dry / -", value=0.02, format="%0.3f")
        C2H6 = st.number_input(label="y_C2H6_dry / -", value=0.01, format="%0.3f")
        
        mixture_dict = dict([('N2', N2*(1-H2O)),
                            ('O2', 0),
                            ('H2O', H2O),
                            ('H2', H2*(1-H2O)),
                            ('CH4', CH4*(1-H2O)),
                            ('CO2', CO2*(1-H2O)),
                            ('carbon monoxide', CO*(1-H2O)),
                            ('C2H4', C2H4*(1-H2O)),
                            ('C2H6', C2H6*(1-H2O)),
                            ])
        fluid = FluidState(mixture_dict, T=Tf , p=pf)
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


tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs(["Grace", "Summary", "Pressure drop", "Scaling", "Particle boundaries", "Diagrams", "DFB design"])

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
    st.write(f"Gas medium: {fluid.subst}")
    style_dict1 = {"FB":"{:.2e}"}
    df = FB.createTable()
    st.table(df.style.format(style_dict1))
    # download button
    output = BytesIO()
    writer = pd.ExcelWriter(output, engine = 'xlsxwriter')
    df.to_excel(writer, sheet_name="Fluidized bed", float_format="%.5f", startrow=0)
    writer.close()
    filename = f'FluidizedBedNumbers_{af.str_date_time()}'
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
        style_dict2 = {"hot":"{:.2e}", "full Glicksman":"{:.2e}", "simple Glicksman":"{:.2e}", "simple Glicksman (low Re)":"{:.2e}", "used":"{:.2e}",
                       "f. Gl. ratio":"{:.2f}", "s. Gl. ratio":"{:.2f}", "s. Gl. ratio (low Re)":"{:.2f}", "used ratio":"{:.2f}"}
        st.dataframe(df_cm.style.format(style_dict2))
        
        # download button
        output = BytesIO()
        writer = pd.ExcelWriter(output, engine = 'xlsxwriter')
        df_cm.to_excel(writer, sheet_name="Scaled FBs", float_format="%.5f", startrow=0)
        writer.close()
        filename = f'ScaledFBs_{af.str_date_time()}'
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
        
with tab7:
    known_fluids = ["Air", "H2O", "CO2"]
    fluid_parameters = ["U_to_Umf / -", "U_to_Ut / -", "U_to_Use / -", "U_to_Uc / -",
                        "U / (m/s)", "V_dot / (m^3/h)", "Vn_dot / (Nm^3/h)", "m_dot / (kg/h)"]

    st.header("Settings")
    st.markdown("**Bed materials**")
    placeholder_beds = st.empty()
    st.markdown("**Fluids**")
    st.markdown(f"Known fluids: {known_fluids}")
    st.markdown("Additional fluids:")
    placeholder_fluids = st.empty()
    st.markdown("**Reactor zones**")
    placeholder_zones = st.empty()
    col1_tab7, col2_tab7 = st.columns(2)
    with col1_tab7:
        placeholder_export_settings = st.empty()
        placeholder_import_ipse_data = st.empty()
    with col2_tab7:
        placeholder_import_settings = st.empty()
    
        
            
    design_settings_ready = st.button("Calculate", key = "design_settings_ready")
    
    
    txt_ipse = placeholder_import_ipse_data.file_uploader("Import IPSE data (.txt-export file)", key="upload_ipse")
    if txt_ipse is not None:
        df_ipse = af.read_ipse_txt(txt_ipse)
        ipse_items = [i.split('.', 1)[0] for i in list(df_ipse.index)]
        ipse_items = list(dict.fromkeys(ipse_items))
        st.dataframe(df_ipse)
        st.write(ipse_items)
        print(type(ipse_items))

    import_xlsx = placeholder_import_settings.file_uploader("Import settings (.xlsx)", key="upload_settings")
    if import_xlsx is not None and txt_ipse is None:
        df_beds = pd.read_excel(import_xlsx, sheet_name="bed_materials")
        df_fluids = pd.read_excel(import_xlsx, sheet_name="fluids")
        df_zones = pd.read_excel(import_xlsx, sheet_name="FB_zones")
        edited_df_beds, edited_df_fluids, edited_df_zones = load_dfb_settings(df_beds, df_fluids, df_zones, None)
    elif import_xlsx is not None and txt_ipse is not None:
        df_beds = pd.read_excel(import_xlsx, sheet_name="bed_materials")
        df_fluids = pd.read_excel(import_xlsx, sheet_name="fluids")
        df_zones = pd.read_excel(import_xlsx, sheet_name="FB_zones")
        edited_df_beds, edited_df_fluids, edited_df_zones = load_dfb_settings(df_beds, df_fluids, df_zones, ipse_items)
    elif import_xlsx is None and txt_ipse is not None:
        df_beds = pd.DataFrame(
        [
           {"name": "bed_gr", "density / kg/m^3": 2800, "d_sv / 10^-6m": 300},
           {"name": "bed_cr", "density / kg/m^3": 2800, "d_sv / 10^-6m": 300},
        ]
        )
        style_dict2 = {"density / kg/m^3":"{:.2e}"}
        
        df_beds.style.format(precision=2, thousands="", decimal=".")

        
        edited_df_beds = placeholder_beds.data_editor(df_beds,
                                    num_rows="dynamic")
        
        df_fluids = pd.DataFrame(
        [
           {"name": "my personal fluid", "IPSE stream": ipse_items[0], "IPSE object": ipse_items[0], "density / kg/m^3": 1.2, "kin. viscosity / mm^2/s": 12.0},
        ]
        )
        
        df_fluids.style.format(precision=2, thousands="", decimal=".")


        edited_df_fluids = placeholder_fluids.data_editor(df_fluids,
                                                          column_config={       
                                                              "IPSE stream": st.column_config.SelectboxColumn(required=False,
                                                                                                              options=ipse_items),
                                                              "IPSE object": st.column_config.SelectboxColumn(required=False,
                                                                                                              options=ipse_items)
                                                              },
                                                              num_rows="dynamic")
        
        

        try:
            df_zones = pd.DataFrame(
            [
             {"name": "GR inlet",
              "bed material": edited_df_beds.iloc[0, 0],
              "fluid": known_fluids[0],
              "T / °C": 400,
              "p / bar": 1.01325,
              "shape": "circular",
              "D / mm": 80,
              "L / mm": 150,
              "W / mm": 100,
              "Fluid. parameter": fluid_parameters[0],
              "Fluid. value": 8},
             {"name": "GR BFB",
              "bed material": edited_df_beds.iloc[0, 0],
              "fluid": known_fluids[0],
              "T / °C": 400,
              "p / bar": 1.01325,
              "shape": "circular",
              "D / mm": 80,
              "L / mm": 150,
              "W / mm": 100,
              "Fluid. parameter": fluid_parameters[0],
              "Fluid. value": 10},
            ]
            )
            
            df_zones.style.format(precision=2, thousands="", decimal=".")        
            
            
            edited_df_zones = placeholder_zones.data_editor(df_zones,
                                        column_config={
                                            "bed material": st.column_config.SelectboxColumn(required=True,
                                            options=list(edited_df_beds["name"])),
                                            "fluid": st.column_config.SelectboxColumn(required=True,
                                            options=known_fluids+list(edited_df_fluids["name"])),
                                            "shape": st.column_config.SelectboxColumn(required=True,
                                            options=["circular","rectangular"]),
                                            "Fluid. parameter": st.column_config.SelectboxColumn(required=True,
                                            options=fluid_parameters)
                                            },
                                        num_rows="dynamic")
        except:
            st.markdown("Bed materials table and fluids table each must contain minimum one entry!")
    else:
        df_beds = pd.DataFrame(
        [
           {"name": "bed_gr", "density / kg/m^3": 2800, "d_sv / 10^-6m": 300},
           {"name": "bed_cr", "density / kg/m^3": 2800, "d_sv / 10^-6m": 300},
        ]
        )
        style_dict2 = {"density / kg/m^3":"{:.2e}"}
        
        df_beds.style.format(precision=2, thousands="", decimal=".")

        
        edited_df_beds = placeholder_beds.data_editor(df_beds,
                                    num_rows="dynamic")
        
        df_fluids = pd.DataFrame(
        [
           {"name": "my personal fluid", "IPSE stream": None, "IPSE object": None, "density / kg/m^3": 1.2, "kin. viscosity / mm^2/s": 12.0},
        ]
        )
        
        df_fluids.style.format(precision=2, thousands="", decimal=".")


        edited_df_fluids = placeholder_fluids.data_editor(df_fluids,
                                                          column_config={       
                                                              "IPSE stream": st.column_config.SelectboxColumn(required=False,
                                                                                                              options=[]),
                                                              "IPSE object": st.column_config.SelectboxColumn(required=False,
                                                                                                              options=[])
                                                              },
                                                              num_rows="dynamic")
        
        

        try:
            df_zones = pd.DataFrame(
            [
             {"name": "GR inlet",
              "bed material": edited_df_beds.iloc[0, 0],
              "fluid": known_fluids[0],
              "T / °C": 400,
              "p / bar": 1.01325,
              "shape": "circular",
              "D / mm": 80,
              "L / mm": 150,
              "W / mm": 100,
              "Fluid. parameter": fluid_parameters[0],
              "Fluid. value": 8},
             {"name": "GR BFB",
              "bed material": edited_df_beds.iloc[0, 0],
              "fluid": known_fluids[0],
              "T / °C": 400,
              "p / bar": 1.01325,
              "shape": "circular",
              "D / mm": 80,
              "L / mm": 150,
              "W / mm": 100,
              "Fluid. parameter": fluid_parameters[0],
              "Fluid. value": 10},
            ]
            )
            
            df_zones.style.format(precision=2, thousands="", decimal=".")        
            
            
            edited_df_zones = placeholder_zones.data_editor(df_zones,
                                        column_config={
                                            "bed material": st.column_config.SelectboxColumn(required=True,
                                            options=list(edited_df_beds["name"])),
                                            "fluid": st.column_config.SelectboxColumn(required=True,
                                            options=known_fluids+list(edited_df_fluids["name"])),
                                            "shape": st.column_config.SelectboxColumn(required=True,
                                            options=["circular","rectangular"]),
                                            "Fluid. parameter": st.column_config.SelectboxColumn(required=True,
                                            options=fluid_parameters)
                                            },
                                        num_rows="dynamic")
        except:
            st.markdown("Bed materials table and fluids table each must contain minimum one entry!")
     
            
    placeholder_export_settings.download_button(
        label="Export settings",
        data=export_dfb_design_settings(edited_df_beds, edited_df_fluids, edited_df_zones),
        file_name= f'dfb_design_settings_exp_{af.str_date_time()}.xlsx'
        )

        

    
    if design_settings_ready:
        st.header("Results")
        df_results = calculate_dfb_design(edited_df_beds, edited_df_fluids, edited_df_zones)
        # st.session_state.df_dfb_design_results = df_results
    
        st.download_button(
            label="Export results",
            data=export_dfb_design_results(edited_df_beds, edited_df_fluids, edited_df_zones, df_results),
            file_name= f'dfb_design_results_exp_{af.str_date_time()}.xlsx'
            )

            
    with st.expander("References"):
        st.write('''
            Diagram taken from: Schmid, J. C. (2014). Development of a novel dual fluidized bed gasification system for increased fuel flexibility [Dissertation, Technische Universität Wien]. reposiTUm. https://doi.org/10.34726/hss.2014.25397
            ''')
            
        

        
    

    
        

        
        
        

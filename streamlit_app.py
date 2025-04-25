#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 15:32:26 2025

@author: kristianfoerster
"""

import streamlit as st
import matplotlib.pyplot as plt
import numpy as np

from greenroof import GreenRoof

st.set_page_config(page_title="Green Roof Simulation", layout="wide")
st.title("Green Roof Model – Interactive Simulation")

# Input widgets
st.sidebar.header("Simulation Parameters")
k = st.sidebar.number_input("Hydraulic conductivity k [m/d]", value=1600)
b = st.sidebar.number_input("Brooks-Corey parameter b", value=5.5)
n = st.sidebar.number_input("Manning’s roughness n", value=0.08)
length = st.sidebar.slider("Roof length [m]", 5, 20, 10)
slope = st.sidebar.slider("Slope [%]", 0.0, 3.0, 1.0, format='%.1f')
d_init = st.sidebar.number_input("Delta in matrix potential (initial condition) [m]", value=0.15)
rainfall = st.sidebar.number_input("Rainfall total [mm]", value=27)
duration = st.sidebar.slider("Rainfall duration [min]", 1, 30, 15, format='%i')




# Run the model
if st.button("Run Simulation"):
    st.subheader("Simulation Results")
    

    try:
        # Initialize and run the model
        rheight=length*slope/100.
        model = GreenRoof(ksat=k, porosity=0.5, b=b, mannings_n=n, 
                          length=length, rheight=rheight, d_init_pot=d_init,
                          duration=60)
        
        model.set_design_rain(rain_duration=duration, rain_amount=rainfall)
        
        results = model.run()

        outflow = results['Darcy flow'].resample('60S', label='right', closed='right').sum()
        outflow_s = results['Surface runoff'].resample('60S', label='right', closed='right').sum()
        outflow_total = outflow + outflow_s # darcy + surface


        # Example: Plot outflow or another relevant variable
        fig, ax = plt.subplots()
        ax.plot(outflow_total)
        ax.set_xlabel("Time")
        ax.set_ylabel("Flow [l/min]")
        ax.set_title("Green Roof Outflow")
        ax.legend()
        st.pyplot(fig)

        #st.dataframe(results.head())

    except Exception as e:
        st.error(f"Simulation failed: {e}")
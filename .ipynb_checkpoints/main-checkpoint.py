import streamlit as st
import matplotlib.pyplot as plt
from estimator_model import run_model

st.set_page_config(layout="wide")

st.title("Microgrid Policy Dashboard")

# ========================
# Sidebar Inputs
# ========================
st.sidebar.header("Policy Parameters")

tariff = st.sidebar.slider("Tariff ($/kWh)", 0.05, 0.50, 0.20)
no_of_households = st.sidebar.slider("Subsidy (%)", 50.0, 1000, 50)
stoptime = st.sidebar.slider("Dropout Rate", 20, 100, 180)
#discount = st.sidebar.slider("Discount Rate", 0.0, 0.1, 0.05)

# ========================
# Run Model
# ========================
if st.button("Run Simulation"):

    results = run_model(
        tariff=tariff,
        no_of_households=no_of_households,
        stoptime=stoptime
    )
    times = sorted(results.keys())
    
    # Plot 1
    fig1, ax1 = plt.subplots()
    ax1.plot(times, results["DSWF"])
    ax1.set_title("DSWF")
    st.pyplot(fig1)

    # Plot 2
    fig2, ax2 = plt.subplots()
    ax2.plot(times, results["avg_income"])
    ax2.set_title("Average Income")
    st.pyplot(fig2)

    # Plot 3
    fig3, ax3 = plt.subplots()
    ax3.plot(times, results["connected_count"])
    ax3.set_title("Social Welfare")
    st.pyplot(fig3)

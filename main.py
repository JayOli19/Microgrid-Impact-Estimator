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
no_of_households = st.sidebar.slider("no_of_households", 100, 1000, 500)
stoptime = st.sidebar.slider("stoptime", 20, 100, 180)
#discount = st.sidebar.slider("Discount Rate", 0.0, 0.1, 0.05)

# ========================
# Run Model
# ========================
if st.button("Run Simulation"):

    results = run_model(
        tariff_value=tariff,
        no_of_households_value=no_of_households,
        stoptime=stoptime
    )
    times = sorted(results.keys())
    #st.write("Result keys:", list(results[times[0]].keys()))  # Debugging line to check keys in results
    connected_series = [results[t].get("connected_count", 0) for t in times]
    appliance_demand_series = [results[t].get("appliance_demand", 0) for t in times]
    total_cases_series = [results[t].get("total_cases", 0) for t in times]
    average_income_series = [results[t].get("avg_income", 0) for t in times]
    saifi_series = [results[t].get("saifi", 0) for t in times]
    Employment = [results[t].get("Employment", 0) for t in times]
    jobs_available = [results[t].get("Jobs Available", 0) for t in times]
    Total_Demand = [results[t].get("Demand", 0) for t in times]
    income = [results[t].get("Income", 0) for t in times]
    SWF = [results[t].get("SWF", 0) for t in times]
    DSWF = [results[t].get("DSWF", 0) for t in times]
    co2 = [results[t].get("co2_emissions", 0) for t in times]
    average_schooling = [results[t].get("average_schooling", 0) for t in times]


    st.write("Social Welfare Function:", DSWF[-1]) 
    # Plot 1
    fig1, ax1 = plt.subplots()
    ax1.plot(times, DSWF)
    ax1.set_title("DSWF")
    st.pyplot(fig1)

    # Plot 2
    fig2, ax2 = plt.subplots()
    ax2.plot(times, average_income_series)
    ax2.set_title("Average Income")
    st.pyplot(fig2)

    # Plot 3
    fig3, ax3 = plt.subplots()
    ax3.plot(times, connected_series)
    ax3.set_title("Connected Households")
    st.pyplot(fig3)

    # Plot 4
    fig4, ax4 = plt.subplots()
    ax4.plot(times, jobs_available)
    ax4.plot(times, Employment)
    ax4.set_title("Jobs Available")
    st.pyplot(fig4)

     # Plot 5
    fig5, ax5 = plt.subplots()
    ax5.plot(times, average_schooling)
    ax5.set_title("Average Schooling")
    st.pyplot(fig5)
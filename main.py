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
no_of_households = st.sidebar.number_input("Number of Households", min_value=1, value=100)
stoptime = st.sidebar.slider("stoptime", 20, 180, 120)
farm_work_shift = st.sidebar.slider("Farm Work Shift", 0.0, 1.0, 0.5)
microgrid_capacity = st.sidebar.number_input("Microgrid Capacity (kW)", min_value=1, value=1000000)

# ========================
# Run Model
# ========================
if st.button("Run Simulation"):

    results = run_model(
        tariff_value=tariff,
        no_of_households_value=no_of_households,
        stoptime=stoptime,
        farm_work_shift=farm_work_shift,
        microgrid_capacity=microgrid_capacity,
        failure_rate=0.01,
        no_of_components=2.0,
        baseline_household_demand=200.0,
        air_change_rate_daily=24.0,
        kitchen_volume=30.0,
        outdoor_concentration=1e-6,
        primary_completion=0.9,
        lower_secondary_completion=0.5,
        upper_secondary_completion=0.3,
        employment_rate_baseline=0.6,
        dropout_rate_baseline=0.11,
        hourly_farm_wage=2.0,
        hourly_non_farm_wage=3.0,
        job_separation_rate=0.05,
        baseline_monthly_wage_men=100.0,
        baseline_monthly_wage_women=80.0,
        baseline_schooling=5.0,
        electrification_effect_men=1.1,
        electrification_effect_women=0.9
    )
    times = sorted(results.keys())
    
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
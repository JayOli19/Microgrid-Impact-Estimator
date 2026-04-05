from os import times

import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from estimator_model import run_impact_estimator_model

st.set_page_config(layout="wide")
st.title("Microgrid Policy Dashboard")

if "scenarios" not in st.session_state:
    st.session_state.scenarios = {}

st.sidebar.header("Scenario Settings")

scenario_name = st.sidebar.text_input("Scenario Name", "Scenario 1")# Name of the scenario being simulated
runs = st.sidebar.slider("Runs per Scenario", 1, 100, 30)#The number of times each scenario is run before averaging the simulation results
tariff = st.sidebar.slider("Tariff ($/kWh)", 0.05, 0.50, 0.20)#The price per kilowatt-hour (kWh) that households pay for electricity.
subsidy_rate = st.sidebar.slider("Subsidy Rate", 0.0, 1.0, 0.25)#The percentage of electric appliance cost (electric stoves and lighting) that is subsidized, where values closer to 1 indicate a higher subsidy and values closer to 0 indicate a lower subsidy.
no_of_households = st.sidebar.number_input("Households", 1, value=350)#Total number of households within the community that the microgrid can serve
stoptime = st.sidebar.slider("Simulation Time (months)", 20, 420, 120)#The total duration of the simulation in months, where a longer time horizon allows for observing long-term impacts but may also introduce more uncertainty.
farm_work_shift = st.sidebar.slider("Farm Work Shift (for Women)", 0.0, 1.0, 0.5)#The proportion of women who use time saved from electrification for farm work, where values closer to 1 indicate more women move to farm work and values closer to 0 indicate more women move towards non-farm labour.
electrification_effect_women = st.sidebar.slider("Electrification Effect (for Women)", 0.8, 1.2, 1.0)#The multiplier effect of electrification on women's wages, where values greater than 1 indicate a positive impact and values less than 1 indicate a negative impact.
microgrid_capacity = st.sidebar.number_input("Microgrid Capacity (kW)", 1, value=1000)#The maximum capacity of the microgrid in kilowatts (kW)
dropout_reduction_rate = st.sidebar.slider("Dropout Reduction (%)", 0.0, 0.5, 0.25)#The percentage reduction in dropout rate due to electrification

def run_scenario(runs):#Runs the impact estimator model for a given scenario and number of runs, returning a list of results for each run that can be averaged later to account for variability in the simulation outcomes.

    results_list = []

    for _ in range(runs):

        results = run_impact_estimator_model( #Results from the impact estimator model, which simulates the technical, social, and environmental effects of microgrid electrification.
            tariff_value=tariff,
            subsidy_rate=subsidy_rate,
            no_of_households_value=no_of_households,
            stoptime=stoptime,
            farm_work_shift=farm_work_shift,
            microgrid_capacity=microgrid_capacity,
            failure_rate=0.01, #Failure rate of the microgrid, where a higher value indicates more downtime and a lower value indicates better reliability.
            no_of_components=2.0,#The number of components in the microgrid, which can affect the overall reliability and maintenance requirements of the system.
            baseline_household_demand=200.0,#The average electricity demand per household in kilowatt-hours (kWh) before electrification
            air_change_rate_daily=24.0,#The number of times the air within the kitchen of a household is replaced with outdoor air in a day, which can affect indoor air quality and the concentration of pollutants.
            kitchen_volume=30.0,#Volume of the kitchen in cubic meters.
            outdoor_concentration=1e-6, #The concentration of pollutants in the outdoor air in grams per cubic meter (g/m³)
            primary_completion=0.9, #The school completion rate for primary education
            lower_secondary_completion=0.5,#The school completion rate for lower secondary education
            upper_secondary_completion=0.3,#The school completion rate for upper secondary education
            employment_rate_baseline=0.6,#The baseline employment rate in the community before electrification
            dropout_rate_baseline=0.11,#The baseline dropout rate from school before electrification, where a higher value indicates more students leaving school and a lower value indicates better retention.
            hourly_farm_wage=2.0,#The average hourly wage for farm work in dollars
            hourly_non_farm_wage=3.0,#The average hourly wage for non-farm work in dollars
            job_separation_rate=0.05,#The rate at which individuals lose their jobs
            baseline_monthly_wage_men=100.0,#The average monthly wage for men before electrification in dollars
            baseline_monthly_wage_women=80.0,#The average monthly wage for women before electrification in dollars
            baseline_schooling=5.0,#The average years of schooling before electrification
            electrification_effect_men=1.1,#The multiplier effect of electrification on men's wages
            electrification_effect_women=electrification_effect_women,
            dropout_reduction_rate=dropout_reduction_rate
        )

        times = sorted(results.keys())#The time points at which the results are recorded, used to analyse trends and changes over the course of the simulation.

        series = {
            "times": times,
            "connected": [results[t].get("connected_count",0) for t in times],#The number of households connected to the microgrid at each time point.
            "income($)": [results[t].get("avg_income",0) for t in times],#The average income of households at each time point.
            "income_women($)": [results[t].get("average_income_women",0) for t in times],#The average income of women at each time point
            "income_men($)": [results[t].get("average_income_men",0) for t in times],#The average income of men at each time point
            "schooling(years)": [results[t].get("average_schooling",0) for t in times],#The average years of schooling for individuals in the community at each time point.
            "co2": [results[t].get("co2_emissions",0) for t in times],#The amount of CO2 emissions in kilograms at each time point, which can reflect the environmental impact of the microgrid.
            "reliability": [results[t].get("reliability",0) for t in times],#The reliability of the microgrid at each time point, which can indicate the consistency of electricity supply.
            "saidi(hours)": [results[t].get("saidi",0) for t in times],#The System Average Interruption Duration Index (SAIDI) in hours at each time point.
            "impact": [results[t].get("impact_score",0) for t in times],#The overall impact score of the microgrid at each time point, which can indicate the combined economic, social, and environmental effects.
            "welfare": [results[t].get("SWF",0) for t in times],#The Social Welfare Function (SWF) at each time point, which can reflect the overall well-being of the community.
            "dswf": [results[t].get("DSWF",0) for t in times],#The Discounted Welfare Function (DSWF) at each time point, which can indicate the long-term welfare impact of the microgrid, where higher values suggest a more positive long-term outcome.
            "cases": [results[t].get("total_cases",0) for t in times],#The number of cases of a particular health outcome (e.g., IHD cases attributable to PM2.5 exposure) at each time point, which can reflect the health impact of the microgrid.
            "net_microgrid_income": [results[t].get("microgrid_net_income",0) for t in times],#The net income generated by the microgrid at each time point, which can indicate the financial sustainability of the microgrid.
            "microgrid_power_use": [results[t].get("power_use",0) for t in times],#The total power use in kilowatt-hours (kWh) at each time point, which can reflect the energy consumption patterns of the community.
            "appliance_demand": [results[t].get("appliance_demand",0) for t in times],#The total demand for electric appliances at each time point, which can indicate the adoption of electric appliances in the community.
            "no_of_current_working_components": [results[t].get("no_of_current_working_components",0) for t in times],#The number of currently working components in the microgrid at each time point, which can reflect the operational status and reliability of the microgrid.
            "available_jobs": [results[t].get("available_jobs",0) for t in times]#The number of available jobs in the community at each time point, which can indicate the economic impact of the microgrid on employment opportunities.
        }

        results_list.append(series) #A list of dictionaries containing the time series data for each metric across multiple runs of the simulation, which can be averaged to account for variability and provide a more robust analysis of the scenario's outcomes.

    return results_list

def average_runs(results_list):#Averages the results from multiple runs of the simulation for a given scenario, calculating the mean value for each metric across all runs to provide a more stable and representative outcome that accounts for variability in the simulation results.

    avg = {}
    keys = results_list[0].keys()

    for key in keys:

        if key == "times":
            avg[key] = results_list[0]["times"]
        else:
            data = np.array([r[key] for r in results_list])
            avg[key] = np.mean(data, axis=0)

    return avg

if st.button("Run Simulation"):#When the "Run Simulation" button is clicked, the code executes the following steps: it runs the specified number of simulations for the current scenario using the `run_scenario` function, averages the results across all runs using the `average_runs` function, and then stores the averaged results in the session state under the name of the scenario for later comparison and visualization.

    runs_output = run_scenario(runs)
    averaged_results = average_runs(runs_output)

    st.session_state.scenarios[scenario_name] = averaged_results

st.header("Scenario Comparison")

selected_scenarios = st.multiselect(
    "Select Scenarios",
    list(st.session_state.scenarios.keys())
)

metric = st.selectbox(
    "Metric",
    ["connected","net_microgrid_income","microgrid_power_use","appliance_demand","income($)","no_of_current_working_components","income_women($)","income_men($)","schooling(years)","co2","reliability","cases","saidi(hours)","impact","welfare","dswf","available_jobs"]
)

if selected_scenarios:#If the user has selected one or more scenarios from the multiselect dropdown, a line plot is created to compare the selected metric across the chosen scenarios.

    fig, ax = plt.subplots()

    for name in selected_scenarios:

        data = st.session_state.scenarios[name]

        ax.plot(
            data["times"],
            data[metric],
            label=name
        )

    ax.set_title(metric)
    ax.set_xlabel("Time (Months)")
    ax.set_ylabel(metric)
    ax.grid()
    ax.legend()
    st.pyplot(fig)

if selected_scenarios:

    summary = [] 

    for name in selected_scenarios:

        data = st.session_state.scenarios[name]

        summary.append({
            "Scenario": name,
            f"Final {metric}": data[metric][-1] #The final value of the selected metric at the end of the simulation time horizon for each scenario.
        })

    st.table(summary) #Displays a table summarizing the final values of the selected metric for each of the selected scenarios.

if selected_scenarios:

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    metric_left = st.selectbox("Left Axis", ["welfare","cases","connected","net_microgrid_income","income_women($)"])
    metric_right = st.selectbox("Right Axis", ["co2","dswf","reliability","schooling(years)","income_men($)","income($)"])

    for name in selected_scenarios:
        data = st.session_state.scenarios[name]

        ax1.plot(data["times"], data[metric_left], label=f"{name}")
        ax2.plot(data["times"], data[metric_right], linestyle="--")

    ax1.set_ylabel(metric_left)
    ax2.set_ylabel(metric_right)
    ax1.set_xlabel("Time (Months)")
    ax1.grid()
    fig.legend(loc="lower left")
    st.pyplot(fig)

if selected_scenarios:

    summary = [] 

    for name in selected_scenarios:

        data = st.session_state.scenarios[name]

        summary.append({
            "Scenario": name,
            f"Final {metric_left}": data[metric_left][-1], #The final value of the selected metric at the end of the simulation time horizon for each scenario.
            f"Final {metric_right}": data[metric_right][-1] #The final value of the selected metric at the end of the simulation time horizon for each scenario.
        })

    st.table(summary) #Displays a table summarizing the final values of the selected metric for each of the selected scenarios.

if st.button("Clear All Scenarios"):
    st.session_state.scenarios = {}#Clears all stored scenarios from the session state, allowing the user to reset the dashboard and start fresh with new scenarios.
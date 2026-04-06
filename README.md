# Microgrid Project Toolkit – Impact Assessment Model

## Overview
Access to reliable and affordable electricity is essential for social and economic development, yet grid extension is often infeasible in rural and remote communities. This repository contains an open-source impact assessment model designed to support the planning, evaluation, and long-term sustainability of microgrid and mini-grid projects.

The toolkit integrates **system dynamics** and **agent-based modeling** approaches to estimate social, environmental, and technical impacts across the project lifecycle. These would include microgrid adoption, system reliability, education, income, health, fossil fuel use, and CO2 emissions. To measure the impact it includes a multi-criteria impact score and a dicscounted social welfare function (DSWF).

---

## Project Objectives
- Support evidence-based decision-making for microgrid projects
- Capture the interaction between adoption, reliability, and socio-economic outcomes
- Reduce repeated design and planning mistakes through structured data use
- Provide a transparent, extensible modeling framework

---

## Model Structure
The toolkit uses a **hybrid modeling approach**:

- **Agent-Based Model (ABM):**
This is a bottom-up modelling approach that simulates individual agents and their decisions and behaviour.
  - Represents households
  - Models adoption decisions, new approach perchase, income, education, health improvements, and impact score.

- **System Dynamics (SD):**
This is a top-down approach focused on aggregated-level feedback loops to show how the system states change over time.
  - Captures aggregate feedback loops
  - Models reliability, job creation, and discounted social welfare function (DWSF) over time.

Adoption dynamics are based on a **Logit Adoption Model** which accounts for system reliability, household income, fuel cost savings, social influence, and health improvements.

---

## Key Features
- Microgrid adoption modeling
- Reliability indicators (e.g., SAIDI and SAIFI)
- Economic impact estimation (household and business income)
- Health and safety impact representation
- Scenario and sensitivity analysis support
- Fully open-source and extensible


---

## Installation and Requirements
- Python 3.x
- Required libraries listed in `requirements.txt`
- Recommended: virtual environment

Example:
```bash
pip install -r requirements.txt
```
---

## Running the Model

The backend has a single-entry point which is the "run_impact_estimator_model()" function.
```python
results = run_impact_estimator_model(...)
```
It runs the python script “estimator_model.py” which contains the backend impact estimator model and returns a dictionary of estimated social, technical and environmental impacts and indicators for every step.
This can be integrated into a dashboard to analyse and present the estimations.

Due to the probabilistic nature of the model, it is good practice to run the simulation multiple times and take an average.
The model can be used to perform the follwing:
- Comparing tariffs, subsidies, and capacity
- Reliability stress testing
- Monte Carlo averaging
- Sensitivity analysis on key parameters


## Streamlit Dashboard
A prototype dashboard was developed using Streamlit. The dashboard allows users to input policy levers (e.g. tariffs and subsidies), microgrid and community factors (e.g. grid capacity and the number of households), and behavioural responses (e.g. the fraction women that move into farm work). Users can then run various scenarios and compare graphs of the estimated outputs which include the number of households connected, the social, technical, and environmental impact indicators, the composite impact scores, and the DSWF.

The code for the dashboard is in the "main.py" script. To launch the dashboard use the following command:

```bash
streamlit run main.py
```
The dashboard allows users to select that number of times to run the scenario before being averaged by adjusting the “Runs per Scenario” slider.

## Limitations
- Model outputs are scenario-based estimates, not precise forecasts
- The model assumes no prior appliance ownership and completion appliance substitution which may not be realistic
- Tariffs, subsidies, grid reliability, and social influence dynamics are assumed to be homogenous
- Some nonlinear relationships are approximated using lookup tables


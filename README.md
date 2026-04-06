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

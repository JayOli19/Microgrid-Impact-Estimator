from BPTK_Py import Model
from BPTK_Py import sd_functions as sd
from BPTK_Py import Agent
from BPTK_Py import Event
from BPTK_Py import DataCollector
from BPTK_Py import SimultaneousScheduler
import random
import math


DISEASE_PARAMS = {
    "alpha": 0.3,
    "gamma":0.6,
    "delta":0.7,
    "c0":1,
    "baseline_incidence":200/100000
}

APPLIANCE_CATALOG = {
    "electric_cooker": {
        "cost": 500,
        "monthly_kwh": 30,
        "beta_0": -3.5,
        "beta_income": 0.0003,
        "beta_reliability": 5.5,
        "beta_social": 1.0,
        "beta_cost": 0.008,
        "beta_tariff": 3.5
    },
    "tv": {
        "cost": 300,
        "monthly_kwh": 5,
        "beta_0": -2.0,
        "beta_income": 0.0002,
        "beta_reliability": 0.5,
        "beta_social": 3.5,
        "beta_cost": 0.002,
        "beta_tariff": 2.5
    },
    "electric_lighting": {
        "cost": 50,
        "monthly_kwh": 3,
        "beta_0": -4.5,
        "beta_income": 0.0004,
        "beta_reliability": 3.5,
        "beta_social": 3.8,
        "beta_cost": 0.0001,
        "beta_tariff": 2.5
    }
}

FUEL_APPLIANCE_CATALOG = {
    "lpg_stove": {
        "hours_used": 2.5,
        "fuel_used_hour": 0.5,
        "cost_per_unit": 10, 
        "fuel_emission_factor_pm25": 15, #g/L
        "fuel_emission_factor_co2": 2.98,
        "kitchen_fraction": 0.25
    },
    "kerosene_lamp": {
        "hours_used": 5,
        "fuel_used_hour": 0.5,
        "cost_per_unit": 0.5, 
        "fuel_emission_factor_pm25": 10, #g/L
        "fuel_emission_factor_co2": 2.68,
        "kitchen_fraction": 0.3
    }
}

# --- Household Agent ---
class Household(Agent):
    def initialize(self):
        self.agent_type = "household"
        self.state = "unconnected"
        self.members = []
        #if random.random() < 0.4:
        self.members.append({"role":"man","employed":random.random()<0.7,"income":0,"exposure":0,"exposure_fraction":0.450,"age":random.randint(18, 59)})
        num_women = random.randint(1, 2)
        #if random.random() < 0.7:
        for _ in range(num_women):
            self.members.append({"role":"woman","employed":random.random()<0.25,"income":0,"exposure":0,"exposure_fraction":0.628,"age":random.randint(18, 59)})
        num_children = random.randint(0, 3)
        for _ in range(num_children):
            self.members.append({"role":"child","education_state":"","study_time":0,"schooling":0.0,"income":0,"exposure":0,"exposure_fraction":0.742,"gender":random.choice(["male","female"]),"age":random.randint(0, 17), "employed":"False"})
                
        self.fuel_appliances = {}  # name → owned (True/False)
        for name in FUEL_APPLIANCE_CATALOG:
            self.fuel_appliances[name] = True

        self.electric_appliances = {}  # name → owned (True/False)
        for name in APPLIANCE_CATALOG:
            self.electric_appliances[name] = False


    def electrify(self):
        self.state = "connected"

    def time_saved_men(self):
        if self.state == "connected":
            return 1   # hours/day saved
        return 0.0

    def time_saved_women(self):
        if self.state == "connected":
            return 2.5   # hours/day saved
        return 0.0
    

    def set_schooling(self):
        dropout_rate_baseline = 0.15
        for m in self.members:
            if m["age"] < 6:
                m["education_state"] = "not_started"
                m["schooling"] = 0
                continue

            potential_years = min(m["age"] - 6, 12)

            completed_years = 0

            year = 0.0
            while year < potential_years:
                if random.random() < dropout_rate_baseline:
                    break
                completed_years += 1
                year += 1

            m["schooling"] = completed_years

            if completed_years >= 12:
                m["education_state"] = "completed"
            elif completed_years < potential_years:
                m["education_state"] = "dropped_out"
            elif m["age"] < 18:
                m["education_state"] = "enrolled"
            else:
                m["education_state"] = "completed"

    def step(self,available_jobs):
        job_changes = 0.0
        def annual_to_monthly(percentage):
            return ((1+percentage)**(1/12))-1
            
        saved_hours_men = self.time_saved_men() #time saved from electrification that is dedicated to work
        saved_hours_women = self.time_saved_women() #time saved from electrification that is dedicated to work
        work_time_men = saved_hours_men * 0.4
        leisure_time_men = saved_hours_men * 0.5
        study_time_men = saved_hours_men * 0.1
        
        work_time_women = saved_hours_women * 0.4
        leisure_time_women = saved_hours_women * 0.4
        study_time_women = saved_hours_women * 0.2
        dropout_rate_baseline = 0.15
        baseline_monthly_wage_men = 130
        baseline_monthly_wage_women = 85
        baseline_schooling = 6.0
        error = 0.1
        employment_rate_baseline = annual_to_monthly(0.07)#0.07
        employment_rate_women_connected = annual_to_monthly(0.16)#0.16
        hourly_farm_wage = 2
        hourly_non_farm_wage = 7
        job_separation_rate = annual_to_monthly(0.02)#0.02
        electrification_effect_women = 0.9 #some studies suggest reduced female wages in areas of high electrification
        electrification_effect_men = 1.1 #elefrification is linked to increased wages for men
        working_days = 20 # number of working days in a monthly
        farm_work_shift = 0 #for farm work shift scenario. 1 for scenario, 0 otherwise
        non_farm_work_shift = 1 #for non-farm work shift scenario. 1 for scenario, 0 otherwise           

        for m in self.members:
            m["age"] += (1/12)
            if self.state == "connected":
                dropout_rate = dropout_rate_baseline*(1-0.25) #assumes electrification reduces dropout rate by 25%
            else:
                dropout_rate = dropout_rate_baseline

             # Start school
            if m["age"] >= 6 and m["education_state"] == "not_started":
                m["education_state"] = "enrolled"

            # If enrolled → accumulate schooling
            if m["education_state"] == "enrolled":
                # Dropout event
                if random.random() < dropout_rate and m["age"]%1==0: # dropout can only occur at the end of a year of schooling
                    m["education_state"] = "dropped_out"
                else:
                    m["schooling"] += (1/12)
                    # Completion condition
                    if m["schooling"] >= 12:
                        m["education_state"] = "completed"

            # Force exit at leaving age
            if m["age"] >= 18 and m["education_state"] == "enrolled":
                m["education_state"] = "completed"
                
       #     if m["role"]=="child" and m["employed"]==False and m["in_school"]:
       #         if random.random()<dropout_rate:
       #             m["in_school"] = False
       #         else:
       #             m["schooling"] += 1/12
               

         #employment decision 
        for m in self.members:
           
            if available_jobs <= 0:
                break
            employment_probability = 0.0
            if not m["employed"] and m["age"]<60:
                if m["role"] == "woman":
                    employment_probability = employment_rate_women_connected if self.state == "connected" else employment_rate_baseline
                elif m["role"] == "man":
                    employment_probability = employment_rate_baseline
                elif m["role"]=="child" and m["gender"] == "male" and m["age"]>=18:
                    employment_probability = employment_rate_baseline
                elif m["role"]=="child" and m["gender"] == "female" and m["age"]>=18:
                    employment_probability = employment_rate_women_connected if self.state == "connected" else employment_rate_baseline
                if random.random() < employment_probability: #probability of getting employed
                    m["employed"] = True
                    available_jobs -= 1
                    job_changes += 1


        for m in self.members:
            if m["employed"]:
                if m["role"] == "man":
                    if self.state == "connected":
                        m["income"] = ((baseline_monthly_wage_men*electrification_effect_men) + (work_time_men*working_days*hourly_non_farm_wage))
                    else:
                        m["income"] = baseline_monthly_wage_men
                elif m["role"] == "woman":
                    if self.state == "connected":
                        m["income"] = ((baseline_monthly_wage_women*electrification_effect_women) + ((work_time_women*working_days*hourly_farm_wage)*farm_work_shift) + (work_time_women*working_days*hourly_non_farm_wage)*non_farm_work_shift)
                    else:
                        m["income"] = baseline_monthly_wage_women
                elif m["role"]=="child" and m["gender"] == "male":
                    if self.state == "connected":
                        #if (m["schooling"]-baseline_schooling) > 0:
                        m["income"] = ((baseline_monthly_wage_men*electrification_effect_men) + (work_time_men*working_days*hourly_non_farm_wage))*math.exp((0.096*(m["schooling"]-baseline_schooling))+(0.0169*(m["age"]-m["schooling"]-6))+(-0.0003*(m["age"]-m["schooling"]-6)**2))#-baseline_schooling)# uses part of Mincer equation to view income effect of schooling
                        #else:
                         #   m["income"] = ((baseline_monthly_wage_men*electrification_effect_men) + (work_time_men*working_days*hourly_non_farm_wage))# uses part of Mincer equation to view income effect of schooling
                    else:
                        m["income"] = baseline_monthly_wage_men
                elif m["role"]=="child" and m["gender"] == "female":
                    if self.state == "connected":
                        #if (m["schooling"]-baseline_schooling) > 0:
                        m["income"] = (((baseline_monthly_wage_women*electrification_effect_women) + ((work_time_women*working_days*hourly_farm_wage)*farm_work_shift) + (work_time_women*working_days*hourly_non_farm_wage)*non_farm_work_shift))*math.exp((0.115*(m["schooling"]-baseline_schooling))+(0.0169*(m["age"]-m["schooling"]-6))+(-0.0003*(m["age"]-m["schooling"]-6)**2))#-baseline_schooling)
                        #else:
                        #    m["income"] = (((baseline_monthly_wage_women*electrification_effect_women) + ((work_time_women*working_days*hourly_farm_wage)*farm_work_shift) + (work_time_women*working_days*hourly_non_farm_wage)*non_farm_work_shift))
                    else:
                        m["income"] = baseline_monthly_wage_women
            if m["age"]>= 60:
                m["employed"] = False
                m["income"] = 0.0


        for m in self.members:
            if m["employed"] and random.random()<job_separation_rate:
                m["employed"] = False
                m["income"] = 0
                available_jobs += 1
                job_changes -= 1

        
        return job_changes          


    def total_income(self):
        return sum(m.get("income",0) for m in self.members)

    def total_study_time(self):
        return sum(m.get("study_time",0) for m in self.members if m["role"]=="child")
    
    def total_schooling(self):
        return sum(m.get("schooling",0) for m in self.members if m["role"]=="child")
    
    def total_children(self):
        return sum(1 for m in self.members if m["role"]=="child")
        

    def consider_appliance_adoption(self, reliability, social_influence,cost_per_kwh):
        income = self.total_income()

        for name, a in APPLIANCE_CATALOG.items():

            if self.electric_appliances[name]:
                continue  # already owned

            utility = (
                a["beta_0"]
                + a["beta_income"] * income
                + a["beta_reliability"] * reliability
                + a["beta_social"] * social_influence
                - a["beta_cost"] * a["cost"]
                - a["beta_tariff"] * cost_per_kwh
            )

            adoption_prob = 1 / (1 + math.exp(-utility))

            if random.random() < adoption_prob:
                if name == "electric_cooker":
                    self.electric_appliances[name] = True
                    self.fuel_appliances["lpg_stove"] = False
                elif name == "electric_lighting":
                    self.electric_appliances[name] = True
                    self.fuel_appliances["kerosene_lamp"] = False
                elif name == "tv":
                    self.electric_appliances[name] = True
        
        

    def fossil_fuel_use(self):
        total_fuel_use = 0.0
        for name, a in FUEL_APPLIANCE_CATALOG.items():

            if self.fuel_appliances[name]:
                fuel_use = a["hours_used"]*a["fuel_used_hour"] * 30
            else:
                fuel_use = 0.0
            total_fuel_use += fuel_use
        return total_fuel_use

    def energy_cost(self,cost_per_kwh):
        total_fuel_cost = 0.0
        for name, a in FUEL_APPLIANCE_CATALOG.items():
            if self.fuel_appliances[name]:
                fuel_cost = a["hours_used"]*a["fuel_used_hour"]*a["cost_per_unit"] * 30
            else:
                fuel_cost = 0.0
                
            total_fuel_cost += fuel_cost
            
        total_electric_appliance_cost = 0.0
        for name, a in APPLIANCE_CATALOG.items():
            if self.electric_appliances[name]:
                electric_appliance_cost = a["monthly_kwh"] * cost_per_kwh
            else:
                electric_appliance_cost = 0.0
                
            total_electric_appliance_cost += electric_appliance_cost
        total_energy_cost = total_fuel_cost + total_electric_appliance_cost
        return total_energy_cost

    def baseline_fuel_cost(self):
        baseline_total_fuel_cost = 0.0
        for name, a in FUEL_APPLIANCE_CATALOG.items():
            baseline_fuel_cost = a["hours_used"]*a["fuel_used_hour"]*a["cost_per_unit"] * 30
            baseline_total_fuel_cost += baseline_fuel_cost
        return baseline_total_fuel_cost
        
    def co2_emissions(self):
        total_co2_emissions = 0.0
        for name, a in FUEL_APPLIANCE_CATALOG.items():
            if self.fuel_appliances[name]:
                co2_emissions = a["hours_used"]*a["fuel_used_hour"]*a["fuel_emission_factor_co2"] * 30
            else:
                co2_emissions = 0.0
        
            total_co2_emissions += co2_emissions
        return total_co2_emissions

    def pm25_concentration(self):
        air_change_rate_daily = 24.0 #air changes per day in the kitchen
        kitchen_volume = 20.0 # m^3
        outdoor_concentration = 1e-6

        weighted_fraction = 0.0
        total_emissions = 0.0

        for name, a in FUEL_APPLIANCE_CATALOG.items():
            if self.fuel_appliances[name]:
                emissions = a["hours_used"] * a["fuel_used_hour"] * a["fuel_emission_factor_pm25"]
                total_emissions += emissions
                weighted_fraction += emissions * a["kitchen_fraction"]

        if total_emissions > 0:
            frac = weighted_fraction / total_emissions
        else:
            frac = 0.0

        pm25_concentration = (total_emissions * frac) / (air_change_rate_daily * kitchen_volume) + outdoor_concentration
        
        return pm25_concentration

    def pm25_concentration_baseline(self):
        air_change_rate_daily = 24.0 #air changes per day in the kitchen
        kitchen_volume = 20.0 # m^3
        outdoor_concentration = 1e-6

        weighted_fraction = 0.0
        total_emissions_baseline = 0.0
        for name, a in FUEL_APPLIANCE_CATALOG.items():
            #if self.fuel_appliances[name]:
            emissions = a["hours_used"] * a["fuel_used_hour"] * a["fuel_emission_factor_pm25"]
            total_emissions_baseline += emissions
            weighted_fraction += emissions * a["kitchen_fraction"]
        if total_emissions_baseline > 0:
            frac = weighted_fraction / total_emissions_baseline
        else:
            frac = 0.0
        pm25_concentration_baseline = ((total_emissions_baseline * frac) / (air_change_rate_daily * kitchen_volume)) + outdoor_concentration
        return pm25_concentration_baseline

    def baseline_cases(self):
        concentration_baseline = self.pm25_concentration_baseline()
        alpha = DISEASE_PARAMS["alpha"]
        gamma = DISEASE_PARAMS["gamma"]
        delta = DISEASE_PARAMS["delta"]
        c0 = DISEASE_PARAMS["c0"]
        baseline_incidence = DISEASE_PARAMS["baseline_incidence"] # number of cases per month
        total_cases_baseline = 0.0
        for member in self.members:
            baseline_exposure = member["exposure_fraction"]*concentration_baseline*1e6 # the individual exposure is a fraction of the concentration in the house. Scaled from g/m3 to ug/m3
            excess = max(0.0, baseline_exposure - c0)
            RR_baseline = 1 + alpha * (1 - math.exp(-gamma * (excess ** delta))) # Relative risk
            PAF_baseline = (RR_baseline - 1) / RR_baseline # Attributable fraction
            attributable_cases_baseline = PAF_baseline * baseline_incidence * len(self.members)
            total_cases_baseline += attributable_cases_baseline
        return total_cases_baseline
        
    def cases(self):
        concentration = self.pm25_concentration()
        alpha = DISEASE_PARAMS["alpha"]
        gamma = DISEASE_PARAMS["gamma"]
        delta = DISEASE_PARAMS["delta"]
        c0 = DISEASE_PARAMS["c0"]
        baseline_incidence = DISEASE_PARAMS["baseline_incidence"] # number of cases per month
        total_cases = 0.0
        for member in self.members:
            member["exposure"] = member["exposure_fraction"]*concentration*1e6 # the individual exposure is a fraction of the concentration in the house. Scaled from g/m3 to ug/m3
            excess = max(0.0, member["exposure"] - c0)
            RR = 1 + alpha * (1 - math.exp(-gamma * (excess ** delta)))# Relative risk
            PAF = (RR - 1) / RR # Attributable fraction
            attributable_cases = PAF * baseline_incidence * len(self.members)
            total_cases += attributable_cases

        return total_cases
         
                
    def appliance_demand(self):
        demand = 0.0
        for name, owned in self.electric_appliances.items():
            if owned:
                demand += APPLIANCE_CATALOG[name]["monthly_kwh"]
        return demand


    def social_welfare_function(self, reliability, social_influence,cost_per_kwh,eta):
        #alpha_income = 0.4;
        #alpha_reliability = 0.3;
        #alpha_cost = 0.2;
        #alpha_health = 0.1
        #eta = 0.95
        income = self.total_income()
        #baseline_cost = self.baseline_fuel_cost()
        #energy_cost = self.energy_cost(cost_per_kwh)
        #cases_baseline = self.baseline_cases()
        #cases = self.cases()
        #health_improvement = cases_baseline - cases
        #cost_savings = 1 - (energy_cost/baseline_cost)
        
        household_utility = (
                #math.log2(1 if income<=0 else income)
               # + alpha_reliability * reliability
                 #health_improvement
               # + alpha_cost * cost_savings
            (income**(1-eta))/(1-eta)
            )
        
        return household_utility

    def consider_microgrid_adoption(self, reliability, social_influence,cost_per_kwh):
        alpha_income = 0.16;
        alpha_reliability = 0.5;
        alpha_health = 0.1
        alpha_social = 0.2
        alpha_savings = 0.4
        alpha_tariff = 0.5
        aversion_to_adoption = -4.5
        income = self.total_income()
        baseline_cost = self.baseline_fuel_cost()
        energy_cost = self.energy_cost(cost_per_kwh)
        cases_baseline = self.baseline_cases()
        cases = self.cases()
        health_improvement = (cases_baseline - cases)/cases_baseline
        cost_savings = (baseline_cost - energy_cost)/baseline_cost

        utility = (
                aversion_to_adoption
                + alpha_income * math.log2(income if income>0 else 1)
                + alpha_reliability * reliability
                + alpha_social * social_influence
                - alpha_tariff * cost_per_kwh
                + alpha_health * health_improvement
                + alpha_savings * cost_savings
            )
        
        microgrid_adoption_prob = 1 / (1 + math.exp(-utility))
        
        if random.random() < microgrid_adoption_prob:
            self.electrify()
            


# --- System Dynamics Model ---
class ElectrificationSD:
    def __init__(self, model):
        self.model = model

        # Stocks
        self.households_not_connected = model.stock("households_not_connected")
        self.net_microgrid_income = model.stock("net_microgrid_income")
        self.no_of_failures = model.stock("no_of_failures")
        self.no_of_components = model.stock("no_of_components")
        self.no_of_current_repairs = model.stock("no_of_current_repairs")
        self.available_jobs = model.stock("available_jobs")
        self.discounted_social_welfare_function = model.stock("discounted_social_welfare_function")
        self.employment_change = model.flow("employment_change")#################
        
        # Flows
        self.household_microgrid_adoption = model.flow("household_microgrid_adoption")
        self.microgrid_income = model.flow("microgrid_income")
        self.microgrid_expenditures = model.flow("microgrid_expenditures")
        self.microgrid_income_flow = model.flow("microgrid_income_flow")
        self.failure_rate = model.flow("failure_rate")
        self.fixture_rate = model.flow("fixture_rate")
        self.job_creation = model.flow("job_creation")
        self.social_welfare_function = model.flow("social_welfare_function")
        self.carbon_cost = model.flow("carbon_cost")

        # Converters
        self.households_connected = model.converter("households_connected")
        self.adopting_households = model.converter("adopting_households")
        self.power_use = model.converter("power_use")
        self.failure_rate_multiplier =  model.converter("failure_rate_multiplier")
        self.mttr = model.converter("mttr")
        self.saifi = model.converter("saifi")
        self.saidi = model.converter("saidi")
        self.no_of_customers = model.converter("no_of_customers")
        self.downtime = model.converter("downtime")
        self.reliability = model.converter("reliability")
        self.social_influence_microgrid = model.converter("social_influence_microgrid")
        self.baseline_demand = model.converter("baseline_demand")
        self.total_demand = model.converter("total_demand")
        self.cost_savings = model.converter("cost_savings")
        self.total_energy_cost = model.converter("energy_cost")
        self.baseline_fuel_cost = model.converter("baseline_fuel_cost")
        self.co2_emissions = model.converter("co2_emissions")
        self.total_cases = model.converter("total_cases")        
        self.baseline_cases = model.converter("baseline_cases")
        self.health_improvements = model.converter("health_improvements")

       # self.human_capital = model.converter("human_capital")
        self.employment = model.converter("employment")
        self.business_income = model.converter("business_income")
        self.capital_investment = model.converter("capital_investment")
        
        self.tariff = model.converter("tariff")
        self.subsidy = model.converter("subsidy")

        # Aggregate social/economic converters
        self.avg_income = model.converter("avg_income")
        self.avg_study_time = model.converter("avg_study_time")
        self.appliance_demand = model.converter("appliance_demand")
        self.utility_sum = model.converter("utility_sum")
        self.changes_in_jobs = model.flow("changes_in_jobs")

        #Constants
        self.baseline_household_demand = model.constant("baseline_household_demand")
        self.no_of_households = model.constant("no_of_households")
        self.microgrid_capacity = model.constant("microgrid_capacity")
        self.initial_failure_rate = model.constant("initial_failure_rate")
        self.cost_per_kwh = model.constant("cost_per_kwh")
        self.cost_per_repair = model.constant("cost_per_repair")
        self.operating_expenditures = model.constant("operating_expenditures")
        self.no_of_components_initial_value = model.constant("no_of_components_initial_value")
        self.attrition_rate = model.constant("attrition_rate")
        self.mpc = model.constant("mpc")
        self.local_spending_fraction = model.constant("local_spending_fraction")
        self.investment_rate = model.constant("investment_rate")
        self.job_creation_efficiency = model.constant("job_creation_efficiency")
        self.initial_jobs = model.constant("initial_jobs")
        self.social_cost_of_carbon = model.constant("social_cost_of_carbon")
        self.eta = model.constant("eta")

        # Equations
        self.social_influence_microgrid.equation = self.households_connected/self.no_of_households
        #self.household_microgrid_adoption.equation = self.households_not_connected * ((1e-2 * self.reliability) + (0.1e-2*self.social_influence_microgrid) + (3e-3*self.cost_savings)+(5e-2*self.health_improvements))#
        self.household_microgrid_adoption.equation = self.households_not_connected - self.households_connected
        #self.households_connected.equation = self.household_microgrid_adoption
        self.households_not_connected.equation = -self.household_microgrid_adoption
        self.baseline_demand.equation = ((self.baseline_household_demand) * self.households_connected)
        self.total_demand.equation = self.baseline_demand+self.appliance_demand
        self.power_use.equation = sd.min(self.total_demand,self.microgrid_capacity)
        self.microgrid_income_flow.equation = self.power_use * self.cost_per_kwh
        self.microgrid_expenditures.equation = ((self.failure_rate - self.fixture_rate)*self.cost_per_repair) + self.operating_expenditures
        self.net_microgrid_income.equation = self.microgrid_income_flow - self.microgrid_expenditures
        self.mttr.equation = sd.lookup(
            self.net_microgrid_income,
            [(0, 48), (30000000, 42), (70000000, 36), (100000000, 30)]
        )

        self.failure_rate_multiplier.equation = sd.lookup(
            self.power_use,
            [(90000, 1), (100000, 3), (115000, 8), (130000, 12)]
        )
        self.failure_rate.equation = self.initial_failure_rate*self.no_of_components*self.failure_rate_multiplier
        self.no_of_failures.equation = self.failure_rate
        self.fixture_rate.equation = self.no_of_current_repairs/self.mttr
        #number of components stock
        self.no_of_components.equation = self.fixture_rate - self.failure_rate
        #number of components being repaired
        self.no_of_current_repairs.equation = self.failure_rate - self.fixture_rate
        # Failure rate depends on maintenance delay
        self.downtime.equation = self.mttr * self.no_of_failures
        # SAIDI (interruptions per customer)
        self.saidi.equation = self.downtime / sd.max(1, self.households_connected)
        # SAIFI (interruptions per customer)
        self.saifi.equation = self.no_of_failures / sd.max(1, self.households_connected)
        # Reliability improves as SAIDI decreases
        self.reliability.equation = 1/(1+5.83*self.saidi)

        self.cost_savings.equation = 1 - (self.total_energy_cost/self.baseline_fuel_cost)
        self.health_improvements.equation = 1 - (self.total_cases/self.baseline_cases)

        # R4: Income → spending → investment → jobs
        self.business_income.equation = self.mpc * self.local_spending_fraction * self.avg_income
        self.capital_investment.equation = self.investment_rate * self.business_income
        self.job_creation.equation = self.job_creation_efficiency * self.capital_investment
        self.employment_change.equation = self.changes_in_jobs
        self.available_jobs.equation = self.job_creation + self.employment_change

        self.social_welfare_function.equation = self.utility_sum*sd.exp(-4.074e-3*sd.time())
        self.carbon_cost.equation = self.co2_emissions*self.social_cost_of_carbon*0.00110231*sd.exp(-4.074e-3*sd.time())
        self.discounted_social_welfare_function.equation = (self.social_welfare_function - self.carbon_cost)#*sd.exp(-4.074e-3*sd.time())
               

        # Initial values
        self.no_of_households.equation = 1000.0
        self.households_connected.initial_value = 0.0
        self.households_not_connected.initial_value = self.no_of_households
        self.microgrid_capacity.equation = 1000000.0 #KW
        self.cost_per_kwh.equation = 0.10
        self.cost_per_repair.equation = 50000.0
        self.operating_expenditures.equation = 100000.0
        self.initial_failure_rate.equation = 0.1
        self.no_of_failures.initial_value = 0.0
        self.no_of_components_initial_value.equation = 2.0
        self.no_of_components.initial_value = self.no_of_components_initial_value
        self.no_of_current_repairs.initial_value = 0.0
        self.baseline_household_demand.equation = 200.0 #in kWh
        self.available_jobs.initial_value = self.initial_jobs
        self.initial_jobs.equation = 75
        self.mpc.equation = 0.75
        self.local_spending_fraction.equation = 0.6
        self.investment_rate.equation = 0.25
        self.job_creation_efficiency.equation = 4.074e-3
        self.social_cost_of_carbon.equation = 130
        #self.eta.equation = 0.95
        #self.attrition_rate.equation = 0.02


# --- Hybrid Model ---
class ElectrificationHybrid(Model):
    def instantiate_model(self):
        super().instantiate_model()
        self.register_agent_factory("household", lambda agent_id, model, properties: Household(agent_id, model, properties))
        self.sd_model = ElectrificationSD(self)
        self.custom_stats = {}

    def configure(self, config):
        super().configure(config)
        #self.sd_model.reliability.equation = self.reliability
        self.sd_model.tariff.equation = self.tariff
        self.sd_model.subsidy.equation = self.subsidy

        self.sd_model.no_of_households.equation = self.no_of_households
        self.sd_model.microgrid_capacity.equation = self.microgrid_capacity
        self.sd_model.cost_per_kwh.equation = self.cost_per_kwh
        self.sd_model.initial_failure_rate.equation = self.initial_failure_rate
        self.sd_model.no_of_components_initial_value.equation = self.no_of_components_initial_value
        self.sd_model.baseline_household_demand.equation = self.baseline_household_demand
        self.sd_model.cost_per_repair.equation = self.cost_per_repair
        self.sd_model.operating_expenditures.equation = self.operating_expenditures
        self.sd_model.initial_jobs.equation = self.initial_jobs
        self.sd_model.social_cost_of_carbon.equation = self.social_cost_of_carbon
        self.sd_model.eta.equation = self.eta
        


    def begin_round(self, time, sim_round, step):
        # Update agent states
        # begin_round
        available_jobs = int(self.sd_model.available_jobs(time))
        cost_per_kwh = self.cost_per_kwh
        eta = self.eta
        reliability = self.sd_model.reliability(time)
        social_influence = self.sd_model.social_influence_microgrid(time)
        total_job_changes = 0.0
        if time == self.starttime:
            for h in self.agents:
                h.set_schooling()
            baseline_cases = sum(h.baseline_cases() for h in self.agents)
            self.sd_model.baseline_cases.equation = baseline_cases
        for h in self.agents:
            if h.state != "connected":
                h.consider_microgrid_adoption(reliability, social_influence,cost_per_kwh)
            job_changes = h.step(available_jobs)
            available_jobs += job_changes
            total_job_changes += job_changes
            if h.state == "connected":
                h.consider_appliance_adoption(reliability, social_influence,cost_per_kwh)
                

        # Aggregate outcomes from agents
        connected_count = sum(1 for h in self.agents if h.state == "connected")
        avg_income = sum(h.total_income() for h in self.agents) / len(self.agents)
        avg_study = sum(h.total_study_time() for h in self.agents) / len(self.agents)
        appliance_demand = sum(h.appliance_demand() for h in self.agents)
        baseline_fuel_cost = sum(h.baseline_fuel_cost () for h in self.agents)
        social_welfare_function = sum(h.social_welfare_function (reliability, social_influence,cost_per_kwh,eta) for h in self.agents)
        
        
        
        total_energy_cost = sum(h.energy_cost(cost_per_kwh) for h in self.agents)
        baseline_fuel_cost = sum(h.baseline_fuel_cost() for h in self.agents)
        co2_emissions = sum(h.co2_emissions() for h in self.agents)
        total_cases = sum(h.cases() for h in self.agents)
        average_schooling = sum(h.total_schooling() for h in self.agents)/sum(h.total_children() for h in self.agents) if sum(h.total_children() for h in self.agents) > 0 else 0
        
        # Feed back into SD converters
        self.sd_model.households_connected.equation = connected_count
        self.sd_model.avg_income.equation = avg_income
        self.sd_model.avg_study_time.equation = avg_study
        self.sd_model.appliance_demand.equation = appliance_demand
        self.sd_model.total_energy_cost.equation = total_energy_cost
        self.sd_model.baseline_fuel_cost.equation = baseline_fuel_cost
        self.sd_model.total_cases.equation = total_cases
        self.sd_model.utility_sum.equation = social_welfare_function
        self.sd_model.co2_emissions.equation = co2_emissions
        self.sd_model.changes_in_jobs.equation = total_job_changes
        employment = sum(
        1 for h in self.agents
        for m in h.members
        if m["employed"]
        )

        self.sd_model.employment.equation = employment


        if time not in self.custom_stats:
            self.custom_stats[time] = {}

        self.custom_stats[time]["connected_count"] = connected_count
        self.custom_stats[time]["appliance_demand"] = appliance_demand
        self.custom_stats[time]["total_cases"] = total_cases
        self.custom_stats[time]["avg_income"] = avg_income
        self.custom_stats[time]["saifi"] = self.sd_model.saifi(time)
        self.custom_stats[time]["Employment"] = employment
        self.custom_stats[time]["Jobs Available"] = available_jobs
        self.custom_stats[time]["Demand"] = self.sd_model.total_demand(time)
        self.custom_stats[time]["Income"] = self.sd_model.net_microgrid_income(time)
        self.custom_stats[time]["SWF"] = social_welfare_function
        self.custom_stats[time]["DSWF"] = self.sd_model.discounted_social_welfare_function(time)
        self.custom_stats[time]["co2_emissions"] = co2_emissions
        self.custom_stats[time]["average_schooling"] = average_schooling
        

        
# --- Simulation Setup ---
def run_model(tariff_value,no_of_households_value,stoptime):
    electrification_hybrid = ElectrificationHybrid(
        1, stoptime, dt=1, name="Electrification Hybrid",
        scheduler=SimultaneousScheduler(),
        #data_collector=DataCollector()

    )

    electrification_hybrid.instantiate_model()

    electrification_config = {
        "runspecs": {
            "starttime": 1,
            "stoptime": stoptime,
            "dt": 1.0
        },
        "properties": {
            "tariff": {"type": "Double", "value": tariff_value},
            "subsidy": {"type": "Double", "value": 0.0},
            "no_of_households": {"type": "Double", "value": no_of_households_value},
            "microgrid_capacity": {"type": "Double", "value": 1000000.0},#KW
            "cost_per_kwh": {"type": "Double", "value": tariff_value},
            "initial_failure_rate": {"type": "Double", "value": 0.0147},
            "no_of_components_initial_value": {"type": "Double", "value": 2.0},
            "baseline_household_demand": {"type": "Double", "value": 200.0},#KW
            "cost_per_repair": {"type": "Double", "value": 4000.0},
            "operating_expenditures": {"type": "Double", "value": 3000.0},
            "initial_jobs": {"type": "Double", "value": 75.0},
            "social_cost_of_carbon":{"type": "Double", "value": 130.0},
            "eta":{"type": "Double", "value": 0.4}
        },
        "agents": [
            {"name": "household", "count": no_of_households_value}
        ]
    }
   
  
    electrification_hybrid.configure(electrification_config)
    
    #custom = electrification_hybrid.custom_stats
    electrification_hybrid.run()
    return electrification_hybrid.custom_stats.copy()
from BPTK_Py import Model
from BPTK_Py import sd_functions as sd
from BPTK_Py import Agent
from BPTK_Py import Event
from BPTK_Py import DataCollector
from BPTK_Py import SimultaneousScheduler
import random
import math

beta_mincer_schooling_men = 0.096# The rate of return to an additional year of schooling for men, based on Mincerian wage regression estimates
beta_mincer_schooling_women = 0.115# The rate of return to an additional year of schooling for women, based on Mincerian wage regression estimates

DISEASE_PARAMS = {#Parameters for the integrated exposure-response function, which models the relationship between PM2.5 exposure and the incidence of ischemic heart disease (IHD).
    "alpha": 1.217816, #The maximum relative risk increase for IHD associated with PM2.5 exposure, representing the upper limit of the risk curve as exposure increases.
    "gamma":0.028152, #The shape parameter for the exposure-response function, determining the steepness of the curve.
    "delta":0.786667 , #The scale parameter for the exposure-response function, affecting the location of the curve.
    "c0":6.533716 , #The counterfactual concentration of PM2.5 below which no excess risk of IHD is assumed, representing a threshold for health effects.
    "baseline_incidence":200/100000 #The baseline incidence rate of IHD in the population.
}

APPLIANCE_CATALOG = {#A catalog of electric appliances that households may adopt, including their costs, energy consumption, and parameters for the adoption utility function.
    "electric_cooker": {
        "cost": 500,
        "monthly_kwh": 30,
        "beta_0": -4.5,
        "beta_income": 0.0002,
        "beta_reliability": 0.5,
        "beta_social": 1.0,
        "beta_cost": 0.0025,
        "beta_tariff": 3.5
    },
    "tv": {
        "cost": 300,
        "monthly_kwh": 7.5,
        "beta_0": -2.0,
        "beta_income": 0.0002,
        "beta_reliability": 0.5,
        "beta_social": 1.0,
        "beta_cost": 0.003,
        "beta_tariff": 2.5
    },
    "electric_lighting": {
        "cost": 50,
        "monthly_kwh": 4,
        "beta_0": -1.8,
        "beta_income": 0.0004,
        "beta_reliability": 0.5,
        "beta_social": 1.0,
        "beta_cost": 0.02,
        "beta_tariff": 2.0
    }
}

FUEL_APPLIANCE_CATALOG = {#A catalog of fuel-based appliances that households may use, including their usage patterns, costs, and emission factors for PM2.5 and CO2.
    "lpg_stove": {
        "hours_used": 1.0, # hours per day
        "fuel_used_hour": 0.2, # kilograms per hour
        "cost_per_unit": 1.568, # cost per kilogram
        "fuel_emission_factor_pm25": 0.45, #g/kg
        "fuel_emission_factor_co2": 2.98,
        "kitchen_fraction": 0.25
    },
    "kerosene_lamp": {
        "hours_used": 4.0,
        "fuel_used_hour": 0.06,
        "cost_per_unit": 1.05, 
        "fuel_emission_factor_pm25": 10, #g/L
        "fuel_emission_factor_co2": 2.68,
        "kitchen_fraction": 0.3
    }
}

# --- Household Agent ---
class Household(Agent): #The Household class represents a household agent in the model, encapsulating its characteristics, behaviours, and interactions with the environment. It includes methods for initializing household members, making decisions about electrification and appliance adoption, calculating energy use and emissions, and assessing health impacts and social welfare.
    def initialize(self): #Initializes the household agent by setting up its members, appliances, and other attributes.
        self.agent_type = "household"
        self.state = "unconnected" #The initial state of the household is set to "unconnected".
        self.members = []
        self.members.append({"role":"man","employed":random.random()<0.8,"income":0,"exposure":0,"exposure_fraction":0.450,"age":random.randint(18, 59),"gender":"male"})
        num_women = random.randint(1, 2)
        for _ in range(num_women):
            self.members.append({"role":"woman","employed":random.random()<0.4,"income":0,"exposure":0,"exposure_fraction":0.628,"age":random.randint(18, 59),"gender":"female"})
        num_children = random.randint(0, 3)
        for _ in range(num_children):
            self.members.append({"role":"child","education_state":"","schooling":0.0,"income":0,"exposure":0,"exposure_fraction":0.742,"gender":random.choice(["male","female"]),"age":random.randint(0, 17), "employed":"False"})
                
        self.fuel_appliances = {} #A dictionary to track the ownership of fuel-based appliances in the household, initialized based on the FUEL_APPLIANCE_CATALOG.
        for name in FUEL_APPLIANCE_CATALOG:
            self.fuel_appliances[name] = True

        self.electric_appliances = {}  #A dictionary to track the ownership of electric appliances in the household, initialized to False for all appliances in the APPLIANCE_CATALOG since the household starts as unconnected.
        for name in APPLIANCE_CATALOG:
            self.electric_appliances[name] = False

    def birth_process(self, baseline_fertility_rate, fertility_reduction_rate):


        # Convert annual fertility to monthly
        monthly_fertility = 1 - (1 - baseline_fertility_rate) ** (1 / 12)

        for m in self.members:
            if m["gender"] == "female" and 18 <= m["age"] <= 45:
                fertility = monthly_fertility

                # Electrification reduces fertility
                if self.state == "connected":
                    fertility *= (1 - fertility_reduction_rate)

                if random.random() < fertility:
                    self.members.append({"role":"child","education_state":"","schooling":0.0,"income":0,"exposure":0,"exposure_fraction":0.742,"gender":random.choice(["male","female"]),"age":random.randint(0, 17), "employed":"False"})

    def death_process(self,baseline_mortality_rate,child_multiplier,elderly_multiplier):
    
        survivors = []

        for m in self.members:
            mortality = baseline_mortality_rate

            if m["age"] < 5:
                mortality *= child_multiplier
            elif m["age"] >= 60:
                mortality *= elderly_multiplier

            # Convert annual mortality to monthly
            monthly_mortality = 1 - (1 - mortality) ** (1 / 12)

            if random.random() > monthly_mortality:
                survivors.append(m)

        self.members = survivors

    def electrify(self):#Changes the state of the household to "connected".
        self.state = "connected"

    def time_saved_men(self): # Hours per day saved for men
        if self.state == "connected":
            return 1   
        return 0.0

    def time_saved_women(self): # Hours per day saved for women
        if self.state == "connected":
            return 2.5   
        return 0.0
    

    def set_schooling(self,primary_completion,lower_secondary_completion,upper_secondary_completion):# Initializes the schooling levels of child members based on their age and the probabilities of completing different education levels.
    

      

        for m in self.members:

            if m["role"] != "child":#Sets schooling for children only
                continue

            if m["age"] < 6:
            # too young for school
                m["education_state"] = "not_started"
                m["schooling"] = 0
                continue

            if m["age"] < 12:
            # still in primary
                m["education_state"] = "enrolled"
                m["schooling"] = m["age"] - 6
                continue

        # age >= 12
            if random.random() < primary_completion:
                years = 6
            else:
                years = random.randint(0, 5)
                m["education_state"] = "dropped_out"

            if m["age"] >= 16:
                if random.random() < lower_secondary_completion:
                    years = 10
                else:
                    years = random.randint(6, 9)
                    m["education_state"] = "dropped_out"

            if m["age"] >= 18:
                if random.random() < upper_secondary_completion:
                    years = 12
                    m["education_state"] = "completed"
                else:
                    years = random.randint(10, 11)
                    m["education_state"] = "dropped_out"

            m["schooling"] = years

    def income_and_employment(self,available_jobs,employment_rate_baseline,hourly_farm_wage,hourly_non_farm_wage,job_separation_rate, dropout_rate_baseline,baseline_monthly_wage_men,baseline_monthly_wage_women,baseline_schooling,electrification_effect_men,electrification_effect_women,farm_work_shift,dropout_reduction_rate):
        job_changes = 0.0#A counter to track the net change in employment status among household members during this step, accounting for both new job acquisitions and job separations.
        def annual_to_monthly(percentage): #Converts an annual percentage rate to a monthly percentage rate.
            return ((1+percentage)**(1/12))-1
            
        saved_hours_men = self.time_saved_men() #time saved from electrification that is dedicated to work
        saved_hours_women = self.time_saved_women() #time saved from electrification that is dedicated to work
        work_time_men = saved_hours_men * 0.4
        work_time_women = saved_hours_women * 0.4
        error_term = random.uniform(-0.0008, 0.0008) # random error term to reflect unobserved factors affecting employment
        employment_rate_baseline = annual_to_monthly(employment_rate_baseline)
        employment_rate_women_connected = annual_to_monthly(employment_rate_baseline+0.09)#Assumes a 9 percentage point increase in employment rate for women due to electrification
       
        job_separation = annual_to_monthly(job_separation_rate)
        
        working_days = 20 # number of working days in a month        

        for m in self.members:
            m["age"] += (1/12)
            if self.total_income() > 350.0: #Assumes that if the household income exceeds $350, the dropout rate for children decreases by a certain percentage due to improved financial stability.
                dropout_rate_baseline *= 0.8
            if self.state == "connected":
                dropout_rate = annual_to_monthly(dropout_rate_baseline*(1-dropout_reduction_rate)) #assumes electrification reduces dropout rate
            else:
                dropout_rate = annual_to_monthly(dropout_rate_baseline)
            if m["role"]=="child":
                # Start school
                if m["age"] >= 6 and m["education_state"] == "not_started":
                    m["education_state"] = "enrolled"

                 # If enrolled, accumulate schooling
                if m["education_state"] == "enrolled":
                    # Dropout event
                    if random.random() < dropout_rate:
                        m["education_state"] = "dropped_out"
                    else:
                        m["schooling"] += (1/12)
                        # Completion condition
                        if m["schooling"] >= 12:
                            m["education_state"] = "completed"
                    # Force exit at leaving age
                    if m["age"] >= 18 and m["education_state"] == "enrolled":
                        m["education_state"] = "completed"
                

         #employment decision 
        for m in self.members:
           
            if available_jobs <= 0:
                break #No more jobs available in the market, so no one else can get employed
            employment_probability = 0.0
            if not m["employed"] and m["age"]<60:
                if m["gender"] == "female" and m["age"] >= 18:
                    employment_probability = employment_rate_women_connected if self.state == "connected" else employment_rate_baseline #Assumes that electrification increases the employment probability for women
                elif m["gender"] == "male" and m["age"] >= 18:
                    employment_probability = employment_rate_baseline
                
                if random.random() < employment_probability: # Probability of getting employed
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
                        if random.random() < farm_work_shift:
                            m["income"] = ((baseline_monthly_wage_women*electrification_effect_women) + (work_time_women*working_days*hourly_farm_wage))
                        else:
                            m["income"] = ((baseline_monthly_wage_women*electrification_effect_women) + (work_time_women*working_days*hourly_non_farm_wage))
                    else:
                        m["income"] = baseline_monthly_wage_women
                elif m["role"]=="child" and m["gender"] == "male":
                    if self.state == "connected":
                        m["income"] = ((baseline_monthly_wage_men*electrification_effect_men) + (work_time_men*working_days*hourly_non_farm_wage))*math.exp((beta_mincer_schooling_men*(m["schooling"]-baseline_schooling))+(0.0219*(m["age"]-m["schooling"]-6))+(-0.0003*(m["age"]-m["schooling"]-6)**2)+error_term)#-baseline_schooling)# uses part of Mincer equation to view income effect of schooling

                    else:
                        m["income"] = baseline_monthly_wage_men
                elif m["role"]=="child" and m["gender"] == "female":
                    if self.state == "connected":
                        if random.random()< farm_work_shift:
                            m["income"] = ((baseline_monthly_wage_women*electrification_effect_women) + (work_time_women*working_days*hourly_farm_wage))*math.exp((beta_mincer_schooling_women*(m["schooling"]-baseline_schooling))+(0.0219*(m["age"]-m["schooling"]-6))+(-0.0003*(m["age"]-m["schooling"]-6)**2)+error_term)
                        else:
                            m["income"] = ((baseline_monthly_wage_women*electrification_effect_women) + (work_time_women*working_days*hourly_non_farm_wage))*math.exp((beta_mincer_schooling_women*(m["schooling"]-baseline_schooling))+(0.0219*(m["age"]-m["schooling"]-6))+(-0.0003*(m["age"]-m["schooling"]-6)**2)+error_term)
                    else:
                        m["income"] = baseline_monthly_wage_women
            if m["age"]>= 60: #Assumes that everyone retires at age 60
                m["employed"] = False
                m["income"] = 0.0

        # Job separation process
        for m in self.members:
            if m["employed"] and random.random()<job_separation:
                m["employed"] = False
                m["income"] = 0
                available_jobs += 1
                job_changes -= 1

        
        return job_changes          


    def total_income(self): #Calculates the total monthly income of the household by summing the incomes of all members.
        return sum(m.get("income",0) for m in self.members)
    
    def income_women(self): #Calculates the total monthly income of the household from female members.
        return sum(m.get("income",0) for m in self.members if m["gender"] == "female")

    def income_men(self): #Calculates the total monthly income of the household from male members.
        return sum(m.get("income",0) for m in self.members if m["gender"] == "male")

    def initial_income(self): #Calculates the initial total monthly income of the household.
        return sum(m.get("income",0) for m in self.members)
    
    def total_women(self): #Calculates the total number of adult female members in the household.
        return sum(1 for m in self.members if m["gender"] == "female" and m["age"] >= 18)
    
    def total_men(self): #Calculates the total number of adult male members in the household.
        return sum(1 for m in self.members if m["gender"] == "male" and m["age"] >= 18)
    
    def total_schooling(self): #Calculates the total years of schooling for all child members over the age of 18.
        return sum(m.get("schooling",0) for m in self.members if m["role"]=="child" and m["age"] >= 18)
    
    def total_adult_children(self): #Calculates the total number of child members in the household who are over the age of 18.
        return sum(1 for m in self.members if m["role"]=="child" and m["age"] >= 18)

    def consider_appliance_adoption(self, reliability, social_influence,cost_per_kwh,subsidy_rate):
        income = self.total_income()
        k = 1.5  # sensitivity parameter for adoption probability
        for name, a in APPLIANCE_CATALOG.items():

            if self.electric_appliances[name]:
                continue  # already owned
            if name != "tv": #Assumes that the electric cooker and electric lighting are subsidized, while the TV is not subsidized.
                appliance_utility = (
                     a["beta_income"] * income
                    + a["beta_reliability"] * reliability
                    + a["beta_social"] * social_influence
                    - a["beta_cost"] * (a["cost"]*(1-subsidy_rate)) # upfront cost of the appliance, adjusted for subsidy
                    - a["beta_tariff"] * cost_per_kwh
                )
            else:  # assumes TV is not subsidized
                appliance_utility = ( 
                     a["beta_income"] * income
                    + a["beta_reliability"] * reliability
                    + a["beta_social"] * social_influence
                    - a["beta_cost"] * (a["cost"])
                    - a["beta_tariff"] * cost_per_kwh
                )
            adoption_prob = 1 / (1 + math.exp(-(a["beta_0"] + k*appliance_utility))) #Calculates the probability of adopting the appliance using a logistic function based on the calculated utility.

            if random.random() < adoption_prob:
                if name == "electric_cooker":
                    self.electric_appliances[name] = True
                    self.fuel_appliances["lpg_stove"] = False #Assumes that adopting an electric cooker leads to discontinuing the use of an LPG stove.
                elif name == "electric_lighting":
                    self.electric_appliances[name] = True
                    self.fuel_appliances["kerosene_lamp"] = False #Assumes that adopting electric lighting leads to discontinuing the use of a kerosene lamp.
                elif name == "tv":
                    self.electric_appliances[name] = True
        
        

    def fossil_fuel_use(self): #Calculates the total monthly fuel use of the household based on the usage of fuel-based appliances, accounting for the hours used, fuel consumption rates, and whether the appliances are still in use.
        total_fuel_use = 0.0
        for name, a in FUEL_APPLIANCE_CATALOG.items():

            if self.fuel_appliances[name]: # If the household is still using the fuel appliance, calculate the monthly fuel use.
                fuel_use = a["hours_used"]*a["fuel_used_hour"] * 30 # Assumes that the usage pattern (hours used per day and fuel consumption rate) applies to every day of the month, multiplying by 30 to get monthly fuel use.
            else:
                fuel_use = 0.0
            total_fuel_use += fuel_use
        return total_fuel_use

    def energy_cost(self,cost_per_kwh): #Calculates the total monthly energy cost for the household based on the usage of electric and fuel-based appliances.
        total_fuel_cost = 0.0
        for name, a in FUEL_APPLIANCE_CATALOG.items():
            if self.fuel_appliances[name]:
                fuel_cost = a["hours_used"]*a["fuel_used_hour"]*a["cost_per_unit"] * 30 # Assumes that the usage pattern (hours used per day and fuel consumption rate) applies to every day of the month, multiplying by 30 to get monthly fuel cost.
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

    def baseline_fuel_cost(self): #Calculates the baseline monthly fuel cost for the household, assuming that all fuel-based appliances are in use and no electric appliances have been adopted.
        baseline_total_fuel_cost = 0.0
        for name, a in FUEL_APPLIANCE_CATALOG.items():
            baseline_fuel_cost = a["hours_used"]*a["fuel_used_hour"]*a["cost_per_unit"] * 30
            baseline_total_fuel_cost += baseline_fuel_cost
        return baseline_total_fuel_cost
        
    def co2_emissions(self):#Calculates the total monthly CO2 emissions for the household based on the usage of fuel-based appliances.
        total_co2_emissions = 0.0
        for name, a in FUEL_APPLIANCE_CATALOG.items():
            if self.fuel_appliances[name]:# If the household is still using the fuel appliance, calculate the monthly CO2 emissions.
                co2_emissions = a["hours_used"]*a["fuel_used_hour"]*a["fuel_emission_factor_co2"] * 30
            else:
                co2_emissions = 0.0
        
            total_co2_emissions += co2_emissions
        return total_co2_emissions
    
    def co2_emissions_initial(self):#Calculates the initial total monthly CO2 emissions for the household, assuming that all fuel-based appliances are in use and no electric appliances have been adopted.
        total_co2_emissions = 0.0
        for name, a in FUEL_APPLIANCE_CATALOG.items():
            co2_emissions = a["hours_used"]*a["fuel_used_hour"]*a["fuel_emission_factor_co2"] * 30
            total_co2_emissions += co2_emissions
        return total_co2_emissions

    def pm25_concentration(self,air_change_rate_daily,kitchen_volume,outdoor_concentration):#Calculates the indoor PM2.5 concentration in the household.

        

        weighted_fraction = 0.0 # A variable to accumulate the weighted fraction of emissions that contribute to indoor PM2.5 concentration, accounting for the kitchen fraction for each fuel appliance.
        total_emissions = 0.0

        for name, a in FUEL_APPLIANCE_CATALOG.items():
            if self.fuel_appliances[name]:
                emissions = a["hours_used"] * a["fuel_used_hour"] * a["fuel_emission_factor_pm25"]
                total_emissions += emissions
                weighted_fraction += emissions * a["kitchen_fraction"]

        if total_emissions > 0:
            frac = weighted_fraction / total_emissions #Calculates the fraction of total emissions that contribute to indoor PM2.5 concentration.
        else:
            frac = 0.0

        pm25_concentration = (total_emissions * frac) / (air_change_rate_daily * kitchen_volume) + outdoor_concentration
        
        return pm25_concentration

    def pm25_concentration_baseline(self,air_change_rate_daily,kitchen_volume,outdoor_concentration):#Calculates the baseline indoor PM2.5 concentration in the household, assuming that all fuel-based appliances are in use and no electric appliances have been adopted. This serves as a reference point for assessing the health impacts of electrification and appliance adoption.

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

    def baseline_cases(self,air_change_rate_daily,kitchen_volume,outdoor_concentration): #Calculates the baseline number of ischemic heart disease (IHD) cases attributable to PM2.5 exposure in the household.
        concentration_baseline = self.pm25_concentration_baseline(air_change_rate_daily,kitchen_volume,outdoor_concentration)
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
        
    def cases(self,air_change_rate_daily,kitchen_volume,outdoor_concentration):#Calculates the number of ischemic heart disease (IHD) cases attributable to PM2.5 exposure in the household under the current conditions.
        concentration = self.pm25_concentration(air_change_rate_daily,kitchen_volume,outdoor_concentration)
        alpha = DISEASE_PARAMS["alpha"]
        gamma = DISEASE_PARAMS["gamma"]
        delta = DISEASE_PARAMS["delta"]
        c0 = DISEASE_PARAMS["c0"]
        baseline_incidence = DISEASE_PARAMS["baseline_incidence"] # number of cases per month
        total_cases = 0.0
        for member in self.members:
            member["exposure"] = member["exposure_fraction"]*concentration*1e6 # the individual exposure is a fraction of the concentration in the house. Scaled from g/m3 to ug/m3
            excess = max(0.0, member["exposure"] - c0)#Calculates the excess exposure above the counterfactual concentration (c0) for the member, which is used to determine the relative risk of IHD associated with their PM2.5 exposure.
            RR = 1 + alpha * (1 - math.exp(-gamma * (excess ** delta)))# Relative risk
            PAF = (RR - 1) / RR # Population Attributable fraction
            attributable_cases = PAF * baseline_incidence * len(self.members)
            total_cases += attributable_cases

        return total_cases
         
                
    def appliance_demand(self): #Calculates the total monthly electricity demand of the household based on the electric appliances that have been adopted, summing their monthly kWh consumption.
        demand = 0.0
        for name, owned in self.electric_appliances.items():
            if owned:
                demand += APPLIANCE_CATALOG[name]["monthly_kwh"]
        return demand


    def impact_score(self,air_change_rate_daily, kitchen_volume, outdoor_concentration,baseline_schooling,average_income_men): #Calculates an overall impact score for the household based on improvements in income, health, CO2 emissions, and cost savings compared to the baseline scenario.
        cases_baseline = self.baseline_cases(air_change_rate_daily, kitchen_volume, outdoor_concentration)
        cases = self.cases(air_change_rate_daily, kitchen_volume, outdoor_concentration)
        health_improvement = (cases_baseline - cases)/cases_baseline if cases_baseline > 0 else 0 #Calculates the improvement in health outcomes as the percentage reduction in IHD cases attributable to PM2.5 exposure compared to the baseline scenario.
        initial_co2_emissions = self.co2_emissions_initial()
        current_co2_emissions = self.co2_emissions()
        co2_reduction = (initial_co2_emissions - current_co2_emissions)/initial_co2_emissions if initial_co2_emissions > 0 else 0
        average_schooling = self.total_schooling() / self.total_adult_children() if self.total_adult_children() > 0 else 0
        education_improvement = (average_schooling - baseline_schooling)/(12-baseline_schooling) if baseline_schooling > 0 else 0
        average_household_income_women = self.income_women()/self.total_women() if self.total_women() > 0 else 0
        average_household_income_men = self.income_men()/self.total_men() if self.total_men() > 0 else 0
        if average_schooling < baseline_schooling:
            education_improvement = 0.0
        if self.total_men() == 0: # If there are no men in the household,compare the females' income to the overall average income of men in the community
            gender_wage_gap = (1 - ((average_income_men - average_household_income_women) / (average_income_men))) if average_household_income_men > 0 else 0
        else:
            gender_wage_gap = (1 - ((average_household_income_men - average_household_income_women) / (average_household_income_men))) if average_household_income_men > 0 else 0
        impact_score =  (health_improvement * 0.3) + (co2_reduction * 0.3)  + (education_improvement * 0.2) + (gender_wage_gap*0.2) # Calculates the overall impact score as a weighted sum of the improvements in health, CO2 emissions, and education, gender wage gap with weights of 0.3, 0.3, and 0.2, 0.2 respectively.
        return impact_score
    
    def utility_function(self,eta,air_change_rate_daily, kitchen_volume, outdoor_concentration, baseline_schooling,average_income_men): #Calculates the utility for the household based on its total income and a given inequality aversion parameter (eta), using a constant relative risk aversion (CRRA) utility function.
        income = self.total_income()
        if eta != 1:
            income_utility = (income**(1-eta))/(1-eta)
        else:
            income_utility = math.log(income+1)
        impact_score = self.impact_score(air_change_rate_daily, kitchen_volume, outdoor_concentration, baseline_schooling,average_income_men)
        utility_score = impact_score * income_utility
        return utility_score

    def consider_microgrid_adoption(self, reliability, social_influence,cost_per_kwh,air_change_rate_daily, kitchen_volume, outdoor_concentration):
        alpha_income = 0.001 #The weight assigned to income in the utility function.
        alpha_reliability = 0.5#The weight assigned to the reliability of the electricity supply in the utility function.
        alpha_health = 0.1 #The weight assigned to health improvements from reduced PM2.5 exposure in the utility function.
        alpha_social = 0.2 #The weight assigned to social influence from neighbours who have already connected to the microgrid in the utility function.
        alpha_savings = 1.0 #The weight assigned to cost savings from switching to electricity compared to fuel costs in the utility function.
        aversion_to_adoption = -3.5 #A baseline aversion to adopting the microgrid, representing factors such as uncertainty, risk aversion, or attachment to traditional energy sources that may make households hesitant to connect.
        k = 0.8  # sensitivity parameter for adoption probability
        reliability_penalty = -1.4 if reliability < 0.8 else 0.0 #A penalty applied to the utility if the reliability of the electricity supply is below a certain threshold, reflecting the negative impact of unreliable electricity on the household's willingness to adopt the microgrid.
        income = self.total_income()
        cases_baseline = self.baseline_cases(air_change_rate_daily, kitchen_volume, outdoor_concentration)
        cases = self.cases(air_change_rate_daily, kitchen_volume, outdoor_concentration)
        health_improvement = (cases_baseline - cases)/cases_baseline if cases_baseline > 0 else 0
        cost_savings = (self.baseline_fuel_cost() - self.energy_cost(cost_per_kwh))/(self.baseline_fuel_cost()) if self.baseline_fuel_cost() > 0 else 0
        price_penalty = 0.0
        if cost_savings < 0:
            cost_savings = 0.0
            price_penalty = -1.0 #A penalty applied to the utility if the cost savings from switching to electricity are negative, reflecting the negative impact of higher costs on the household's willingness to adopt the microgrid.
    
        microgrid_utility = (
                alpha_income * income
                + alpha_reliability * reliability
                + alpha_social * social_influence
                + alpha_health * health_improvement
                + alpha_savings * cost_savings
                + reliability_penalty
                + price_penalty
            )
        
        microgrid_adoption_prob = 1 / (1 + math.exp(-(aversion_to_adoption + k*microgrid_utility)))
        
        if random.random() < microgrid_adoption_prob:
            self.electrify()            


# --- System Dynamics Model ---
class ElectrificationSD:
    def __init__(self, model): #The constructor for the ElectrificationSD class initializes the stocks, flows, converters, and constants of the SD model.
        self.model = model

        # Initialization of Stocks
        self.net_microgrid_income = model.stock("net_microgrid_income")
        self.no_of_failures = model.stock("no_of_failures")
        self.no_of_working_components = model.stock("no_of_working_components")
        self.no_of_current_repairs = model.stock("no_of_current_repairs")
        self.available_jobs = model.stock("available_jobs")
        self.discounted_social_welfare_function = model.stock("discounted_social_welfare_function")
        
        # Initialization of Flows
        self.microgrid_income = model.flow("microgrid_income")
        self.microgrid_expenditures = model.flow("microgrid_expenditures")
        self.microgrid_income_flow = model.flow("microgrid_income_flow")
        self.failure_rate = model.flow("failure_rate")
        #self.fixture_rate = model.flow("fixture_rate")
        self.job_creation = model.flow("job_creation")
        self.social_welfare_function = model.flow("social_welfare_function")
        self.employment_change = model.flow("employment_change")
        self.current_working_components = model.flow("current_working_components")
        self.previous_working_components = model.flow("previous_working_components")

        # Initialization of Converters
        self.households_connected = model.converter("households_connected")
        self.adopting_households = model.converter("adopting_households")
        self.power_use = model.converter("power_use")
        self.failure_rate_multiplier =  model.converter("failure_rate_multiplier")
        self.mttr = model.converter("mttr")
        self.saifi = model.converter("saifi")
        self.saidi = model.converter("saidi")
        self.availability = model.converter("availability")
        self.downtime = model.converter("downtime")
        self.reliability = model.converter("reliability")
        self.baseline_demand = model.converter("baseline_demand")
        self.total_demand = model.converter("total_demand")

        self.business_income = model.converter("business_income")
        self.capital_investment = model.converter("capital_investment")
        self.load_ratio = model.converter("load_ratio")
        self.income_ratio = model.converter("income_ratio")


        # Aggregate social/economic flows and converters
        self.changes_in_jobs = model.converter("changes_in_jobs")
        self.avg_income = model.converter("avg_income")
        self.appliance_demand = model.converter("appliance_demand")
        self.utility_sum = model.converter("utility_sum")
        

        # Initialization of Constants
        self.baseline_household_demand = model.constant("baseline_household_demand")
        self.no_of_households = model.constant("no_of_households")
        self.microgrid_capacity = model.constant("microgrid_capacity")
        self.initial_failure_rate = model.constant("initial_failure_rate")
        self.cost_per_kwh = model.constant("cost_per_kwh")
        self.cost_per_repair = model.constant("cost_per_repair")
        self.operating_expenditures = model.constant("operating_expenditures")
        self.no_of_components = model.constant("no_of_components")
        self.attrition_rate = model.constant("attrition_rate")
        self.mpc = model.constant("mpc") #The marginal propensity to consume (MPC) represents the proportion of additional income that households are expected to spend on consumption.
        self.local_spending_fraction = model.constant("local_spending_fraction") #The local spending fraction represents the proportion of household spending that is expected to be directed towards local businesses and services.
        self.investment_rate = model.constant("investment_rate")#The investment rate represents the proportion of business income that is expected to be reinvested in the local economy, contributing to economic growth and job creation.
        self.job_creation_efficiency = model.constant("job_creation_efficiency")
        self.initial_jobs = model.constant("initial_jobs")
        self.eta = model.constant("eta")
        self.discount_rate = model.constant("discount_rate")
        self.capacity_factor = model.constant("capacity_factor")

        # Equations
        self.baseline_demand.equation = ((self.baseline_household_demand) * self.households_connected)
        self.total_demand.equation = (self.baseline_demand+self.appliance_demand)
        self.power_use.equation = sd.min(self.total_demand,(self.microgrid_capacity*720.0*self.capacity_factor))
        self.microgrid_income_flow.equation = self.power_use * self.cost_per_kwh
        self.microgrid_expenditures.equation = (self.no_of_failures*self.cost_per_repair) + self.operating_expenditures
        self.net_microgrid_income.equation = self.microgrid_income_flow - self.microgrid_expenditures
        self.load_ratio.equation = self.total_demand / sd.max(1, self.microgrid_capacity*720.0*self.capacity_factor)
        self.income_ratio.equation = self.microgrid_income_flow / self.operating_expenditures
        self.mttr.equation = sd.lookup( #The mean time to repair (MTTR) is determined by looking up the income ratio in a predefined table.
            self.income_ratio,
            [[0.2, 336], [0.4, 240], [0.6, 168], [0.8,120], [1.0, 96], [1.5, 72],[2.0, 60],[3.0, 48]]
        )

        self.failure_rate_multiplier.equation = sd.lookup(#The failure rate multiplier is determined by looking up the load ratio in a predefined table.
            self.load_ratio,
            [[0.3, 0.6], [0.5, 0.8], [0.7, 1], [0.8, 1.2],[0.9, 1.6],[1.0, 2.5],[1.1, 4], [1.2, 7]]
        )
        self.failure_rate.equation = self.initial_failure_rate * self.no_of_working_components * self.failure_rate_multiplier
        self.no_of_failures.equation = self.failure_rate
        self.availability.equation = (1/self.failure_rate)/sd.max(1e-3, (1/self.failure_rate)+(self.mttr))
        #number of working components
        self.current_working_components.equation = self.no_of_components*self.availability
        self.previous_working_components.equation = self.no_of_working_components
        self.no_of_working_components.equation = self.current_working_components - self.previous_working_components
        # Failure rate depends on maintenance delay
        self.downtime.equation = self.mttr * self.no_of_failures
        # SAIDI (interruptions per customer)
        self.saidi.equation = self.downtime / sd.max(1, self.households_connected)
        # SAIFI (interruptions per customer)
        self.saifi.equation = self.no_of_failures / sd.max(1, self.households_connected)
        # Reliability improves as SAIDI decreases
        self.reliability.equation = 1/(1+0.25*self.saidi)

        # Impact of Income on spending, business investment and job creation
        self.business_income.equation = self.mpc * self.local_spending_fraction * self.avg_income
        self.capital_investment.equation = self.investment_rate * self.business_income
        self.job_creation.equation = self.job_creation_efficiency * self.capital_investment
        self.employment_change.equation = self.changes_in_jobs
        self.available_jobs.equation = self.job_creation + self.employment_change
        # Utility function and social welfare
        self.social_welfare_function.equation = self.utility_sum*sd.exp(-self.discount_rate*sd.time())
        self.discounted_social_welfare_function.equation = self.social_welfare_function
               

        # Initial values
        self.households_connected.initial_value = 0.0
        self.operating_expenditures.equation = 5000.0
        self.initial_failure_rate.equation = 0.1
        self.no_of_failures.initial_value = 0.0
        self.no_of_components.equation = 2.0
        self.no_of_working_components.initial_value = self.no_of_components
        self.no_of_current_repairs.initial_value = 0.0
        self.available_jobs.initial_value = self.initial_jobs
        self.initial_jobs.equation = 75.0
        self.mpc.equation = 0.7
        self.local_spending_fraction.equation = 0.6
        self.investment_rate.equation = 0.1
        self.job_creation_efficiency.equation = 4.074e-3
        self.discount_rate.equation = 4.074e-3
        self.capacity_factor.equation = 0.8


# --- Hybrid Model ---
class ElectrificationHybrid(Model): #The constructor for the ElectrificationHybrid class initializes the agent factory for creating household agents and instantiates the system dynamics model, as well as a dictionary to store custom statistics.
    def instantiate_model(self): #The instantiate_model method is called to set up the model, where it registers the agent factory for creating household agents and initializes the system dynamics model.
        super().instantiate_model()
        self.register_agent_factory("household", lambda agent_id, model, properties: Household(agent_id, model, properties))
        self.sd_model = ElectrificationSD(self)
        self.custom_stats = {}

    def configure(self, config):# The configure method sets up the equations in the system dynamics model based on the parameters defined in the configuration. It calls the parent class's configure method and then assigns the appropriate equations to each variable in the SD model using the parameters from the configuration.
        super().configure(config)
        self.sd_model.no_of_households.equation = self.no_of_households
        self.sd_model.microgrid_capacity.equation = self.microgrid_capacity
        self.sd_model.cost_per_kwh.equation = self.cost_per_kwh
        self.sd_model.initial_failure_rate.equation = self.initial_failure_rate
        self.sd_model.no_of_components.equation = self.no_of_components
        self.sd_model.baseline_household_demand.equation = self.baseline_household_demand
        self.sd_model.cost_per_repair.equation = self.cost_per_repair
        self.sd_model.operating_expenditures.equation = self.operating_expenditures
        self.sd_model.initial_jobs.equation = self.initial_jobs
        self.sd_model.eta.equation = self.eta


    def begin_round(self, time, sim_round, income_and_employment): #The begin_round method is called at the beginning of each simulation round to update the states of the agents and aggregate outcomes.
        # Update agent states
        # begin_round
        available_jobs = int(self.sd_model.available_jobs(time))
        reliability = self.sd_model.reliability(time)
        social_influence = self.sd_model.households_connected(time)/self.no_of_households if self.no_of_households > 0 else 0
        cost_per_kwh = self.cost_per_kwh
        eta = self.eta
        subsidy_rate = self.subsidy_rate
        primary_completion = self.primary_completion
        lower_secondary_completion = self.lower_secondary_completion
        upper_secondary_completion = self.upper_secondary_completion
        employment_rate_baseline = self.employment_rate_baseline
        dropout_rate_baseline = self.dropout_rate_baseline
        hourly_farm_wage = self.hourly_farm_wage
        hourly_non_farm_wage = self.hourly_non_farm_wage
        job_separation_rate = self.job_separation_rate
        baseline_monthly_wage_men = self.baseline_monthly_wage_men
        baseline_monthly_wage_women = self.baseline_monthly_wage_women
        baseline_schooling = self.baseline_schooling
        electrification_effect_men = self.electrification_effect_men
        electrification_effect_women = self.electrification_effect_women
        farm_work_shift = self.farm_work_shift
        air_change_rate_daily = self.air_change_rate_daily
        kitchen_volume = self.kitchen_volume
        outdoor_concentration = self.outdoor_concentration
        dropout_reduction_rate = self.dropout_reduction_rate
        # Demographic parameters
        self.baseline_fertility_rate = 0.08      # annual probability per adult woman
        self.fertility_reduction_rate = 0.3      # reduction when electrified

        self.baseline_mortality_rate = 0.01      # annual base mortality
        self.child_mortality_multiplier = 2.0
        self.elderly_mortality_multiplier = 3.0
        total_job_changes = 0.0
        if time == self.starttime: #At the start of the simulation, it initializes the schooling, income, and CO2 emissions for each household agent, and calculates the baseline number of IHD cases attributable to PM2.5 exposure.
            for h in self.agents:
                h.set_schooling(primary_completion,lower_secondary_completion,upper_secondary_completion)
                h.initial_income() #initialize income for the first round
                h.co2_emissions_initial() #initialize CO2 emissions for the first round
                
            baseline_cases = sum(h.baseline_cases(air_change_rate_daily,kitchen_volume,outdoor_concentration) for h in self.agents)
        for h in self.agents:
            if h.state != "connected" and time%3==0: #household considers microgrid adoption every 3 months
                h.consider_microgrid_adoption(reliability, social_influence,cost_per_kwh,air_change_rate_daily, kitchen_volume, outdoor_concentration)
            job_changes = h.income_and_employment(available_jobs,employment_rate_baseline,hourly_farm_wage,hourly_non_farm_wage,job_separation_rate, dropout_rate_baseline,baseline_monthly_wage_men,baseline_monthly_wage_women,baseline_schooling,electrification_effect_men,electrification_effect_women,farm_work_shift,dropout_reduction_rate)
            available_jobs += job_changes
            total_job_changes += job_changes#The total change in employment across all households is accumulated to update the available jobs in the system dynamics model.
            #for h in self.agents:

            h.birth_process(
            baseline_fertility_rate=self.baseline_fertility_rate,
            fertility_reduction_rate=self.fertility_reduction_rate
            )

            h.death_process(
            baseline_mortality_rate=self.baseline_mortality_rate,
            child_multiplier=self.child_mortality_multiplier,
            elderly_multiplier=self.elderly_mortality_multiplier
            )
            if h.state == "connected":
                h.consider_appliance_adoption(reliability, social_influence,cost_per_kwh,subsidy_rate)#If the household is connected to the microgrid, it considers adopting electric appliances.
                
    
        # Aggregate outcomes from agents
        connected_count = sum(1 for h in self.agents if h.state == "connected")
        avg_income = sum(h.total_income() for h in self.agents) / len(self.agents) if len(self.agents) > 0 else 0
        appliance_demand = sum(h.appliance_demand() for h in self.agents)
        co2_emissions = sum(h.co2_emissions() for h in self.agents)
        total_cases = sum(h.cases(air_change_rate_daily,kitchen_volume,outdoor_concentration) for h in self.agents)
        average_income_women = sum(h.income_women() for h in self.agents)/sum(h.total_women() for h in self.agents) if sum(h.total_women() for h in self.agents) > 0 else 0
        average_income_men = sum(h.income_men() for h in self.agents)/sum(h.total_men() for h in self.agents) if sum(h.total_men() for h in self.agents) > 0 else 0
        average_schooling = sum(h.total_schooling() for h in self.agents)/sum(h.total_adult_children() for h in self.agents) if sum(h.total_adult_children() for h in self.agents) > 0 else 0
        impact_score = sum(h.impact_score(air_change_rate_daily, kitchen_volume, outdoor_concentration,baseline_schooling,average_income_men) for h in self.agents)/len(self.agents) if len(self.agents) > 0 else 0
        social_welfare_function = sum(h.utility_function (eta,air_change_rate_daily, kitchen_volume, outdoor_concentration, baseline_schooling,average_income_men) for h in self.agents)
        
        # Feed back into SD converters
        self.sd_model.households_connected.equation = connected_count
        self.sd_model.avg_income.equation = avg_income
        self.sd_model.appliance_demand.equation = appliance_demand
        self.sd_model.utility_sum.equation = social_welfare_function
        self.sd_model.changes_in_jobs.equation = total_job_changes
        employment = sum(
        1 for h in self.agents
        for m in h.members
        if m["employed"]
        )



        if time not in self.custom_stats:
            self.custom_stats[time] = {}

        self.custom_stats[time]["connected_count"] = connected_count
        self.custom_stats[time]["not_connected_count"] = self.no_of_households - connected_count
        self.custom_stats[time]["appliance_demand"] = appliance_demand
        self.custom_stats[time]["total_cases"] = total_cases
        self.custom_stats[time]["avg_income"] = avg_income
        self.custom_stats[time]["saifi"] = self.sd_model.saifi(time)
        self.custom_stats[time]["saidi"] = self.sd_model.saidi(time)
        self.custom_stats[time]["reliability"] = reliability
        self.custom_stats[time]["Employment"] = employment
        self.custom_stats[time]["Jobs Available"] = available_jobs
        self.custom_stats[time]["Demand"] = self.sd_model.total_demand(time)
        self.custom_stats[time]["Income"] = self.sd_model.net_microgrid_income(time)
        self.custom_stats[time]["SWF"] = social_welfare_function
        self.custom_stats[time]["DSWF"] = self.sd_model.discounted_social_welfare_function(time)
        self.custom_stats[time]["co2_emissions"] = co2_emissions
        self.custom_stats[time]["average_schooling"] = average_schooling
        self.custom_stats[time]["average_income_women"] = average_income_women
        self.custom_stats[time]["average_income_men"] = average_income_men
        self.custom_stats[time]["impact_score"] = impact_score
        self.custom_stats[time]["microgrid_net_income"] = self.sd_model.net_microgrid_income(time)
        self.custom_stats[time]["power_use"] = self.sd_model.power_use(time)
        self.custom_stats[time]["no_of_current_working_components"] = self.sd_model.no_of_working_components(time)
        self.custom_stats[time]["available_jobs"] = self.sd_model.available_jobs(time)

        
# --- Simulation Setup ---
def run_impact_estimator_model(tariff_value,subsidy_rate,no_of_households_value,stoptime,microgrid_capacity,failure_rate,no_of_components,farm_work_shift,baseline_household_demand,air_change_rate_daily,kitchen_volume,outdoor_concentration,primary_completion,lower_secondary_completion,upper_secondary_completion,employment_rate_baseline,dropout_rate_baseline,hourly_farm_wage,hourly_non_farm_wage,job_separation_rate,baseline_monthly_wage_men,baseline_monthly_wage_women,baseline_schooling,electrification_effect_men,electrification_effect_women,dropout_reduction_rate):
    electrification_hybrid = ElectrificationHybrid(
        1, stoptime, dt=1, name="Electrification Hybrid",
        scheduler=SimultaneousScheduler(),

    )

    electrification_hybrid.instantiate_model()

    electrification_config = {
        "runspecs": {
            "starttime": 1,
            "stoptime": stoptime,
            "dt": 1.0
        },
        "properties": {
            "subsidy_rate": {"type": "Double", "value": subsidy_rate},
            "no_of_households": {"type": "Double", "value": no_of_households_value},
            "microgrid_capacity": {"type": "Double", "value": microgrid_capacity},#KW
            "cost_per_kwh": {"type": "Double", "value": tariff_value},
            "initial_failure_rate": {"type": "Double", "value": failure_rate},
            "no_of_components": {"type": "Double", "value": no_of_components},
            "baseline_household_demand": {"type": "Double", "value": baseline_household_demand},#KW
            "cost_per_repair": {"type": "Double", "value": 1000.0},
            "operating_expenditures": {"type": "Double", "value": 5000.0},
            "initial_jobs": {"type": "Double", "value": 75.0},
            "eta":{"type": "Double", "value": 1.0},
            "primary_completion":{"type": "Double", "value": primary_completion},
            "lower_secondary_completion":{"type": "Double", "value": lower_secondary_completion},
            "upper_secondary_completion":{"type": "Double", "value": upper_secondary_completion},
            "employment_rate_baseline":{"type": "Double", "value": employment_rate_baseline},
            "dropout_rate_baseline":{"type": "Double", "value": dropout_rate_baseline},
            "hourly_farm_wage":{"type": "Double", "value": hourly_farm_wage},
            "hourly_non_farm_wage":{"type": "Double", "value": hourly_non_farm_wage},
            "job_separation_rate":{"type": "Double", "value": job_separation_rate},
            "baseline_monthly_wage_men":{"type": "Double", "value": baseline_monthly_wage_men},
            "baseline_monthly_wage_women":{"type": "Double", "value": baseline_monthly_wage_women},
            "baseline_schooling":{"type": "Double", "value": baseline_schooling},
            "electrification_effect_men":{"type": "Double", "value": electrification_effect_men},
            "electrification_effect_women":{"type": "Double", "value": electrification_effect_women},
            "farm_work_shift":{"type": "Double", "value": farm_work_shift},
            "air_change_rate_daily":{"type": "Double", "value": air_change_rate_daily},
            "kitchen_volume":{"type": "Double", "value": kitchen_volume},
            "outdoor_concentration":{"type": "Double", "value": outdoor_concentration},
            "dropout_reduction_rate":{"type": "Double", "value": dropout_reduction_rate}
        },
        "agents": [
            {"name": "household", "count": no_of_households_value}
        ]
    }
   
  
    electrification_hybrid.configure(electrification_config)
    electrification_hybrid.run()
    return electrification_hybrid.custom_stats.copy()#The function run_impact_estimator_model sets up and runs the electrification hybrid model with the specified parameters, and returns the custom statistics collected during the simulation.
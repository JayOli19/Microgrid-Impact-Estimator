import pytest
from BPTK_Py import Model
from impact_estimator_model import Household
from impact_estimator_model import run_impact_estimator_model


class TestModel(Model):
    def __init__(self):
        super().__init__(starttime=0, stoptime=1)


def test_household_initialization():
    model = TestModel()

    h = Household(
        agent_id=1,
        model=model,
        properties={"electrified": False}
    )
    h.initialize()
    assert h.state == "unconnected"
    assert len(h.members) >= 1
    assert all("role" in m for m in h.members)
    assert all(name in h.fuel_appliances for name in ["lpg_stove", "kerosene_lamp"])


def test_electrification_changes_state():
    model = TestModel()

    h = Household(
        agent_id=1,
        model=model,
        properties={"electrified": False}
    )
    h.initialize()

    h.electrify()
    
    assert h.state == "connected"

def test_pm25_reduction_after_electric_cooking():
    model = TestModel()
    h = Household(
        agent_id=1,
        model=model,
        properties={"electrified": False}
    )
    h.initialize()

    baseline = h.pm25_concentration_baseline(
        air_change_rate_daily=24,
        kitchen_volume=30,
        outdoor_concentration=1e-6
    )

    h.fuel_appliances["lpg_stove"] = False

    after = h.pm25_concentration(
        air_change_rate_daily=24,
        kitchen_volume=30,
        outdoor_concentration=1e-6
    )

    assert after < baseline

def test_appliance_demand_zero_initially():
    model = TestModel()
    h = Household(
        agent_id=1,
        model=model,
        properties={"electrified": False}
    )
    h.initialize()

    assert h.appliance_demand() == 0.0
    

def test_impact_score_bounds():
    model = TestModel()
    h = Household(
        agent_id=1,
        model=model,
        properties={"electrified": False}
    )
    h.initialize()
    h.electrify()

    score = h.impact_score(
        air_change_rate_daily=24,
        kitchen_volume=30,
        outdoor_concentration=1e-6,
        baseline_schooling=5,
        average_income_men=100
    )

    assert 0.0 <= score <= 1.0

def test_baseline_cases_remains_constant():
    model = TestModel()
    h = Household(
        agent_id=1,
        model=model,
        properties={"electrified": False}
    )
    h.initialize()

    baseline = h.pm25_concentration_baseline(
        air_change_rate_daily=24,
        kitchen_volume=30,
        outdoor_concentration=1e-6
    )
    h.fuel_appliances["lpg_stove"] = False
    h.fuel_appliances["kerosene_lamp"] = False
    new_baseline = h.pm25_concentration_baseline(
        air_change_rate_daily=24,
        kitchen_volume=30,
        outdoor_concentration=1e-6
    )
    assert baseline == new_baseline

def test_energy_cost_calculation():
    model = TestModel()
    h = Household(
        agent_id=1,
        model=model,
        properties={"electrified": False}
    )
    h.initialize()
    h.fuel_appliances["lpg_stove"] = False
    h.fuel_appliances["kerosene_lamp"] = False
    h.electric_appliances["electric_cooker"] = True
    h.electric_appliances["electric_lighting"] = True
    h.electric_appliances["tv"] = True
    cost = h.energy_cost(cost_per_kwh=0.2)

    assert cost == 0.2 * (30.0 + 7.5 + 4.0)

def test_fossil_fuel_use_calculation():
    model = TestModel()
    h = Household(
        agent_id=1,
        model=model,
        properties={"electrified": False}
    )
    h.initialize()
    h.fuel_appliances["lpg_stove"] = True
    h.fuel_appliances["kerosene_lamp"] = True
    fuel_use = h.fossil_fuel_use()

    assert fuel_use == (0.2*1*30)+(4*0.06*30)

def test_connected_households_never_exceed_total():
    results = run_impact_estimator_model(
        tariff_value=0.2,
        subsidy_rate=0.3,
        no_of_households_value=40,
        stoptime=24,
        microgrid_capacity=500,
        failure_rate=0.05,
        no_of_components=2,
        farm_work_shift=0.5,
        baseline_household_demand=200,
        air_change_rate_daily=24,
        kitchen_volume=30,
        outdoor_concentration=1e-6,
        primary_completion=0.9,
        lower_secondary_completion=0.6,
        upper_secondary_completion=0.3,
        employment_rate_baseline=0.6,
        dropout_rate_baseline=0.1,
        hourly_farm_wage=2.0,
        hourly_non_farm_wage=3.0,
        job_separation_rate=0.05,
        baseline_monthly_wage_men=100,
        baseline_monthly_wage_women=80,
        baseline_schooling=5,
        electrification_effect_men=1.1,
        electrification_effect_women=1.15,
        dropout_reduction_rate=0.3
    )
    times = sorted(results.keys())
    connected_count = [results[t].get("connected_count",0) for t in times]
    assert len(connected_count) == 24
    assert all(0 <= count <= 40 for count in connected_count)

def test_number_of_households_is_constant():
    results = run_impact_estimator_model(
        tariff_value=0.2,
        subsidy_rate=0.3,
        no_of_households_value=40,
        stoptime=24,
        microgrid_capacity=500,
        failure_rate=0.05,
        no_of_components=2,
        farm_work_shift=0.5,
        baseline_household_demand=200,
        air_change_rate_daily=24,
        kitchen_volume=30,
        outdoor_concentration=1e-6,
        primary_completion=0.9,
        lower_secondary_completion=0.6,
        upper_secondary_completion=0.3,
        employment_rate_baseline=0.6,
        dropout_rate_baseline=0.1,
        hourly_farm_wage=2.0,
        hourly_non_farm_wage=3.0,
        job_separation_rate=0.05,
        baseline_monthly_wage_men=100,
        baseline_monthly_wage_women=80,
        baseline_schooling=5,
        electrification_effect_men=1.1,
        electrification_effect_women=1.15,
        dropout_reduction_rate=0.3
    )
    times = sorted(results.keys())
    not_connected_count = [results[t].get("not_connected_count",0) for t in times]
    connected_count = [results[t].get("connected_count",0) for t in times]
    print(not_connected_count)
    print(connected_count)
    #assert len(not_connected_count) == 24
    i=0
    while i < (len(not_connected_count)):
        assert not_connected_count[i]+connected_count[i] == 40#Check that the sum of connected and not connected households equals the total number of households (40) at each time step
        i += 1
    
def test_no_of_working_components_bounds():
    results = run_impact_estimator_model(
        tariff_value=0.2,
        subsidy_rate=0.3,
        no_of_households_value=40,
        stoptime=24,
        microgrid_capacity=500,
        failure_rate=0.05,
        no_of_components=2,
        farm_work_shift=0.5,
        baseline_household_demand=200,
        air_change_rate_daily=24,
        kitchen_volume=30,
        outdoor_concentration=1e-6,
        primary_completion=0.9,
        lower_secondary_completion=0.6,
        upper_secondary_completion=0.3,
        employment_rate_baseline=0.6,
        dropout_rate_baseline=0.1,
        hourly_farm_wage=2.0,
        hourly_non_farm_wage=3.0,
        job_separation_rate=0.05,
        baseline_monthly_wage_men=100,
        baseline_monthly_wage_women=80,
        baseline_schooling=5,
        electrification_effect_men=1.1,
        electrification_effect_women=1.15,
        dropout_reduction_rate=0.3
    )
    times = sorted(results.keys())
    no_of_components = [results[t].get("no_of_current_working_components",0) for t in times]
    i=0
    while i < (len(no_of_components)):
        print(no_of_components[i])
        assert 0<= no_of_components[i]<=2
        i += 1
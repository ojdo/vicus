import coopr.environ
import pandas as pd
import vicus
from coopr.opt.base import SolverFactory

# INIT
filename = 'data-example.xlsx'
(offset, length) = (1500, 7*24) # timestep selection
timesteps = range(offset, offset+length+1)


# SCENARIOS

def scenario_base(data):
    # don't change input data
    pass

def scenario_storage_expensive(data):
    # increse storage price by +200%
    sto = data['storage']
    sto.loc[('bat', 'ElecDC'), 'inv-cost-p'] *= 3.0
    sto.loc[('pst', 'ElecAC'), 'inv-cost-p'] *= 3.0
    sto.loc[('bat', 'ElecDC'), 'inv-cost-c'] *= 3.0
    sto.loc[('pst', 'ElecAC'), 'inv-cost-c'] *= 3.0


def scenario_co2_limit(data):
    # change CO2 limit
    co = data['commodity']
    co.loc[('CO2', 'Env'), 'max'] *= 0.05
    

def scenario_co2_limit_and_storage_expensive(data):
    # combine two scenarios
    scenario_storage_expensive(data)
    scenario_co2_limit(data)

scenarios = [
    scenario_base,
    scenario_storage_expensive,
    scenario_co2_limit,
    scenario_co2_limit_and_storage_expensive]


# MAIN 

for scenario in scenarios:
    # scenario name, read and modify data for scenario
    sce = scenario.__name__
    data = vicus.read_excel(filename)
    scenario(data)

    # create model, solve it, read results
    model = vicus.create_model(data, timesteps)
    prob = model.create()
    optim = SolverFactory('glpk')  # cplex, glpk, gurobi, ...
    result = optim.solve(prob, tee=True)
    prob.load(result)
    
    # write report to spreadsheet
    vicus.report(prob, '{}.xlsx'.format(sce), ['ElecAC', 'ElecDC'])
    
    # example: customise color of commodity
    vicus.COLOURS['ElecDC'] = (0, 121, 239)
    
    # timeseries plot for each commodity
    for co in ['ElecAC', 'ElecDC']:
        fig = vicus.plot(prob, co=co)
        
        for ext in ['png', 'pdf']:
            fig_filename = '{}-{}.{}'.format(sce, co, ext)
            fig.savefig(fig_filename, bbox_inches='tight') 



import coopr.environ
import pandas as pd
import vicus
from coopr.opt.base import SolverFactory

# INIT
filename = 'data-example.xlsx'
(offset, length) = (1500, 7*24)
timesteps = range(offset, offset+length+1) # timestep selection

# create model
# model instance
# optimizer to calculate result
model = vicus.create_model(filename, timesteps)
prob = model.create()
optim = SolverFactory('glpk')
result = optim.solve(prob, tee=True)
prob.load(result)

# plot
for co in ['ElecAC', 'ElecDC']:
    fig = vicus.plot(prob, co=co)
    fig.savefig('plot-{}.png'.format(co), bbox_inches='tight') 

# report
vicus.report(prob, 'report.xlsx', ['ElecAC', 'ElecDC'])


"""VICUS: A linear optimisation model for localised energy systems

VICUS is a linear optimization model for a localized energy system. It
minimises total cost for providing energy in form of desired commodities
(usually electricity). The model contains commodities (electricity, fossil
fuels, renewable energy sources, greenhouse gases), processes that convert one
commodity to another (while emitting greenhouse gases as a secondary output),
and storage for saving/retrieving commodities.

It operates on a time-discrete basis. The word "localised" means that
all process and storage entities are connected to a single virtual node 
without transmission losses. All physical quantities in PICUS (except 
for greenhouse gases) represent energy contents of commodities in the 
unit kWh.

Commodities can have one of four types:
  - Stock: Buyable at any time for a given price. Supply can be limited
    per timestep or for a whole year. Examples are coal, gas, uranium
    or biomass.
  - SupIm: Supply intermittent stands for fluctuating resources like
    solar radiation and wind energy, which are available according to 
    a timeseries of values, which could be derived from weather data.
  - Demand: These commodities have a timeseries for the requirement
    associated and must be provided by output from other process or 
    from storage. Usually, there is only one demand commodity called 
    electricity (abbreviated to Elec or ElecAC), but
  - Env: The special commodity CO2 is of this type and represents the
    amount (in tons) of greenhouse gas emissions from processes. Its
    total amount can be limited, to investigate the effect of policies
    on the.
    
    
Process entities are defined over the tuple
    (process name, input commodity, output commodity)
and thus represent the conversion of energy stored in form of input
commodity to the output commodity. For example, the process
    (gas turbine, natural gas, electricity)
represents a natural gas fired power plant operating with a gas turbine.

Storage entities are defined over the tuple
    (storage name, stored commodity)
representing a label for the storage and the commodity it can save and
later retrieve. For example, the storage
    (pump storage, electricity)
represents pumped storage hydro plant that uses electricity to pump
water from one reservoir to a higher reservoir, converting electrical
energy to potential energy. To retrieve electricity, water from the 
higher reservoir is channeled through a generator.

"""
import coopr.pyomo as pyomo
import pandas as pd
from datetime import datetime
from operator import itemgetter

COLOURS = {
    'Biomass': (0, 122, 55),
    'Coal': (100, 100, 100),
    'Demand': (25, 25, 25),
    'Diesel': (116, 66, 65),
    'Gas': (237, 227, 0),   
    'ElecAC': (0, 101, 189),
    'ElecDC': (0, 101, 189),
    'Heat': (230, 112, 36), 
    'Hydro': (198, 188, 240),
    'Import': (128, 128, 200), 
    'Lignite': (116, 66, 65), 
    'Oil': (116, 66, 65), 
    'Overproduction': (190, 0, 99),
    'Slack': (163, 74, 130), 
    'Solar': (243, 174, 0), 
    'Storage': (60, 36, 154),
    'Wind': (122, 179, 225), 
    'Stock': (222, 222, 222),
    'Decoration': (128, 128, 128)}

def read_excel(filename):
    """Read Excel input file and prepare VICUS input dict.

    Reads an Excel spreadsheet that adheres to the structure shown in
    data-example.xlsx, returning a dict of DataFrames, one for each sheet.
    The attribute 'annuity-factor' is derived here from the columns 'wacc'
    and 'depreciation' for 'Process' and 'Storage'.

    Args:
        filename: filename to an Excel spreadsheet with the required sheets
            'Commodity', 'Process', 'Storage', 'Demand' and 'SupIm'.

    Returns:
        a dict of 5 DataFrames
        
    Example:
        >>> data = read_excel('data-example.xlsx')
        >>> data['commodity'].loc[('CO2', 'Env'), 'max']
        25000.0
    """
    with pd.ExcelFile(filename) as xls:
        commodity = xls.parse('Commodity', index_col=['Co', 'Type'])
        process = xls.parse('Process', index_col=['Pro', 'CoIn', 'CoOut'])
        storage = xls.parse('Storage', index_col=['Sto', 'Co'])
        demand = xls.parse('Demand', index_col=['t']) 
        supim = xls.parse('SupIm', index_col=['t'])
    
    # derive annuity factor for process and storage
    process['annuity_factor'] = annuity_factor(
        process['depreciation'], process['wacc'])
    storage['annuity_factor'] = annuity_factor(
        storage['depreciation'], storage['wacc'])

    data = {
        'commodity': commodity,
        'process': process,
        'storage': storage,
        'demand': demand,
        'supim': supim}

    # sort nested indexes to make direct assignments work, cf
    # http://pandas.pydata.org/pandas-docs/stable/indexing.html#the-need-for-sortedness-with-multiindex
    for key in data:
        if isinstance(data[key].index, pd.core.index.MultiIndex):
            data[key].sortlevel(inplace=True)
    return data


def create_model(data, timesteps):
    """ Create a VICUS model object from input data.
    
    Creates and returns a Pyomo ConcreteModel object, given a model
    input file (supported formats: Excel spreadsheet [planned: SQLite DB])
    
    Args:
        data: input dict with fields for Commodities, Processes, 
            Storage and timeseries for Demand and SupIm
        timesteps: numpy array of timestep labels, matching the ones 
            used in the Demand and SupIm timeseries
            
    Returns:
        A coopr.pyomo ConcreteModel object, ready to be instantiated and
        solved.
                  
    """
    m = pyomo.ConcreteModel()
    m.name = 'VICUS'
    m.created = datetime.now().strftime('%Y%m%dT%H%M%S',)
    
    # Preparations
    # ============
    # Data import. Syntax to access a value within equation definitions looks
    # like this:
    #
    #     m.process.loc[sit, pro, coin, cout][attribute]
    #
    get_inputs = itemgetter(
        "commodity", "process", "storage", "demand", "supim")
    (m.commodity, m.process, m.storage, m.demand, m.supim) = get_inputs(data)
      
    # Sets
    # ====
    # Syntax: m.{name} = Set({domain}, initialize={values})
    # where name: set name
    #       domain: set domain for tuple sets, a cartesian set product
    #       values: set values, a list or array of element tuples
    m.t = pyomo.Set(ordered=True, initialize=timesteps)
    m.tm = pyomo.Set(within=m.t, initialize=timesteps[1:])
    m.co = pyomo.Set(initialize=m.commodity.index.levels[0])
    m.coin = pyomo.Set(within=m.co,initialize=m.process.index.levels[1])
    m.cout = pyomo.Set(within=m.co,initialize=m.process.index.levels[2])
    m.co_type = pyomo.Set(initialize=m.commodity.index.levels[1])
    m.pro = pyomo.Set(initialize=m.process.index.levels[0])
    m.sto = pyomo.Set(initialize=m.storage.index.levels[0])
    m.cost_type = pyomo.Set(initialize=['Inv', 'Fix', 'Var', 'Fuel'])

    # sets of existing tuples:
    # co_tuples = [('Coal', 'Stock'), ('Wind', 'SupIm'), ('ElecAC', 'Demand')...]
    # pro_tuples = [('pp', 'Coal', 'ElecAC'), ('wt', 'Wind', 'ElecAC')...]
    # sto_tuples = [('bat', 'ElecDC'), ('pst', 'ElecAC')...]
    m.co_tuples = pyomo.Set(within=m.co*m.co_type, initialize=m.commodity.index)  
    m.pro_tuples = pyomo.Set(within=m.pro*m.coin*m.cout, initialize=m.process.index)
    m.sto_tuples = pyomo.Set(within=m.sto*m.co, initialize=m.storage.index)
    
    # subsets of commodities by type
    # for shorter equations that apply to only one commodity type
    m.co_supim = pyomo.Set(within=m.co, initialize=
                          (c[0] for c in m.co_tuples if c[1] == 'SupIm'))
    m.co_stock = pyomo.Set(within=m.co, initialize=
                          (c[0] for c in m.co_tuples if c[1] == 'Stock'))
    m.co_demand = pyomo.Set(within=m.co, initialize=
                           (c[0] for c in m.co_tuples if c[1] == 'Demand'))
    
    # Parameters
    # ==========
    # for model entities (commodity, process, storage) no Pyomo params
    # are needed, just use the DataFrames m.commodity, m.process and
    # m.storage directly.
    # Syntax: m.{name} = Param({domain}, initialize={values})
    # where name: param name
    #       domain: one or multiple model sets; empty for scalar parameters
    #       values: dict of values, addressed by elements of domain sets
    m.weight = pyomo.Param(initialize=float(8760) / len(m.t))
 
    # Variables
    # =========
    # listed alphabetically
    # Syntax: m.{name} = Var({domain}, within={range})
    # where name: variable name
    #       domain: variable domain, consisting of one or multiple sets
    #       range: variable values, like Binary, Integer, NonNegativeReals
    m.cap_pro = pyomo.Var(
        m.pro_tuples, within=pyomo.NonNegativeReals,
        doc='Total process capacity (kW)')
    m.cap_pro_new = pyomo.Var(
        m.pro_tuples, within=pyomo.NonNegativeReals,
        doc='New process capacity (kW)')
    m.cap_sto_c = pyomo.Var(
        m.sto_tuples, within=pyomo.NonNegativeReals,
        doc='Total storage size (kWh)')
    m.cap_sto_c_new = pyomo.Var(
        m.sto_tuples, within=pyomo.NonNegativeReals,
        doc='New storage capacity (kWh)')
    m.cap_sto_p = pyomo.Var(
        m.sto_tuples, within=pyomo.NonNegativeReals,
        doc='Total storage power (kW)')
    m.cap_sto_p_new = pyomo.Var(
        m.sto_tuples, within=pyomo.NonNegativeReals,
        doc='New storage power (kW)')
    m.co2_pro_out = pyomo.Var(
        m.tm, m.pro_tuples, within=pyomo.NonNegativeReals,
        doc='CO2 emissions from processes (kg) per timestep')
    m.costs = pyomo.Var(
        m.cost_type, within=pyomo.NonNegativeReals,
        doc='Costs by type (EUR/a)')
    m.e_co_stock = pyomo.Var(
        m.tm, m.co_stock, within=pyomo.NonNegativeReals,
        doc='Source power flow from stock commodities (kW) per timestep')
    m.e_pro_in = pyomo.Var(
        m.tm, m.pro_tuples, within=pyomo.NonNegativeReals,
        doc='Power flow into process (kW) per timestep')
    m.e_pro_out = pyomo.Var(
        m.tm, m.pro_tuples, within=pyomo.NonNegativeReals,
        doc='Power flow out of process (kW) per timestep')
    m.e_sto_in = pyomo.Var(
        m.tm, m.sto_tuples, within=pyomo.NonNegativeReals,
        doc='Power flow into storage (kW) per timestep')
    m.e_sto_out = pyomo.Var(
        m.tm, m.sto_tuples, within=pyomo.NonNegativeReals,
        doc='Power flow out of storage (kW) per timestep')
    m.e_sto_con = pyomo.Var(
        m.t, m.sto_tuples, within=pyomo.NonNegativeReals,
        doc='Energy content of storage (kWh) in timestep')
    
    # Equation definition
    # ===================
    # listed by topic. All equations except the Objective function are 
    # of type Constraint, although there are two semantics for those, 
    # indicated by the name prefix (def, res).
    #  - def: definition, usually equations, defining variable values
    #  - res: restriction, usually inequalities, limiting variable values
    # topics
    #  - commodity
    #  - process
    #  - storage
    #  - emissions
    #  - costs
    
    # commodity
    def res_demand_rule(m, tm, co, co_type):
        if co not in m.co_demand:
            return pyomo.Constraint.Skip
        else:
            provided_energy = - commodity_balance(m, tm, co)
            return provided_energy >= \
                   m.demand.loc[tm][co] * \
                   m.commodity.loc[co, co_type]['peak']
    
    def def_e_co_stock_rule(m, tm, co, co_type):
        if co not in m.co_stock:
            return pyomo.Constraint.Skip
        else:
            return m.e_co_stock[tm, co] == commodity_balance(m, tm, co)
    
    def res_stock_hour_rule(m, tm, co, co_type):
        if co not in m.co_stock:
            return pyomo.Constraint.Skip
        else:          
            return m.e_co_stock[tm, co] <= \
                   m.commodity.loc[co, co_type]['maxperhour']
        
    def res_stock_total_rule(m, co, co_type):
        if co not in m.co_stock:
            return pyomo.Constraint.Skip
        else:
            # calculate total consumption of commodity co
            total_consumption = 0
            for tm in m.tm:
                total_consumption += m.e_co_stock[tm, co] * m.weight
            return total_consumption <= m.commodity.loc[co, co_type]['max']
            
    # process
    def def_process_capacity_rule(m, tm, pro, coin, cout):
        return m.cap_pro[pro,coin,cout] == \
               m.cap_pro_new[pro,coin,cout] + \
               m.process.loc[pro,coin,cout]['inst-cap']
                                      
    def def_process_output_rule(m, tm, pro, coin, cout):
        return m.e_pro_out[tm, pro, coin, cout] == \
               m.e_pro_in[tm, pro, coin, cout] * \
               m.process.loc[pro, coin, cout]['eff']

    def def_intermittent_supply_rule(m, tm, pro, coin, cout):
        if coin in m.co_supim:
            return m.e_pro_in[tm, pro, coin, cout] == \
                   m.cap_pro[pro, coin, cout] * m.supim.loc[tm][coin]
        else:
            return pyomo.Constraint.Skip
        
    def def_co2_emissions_rule(m, tm, pro, coin, cout):
        return m.co2_pro_out[tm, pro, coin, cout] == \
               m.e_pro_in[tm, pro, coin, cout] * \
               m.process.loc[pro, coin, cout]['co2'] * \
               m.weight
        
    def res_process_output_by_capacity_rule(m, tm, pro, coin, cout):
        return m.e_pro_out[tm, pro, coin, cout] <= m.cap_pro[pro, coin, cout]
        
    def res_process_capacity_rule(m, pro, coin, cout):
        return (m.process.loc[pro, coin, cout]['cap-lo'],
                m.cap_pro[pro, coin, cout],
                m.process.loc[pro, coin, cout]['cap-up'])
    
    # storage
    def def_storage_state_rule(m, t, sto, co):
        return m.e_sto_con[t, sto, co] == \
               m.e_sto_con[t-1, sto, co] + \
               m.e_sto_in[t, sto, co] * m.storage.loc[sto, co]['eff-in'] - \
               m.e_sto_out[t, sto, co] / m.storage.loc[sto, co]['eff-out']
    
    def def_storage_power_rule(m, sto, co):
        return m.cap_sto_p[sto, co] == \
               m.cap_sto_p_new[sto, co] + \
               m.storage.loc[sto, co]['inst-cap-p']

    def def_storage_capacity_rule(m, sto, co):
        return m.cap_sto_c[sto, co] == \
               m.cap_sto_c_new[sto, co] + \
               m.storage.loc[sto, co]['inst-cap-p']

    def res_storage_input_by_power_rule(m, t, sto, co):
        return m.e_sto_in[t, sto, co] <= m.cap_sto_p[sto, co]
        
    def res_storage_output_by_power_rule(m, t, sto, co):
        return m.e_sto_out[t, sto, co] <= m.cap_sto_p[sto, co]
        
    def res_storage_state_by_capacity_rule(m, t, sto, co):
        return m.e_sto_con[t, sto, co] <= m.cap_sto_c[sto, co]
               
    def res_storage_power_rule(m, sto, co):
        return (m.storage.loc[sto, co]['cap-lo-p'],
                m.cap_sto_p[sto, co],
                m.storage.loc[sto, co]['cap-up-p'])
        
    def res_storage_capacity_rule(m, sto, co):
        return (m.storage.loc[sto, co]['cap-lo-c'],
                m.cap_sto_c[sto, co],
                m.storage.loc[sto, co]['cap-up-c'])
                
    def res_initial_and_final_storage_state_rule(m, t, sto, co):
        if t == m.t[1]: # first timestep (Pyomo uses 1-based indexing)
            return m.e_sto_con[t, sto, co] == \
                   m.cap_sto_c[sto, co] * \
                   m.storage.loc[sto, co]['init']
        elif t == m.t[-1]: # last timestep
            return m.e_sto_con[t, sto, co] >= \
                   m.cap_sto_c[sto, co] * \
                   m.storage.loc[sto, co]['init']
        else:
            return pyomo.Constraint.Skip
    
    # emissions
    def res_co2_emission_rule(m):
        return pyomo.summation(m.co2_pro_out) <= \
               m.commodity.loc['CO2','Env']['max']
    
    # costs
    def def_costs_rule(m, cost_type):
        """ Calculate total costs by cost type.
        
        Sums up process activity and capacity expansions
        and sums them in the cost types that are specified in the set
        m.cost_type. To change or add cost types, add/change entries 
        there and modify the if/elif cases in this function accordingly.
        
        Cost types are
          - Investment costs for process power, storage power and 
            storage capacity, annualized.
          - Fixed costs for process & storage power, storage capacity.
          - Variable costs for process and storage activity.
          - Fuel costs for purchased stock commodities.
        """
        if cost_type == 'Inv':
            return m.costs['Inv'] == \
                sum(m.cap_pro_new[p] * 
                    m.process.loc[p]['inv-cost'] * 
                    m.process.loc[p]['annuity_factor'] 
                    for p in m.pro_tuples) + \
                sum(m.cap_sto_p_new[s] * 
                    m.storage.loc[s]['inv-cost-p'] * 
                    m.storage.loc[s]['annuity_factor'] +
                    m.cap_sto_c_new[s] * 
                    m.storage.loc[s]['inv-cost-c'] * 
                    m.storage.loc[s]['annuity_factor'] 
                    for s in m.sto_tuples)
                
        elif cost_type == 'Fix':
            return m.costs['Fix'] == \
                sum(m.cap_pro[p] * m.process.loc[p]['fix-cost'] 
                    for p in m.pro_tuples) + \
                sum(m.cap_sto_p[s] * m.storage.loc[s]['fix-cost-p'] +
                    m.cap_sto_c[s] * m.storage.loc[s]['fix-cost-c'] 
                    for s in m.sto_tuples)
                
        elif cost_type == 'Var':
            return m.costs['Var'] == \
                sum(m.e_pro_out[(tm,) + p] * 
                    m.process.loc[p]['var-cost'] * 
                    m.weight 
                    for tm in m.tm for p in m.pro_tuples) + \
                sum(m.e_sto_con[(tm,) + s] * 
                    m.storage.loc[s]['var-cost-c'] * m.weight +
                    (m.e_sto_in[(tm,) + s] + m.e_sto_out[(tm,) + s]) * 
                    m.storage.loc[s]['var-cost-p'] * m.weight
                    for tm in m.tm for s in m.sto_tuples)
            
        elif cost_type == 'Fuel':
            return m.costs['Fuel'] == \
                sum(m.e_co_stock[(tm,c[0])] * 
                    m.commodity.loc[c]['price'] * 
                    m.weight 
                    for tm in m.tm for c in m.co_tuples 
                    if c[0] in m.co_stock)
            
        else:
            raise NotImplementedError("Unknown cost type!")
            
    def obj_rule(m):
        """ Return sum of total costs over all cost types.
        
        Simply calculates the sum of m.costs over all m.cost_types.
        """
        return pyomo.summation(m.costs)
    
    
    # Equation declaration
    # ====================
    # declarations connect rule functions to the model, specifying
    # the model sets for which the constraints are enforced. 
    
    # commodity
    m.res_demand = pyomo.Constraint(m.tm, m.co_tuples,
        doc='storage + process balance >= demand')
    m.def_e_co_stock = pyomo.Constraint(m.tm, m.co_tuples,
        doc='commodity source term = hourly commodity consumption')
    m.res_stock_hour = pyomo.Constraint(m.tm, m.co_tuples,
        doc='hourly commodity source term <= commodity.maxperhour')
    m.res_stock_total = pyomo.Constraint(m.co_tuples,
        doc='total commodity source term <= commodity.max')
    
    # process
    m.def_process_capacity = pyomo.Constraint(m.tm, m.pro_tuples,
        doc='total process capacity = inst-cap + new capacity')
    m.def_process_output = pyomo.Constraint(m.tm, m.pro_tuples,
        doc='process output = process input * efficiency')
    m.def_intermittent_supply = pyomo.Constraint(m.tm, m.pro_tuples,
        doc='process output = process capacity * supim timeseries')
    m.def_co2_emissions = pyomo.Constraint(m.tm, m.pro_tuples,
        doc='process co2 output = process input * process.co2 * weight')
    m.res_process_output_by_capacity = pyomo.Constraint(m.tm, m.pro_tuples,
        doc='process output <= process capacity')
    m.res_process_capacity = pyomo.Constraint(m.pro_tuples,
        doc='process.cap-lo <= process capacity <= process.cap-up')
    
    # storage 
    m.def_storage_state = pyomo.Constraint(m.tm, m.sto_tuples,
        doc='storage[t] = storage[t-1] + input - output')
    m.def_storage_power = pyomo.Constraint(m.sto_tuples,
        doc='storage power = inst-cap + new power')
    m.def_storage_capacity = pyomo.Constraint(m.sto_tuples,
        doc='storage capacity = inst-cap + new capacity')
    m.res_storage_input_by_power = pyomo.Constraint(m.tm, m.sto_tuples,
        doc='storage input <= storage power')
    m.res_storage_output_by_power = pyomo.Constraint(m.tm, m.sto_tuples,
        doc='storage output <= storage power')
    m.res_storage_state_by_capacity = pyomo.Constraint(m.t, m.sto_tuples,
        doc='storage content <= storage capacity')
    m.res_storage_power = pyomo.Constraint(m.sto_tuples,
        doc='storage.cap-lo <= storage power <= storage.cap-up')
    m.res_storage_capacity = pyomo.Constraint(m.sto_tuples,
        doc='storage.cap-lo <= storage capacity <= storage.cap-up')
    m.res_initial_and_final_storage_state = pyomo.Constraint(m.t, m.sto_tuples,
        doc='storage content initial == and final >= storage.init * capacity')
    
    # emissions
    m.res_co2_emission = pyomo.Constraint(
        doc='total co2 emissions <= commodity.co2.max')

    # costs
    m.def_costs = pyomo.Constraint(m.cost_type,
        doc='main cost function by cost type')
    m.obj = pyomo.Objective(
        sense=pyomo.minimize,
        doc='cost = sum of all cost types')    
        
    return m

          
def annuity_factor(n, i):
    """ Calculate annuity factor from depreciation and interest.
    
    Evaluates the annuity factor formula for depreciation duration
    and interest rate. Works also well for equally sized numpy arrays 
    of values for n and i.
    
    Args:
        n: depreciation period (years)
        i: interest rate (percent, e.g. 0.06 means 6 %)
        
    Returns:
        Value of the annuity factor formula
            (1+i)**n * i / ((1+i)**n - 1)

    """
    return (1+i)**n * i / ((1+i)**n - 1)


def commodity_balance(m, tm, co):
    """ Calculate commodity balance at given timestep.
    
    For a given commodity co and timestep tm, calculate the balance of
    consumed (to process/storage, counts positive) and provided (from
    process/storage, counts negative) energy. Used as helper function 
    in create_model for constraints on demand and stock commodities.
    
    Args:
        m: the model object
        tm: the timestep
        co: the commodity
        
    Returns
        Pyomo sum expression for net value of consumed (positive) or 
        provided (negative) energy.

    """
    balance = 0
    for p in m.pro_tuples:
        if p[1] == co:
            # usage as input for process increases balance
            balance += m.e_pro_in[(tm,)+p]
        if p[2] == co:
            # output from processes decreases balance
            balance -= m.e_pro_out[(tm,)+p]
    for s in m.sto_tuples:
        # usage as input for storage increases consumption
        # output from storage decreases consumption
        if s[1] == co:
            balance += m.e_sto_in[(tm,)+s]
            balance -= m.e_sto_out[(tm,)+s]
    return balance


def get_entity(instance, name):
    """ Return a DataFrame for an entity in model instance.

    Args:
        instance: a Pyomo ConcreteModel instance
        name: name of a Set, Param, Var, Constraint or Objective

    Returns:
        a single-columned Pandas DataFrame with domain as index
    """

    # retrieve entity, its type and its onset names
    entity = instance.__getattribute__(name)
    labels = get_onset_names(entity)

    # extract values
    if isinstance(entity, pyomo.Set):
        # Pyomo sets don't have values, only elements
        results = pd.DataFrame([(v, 1) for v in entity.value])

        # for unconstrained sets, the column label is identical to their index
        # hence, make index equal to entity name and append underscore to name
        # (=the later column title) to preserve identical index names for both
        # unconstrained supersets
        if not labels:
            labels = [name]
            name = name+'_'

    elif isinstance(entity, pyomo.Param):
        if entity.dim() > 1:
            results = pd.DataFrame([v[0]+(v[1],) for v in entity.iteritems()])
        else:
            results = pd.DataFrame(entity.iteritems())
    else:
        # create DataFrame
        if entity.dim() > 1:
            # concatenate index tuples with value if entity has
            # multidimensional indices v[0]
            results = pd.DataFrame(
                [v[0]+(v[1].value,) for v in entity.iteritems()])
        else:
            # otherwise, create tuple from scalar index v[0]
            results = pd.DataFrame(
                [(v[0], v[1].value) for v in entity.iteritems()])

    # check for duplicate onset names and append one to several "_" to make
    # them unique, e.g. ['sit', 'sit', 'com'] becomes ['sit', 'sit_', 'com']
    for k, label in enumerate(labels):
        if label in labels[:k]:
            labels[k] = labels[k] + "_"

    # name columns according to labels + entity name
    results.columns = labels + [name]
    results.set_index(labels, inplace=True)

    return results


def get_entities(instance, names):
    """ Return one DataFrame with entities in columns and a common index.

    Works only on entities that share a common domain (set or set_tuple), which
    is used as index of the returned DataFrame.

    Args:
        instance: a Pyomo ConcreteModel instance
        names: list of entity names (as returned by list_entities)

    Returns:
        a Pandas DataFrame with entities as columns and domains as index
    """

    df = pd.DataFrame()
    for name in names:
        other = get_entity(instance, name)

        if df.empty:
            df = other
        else:
            index_names_before = df.index.names

            df = df.join(other, how='outer')

            if index_names_before != df.index.names:
                df.index.names = index_names_before

    return df


def list_entities(instance, entity_type):
    """ Return list of sets, params, variables, constraints or objectives

    Args:
        instance: a Pyomo ConcreteModel object
        entity_type: "set", "par", "var", "con" or "obj"

    Returns:
        DataFrame of entities

    Example:
        >>> data = read_excel('data-example.xlsx')
        >>> model = create_model(data, range(1,25))
        >>> list_entities(model, 'obj')  #doctest: +NORMALIZE_WHITESPACE
                                         Description Domain
        Name
        obj   minimize(cost = sum of all cost types)     []
        [1 rows x 2 columns]

    """

    iter_entities = instance.__dict__.iteritems()

    if entity_type == 'set':
        entities = sorted(
            (x, y.doc, get_onset_names(y)) for (x, y) in iter_entities
            if isinstance(y, pyomo.Set) and not y.virtual)

    elif entity_type == 'par':
        entities = sorted(
            (x, y.doc, get_onset_names(y)) for (x, y) in iter_entities
            if isinstance(y, pyomo.Param))

    elif entity_type == 'var':
        entities = sorted(
            (x, y.doc, get_onset_names(y)) for (x, y) in iter_entities
            if isinstance(y, pyomo.Var))

    elif entity_type == 'con':
        entities = sorted(
            (x, y.doc, get_onset_names(y)) for (x, y) in iter_entities
            if isinstance(y, pyomo.Constraint))

    elif entity_type == 'obj':
        entities = sorted(
            (x, y.doc, get_onset_names(y)) for (x, y) in iter_entities
            if isinstance(y, pyomo.Objective))

    else:
        raise ValueError("Unknown parameter entity_type")

    entities = pd.DataFrame(entities,
                            columns=['Name', 'Description', 'Domain'])
    entities.set_index('Name', inplace=True)
    return entities


def get_onset_names(entity):
    """
        Example:
            >>> data = read_excel('data-example.xlsx')
            >>> model = create_model(data, range(1,25))
            >>> get_onset_names(model.e_co_stock)
            ['t', 'sit', 'com', 'com_type']
    """
    # get column titles for entities from domain set names
    labels = []

    if isinstance(entity, pyomo.Set):
        if entity.dimen > 1:
            # N-dimensional set tuples, possibly with nested set tuples within
            if entity.domain:
                domains = entity.domain.set_tuple
            else:
                domains = entity.set_tuple

            for domain_set in domains:
                labels.extend(get_onset_names(domain_set))

        elif entity.dimen == 1:
            if entity.domain:
                # 1D subset; add domain name
                labels.append(entity.domain.name)
            else:
                # unrestricted set; add entity name
                labels.append(entity.name)
        else:
            # no domain, so no labels needed
            pass

    elif isinstance(entity, (pyomo.Param, pyomo.Var, pyomo.Constraint,
                    pyomo.Objective)):
        if entity.dim() > 0 and entity._index:
            labels = get_onset_names(entity._index)
        else:
            # zero dimensions, so no onset labels
            pass

    else:
        raise ValueError("Unknown entity type!")

    return labels


def get_constants(instance):
    """Return summary DataFrames for important variables
    
    Usage:
        costs, cpro, csto, co2 = get_constants(instance)
    
    Args:
        instance: a picus model instance
        
    Returns:
        costs, cpro, csto, co2)
    """
    costs = get_entity(instance, 'costs')
    cpro = get_entities(instance, ['cap_pro', 'cap_pro_new'])
    csto = get_entities(instance, ['cap_sto_c', 'cap_sto_c_new',
                                        'cap_sto_p', 'cap_sto_p_new'])

    # co2 timeseries
    co2 = get_entity(instance, 'co2_pro_out')
    co2 = co2.unstack(0).sum(1) # sum co2 emissions over timesteps
    
    # better labels
    cpro.columns = ['Total', 'New']
    csto.columns = ['C Total', 'C New', 'P Total', 'P New']
    co2.name = 'CO2'
    
    return costs, cpro, csto, co2

def get_timeseries(instance, co, timesteps=None):
    """Return DataFrames of all timeseries referring to given commodity
    
    Usage:
        created, consumed, storage = get_timeseries(instance, co)
    
    Args:
        instance: a picus model instance
        co: a commodity
        timesteps: optional list of timesteps, defaults to modelled timesteps
        
    Returns:
        created: timeseries of commodity creation, including stock source
        consumed: timeseries of commodity consumption, including demand
        storage: timeseries of commodity storage (level, stored, retrieved)
    """
    if timesteps is None:
        # default to all simulated timesteps
        timesteps = sorted(get_entity(instance, 'tm').index)
        
    # DEMAND
    # default to zeros if commodity has no demand, otherwise recreate timeseries
    if co not in instance.co_demand:
        demand = pd.Series(0, index=timesteps)
    else:
        demand = instance.demand.loc[timesteps][co] * \
                 instance.commodity.loc[co, 'Demand']['peak']
    demand.name = 'Demand'
    
    # STOCK
    eco = get_entity(instance, 'e_co_stock')['e_co_stock'].unstack()
    try:
        stock = eco.loc[timesteps][co]
    except KeyError:
        stock = pd.Series(0, index=timesteps)
    stock.name = 'Stock'

    # PROCESS
    # group process energies by input/output commodity
    # select all entries of created and consumed desired commodity co
    # and slice to the desired timesteps
    epro = get_entities(instance, ['e_pro_in', 'e_pro_out'])
    epro.index.names = ['tm', 'pro', 'coin', 'cout']
    epro = epro.groupby(level=['tm','coin','cout']).sum()
    try:
        created = epro.xs(co, level='cout')['e_pro_out'].unstack()
        created = created.loc[timesteps]
    except KeyError:
        created = pd.DataFrame(index=timesteps)
        
    try:
        consumed = epro.xs(co, level='coin')['e_pro_in'].unstack()
        consumed = consumed.loc[timesteps]
    except KeyError:
        consumed = pd.DataFrame(index=timesteps)
    
    # remove Slack if zero, keep else
    if 'Slack' in created.columns and not created['Slack'].any():
        created.pop('Slack')
    
    # STORAGE
    # group storage energies by commodity
    # select all entries with desired commodity co
    esto = get_entities(instance, ['e_sto_con', 'e_sto_in', 'e_sto_out'])
    esto = esto.groupby(level=['t','co']).sum()
    try:
        storage = esto.xs(co, level='co')
        storage = storage.loc[timesteps]
        storage.columns=['Level', 'Stored', 'Retrieved']
    except KeyError:
        storage = pd.DataFrame(0, index=timesteps,
                               columns=['Level', 'Stored', 'Retrieved'])
 
    # show stock as created
    created = created.join(stock)
    
    # show demand as consumed
    consumed = consumed.join(demand)
    
    return created, consumed, storage

def report(instance, filename, commodities):
    """Write result summary to a spreadsheet file
    
    Args:
        instance: a picus model instance
        filename: Excel spreadsheet filename, will be overwritten if exists
        commodities: list of commodities for which to create timeseries sheets
        
    Returns:
        Nothing
    """
    # get the data
    costs, cpro, csto, co2 = get_constants(instance)
    
    # write to Excel
    writer = pd.ExcelWriter(filename)

    # write to excel
    costs.to_excel(writer, 'Costs')
    co2.to_frame('CO2').to_excel(writer, 'CO2')
    cpro.to_excel(writer, 'Process caps')
    csto.to_excel(writer, 'Storage caps')
    energies = []
    timeseries = {}

    # timeseries 
    for co in commodities:
        created, consumed, storage = get_timeseries(instance, co)
        
        overprod = pd.DataFrame(columns=['Overproduction'], data=
                   created.sum(axis=1) - consumed.sum(axis=1) + 
                   storage['Retrieved'] - storage['Stored'])

        tableau = pd.concat([created, consumed, storage, overprod], axis=1,
                            keys=['Created', 'Consumed', 'Storage', 'Balance'])
        timeseries[co] = tableau.copy()
        
        # timeseries sums
        sums = pd.concat([created.sum(), 
                          consumed.sum(), 
                          storage.sum().drop('Level'), 
                          overprod.sum()], axis=0, 
                         keys=['Created', 'Consumed', 'Storage', 'Balance'])
        energies.append(sums.to_frame(co))
        
    # concatenate the timeseries sums
    energy = pd.concat(energies, axis=1).fillna(0)
    energy.to_excel(writer, 'Energy sums')
    
    for co in commodities:
        timeseries[co].to_excel(writer, "{} timeseries".format(co))
    
    writer.save()

def plot(instance, co, timesteps=None):
    """Stacked timeseries of commodity balance 

    Creates a stackplot of the energy balance of a given commodity, together
    with stored energy in a second subplot. 

    Args:
        instance: a picus model instance
        co: (output) commodity to plot
        timesteps: optional list of modelled timesteps to plot 
                   (e.g. range(1,197)), defaults to set instance.tm
        
    Returns:
        fig: figure handle
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    
    if timesteps is None:
        # default to all simulated timesteps
        timesteps = sorted(get_entity(instance, 'tm').index)

    # FIGURE
    fig = plt.figure(figsize=(16,8))
    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[2, 1])
    
    created, consumed, storage = get_timeseries(instance, co, timesteps)

    # move retrieved/stored storage timeseries to created/consumed and
    # rename storage columns back to 'storage' for color mapping
    created = created.join(storage['Retrieved'])
    consumed = consumed.join(storage['Stored'])
    created.rename(columns={'Retrieved':'Storage'}, inplace=True)
    consumed.rename(columns={'Stored':'Storage'}, inplace=True)
    
    # only keep storage content in storage timeseries
    storage = storage['Level']
    
    # move demand to its own plot
    demand = consumed.pop('Demand')    
    
    # remove all columns from created which are all-zeros in both created and 
    # consumed (except the last one, to prevent a completely empty frame)
    for col in created.columns:
        if not created[col].any() and len(created.columns) > 1:
            if not col in consumed.columns or not consumed[col].any():
                created.pop(col)
    
    # PLOT CREATED
    ax0 = plt.subplot(gs[0])
    sp0 = ax0.stackplot(created.index, created.as_matrix().T, linewidth=0.15)
    
    # Unfortunately, stackplot does not support multi-colored legends by itself.
    # Therefore, a so-called proxy artist - invisible objects that have the 
    # correct color for the legend entry - must be created. Here, Rectangle 
    # objects of size (0,0) are used. The technique is explained at 
    # http://stackoverflow.com/a/22984060/2375855
    proxy_artists = []
    for k, commodity in enumerate(created.columns):
        this_color = to_color(commodity)
            
        sp0[k].set_facecolor(this_color)
        sp0[k].set_edgecolor(to_color('Decoration'))
        
        proxy_artists.append(mpl.patches.Rectangle((0,0), 0,0, 
                                                   facecolor=this_color))
    
    # label
    ax0.set_title('Energy balance of commodity "{}"'.format(co))
    ax0.set_ylabel('Power (kW)')
    
    # legend
    lg = ax0.legend(reversed(proxy_artists), 
                    reversed(tuple(created.columns)),
                    frameon=False,
                    ncol=created.shape[1])
    plt.setp(lg.get_patches(), edgecolor=to_color('Decoration'), linewidth=0.15)
    plt.setp(ax0.get_xticklabels(), visible=False)
    
    # PLOT CONSUMED
    sp00 = ax0.stackplot(consumed.index, -consumed.as_matrix().T, linewidth=0.15)
    
    # color
    for k, commodity in enumerate(consumed.columns):
        this_color = to_color(commodity)
            
        sp00[k].set_facecolor(this_color)
        sp00[k].set_edgecolor(to_color('Decoration'))
    
    # PLOT DEMAND
    dp = ax0.plot(demand.index, demand.values, linewidth=1.2, 
        color=to_color('Demand'))
    
    # PLOT STORAGE
    ax1 = plt.subplot(gs[1], sharex=ax0)
    sp1 = ax1.stackplot(storage.index, storage.values, linewidth=0.15)
    
    # color
    sp1[0].set_facecolor(to_color('Storage'))
    sp1[0].set_edgecolor(to_color('Decoration'))
    
    # labels
    ax1.set_title('Energy storage content of commodity "{}"'.format(co))
    ax1.set_xlabel('Time in year (h)')
    ax1.set_ylabel('Energy (kWh)')
    
    # make xtick distance duration-dependent
    if len(timesteps) > 26*168:
        steps_between_ticks = 168*4
    elif len(timesteps) > 3*168:
        steps_between_ticks = 168
    else:
        steps_between_ticks = 24
    xticks = timesteps[::steps_between_ticks]
    
    # set limits and ticks for both axes
    for ax in [ax0, ax1]:
        #ax.set_axis_bgcolor((0,0,0,0))
        plt.setp(ax.spines.values(), color=to_color('Decoration'))
        ax.set_xlim((timesteps[0], timesteps[-1]))
        ax.set_xticks(xticks)
        ax.xaxis.grid(True, 'major', color=to_color('Decoration'))
        ax.yaxis.grid(True, 'major', color=to_color('Decoration'))
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
    return fig

def to_color(obj=None):
    from random import random
    
    if obj is None:
        obj = random()
    try:
        color = tuple(rgb/255.0 for rgb in COLOURS[obj])
    except KeyError:
        color = "#{:06x}".format(abs(hash(obj)))[:7] # random deterministic color
    return color


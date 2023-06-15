# Python module for the TRNSYS Type calling Python using CFFI
# Data exchange with TRNSYS uses a dictionary, called TRNData in this file (it is the argument of all functions).
# Data for this module will be in a nested dictionary under the module name,
# i.e. if this file is called "MyScript.py", the inputs will be in TRNData["MyScript"]["inputs"]
# for convenience the module name is saved in thisModule
#
# MKu, 2022-02-15

import math
import os

thisModule = os.path.splitext(os.path.basename(__file__))[0]

def calcNu(Re,Pr,charL,l_pipe): # calculate Nusselt numbers
    if Re<2300 and 0.6<Pr<1000:
        nu=lamNuFDev(Re,Pr,charL,l_pipe)
        #print("Re = "+str(round(Re))+" ist laminar!")
    elif Re>pow(10,4) and 0.6<Pr and Pr<1000:
        #print("Re = "+str(round(Re))+" ist turbulent!")
        nu=turbNu(Re,Pr,charL,l_pipe)
    elif 2300<Re and 0.6<Pr<1000: # Übergangsbereich
        #print("Re = "+str(round(Re))+" ist grenzwertig!")
        gamma = (Re-2300)/(pow(10,4)-2300)
        if gamma<0 or gamma>1:
            print("Caution your gamma is out of control")
        #nu = (1-gamma)*lamNuFDev(2300,Pr,charL,l_pipe)+gamma*turbNu(10000,Pr,charL,l_pipe) # laut VDI Wärmeatlas mit therm./strömungstechn. Anlauf
        nu = (1-gamma)*lamNuNDev(2300,Pr,charL,l_pipe)+gamma*turbNu(10000,Pr,charL,l_pipe)
    else:
        print("Please choose a setup with nicer flow properties.")
        # This will cause a crash!
    return nu

def lamNuFDev(Re,Pr,charL,l_pipe): # laminar, flow fully developed, constant wall temperature
    nu = pow((pow(3.66,3)+pow(0.7,3)+pow((1.615*pow((Re*Pr*charL/l_pipe),(-1/3))-0.7),3)),(1/3))
    return nu

def lamNuNDev(Re,Pr,charL,l_pipe): # laminar, fluid + thermal flow not fully developed, constant wall temperature
    nu = pow((pow(lamNuFDev(Re,Pr,charL,l_pipe),3)+pow(0.924*pow(Pr,1/3)*pow(Re*charL/l_pipe,1/2),3)),(1/3))
    return nu

def turbNu(Re,Pr,charL,l_pipe): # turbulent
    xi = pow((1.8*math.log10(Re)-1.5),(-2))
    nu = xi/8*Re*Pr/(1+12.7*math.sqrt(xi)/8*(pow(Pr,(2/3))-1))*(1+pow((charL/l_pipe),(2/3)))
    return nu

# Initialization: function called at TRNSYS initialization
# ----------------------------------------------------------------------------------------------------------------------
def Initialization(TRNData):

    # This model has nothing to initialize
    
    return


# StartTime: function called at TRNSYS starting time (not an actual time step, initial values should be reported)
# ----------------------------------------------------------------------------------------------------------------------
def StartTime(TRNData):

    # Define local short names for convenience (this is optional)
    mflow = TRNData[thisModule]["inputs"][0] # kg/h total fluid mass flow rate
    rho_fluid = TRNData[thisModule]["inputs"][1] # kg/m3 fluid density
    cp_fluid = TRNData[thisModule]["inputs"][2] # kJ/(kg*K) fluid heat capacity
    tc_fluid = TRNData[thisModule]["inputs"][3] # kJ/(h*m*K) fluid thermal conductivity
    visc_dyn_fluid = TRNData[thisModule]["inputs"][4] # kg/(m*h) fluid dynamic viscosity
    di_pipe = TRNData[thisModule]["inputs"][5] # m inner diameter of the pipe
    do_pipe = TRNData[thisModule]["inputs"][6] # m outer diameter of the pipe
    tc_pipe = TRNData[thisModule]["inputs"][7] # kJ/(h*m*K) pipe material thermal conductivity
    l_pipe = TRNData[thisModule]["inputs"][8] # m pipe length per circuit
    res_contact = TRNData[thisModule]["inputs"][9] # m2*K*h/kJ Thermal contact resistance pipe - floor
    nCircuits =  TRNData[thisModule]["inputs"][10] # number of parallel circuits, used to divide mass flow

    # Calculate the outputs
    ri_pipe=di_pipe/2 # m inner radius of the pipe
    d_pipe=(do_pipe-di_pipe)/2 # m thickness of the pipe
    mflow_SI = max(0.0001,mflow/(3600*nCircuits)) # kg/s fluid mass flow rate per circuit
    cp_fluid_SI = cp_fluid*1000 # J/(kg*K) fluid heat capacity
    tc_fluid_SI = tc_fluid/3.6 # W/(m*K) fluid thermal conductivity
    visc_dyn_fluid_SI = visc_dyn_fluid/3600 # kg/(m*s) fluid dynamic viscosity
    tc_pipe_SI = tc_pipe/3.6 # W/(m*K) pipe material thermal conductivity
    res_contact_SI = res_contact*3.6 # m2*K/W Thermal contact resistance pipe - floor
    perimeter = 2*math.pi*ri_pipe # m perimeter of the pipe
    a_crossSec = math.pi*pow(ri_pipe,2) # m2 cross section of the pipe
    ro_pipe = ri_pipe+d_pipe # m outer radius of the pipe
    charL = 4*a_crossSec/perimeter # m characteristic length of the problem
    hx_area = perimeter*l_pipe # shell area of the pipe
    # fluid properties
    visc_kin = visc_dyn_fluid_SI/rho_fluid # m2/s kinematic viscosity of the fluid
    cmin = mflow_SI*cp_fluid_SI  # J/(K*s) minimum heat capacity flow of the two media
    a = tc_fluid_SI/(rho_fluid*cp_fluid_SI) # m2/s thermal diffusivity of the fluid
    # flow characteristics
    vinf_fluid = mflow_SI/(a_crossSec*rho_fluid) #m/s fluid velocity
    Re = vinf_fluid*charL/visc_kin # Reynolds number
    Pr = visc_kin/a # Prandtl number
    nu = calcNu(Re,Pr,charL,l_pipe) # Nusselt number
    alpha = nu*tc_fluid_SI/charL # W/(m2*K) heat transfer coefficient fluid/pipe
    # results
    res_alpha = 1/(alpha*hx_area) # K/W
    res_tc = math.log(ro_pipe/ri_pipe)/(tc_pipe_SI*2*math.pi*l_pipe) # K/W
    res_tot = res_alpha + res_tc + res_contact_SI/hx_area # K/W overall thermal resistance
    k = 1/(res_tot*hx_area) # W/(m2*K) thermal transmittance (U value) related to inner radius
    ntu = k*hx_area/cmin # Number of transfer units
    effectiveness = min(0.9999,max(0.0001,1-math.exp(-ntu))) # heat exchanger effectiveness

    # Set outputs in TRNData
    TRNData[thisModule]["outputs"][0] = effectiveness # heat exchanger effectiveness
    TRNData[thisModule]["outputs"][1] = ntu # number of transfer units
    TRNData[thisModule]["outputs"][2] = alpha # W/(m2*K) heat transfer coefficient fluid-pipe
    TRNData[thisModule]["outputs"][3] = vinf_fluid # m/s fluid velocity far from the wall
    TRNData[thisModule]["outputs"][4] = nu # Nusselt number
    TRNData[thisModule]["outputs"][5] = Pr # Prandtl number
    TRNData[thisModule]["outputs"][6] = Re # Reynolds number

    return


# Iteration: function called at each TRNSYS iteration within a time step
# ----------------------------------------------------------------------------------------------------------------------
def Iteration(TRNData):

     # Define local short names for convenience (this is optional)
    mflow = TRNData[thisModule]["inputs"][0] # kg/h total fluid mass flow rate
    rho_fluid = TRNData[thisModule]["inputs"][1] # kg/m3 fluid density
    cp_fluid = TRNData[thisModule]["inputs"][2] # kJ/(kg*K) fluid heat capacity
    tc_fluid = TRNData[thisModule]["inputs"][3] # kJ/(h*m*K) fluid thermal conductivity
    visc_dyn_fluid = TRNData[thisModule]["inputs"][4] # kg/(m*h) fluid dynamic viscosity
    di_pipe = TRNData[thisModule]["inputs"][5] # m inner diameter of the pipe
    do_pipe = TRNData[thisModule]["inputs"][6] # m outer diameter of the pipe
    tc_pipe = TRNData[thisModule]["inputs"][7] # kJ/(h*m*K) pipe material thermal conductivity
    l_pipe = TRNData[thisModule]["inputs"][8] # m pipe length per circuit
    res_contact = TRNData[thisModule]["inputs"][9] # m2*h*K/kJ Thermal contact resistance pipe - floor
    nCircuits =  TRNData[thisModule]["inputs"][10] # number of parallel circuits, used to divide mass flow

    # Calculate the outputs
    ri_pipe=di_pipe/2 # m inner radius of the pipe
    d_pipe=(do_pipe-di_pipe)/2 # m thickness of the pipe
    ro_pipe = ri_pipe+d_pipe # m outer radius of the pipe
    mflow_SI = max(0.0001,mflow/(3600*nCircuits)) # kg/s fluid mass flow rate per circuit
    cp_fluid_SI = cp_fluid*1000 # J/(kg*K) fluid heat capacity
    tc_fluid_SI = tc_fluid/3.6 # W/(m*K) fluid thermal conductivity
    visc_dyn_fluid_SI = visc_dyn_fluid/3600 # kg/(m*s) fluid dynamic viscosity
    tc_pipe_SI = tc_pipe/3.6 # W/(m*K) pipe material thermal conductivity
    res_contact_SI = res_contact*3.6 # m2*K/W Thermal contact resistance pipe - floor
    perimeter_i = 2*math.pi*ri_pipe # m inner perimeter of the pipe
    a_crossSec = math.pi*pow(ri_pipe,2) # m2 cross section of the pipe
    charL = 4*a_crossSec/perimeter_i # m characteristic length of the problem
    hx_area_i = perimeter_i*l_pipe # shell area of the pipe
    # fluid properties
    visc_kin = visc_dyn_fluid_SI/rho_fluid # m2/s kinematic viscosity of the fluid
    cmin = mflow_SI*cp_fluid_SI # J/(K*s) minimum heat capacity flow of the two media
    a = tc_fluid_SI/(rho_fluid*cp_fluid_SI) # m2/s thermal diffusivity of the fluid
    # flow characteristics
    vinf_fluid = mflow_SI/(a_crossSec*rho_fluid) #m/s fluid velocity
    Re = vinf_fluid*charL/visc_kin # Reynolds number
    Pr = visc_kin/a # Prandtl number
    nu = calcNu(Re,Pr,charL,l_pipe) # Nusselt number
    alpha = nu*tc_fluid_SI/charL # W/(m2*K) heat transfer coefficient fluid/pipe
    # results
    res_alpha = 1/(alpha*hx_area_i) # K/W
    res_tc = math.log(ro_pipe/ri_pipe)/(tc_pipe_SI*2*math.pi*l_pipe) # K/W
    res_tot = res_alpha + res_tc + res_contact_SI/hx_area_i # K/W overall thermal resistance
    k = 1 / (res_tot*hx_area_i) # W/(m2*K) thermal transmittance (U value) related to inner radius
    ntu = k*hx_area_i/cmin # Number of transfer units
    effectiveness = min(0.9999,max(0.0001,1-math.exp(-ntu))) # heat exchanger effectiveness

    # Set outputs in TRNData
    TRNData[thisModule]["outputs"][0] = effectiveness # heat exchanger effectiveness
    TRNData[thisModule]["outputs"][1] = ntu # number of transfer units
    TRNData[thisModule]["outputs"][2] = alpha # W/(m2*K) heat transfer coefficient fluid/pipe
    TRNData[thisModule]["outputs"][3] = vinf_fluid # m/s fluid velocity far from the wall
    TRNData[thisModule]["outputs"][4] = nu # Nusselt number
    TRNData[thisModule]["outputs"][5] = Pr # Prandtl number
    TRNData[thisModule]["outputs"][6] = Re # Reynolds number

    return

# EndOfTimeStep: function called at the end of each time step, after iteration and before moving on to next time step
# ----------------------------------------------------------------------------------------------------------------------
def EndOfTimeStep(TRNData):

    # This model has nothing to do during the end-of-step call
    
    return


# LastCallOfSimulation: function called at the end of the simulation (once) - outputs are meaningless at this call
# ----------------------------------------------------------------------------------------------------------------------
def LastCallOfSimulation(TRNData):

    # NOTE: TRNSYS performs this call AFTER the executable (the online plotter if there is one) is closed. 
    # Python errors in this function will be difficult (or impossible) to diagnose as they will produce no message.
    # A recommended alternative for "end of simulation" actions it to implement them in the EndOfTimeStep() part, 
    # within a condition that the last time step has been reached.
    #
    # Example (to be placed in EndOfTimeStep()):
    #
    # stepNo = TRNData[thisModule]["current time step number"]
    # nSteps = TRNData[thisModule]["total number of time steps"]
    # if stepNo == nSteps-1:     # Remember: TRNSYS steps go from 0 to (number of steps - 1)
    #     do stuff that needs to be done only at the end of simulation

    return
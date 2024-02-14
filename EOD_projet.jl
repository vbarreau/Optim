######################################################################### 
#                       PROJET OPTIM EOD                                #
#                ELOI BUFFET   -  VICTOR BARREAU                        #
#                     MIX 100% RENOUVELABLES                            #
#########################################################################

#packages
using JuMP
#use the solver you want
using HiGHS
#package to read excel files
using XLSX

Tmax = 168 #optimization for 1 year (7*24*52 hours)
data_file = "Donnees_etude_de_cas_ETE305.xlsx"
#data for load and fatal generation
load = XLSX.readdata(data_file, "Consomation", "C4:C171")
wind = XLSX.readdata(data_file, "Consomation", "D4:D171")
solar = XLSX.readdata(data_file, "Consomation", "E4:E171")
hydro_fatal = XLSX.readdata(data_file, "Consomation", "F4:F171")
thermal_fatal = XLSX.readdata(data_file, "Consomation", "G4:G171")

Nsolref = 4
Nwref = 4

NmaxSol  = 20
NmaxW    = 20
# #total of RES
# Pres = wind + solar + hydro_fatal + thermal_fatal

#data for thermal clusters
# Nth = 5 #number of thermal generation units
# names = XLSX.readdata(data_file, "data", "J4:J8")
# dict_th = Dict(i=> names[i] for i in 1:Nth)
# costs_th = XLSX.readdata(data_file, "data", "K4:K8")
# Pmin_th = XLSX.readdata(data_file, "data", "M4:M8") #MW
# Pmax_th = XLSX.readdata(data_file, "data", "L4:L8") #MW
# dmin = XLSX.readdata(data_file, "data", "N4:N8") #hours
# ramp_th = XLSX.readdata(data_file, "data", "O4:O8") #MW

#data for hydro reservoir
# Nhy = 1 #number of hydro generation units
NmaxHy = 5
Pmin_hy = zeros(Nhy)
Pmax_hy = XLSX.readdata(data_file, "data", "R4") *ones(Nhy) #MW
e_hy = XLSX.readdata(data_file, "data", "S4")*ones(Nhy) #MWh
costs_hy = XLSX.readdata(data_file, "data", "Q4")*ones(Nhy) #MWh

#costs
cth = repeat(costs_th', Tmax) #cost of thermal generation €/MWh
chy = repeat(costs_hy', Tmax) #cost of hydro generation €/MWh
cuns = 5000*ones(Tmax) #cost of unsupplied energy €/MWh
cexc = 0*ones(Tmax) #cost of in excess energy €/MWh


#data for STEP/battery
#weekly STEP
Pmax_STEP = XLSX.readdata(data_file, "data", "R5") #MW
rSTEP = XLSX.readdata(data_file, "data", "T5")
stock_step_max = XLSX.readdata(data_file, "data", "S5")

#battery
Pmax_battery = 280 #MW
rbattery = 0.85
d_battery = 2 #hours


#############################
#create the optimization model
#############################
model = Model(HiGHS.Optimizer)

#############################
#define the variables
#############################
# #thermal generation variables
# @variable(model, Pth[1:Tmax,1:Nth] >= 0)
# @variable(model, UCth[1:Tmax,1:Nth], Bin)
# @variable(model, UPth[1:Tmax,1:Nth], Bin)
# @variable(model, DOth[1:Tmax,1:Nth], Bin)

# installation variable 
@variable(model, Nsol>=0)
@variable(model, Nw>=0)
@variable(model, Nhy>=1)

#hydro generation variables
@variable(model, Phy[1:Tmax,1:Nhy] >= 0)   # Pas sur que ça soit linéaire ça !

#unsupplied energy variables
@variable(model, Puns[1:Tmax] >= 0)
#in excess energy variables
@variable(model, Pexc[1:Tmax] >= 0)

#weekly STEP variables
@variable(model, Nstep>=0)
@variable(model, Pcharge_STEP[1:Tmax,1:Nstep] >= 0)
@variable(model, Pdecharge_STEP[1:Tmax,1:Nstep] >= 0)
@variable(model, stock_STEP[1:Tmax,1:Nstep] >= 0)

#battery variables
@variable(model, Nbat>=0)
@variable(model, Pcharge_battery[1:Tmax] >= 0)
@variable(model, Pdecharge_battery[1:Tmax] >= 0)
@variable(model, stock_battery[1:Tmax] >= 0)


# #############################
#define the objective function
#############################
@objective(model, Min, sum(Phy.*chy)+Puns'cuns+Pexc'cexc)

#############################
#define the constraints
#############################
#balance constrain
# @constraint(model, balance[t in 1:Tmax], sum(Pth[t,g] for g in 1:Nth) + sum(Phy[t,h] for h in 1:Nhy) + Pres[t] + Puns[t] - load[t] - Pexc[t] - Pcharge_STEP[t] + Pdecharge_STEP[t] - Pcharge_battery[t] + Pdecharge_battery[t] == 0)
@constraint(model, balance[t in 1:Tmax], sum(Phy[t,h] for h in 1:Nhy) + thermal_fatal[t] + Phy[t].*Nhy + solar[t].*Nsol/Nsolref + wind[t].*Nw/Nwref + Puns[t] - load[t] - Pexc[t] - Pcharge_STEP[t] + Pdecharge_STEP[t] - Pcharge_battery[t] + Pdecharge_battery[t] == 0)

# #thermal unit Pmax constraints
# @constraint(model, max_th[t in 1:Tmax, g in 1:Nth], Pth[t,g] <= Pmax_th[g]*UCth[t,g])

# #thermal unit Pmin constraints
# @constraint(model, min_th[t in 1:Tmax, g in 1:Nth], Pmin_th[g]*UCth[t,g] <= Pth[t,g])

# #thermal unit Dmin constraints
# for g in 1:Nth
#         if (dmin[g] > 1)
#             @constraint(model, [t in 2:Tmax], UCth[t,g]-UCth[t-1,g]==UPth[t,g]-DOth[t,g],  base_name = "fct_th_$g")
#             @constraint(model, [t in 1:Tmax], UPth[t]+DOth[t]<=1,  base_name = "UPDOth_$g")
#             @constraint(model, UPth[1,g]==0,  base_name = "iniUPth_$g")
#             @constraint(model, DOth[1,g]==0,  base_name = "iniDOth_$g")
#             @constraint(model, [t in dmin[g]:Tmax], UCth[t,g] >= sum(UPth[i,g] for i in (t-dmin[g]+1):t),  base_name = "dminUPth_$g")
#             @constraint(model, [t in dmin[g]:Tmax], UCth[t,g] <= 1 - sum(DOth[i,g] for i in (t-dmin[g]+1):t),  base_name = "dminDOth_$g")
#             @constraint(model, [t in 1:dmin[g]-1], UCth[t,g] >= sum(UPth[i,g] for i in 1:t), base_name = "dminUPth_$(g)_init")
#             @constraint(model, [t in 1:dmin[g]-1], UCth[t,g] <= 1-sum(DOth[i,g] for i in 1:t), base_name = "dminDOth_$(g)_init")
#     end
# end

# Ramp constraint
#@constraint(model, rampupth[t in 1:Tmax-1, g in 1:Nth], (Pth[t+1,g] - Pth[t,g]) .- ramp_th[g]<= 0 )
# @constraint(model, rampupdo[t in 1:Tmax-1, g in 1:Nth], -Pth[t+1,g] + Pth[t]<=ramp_th[g] )


# hydro unit constraints
@constraint(model, Nsol)
@constraint(model, bounds_hy[t in 1:Tmax, h in 1:Nhy], Pmin_hy[h] <= Phy[t,h] <= Pmax_hy[h])

# hydro stock constraint
@constraint(model,StockHy[t in 1:Tmax],sum(Phy[1:t]) .- e_hy <= 0)
@constraint(model,MaxPhy[t in 1:Tmax], Phy[t] .- Pmax_hy     <= 0)

# weekly STEP
@constraint(model,Cstep1[t in 1:Tmax-1],stock_STEP[t]+Pcharge_STEP[t]*rSTEP-Pdecharge_STEP[t]== stock_STEP[t+1] )
@constraint(model,Cstep2[t in 1:Tmax],Pcharge_STEP[t]<=Pmax_STEP)
@constraint(model,Cstep3[t in 1:Tmax],Pdecharge_STEP[t]<=Pmax_STEP)
@constraint(model,Cstep4[t in 1:Tmax], stock_STEP[t]<=stock_step_max)

# battery
@constraint(model,bat_bilan[t in 1:Tmax-1], stock_battery[t]+Pcharge_battery[t]*rbattery-Pdecharge_battery[t]/rbattery==stock_battery[t+1])
@constraint(model,bat_emax[t in 1:Tmax],stock_battery[t]<=d_battery*Pmax_battery)
@constraint(model,bat_Pmax1[t in 1:Tmax],Pcharge_battery[t]<=Pmax_battery)
@constraint(model,bat_Pmax2[t in 1:Tmax],Pdecharge_battery[t]<=Pmax_battery)





# no need to print the model when it is too big
# solve the model
optimize!(model)
#------------------------------
# Results
@show termination_status(model)
@show objective_value(model)


#exports results as csv file
hy_gen = value.(Phy)
STEP_charge = value.(Pcharge_STEP)
STEP_decharge = value.(Pdecharge_STEP)
battery_charge = value.(Pcharge_battery)
battery_decharge = value.(Pdecharge_battery)


# new file created
touch("results.csv")

# file handling in write mode
f = open("results.csv", "w")

for name in names
    write(f, "$name ;")
end
write(f, "Hydro ; STEP pompage ; STEP turbinage ; Batterie injection ; Batterie soutirage ; RES ; load ; Net load \n")

for t in 1:Tmax
    for g in 1:Nth
        write(f, "$(th_gen[t,g]) ; ")
    end
    for h in 1:Nhy
        write(f, "$(hy_gen[t,h]) ;")
    end
    write(f, "$(STEP_charge[t]) ; $(STEP_decharge[t]) ;")
    write(f, "$(battery_charge[t]) ; $(battery_decharge[t]) ;")
    write(f, "$(Pres[t]) ;  $(load[t]) ; $(load[t]-Pres[t]) \n")

end

close(f)

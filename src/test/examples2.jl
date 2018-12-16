using Pkg
Pkg.activate(".")
using Todjl

turb=Turbine("turbine")
turb.mode=:design
setattrs(turb;pr=0.1,eta=0.9)
cond=Condensor("cond")
pump=Pump("pump")
setattrs(pump;prise=1.,eta=0.8)
heater=GeneralHeater("gh")
setattrs(heater,temp=300)

c_heater_turb=connect(heater,:out,turb,:in)
setattrs(c_heater_turb;m=1,p=1)
c_pump_heater=connect(pump,:out,heater,:in)
c_cond_pump=connect(cond,:waterout,pump,:in)
c_turb_cond=connect(turb,:out,cond,:steam)

nw=Network()
addcomponents(nw,turb,cond,pump,heater)
addconnections(nw,c_heater_turb,c_pump_heater,c_cond_pump,c_turb_cond)

initconnection(nw)

vec,mat=assemble(nw)

for c in nw.conns
    show(c)
end


@show vec
show(mat)
@show length(vec)
@show size(mat)

inv(mat)




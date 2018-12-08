using Revise
using Todjl

t0=Turbine("hp_throttle")
t1=Turbine("hp1")
t2=Turbine("hp2")
t3=Turbine("ip1")
t4=Turbine("ip2")
t5=Turbine("lp1")
t6=Turbine("lp2")
t7=Turbine("lp3")
t8=Turbine("lp4")
te=Turbine("exhaust")

sp1=Split("hp1after")
sp2=Split("hp2after")
sp3=Split("ip1after")
sp4=Split("ip2after")
sp5=Split("lp1after")
sp6=Split("lp2after")
sp7=Split("lp3after")
sp8=Split("lp4after")

fwh1=FeedwaterHeater("hp1heater")
fwh2=FeedwaterHeater("hp2heater")
fwh3=FeedwaterHeater("hp3heater")
dae=Dearator("dearator")
fwh5=FeedwaterHeater("lp5hater")
fwh6=FeedwaterHeater("lp6heater")
fwh7=FeedwaterHeater("lp7heater")
fwh8=FeedwaterHeater("lp8heater")
cnd=Condensor("cond")

cwp=Pump("cwp")
fwp=Pump("fwp")

sh=GeneralHeater("super heater")
rh=GeneralHeater("reheater")

c_mainsteam=Connection(sh,:out,t0,:in)
c_t0_t1=Connection(t0,:out,t1,:in)
c_t1_sp1=Connection(t1,:out,sp1,:in)
c_sp1_t2=Connection(sp1,:out1,t2,:in)
c_t2_sp2=Connection(t2,:out1,sp2,:in)
c_sp2_rh=Connection(sp2,out1,rh,:in)
c_rh_t3=Connection(rh,:out,t3,:in)
c_t3_sp3=Connection(t3,:out,sp3,:in)
c_sp3_t4=Connection(sp3,out1,t4,:in)
c_t4_sp4=Connection


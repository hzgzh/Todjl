

pt(p::Float64)=CoolProp.PropsSI("T","P",p*1e5,"Q",0.0,"Water")
pth(p::Float64,t::Float64)=CoolProp.PropsSI("H","P",p*1.e5,"T",t+273.15,"Water")/1000.
pts(p::Float64,t::Float64)=CoolProp.PropsSI("S","P",p*1.e5,"T",t+273.15,"Water")/1000.
ptv(p::Float64,t::Float64)=1/CoolProp.PropsSI("D","P",p*1.e5,"T",t+273.15,"Water")
pht(p::Float64,h::Float64)=CoolProp.PropsSI("T","P",p*1.e5,"H",1000h,"Water")-273.15
phs(p::Float64,h::Float64)=CoolProp.PropsSI("S","P",p*1.e5,"H",1000h,"Water")/1000.
phv(p::Float64,h::Float64)=1/CoolProp.PropsSI("D","P",p*1.e5,"H",1000h,"Water")
phq(p::Float64,h::Float64)=CoolProp.PropsSI("Q","P",p*1.e5,"H",1000h,"Water")
psh(p::Float64,s::Float64)=CoolProp.PropsSI("H","P",p*1.e5,"S",1000s,"Water")/1000.
pst(p::Float64,s::Float64)=CoolProp.PropsSI("T","P",p*1.e5,"S",1000s,"Water")-273.15
psv(p::Float64,s::Float64)=1/CoolProp.PropsSI("D","P",p*1.e5,"S",1000s,"Water")
pqh(p::Float64,q::Float64)=CoolProp.PropsSI("H","P",p*1.e5,"Q",q,"Water")/1000.0
pqs(p::Float64,q::Float64)=CoolProp.PropsSI("S","P",p*1.e5,"Q",q,"Water")/1000.0
tp(t::Float64)=CoolProp.PropsSI("P","T",t+273.15,"Q",0,"Water")/1.e5
thp(t::Float64,h::Float64)=CoolProp.PropsSI("P","T",t+273.15,"H",1000h,"Water")/1.e5
tsp(t::Float64,s::Float64)=CoolProp.PropsSI("P","T",t+273.15,"S",1000s,"Water")/1.e5
tqp(t::Float64,q::Float64)=CoolProp.PropsSI("P","T",t+2731.5,"Q",q,"Water")/1.e5
ths(t::Float64,h::Float64)=CoolProp.PropsSI("S","T",t+273.15,"H",1000h,"Water")/1000.
dpdh(p::Float64,h::Float64)=1/dhdp(p,h)
dhdp(p::Float64,t::Float64)=CoolProp.PropsSI("d(Hmass)/dP|T","P",p*1e5,"T",t+273.15,"Water")*100.
dhdt(p::Float64,t::Float64)=CoolProp.PropsSI("d(Hmass)/dT|P","P",p*1e5,"T",t+273.15,"Water")/1000.
sat_dpdt(p::Float64)=CoolProp.PropsSI("d(P)/d(T)|sigma","P",p*1e5,"Q",1,"Water")/1.e5
sat_dtdp(p::Float64)=CoolProp.PropsSI("d(T)/d(P)|sigma","P",p*1e5,"Q",1,"Water")*1.e5
sat_dhdp(p::Float64,q::Float64)=CoolProp.PropsSI("d(Hmass)/d(P)|sigma","P",p*1e5,"Q",0,"Water")*100.
dtdp(p::Float64,h::Float64)=CoolProp.PropsSI("d(T)/d(P)|Hmass","P",p*1e5,"Hmass",h*1e3,"Water")*1.e5
dtdh(p::Float64,h::Float64)=CoolProp.PropsSI("d(T)/d(Hmass)|P","P",p*1e5,"Hmass",h*1e3,"Water")*1.e3

proprange()=(611.0/1.0e5,1.0e4,0.0,2000-273.16)

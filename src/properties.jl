

pt(p::T) where T<:Real=CoolProp.PropsSI("T","P",p*1.e5,"Q",0.0,"Water")
pth(p::S,t::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("H","P",p*1.e5,"T",t+273.15,"Water")/1000.
pts(p::S,t::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("S","P",p*1.e5,"T",t+273.15,"Water")/1000.
ptv(p::S,t::T) where {S<:Real,T<:Real} =1.0/CoolProp.PropsSI("D","P",p*1.e5,"T",t+273.15,"Water")
pht(p::S,h::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("T","P",p*1.e5,"H",1000h,"Water")-273.15
phs(p::S,h::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("S","P",p*1.e5,"H",1000h,"Water")/1000.
phv(p::S,h::T) where {S<:Real,T<:Real} =1.0/CoolProp.PropsSI("D","P",p*1.e5,"H",1000h,"Water")
phq(p::S,h::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("Q","P",p*1.e5,"H",1000h,"Water")
psh(p::S,s::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("H","P",p*1.e5,"S",1000s,"Water")/1000.
pst(p::S,s::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("T","P",p*1.e5,"S",1000s,"Water")-273.15
psv(p::S,s::T) where {S<:Real,T<:Real} =1.0/CoolProp.PropsSI("D","P",p*1.e5,"S",1000s,"Water")
pqh(p::S,q::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("H","P",p*1.e5,"Q",q,"Water")/1000.0
pqs(p::S,q::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("S","P",p*1.e5,"Q",q,"Water")/1000.0
tp(t::T) where T<:Real=CoolProp.PropsSI("P","T",t+273.15,"Q",0,"Water")/1.e5
thp(t::S,h::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("P","T",t+273.15,"H",1000h,"Water")/1.e5
tsp(t::S,s::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("P","T",t+273.15,"S",1000s,"Water")/1.e5
tqp(t::S,q::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("P","T",t+2731.5,"Q",q,"Water")/1.e5
ths(t::S,h::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("S","T",t+273.15,"H",1000h,"Water")/1000.
dpdh(p::S,h::T) where {S<:Real,T<:Real} =1.0/dhdp(p,h)
dhdp(p::S,t::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("d(Hmass)/d(P)|T","P",p*1e5,"T",t+273.15,"Water")*100.
dhdt(p::S,t::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("d(Hmass)/d(T)|P","P",p*1e5,"T",t+273.15,"Water")/1000.
sat_dpdt(p::T) where T<:Real=CoolProp.PropsSI("d(P)/d(T)|sigma","P",p*1e5,"Q",1,"Water")/1.e5
sat_dtdp(p::T) where T<:Real=CoolProp.PropsSI("d(T)/d(P)|sigma","P",p*1e5,"Q",1,"Water")*1.e5
sat_dhdp(p::T,q::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("d(Hmass)/d(P)|sigma","P",p*1e5,"Q",0,"Water")*100.
dtdp(p::S,h::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("d(T)/d(P)|Hmass","P",p*1e5,"Hmass",h*1e3,"Water")*1.e5
dtdh(p::S,h::T) where {S<:Real,T<:Real} =CoolProp.PropsSI("d(T)/d(Hmass)|P","P",p*1e5,"Hmass",h*1e3,"Water")*1.e3

proprange()=(611.0/1.0e5,1.0e4,0.0,2000-273.16)

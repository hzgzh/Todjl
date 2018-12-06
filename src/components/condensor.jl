export Condensor,setattr,addconnection,euqations,jacobi
mutable struct Condensor
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    calcmode::Symbol
    mode::Symbol
end

function setattr(comp::Condensor;kwarg...)
    opts=[:temp]
    for key in keys(kwargs)
        if key in opts
            attrs[key]=kwargs[key]
        end
    end
end

function addconnection(comp::Condensor,c::Connection,id::String)
    port=[:steam,:drain,:waterout]
    if s in port
        comp.conns[s]=c
    else
        println("wrong input")
    end
end

function equations(comp::Condensor)
    si=comp.conns[:steam];di=comp.conns[:drain];wo=comp.conns[:waterout]
    res=[]
    push!(res,si.m.val+di.m.val-wo.m.val)
    push!(res,wo.h.val-pqh(si.p.val,0))
    push!(res,si.p.val-wo.p.val)
end

conns(comp::Condensor)=[:steam,:drain,:waterout]

function jacobi(comp::Condensor,c::Connection)
    si=comp.conns[:steam];di=comp.conns[:drain];wo=comp.conns[:waterout]
    jac=zeros(3,9)
    if c == si
        jac[1,1]=1.0;jac[2,3]=-sat_dhdp(si.p.val,0.0);jac[3,3]=1.0
    end
    if c == di
        jac[1,4]=1.0;
        
    end
    if c == wo
        jac[1,7]=-1.0;jac[2,8]=1.0;jac[3,9]=-1.0
    end
    res
end
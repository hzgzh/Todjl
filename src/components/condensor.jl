@with_kw mutable struct Condensor <: AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    mode=:design
    method=:fixpress

    Condensor(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end

function (cp::Condensor)(;kwargs)

end

function setattr(comp::Condensor;kwarg...)
    opts=[:temp];attrs=comp.attrs
    for key in keys(kwargs)
        if key in opts
            attrs[key]=kwargs[key]
        end
    end
end

inlets(comp::Condensor)=[comp.conns[:steam] comp.conns[:drain]]
outlets(comp::Condensor)=[comp.conns[:waterout]]

function addconnection(comp::Condensor,port::Symbol,c::Connection)
    ports=[:steam,:drain,:waterout]
    if port in ports
        comp.conns[port]=c
    else
        println("wrong input")
    end
end

function equations(comp::Condensor)
    si=comp.conns[:steam];wo=comp.conns[:waterout]
    res=zeros(0)
    if haskey(comp.conns,:drain)
        di=comp.conns[:drain]
        push!(res,si.m.val+di.m.val-wo.m.val)
    else
        push!(res,si.m.val-wo.m.val)
    end
    push!(res,pqh(si.p.val,0)-wo.h.val)
    push!(res,si.p.val-wo.p.val)
    return res
end

conns(comp::Condensor)=[:steam,:drain,:waterout]

function jacobi(comp::Condensor,c::Connection)
    si=comp.conns[:steam];wo=comp.conns[:waterout];
    di=haskey(comp.conns,:drain) ? comp.conns[:drain] : nothing
    
    jac=zeros(3,3)
    if c == si
        jac[1,1]=1.0;jac[2,2]=eps();jac[2,3]=sat_dhdp(si.p.val,0.0);jac[3,3]=1.0
    end
    if c == di
        jac[1,1]=1.0;
    end
    if c == wo
        jac[1,1]=-1.0;jac[2,2]=-1.0;jac[3,3]=-1.0
    end
    return jac
end

export Condensor,equations,jacobi,setattrs,addconnection
@with_kw mutable struct Source <: AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    mode=:design
    method=:no
    Source(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end

function (cp::Source)(;kwargs)
    
end

function setattrs(comp::Source;kwargs)
    opts=[:m,:h,:T];attrs=comp.attrs
    for key in keys(kwargs)
        if key in opts
            attrs[key]=kwargs[key]
        end
    end
    
    if length(attrs)!=3
        print("wrong input attr")
    end
end

inlets(comp::Source)=[];outlets(comp::Source)=[comp.conns[:out]]

function addconnection(comp::Source,port::Symbol,c::Connection)
    ports=[:out]
    if port in ports
        comp.conns[port]=c
    else
        print("wrong connection name")
    end
end

function equations(comp::Source)
     o.m.val_set=attrs[:m];o.h.val_set=attrs[:h];o.p.val_set=attrs[:p]
     o.m.is_set=true;o.h.is_set=true;o.p.is_set=true
    []
end

function jacobi(comp::Source,c::Connection)
    []
end


export Source,equations,jacobi,setattrs,addconnection
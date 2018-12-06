export Source
mutable struct Source <: Component
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
end
 
function setattrs(comp::Source;kwargs)
    opts=[:m,:h,:T]
    for key in keys(kwargs)
        if key in opts
            attrs[key]=kwargs[key]
        end
    end
    
    if length(attrs)!=3
        print("wrong input attr")
    end
end

function addconnection(comp::Source,c::Connection,si::Symbol)
    port=[:out]
    if si in port
        comp.conns[si]=c
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
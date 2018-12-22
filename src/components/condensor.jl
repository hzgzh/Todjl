@with_kw mutable struct Condensor <: AbstractHeater
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,VarProp}
    
    Condensor(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end

designmode(cp::FeedwaterHeater)=[:tl,:pr2]
offdesignmode(cp::FeedwaterHeater)=[]

function additional_equations(cp::Condensor)
    return pqh(hi.p.val,0)-ho.h.val
end

function additional_derivatives(cp::Condensor)
    
    e_der=zeros(1,len,3)
    if haskey(cp.conns,:drain)
        e_der[1,1,3]=dhdp(hi.p.val,0);e_der[1,4,2]=-1
    else
        e_der[1,1,3]=dhdp(hi,p.val,0);e_der[1,3,2]=-1
    end
    return e_der
end

export Condensor,equations,jacobi,setattrs,addconnection
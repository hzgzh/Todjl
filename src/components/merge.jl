@with_kw mutable struct Merge <: AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    mode=:design
    method=:no
    Merge(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end


inlets(cp::Merge)=[cp.conns[:in1] cp.conns[:in2]]
outlets(cp::Merge)=[cp.conns[:out1]]
ports(cp::Merge)=[:in1,:in2,:out]
attrs(cp::Merge)=[]

function equations(cp::Merge)
    in1=cp.conns[:in1];in2=cp.conns[:in2];o=cp.conns[:out]
    vec=[]
    vec+=mass_res(cp)
    res=0.
    for (idx,i) in inlets(cp)
        res+=i.m.val*i.h.val
    end
    for (idx,o) in outlets(cp)
        res-=o.m.val*o.h.val
    end
    vec+=res
    vec+=in1.p.val-o.p.val
    vec+=in2.p.val-o.p.val
    return vec
end

function derivatives(cp::Merge)
    in1=cp.in1;in2=cp.in2;o=cp.out1
    der=mass_deriv(cp)

    e_der=zeros(1,3,3)
    e_der[1,1,1]=in1.h.val;e_der[1,1,2]=in1.m.val
    e_der[1,2,1]=in2.h.val;e_der[1,2,2]=in1.m.val
    e_der[1,3,1]=-o.h.val;e_der[1,3,2]=-o.m.val
    der+=e_der

    p_der=zeros(2,3,3)
    p_der[1,1,3]=1;p_der[1,3,3]=-1
    p_der[2,1,3]=1;p_der[2,3,3]=-1
    der+=p_der
    return der
end

initsource(cp::Merge,c::Connection)=[1.0,500,1.0]
inittarget(cp::Merge,c::Connection)=[1.0,500,1.0]


export Merge,equations,jacobi,setattrs,addconnection







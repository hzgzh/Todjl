@with_kw mutable struct Split <: AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,AbstractVariable}
    mode=:design
    method=:no
end

inlets(cp::Split)=[cp.conns[:in]]
outlets(cp::Split)=[cp.conns[:out1] cp.conns[:out2]]
portnames(cp::Split)=[:in,:out1,:out2]

function equations(cp::Split)
    in=cp.conns[:in];o1=cp.conns[:out1];o2=cp.conns[:out2]
    vec=[]
    vec+=mass_res(cp)
    e_res1=in.h.val-o1.h.val
    vec+=e_res1
    e_res2=in.h.val-o2.h.val
    vec+=e_res2
    p_res1=in.p.val-o1.p.val
    vec+=p_res1
    p_res2=in.p.val-o2.p.val
    vec+=p_res2
    return vecc
end

function derivatives(cp::Split)
    m_der=mass_deriv(cp)
    e_der=zeros(2,3,3)
    e_der[1,1,2]=1;e_der[1,2,2]=-1
    e_der[2,1,2]=1;e_der[2,2,2]=-1
    p_der=zeros(2,3,3)
    p_der[1,1,3]=1;p_der[1,2,3]=-1
    p_der[2,1,3]=1;p_der[2,2,3]=-1
    return vcat(m_der,e_der,p_der)
end

initsource(cp::Split,c::Connection)=[1.0,500.,1.0]
inittarget(cp::Split,c::Connection)=[1.0,500.,1.0]

export Split,equations,jacobi,setattrs,addconnection
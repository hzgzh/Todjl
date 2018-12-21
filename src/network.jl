function addcomponents(nw::Network,comps...)
    for comp in comps
        push!(nw.comps,comp)
    end
end
function addconnections(nw::Network,conns...)
    for conn in conns
        push!(nw.conns,conn)
    end
end

function addbuses(nw::Network,buses...)
    for bus in buses
        push!(nw.buses,bus)
    end
end

function initconnection(nw::Network)
    for c in nw.conns
        ms,hs,ps=initsource(c.source,c)
        mt,ht,pt=inittarget(c.target,c)
        if ps==0.0 && pt==0.0
            c.p.val=1
        elseif ps==0.0
            c.p.val=pt
        elseif pt==0.0
            c.p.val=ps
        else
            c.p.val=(ps+pt)/2.0
        end

        if hs==0.0 && ht==0.0
            c.h.val=1000.0
        elseif hs==0.0
            c.h.val=ht
        elseif ht==0.0
            c.h.val=hs
        else 
            c.h.val=(hs+ht)/2.0
        end
    end
end

function solvecomponent(nw::Network)
    comps=nw.comps;conns=nw.conns
    numvars=3*length(conns)
    mat=zeros(0,numvars)
    vec=zeros(0)
    for cp in comps
        eqs=equations(comp)
        vec=[vec;eqs]
        der=zeros(length(eqs),numvars)
        i=1
        for conn in [inlets(cp) outlets(cp)]
            jac=derivatives(cp)
            pos=find(conns,conn)
            der[1:length(eqs),3*(pos-1)+1:3*pos]=jac[:,i,:]
            i+=1
        end
        mat=vcat(mat,derive)
    end
    (vec,mat)
end

function solve_connection(nw::Network)
    res=[];
    der=zeros(0,3*length(nw.conns))
    i=1
    for conn in nw.conns
        eqs=solve_prop_eq(cp)
        res+=eqs
        der1=zeros(length(eqs),3*length(nw.conns))
        d=solve_prop_der(cp)
        der1[1:length(eqs),3i+1:3(i+1)]=d
        der=vcat(der,der1)

        eqs=solve_prop_ref_eq(cp)
        res+=eqs
        der1=zeros(length(eqs),3*length(nw.conns))
        d=solve_prop_ref_der(cp)
        der1[1:length(eqs),3i-2:3i)]=d[:,1,:]
        pos=find(nw.conns,conn.ref.conn)
        der1[1:length(eqs),3*pos-2,3*pos]=d[:,2,:]
        der=vcat(der,der1)
        i+=1
    end
    return (res,der)
end

function solve_prop_eq(c::Connection)
    res=[]
    if c.m.val_set
        res+=0
    end
    if c.h.val_set
        res+=0
    end
    if c.p.val_set
        res+=0
    end
    if c.T.val_set
        res+=c.T.val-pht(c.p.val,c.h.val)
    end
end

function solve_prop_der(c::Connection)
    der=zeros(0,3)
    if c.m.val_set
        d=zeros(1,3)
        d[1,1]=1
        der=vcat(der,d)
    end
    if c.h.val_set
        d=zeros(1,3)
        d[1,2]=1
        der=vcat(der,d)
    end
    if c.p.val_set
        d=zeros(1,3)
        d[1,3]=1
        der=vcat(der,d)
    end
    if c.T.val_set
        d=zeros(1,3)
        d[1,2]=-dtdh(c.p.val,c.h.val)
        d[1,3]=-dtdp(c.p.val,c.h.val)
        der=vat(der,d)
    end
    return der
end

function solve_prop_ref_eq(c::Connection)
    res=[]
    if c.m.ref_set
        res+=c.m.val-c.m.ref.conn.m.val*c.m.ref.factor+c.m.ref.delta
    end
    if c.h.ref_set
        res+=c.h.val-c.h.ref.conn.h.val*c.h.ref.factor+c.h.ref.delta
    end
    if c.p.ref_set
        res+=c.p.val-c.p.ref.conn.p.val*c.p.ref.factor+c.p.ref.delta
    end
    if c.T.ref_set
        res+=pht(c.p.val,c.p.h.val)-pht(c.p.ref.conn.p.val,c.h.ref.conn.h.val)*c.T.ref.factor+c.T.ref.delta
    end
    return res
end

function solve_prop_ref_der(c::Connection)
    deriv=zeros(0,2,3)
    if c.m.ref_set
        der=zeros(1,2,3)
        der[1,1,1]=1
        der[1,2,1]=-c.m.ref.factor
        deriv=vcat(deriv,der)
    end
    if c.h.ref_set
        der=zeros(1,2,3)
        der[1,1,2]=1
        der[1,2,2]=-c.h.ref.factor
        deriv=vcat(deriv,der)
    end
    if c.p.ref_set
        der=zeros(1,2,3)
        der[1,1,3]=1
        der[1,2,3]=-c.p.ref.factor
        deriv=vcat(deriv,der)
    end
    if c.T.ref_set
        der=zeros(1,2,3)
        der[1,1,1]=dtdp(c.p.val,c.h.val);der[1,1,2]=dtdh(c.p.val,c.h.val)
        der[1,2,1]=-dtdp(c.m.ref.p.val,c.m.ref.h.val);der[1,2,2]=-dtdh(c.m.ref.p.val,c.m.ref.h.val)
        deriv=vcat(deriv,der)
    end
    return deriv
end

function find(conns::Vector{Connection},c:Connection)
    [id for (id,conn) in enumerate(conns) if conn=c]
end


export Network,addcomponents,addconnections,assemble,initconnection
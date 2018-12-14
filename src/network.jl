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

function assemble(nw::Network)
    conns=nw.conns;comps=nw.comps
    numvars=3*length(conns)
    col=1
    vec=zeros(0)
    mat=zeros(0,numvars)
    for comp in comps
        eqs=equations(comp)
        vec=[vec;eqs]
        numofequations=length(eqs)
        derive=zeros(length(eqs),numvars)
        for conn in conns
            jac=jacobi(comp,conn)
            if isempty(jac)
                col+=3
                continue
            else
                derive[1:numofequations,col:col+2]+=jac
                col+=3
            end
        end
        mat=vcat(mat,derive)
        col=1
    end
    (vec,mat)
end


export Network,addcomponents,addconnections,assemble,initconnection
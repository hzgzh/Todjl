function calcturbineparams(p::Vector,h::Vector)
    n=length(p)
    ratio=[];ets=[]
    for i in range(1,n-1)
        push!(ratio,p[i+1]/p[i])
        push!(ets,(h[i]-h[i+1])/(h[i]-psh(p[i+1],phs(p[i],h[i]))))
    end
    (ratio,ets)
end

export calcturbineparams

p=[167,60.21,41.31];h=[3398.6,3130.4,3036.7]
calcturbineparams(p,h)
p=[38.03,23.89,11.74];h=[3534.7,3395.2,3179.1]
calcturbineparams(p,h)
p=[11.55,4.061,2.398,1.371,0.722,0.049];h=[3184.8,2947.4,2845.3,2742.9,2646.3,2327.8]
calcturbineparams(p,h)
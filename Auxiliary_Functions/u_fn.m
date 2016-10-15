function u = u_fn(c,gamma)

if gamma == 1
    
    u = log(c);
    
else
    
u = 1/(1-gamma) * (c.^(1-gamma) - 1);
    
end
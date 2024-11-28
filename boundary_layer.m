function theta = boundary_layer(phi,s)
if s>phi
    theta = 1;
elseif abs(s) < phi
    theta = s/phi;
else
    theta = -1;
end
end

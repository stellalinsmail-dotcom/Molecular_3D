
function EA=newgetEA(dvar,ka)
cb=-0.000122;
EA=0.043844*ka/2*dvar^2*(1+cb*dvar);
end
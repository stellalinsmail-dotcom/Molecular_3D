function ET=newgetET(phi,V1,V2,V3)
ET=0.5*(V1*(1+cosd(phi))+V2*(1-cosd(2*phi)+V3*(1+cosd(3*phi))));
end
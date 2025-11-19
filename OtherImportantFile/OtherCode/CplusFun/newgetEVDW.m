function EVDW=newgetEVDW(R_ij,Rv_ij,epsilon_ij)
EVDW=epsilon_ij*(1.07*Rv_ij/(R_ij+0.07*Rv_ij))^7*(1.12*Rv_ij^7/(R_ij+0.12*Rv_ij^7)-2);
end
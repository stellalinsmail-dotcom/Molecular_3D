function EB=newgetEB(dr_ij,kb)
cs=-2;
EB=143.9525*kb/2*dr_ij^2*(1+cs*dr_ij+7/12*cs^2*dr_ij^2);
end
function EBA=newgetEBA(dr_ij,dr_kj,dvar,kba_ijk,kba_kji)
EBA=2.51210*(kba_ijk*dr_ij+kba_kji*dr_kj)*dvar;
end
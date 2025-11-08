function plegend=PlotScatterline(legend_name,table_xy,sep,line_marker,scat_marker,line_width,marker_size,color)
    hold on
    plot(table_xy(:,1),table_xy(:,2),line_marker,'HandleVisibility','off',LineWidth=line_width,Color=color);
    table_xy_s=table_xy(1:sep:end,:);
    h1 = plot(table_xy_s(:,1),table_xy_s(:,2),scat_marker,'HandleVisibility','off',Color=color);
    set(h1, 'markerfacecolor', get(h1, 'color'),'MarkerSize',marker_size);
    h2=plot(table_xy(end,1),table_xy(end,2),scat_marker,'HandleVisibility','off',Color=color);
    set(h2, 'markerfacecolor', get(h2, 'color'),'MarkerSize',marker_size);
    a_marker=[line_marker,scat_marker];
    plegend=plot(table_xy_s(1,1),table_xy_s(1,2),a_marker,'DisplayName',legend_name,LineWidth=line_width,Color=color);
    set(plegend, 'markerfacecolor', color,'MarkerSize',marker_size*0.75);
end
function  DrawRULine(gca,axes_xcolor,axes_ycolor,lwidth)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明

    hold on
    XL = get(gca,'xlim'); XR = XL(2);
    YL = get(gca,'ylim'); YT = YL(2);
    plot(XL,YT*ones(size(XL)),'-','color',axes_xcolor,'linewidth',lwidth,'HandleVisibility','off')
    plot(XR*ones(size(YL)),YL,'-','color', axes_ycolor,'linewidth',lwidth,'HandleVisibility','off')
end

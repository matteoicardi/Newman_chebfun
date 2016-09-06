
    figure(1)
    plot(ps)
    hold on
    plot(pe)
    figure(2)
    plot(ce)
    hold on
    plot(cs)
    figure(3)
    plot(j(ce,cs,pe,ps)./F./cmax)
    hold on
    drawnow

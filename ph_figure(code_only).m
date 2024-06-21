gif_filename = 'pH_change_animation.gif';

figure(2);
hold on;

fig1 = plot(pH_scale, pH(:,1), 'DisplayName', 'pH change in EKR');
fig2 = plot(pH_scale, pH_B(:,1), 'DisplayName', 'pH change in BKR');

xlabel('Distance (cm)');
ylabel('pH Level');

for h = 1:50401
    set(fig1, 'XData', pH_scale, 'YData', pH(:,h)');
    set(fig2, 'XData', pH_scale, 'YData', pH_B(:,h)');

    h2 = h / 1440;
    title(['pH Change - Day: ', sprintf('%.2f', h2)]);
    
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    if h == 1
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.01);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.01);
    end
    
    pause(0.00001);
end

legend;
hold off;

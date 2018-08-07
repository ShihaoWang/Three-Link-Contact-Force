function Three_Link_Plot(z)

% This function is used to plot the three link robot 

% The input to this function z is assumed to be directly from the ode
% integration function

load('rB_fn.mat');          load('rC_fn.mat');              load('rD_fn.mat');

[m,n] = size(z);

rA = [ 0 0 0]';

Link_Color = [57 132 243] ./ 255;
Pivot_Color = [255 26 9] ./ 255;
for i = 1:m 
    
    state_i = z(i,:);
    Theta = state_i(1);             Alpha = state_i(2);             Beta = state_i(3);
    Thetadot = state_i(4);          Alphadot = state_i(5);          Betadot = state_i(6);
    rB = rB_fn(Theta);              rC = rC_fn(Alpha, Theta);       rD = rD_fn(Alpha, Beta, Theta);
    
    % First is to plot the link image
    plot([rA(1); rB(1)], [rA(2); rB(2)],'LineWidth',2.5,'color',Link_Color);
    hold on
    plot([rB(1); rC(1)], [rB(2); rC(2)],'LineWidth',2.5,'color',Link_Color);
    hold on
    plot([rC(1); rD(1)], [rC(2); rD(2)],'LineWidth',2.5,'color',Link_Color);
    hold on
    
    % Second is to plot the pivot
    plot(rA(1), rA(2),'o',...
        'MarkerSize',5,...
        'MarkerEdgeColor',Pivot_Color,...
        'MarkerFaceColor',Pivot_Color);
    hold on    
    
     plot(rB(1), rB(2),'o',...
        'MarkerSize',5,...
        'MarkerEdgeColor',Pivot_Color,...
        'MarkerFaceColor',Pivot_Color);
    hold on    
    
     plot(rC(1), rC(2),'o',...
        'MarkerSize',5,...
        'MarkerEdgeColor',Pivot_Color,...
        'MarkerFaceColor',Pivot_Color);
    hold on    
    
    plot(rD(1), rD(2),'o',...
        'MarkerSize',5,...
        'MarkerEdgeColor',Pivot_Color,...
        'MarkerFaceColor',Pivot_Color);
    
    xlabel('Horizontal Axis')
    ylabel('Vertical Axis')
    axis equal
    hold off
    pause(0.01)
end

end


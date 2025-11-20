function [state] = TIV(model,V0)


paramsT = model.paramsT;
tspan=0:paramsT.IntT:paramsT.maxT;

%sigma=10;
if paramsT.modelNum ==1
    Para = [paramsT.p(paramsT.modelNum) paramsT.Beta_I(paramsT.modelNum)...
        paramsT.delta(paramsT.modelNum) paramsT.c(paramsT.modelNum)];
    y0 = [paramsT.T0 paramsT.I0 V0];
    
    [~,state]=ode15s(@tiv_simple,tspan,y0,[],Para);
elseif paramsT.modelNum ==2
    Para = [paramsT.p(paramsT.modelNum) paramsT.Beta_I(paramsT.modelNum)...
        paramsT.delta(paramsT.modelNum) paramsT.c(paramsT.modelNum) paramsT.k];
    %y0 = [paramsT.T0 paramsT.L0 paramsT.I0 V0];
    y0 = [paramsT.T0 paramsT.L0 V0 paramsT.I0];

    [~,state]=ode15s(@tiv_delay,tspan,y0,[],Para);
end


end



function f=tiv_simple(~,y,Para)
    f=zeros(3,1);
    f(1)=-Para(2)*y(1)*y(3);
    f(2)=Para(2)*y(1)*y(3)-Para(3)*y(2);
    f(3)=Para(1)*y(2)-Para(4)*y(3);
end

function f=tiv_delay(~,y,Para)
    f = zeros(4,1);
    f(1)=-Para(2)*y(1)*y(3);
    f(2)=Para(2)*y(1)*y(3)-Para(5)*y(2);
    f(3)=Para(1)*y(4)-Para(4)*y(3);
    f(4)=Para(5)*y(2)-Para(3)*y(4);
    %f(1)=-Para(2)*y(1)*y(4);
    %f(2)=Para(2)*y(1)*y(4)-Para(5)*y(2);
    %f(3)=Para(5)*y(2)-Para(3)*y(3);
    %f(4)=Para(1)*y(3)-Para(4)*y(4);
end
% plot TIV model
% V0 = 10.^(0:0.5:2);
% for i = 1:length(V0)
%     State{i,1}=TIV(model,V0(i));
% 
% end
% time = 0:1/24:10;
% for i = 1:length(V0)
%     hline = plot(time,log10(State{i,1}(1:24*10+1,3)),'b-', 'linewidth',1.5);
%     hline.Color = [hline.Color 0.9];
%     hold on
% end

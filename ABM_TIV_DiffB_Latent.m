function [Ytime,states,InfectIdtau,agent]=ABM_TIV_DiffB_Latent(model,seed,type,TIV_VL)

rng(seed)
beta_tau=model.beta_tau;
Get_beta=beta_tau{type};
condition=model.condition;
params = model.params;
I0 = params.I0;
TimeStep = params.IntT;
paramsT=model.paramsT;
tspan =0:paramsT.IntT:paramsT.maxT;
%V0 = 10.^(2*rand(1,params.N));
InitalInfectTime = 0:params.IntT:4;
I0Tau = InitalInfectTime(randi(length(InitalInfectTime),[1,10]));
InfectIdtau = {};
jj=1;
cc=1;
ContactN = poissrnd(params.Arrival*TimeStep,[1,10^6]);
p = rand([1,10^6]);
for i = 1:params.N
    agent(i).I = 0;
    agent(i).S = 1;
    agent(i).R = 0;
    agent(i).InfectT=0;
    agent(i).RecoverT=0;
    agent(i).TIVs = TIV_VL(i,:);
    agent(i).IID = [];
    agent(i).betta=0;
    agent(i).CT = 0;
    agent(i).CTDay = 0;
end

InfectN=1:I0;

for i = 1:length(InfectN)
    agent(i).I=1;
    agent(i).S=0;
    %[stateI] = TIV(model,V0(i));
    RecoverT = recover_fun(agent(i).TIVs,tspan,params.FoC);
    % if type==1 || type == 2
    %     agent(i).betta = Get_beta(agent(i).TIVs(1:single(RecoverT/params.IntT)+1),max(agent(i).TIVs(1:single(RecoverT/params.IntT)+1)));
    % 
    % else
    agent(i).betta = Get_beta(agent(i).TIVs(1:single(RecoverT/params.IntT)+1));
    % end
    agent(i).InfectT=-I0Tau(i);
    agent(i).RecoverT = RecoverT-I0Tau(i);   
    agent(i).IID = 0;
    agent(i).CT = (log10(agent(i).TIVs(1:single(RecoverT/params.IntT)+1))-9.988)/(-0.3152);
    agent(i).CTDay = [ceil(agent(i).InfectT):1:floor(agent(i).RecoverT);...
        agent(i).CT(1:1/params.IntT:(floor(agent(i).RecoverT)-ceil(agent(i).InfectT))/params.IntT+1)];
end

InfectIdtau{1} = InfectN;

states = [sum([agent(:).S]) sum([agent(:).I]) sum([agent(:).R])];
% Record the Susceptiable individual number
InfectN = InfectIdtau{end};
SuscepN=max(InfectIdtau{end})+1:params.N;
% Record the time
Ytime=0;
        
  
for ii = 1:params.T/TimeStep
    tt=ii*TimeStep;
    if ~condition(states(end,:))
        InfectPeople = []; % Record the infected people being recover for each time step
        for j = 1:length(InfectN)
            if tt<=agent(InfectN(j)).RecoverT
                

                if length(SuscepN)>0
                    
                    if jj>=10^5
                        ContactN = poissrnd(params.Arrival*TimeStep,[1,10^5]);                      
                        jj=1;
                    end

                    Contact = binornd(ContactN(jj),length(SuscepN)/(params.N-1));
                    jj = jj +1;
                    %for ii = 1:ContactN
                    %    if rand<length(SuscepN)/params.N
                    %        Contact=Contact+1;
                    %    end
                    %end
                    %Contact = round(ContactN*length(SuscepN)/params.N);
                    if Contact>0
                        SuscepPeople=[];
                        for i = 1:Contact
                            
                            p_ill = agent(InfectN(j)).betta(single((tt-agent(InfectN(j)).InfectT)/params.IntT)+1);
                            if cc>=10^5
                                p = rand([1,10^5]);
                                cc= 1;
                            end
                            cc=cc+1;
        
                            if p(cc)<=p_ill                              
                                agent(SuscepN(i)).InfectT=tt;
                                agent(SuscepN(i)).I=1;                                
                                agent(SuscepN(i)).S=0;     
                                agent(SuscepN(i)).IID = InfectN(j);
                                %[stateI] = TIV(model,V0(SuscepN(i)));                                
                                %agent(SuscepN(i)).TIVs = stateI;                    
                                RecoverT = recover_fun(agent(SuscepN(i)).TIVs,tspan,params.FoC);                               
                                agent(SuscepN(i)).RecoverT=tt+RecoverT;
                                stateVL= agent(SuscepN(i)).TIVs(1:single(RecoverT/params.IntT)+1);
                                % if type ==1 || type == 2
                                %     agent(SuscepN(i)).betta(1:single((RecoverT)/params.IntT)+1) = Get_beta(stateVL,max(stateVL));
                                % else
                                agent(SuscepN(i)).betta(1:single((RecoverT)/params.IntT)+1) = Get_beta(stateVL);
                                %end
                                agent(SuscepN(i)).CT = (log10(agent(SuscepN(i)).TIVs(1:single(RecoverT/params.IntT)+1))-9.988)/(-0.3152);
                                agent(SuscepN(i)).CTDay = [ceil(agent(SuscepN(i)).InfectT):floor(agent(SuscepN(i)).RecoverT);...
                                            agent(SuscepN(i)).CT(1:1/params.IntT:(floor(agent(SuscepN(i)).RecoverT)-ceil(agent(SuscepN(i)).InfectT))/params.IntT+1)];
                                InfectN=[InfectN SuscepN(i)];
                                SuscepPeople = [SuscepPeople i];
                                                                
                            end
                        end
                        SuscepN(SuscepPeople)=[]; % Remove the susceptiable people who have infected
                    end                    
                end
    
            else
                agent(InfectN(j)).I=0;
                agent(InfectN(j)).R=1;
                InfectPeople = [InfectPeople j];
                %agent(InfectN(j)).betta(single((agent(InfectN(j)).RecoverT-agent(InfectN(j)).InfectT)/params.IntT)+1)=0;
                
            end
        
        end
        InfectN(InfectPeople)=[]; % Remove the infected people who have recovered
        InfectIdtau{single(tt/params.IntT)+1}=InfectN;

    end
    states= [states; sum([agent(:).S]) sum([agent(:).I]) sum([agent(:).R])];
    Ytime=[Ytime tt];
end


end

function [RecoverT] = recover_fun(VL_state,tspan,FoC)
    
    NoT = sum(VL_state>=FoC);
    RecoverT = tspan(NoT);

end
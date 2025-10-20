function [SumStats,states] = OutputABM(model,seed,ConFun,V_StateC)
Day = 150;
N =10000;

[Ytime,states,InfectIdtau,agent]=ABM_TIV_DiffB_Latent(model,seed,ConFun,V_StateC);

DI=zeros(1,Day);
SkewCT = zeros(1,Day+1);
MedianCT= zeros(1,Day+1);
MeanCT = zeros(1,Day+1);
sdCT = zeros(1,Day+1);
 for j = 1:Day
    DI(j+1) = states(24*(j-1)+1,1)-states(24*(j)+1,1);   
 end

InfectID = InfectIdtau;
CT_i = [];
for j = 1:Day+1
    if j*24+1<=length(InfectID)
        InfectList = InfectID{(j-1)*24+1};
        CT_i =[];
        for jj = 1:length(InfectList)
            Time = find(agent(InfectList(jj)).CTDay(1,:)==j);
            CT_i = [CT_i agent(InfectList(jj)).CTDay(2,Time)];
        end
        CT_i(find(CT_i>38.033))=[];
        SkewCT(j) = skewness(CT_i);
        MedianCT(j) = median(CT_i);
        MeanCT(j) = mean(CT_i);
        sdCT(j) = std(CT_i);
    else
        SkewCT(j)=0;
        MedianCT(j)=0;
        MeanCT(j)=0;
        sdCT(j) = 0;
    end
end
MeanCT(isnan(MeanCT))=0;
sdCT(isnan(sdCT))=0;
SkewCT(isnan(SkewCT))=0;
SumStats.DI = DI;
SumStats.SkewCT=SkewCT;
SumStats.MedianCT = MedianCT;
SumStats.MeanCT = MeanCT;
SumStats.sdCT = sdCT;
end
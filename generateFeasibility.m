%% Feasibility and stability for Ki and Kp

clear all
[ad,bd,A,B_u,B_d,R,Ts,N,Ndist,K,Ti,Kp,Ki,K_fb,distVec,distVec_ctrl,x0] = readParamsToWorkspace();

n = 70;

KiList = -1*logspace(-3,2,n);
KpList = -1*logspace(-3,2,n);

FeasibleListOsv = zeros(length(KiList)*length(KpList), 4);

infiniteDenominator = zeros(length(KiList),2);

for i = 1:length(KiList)
    Ki = KiList(i);
    
    % Infeasible when a=b(kp-ki)
    infiniteDenominator(i,1) = Ki;
    infiniteDenominator(i,2) = Ki + ad/bd;

    for j = 1:length(KpList)
        Kp = KpList(j);

        K_fb = [-Ki, -Kp];

        FeasibleListOsv((i-1)*length(KpList) + j, 1) = Ki;
        FeasibleListOsv((i-1)*length(KpList) + j, 2) = Kp;

        [Q, Qf] = calculateQdiscForPI(ad,bd,Ki,Kp,R);
        
        % Feasibility (x^TQx>=0). Column 3
        if Q(1,1) > 0 && Q(2,2) >= 0  %min(diag(Q)) < 0
            FeasibleListOsv((i-1)*length(KpList) + j, 3) = 1;
        else
            FeasibleListOsv((i-1)*length(KpList) + j, 3) = 0;
        end

        % Stability abs(eig((A-BK))) <= 1. On column 4
        %z1 = (1 + ad - bd*Kp)/2 + sqrt((-1-ad+bd*Kp)^2/4 - (ad+bd*(Ki-Kp)));
        %z2 = (1 + ad - bd*Kp)/2 - sqrt((-1-ad+bd*Kp)^2/4 - (ad+bd*(Ki-Kp)));
        maxEig = max(abs(eig(A-B_u*K_fb)));
        if maxEig > 1
            FeasibleListOsv((i-1)*length(KpList) + j, 4) = 0;
        else
            if FeasibleListOsv((i-1)*length(KpList) + j, 3) == 0
                FeasibleListOsv((i-1)*length(KpList) + j, 4) = 1;
            else
                FeasibleListOsv((i-1)*length(KpList) + j, 4) = 2;
            end
        end     
    end
end


figure;
scatter(log10(-1*FeasibleListOsv(:, 1)), log10(-1*FeasibleListOsv(:, 2)), [], FeasibleListOsv(:, 4), "filled")
xlabel('log(-Ki)'); ylabel('log(-Kp)'); zlabel('StableAndFeasible');
hold on;
scatter(log10(-1*infiniteDenominator(:,1)),log10(-1*infiniteDenominator(:,2)),'filled')
ylim([-3 2])


% Kp: log10(-1*-0.9)=-0.046
% Ki: log10(-1*-0.0106)=-2
%%
% Save data as txt file
T_cell = table(log10(-1*FeasibleListOsv(:, 1)),log10(-1*FeasibleListOsv(:, 2)),FeasibleListOsv(:, 4), 'VariableNames',["log10(-Ki)","log10(-Kp)","StableAndFeasible"]);
writetable(T_cell,'txtData/StableAndFeasible.txt');

T_cell = table(log10(-1*infiniteDenominator(:,1)),log10(-1*infiniteDenominator(:,2)), 'VariableNames',["log10(-Ki)","log10(-Kp)"]);
writetable(T_cell,'txtData/StableAndFeasible_infiniteDenominator.txt');



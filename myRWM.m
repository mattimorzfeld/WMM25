function [X, D, LogPi, AccRatio] = myRWM(Nsteps,Xo,logpi,H)
nData = sum(H);
n = length(Xo);
X = zeros(n,Nsteps);
D = zeros(nData,Nsteps);
LogPi = zeros(Nsteps,1);
%% initialization
X(:,1) = Xo;
x = Xo;
[logpi_x,dx] = logpi(x);

accMoves = 0;
%% move NSteps
for Step=1:Nsteps-1
    xp = x+sqrt([.001;.001;.001;1;2;10]).*randn(n,1);
    xp(3) = randi(2,1)-1;
    [logpi_xp,dxp] = logpi(xp);
    loga = logpi_xp-logpi_x;
    a = min(1,exp(loga));
    if rand<a
        accMoves = accMoves+1;
        X(:,Step+1) = xp;
        D(:,Step+1) = dx;
        LogPi(Step+1) = logpi_xp;
        dx = dxp;
        x = xp;
        logpi_x = logpi_xp;
    else
        X(:,Step+1) = x;
        D(:,Step+1) = dx;
        LogPi(Step+1) = logpi_x;
    end

    if ~mod(Step,1e2)
        fprintf('Acc. ratio at step %g/%g: %g\r',Step,Nsteps,accMoves/Step)
    end
end
AccRatio = accMoves/Nsteps;
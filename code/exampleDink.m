%%% Example of Dinkelbach's Algorithm
%%% Written by Quang-Doanh Vu
%%% This is a modified version of the one at https://sites.google.com/site/vuquangdoanh/publications?authuser=0
clear 
clc
load ../data/examplechannel
whos
%%%%
%%%% 3 users, 2 small cells, smallcell anttena = 4; macrocell antenna = 8;
%%%% SDPT3+MAXDET
tau = 0; %%Dinkelbach parameter
Time=[];
Sc  = cell(NumUser,1);
LowM  = cell(NumUser,1);
p = sdpvar(NumUser,L);
sumrate = 0;
F0 = []; % constraints that are the same in all iterations
for iUser = 1:NumUser
    Sc{iUser} = sdpvar(User(iUser).Nt,User(iUser).Nt,'hermitian','complex');
    LowM{iUser} = (tril(sdpvar(L,L,'full','complex'),-1));
    sumrate = sumrate + logdet(diag(p(iUser,:)));
    F0 = [F0,Sc{iUser}>=0,p(iUser,:)>=0];
end
solveroptions = sdpsettings('solver','sdpt3','verbose',0);
maxIteration = 100;
tau_seq = zeros(maxIteration,1);
for iIteration =1:maxIteration
    rateE = zeros(1,NumUser);
    BS = struct('Pow',[]); PowTot=[];
    for j=1:M+1
        BS(j).Pow=zeros(K(j));
    end
    F1 = []; % constraints need to be updated in each iteration
    for iUser = 1:NumUser
        F1 = [F1,([eye(L)+User(iUser).H_h*Sc{iUser}*User(iUser).H_h',LowM{iUser}+diag(p(iUser,:));(LowM{iUser}+diag(p(iUser,:)))',diag(p(iUser,:))]>=0)];
        F1=[F1, [p(iUser,1),exp(rmin(iUser)/2);exp(rmin(iUser)/2),p(iUser,2)]>=0];
        BS(1).Pow = BS(1).Pow + User(iUser).A1*User(iUser).V0*Sc{iUser}*(User(iUser).A1*User(iUser).V0)';
        if(iUser>Nm)
            BS(User(iUser).possBS(2)).Pow = BS(User(iUser).possBS(2)).Pow+ User(iUser).A2*User(iUser).V0*Sc{iUser}*(User(iUser).A2*User(iUser).V0)';
        end
    end
    for j=1:M+1
        F1 = [F1, trace(BS(j).Pow) <= Pmax(j)];
        PowTot = blkdiag(PowTot,BS(j).Pow*lambda(j));
        for jj = 1:K(j)
            F1=[F1, BS(j).Pow(jj,jj)<=Pa(j)];
        end
    end
    objective = sumrate - tau*(trace(PowTot)+Pc);
    diagontics = solvesdp([F0,F1],-objective,solveroptions);
    for iUser=1:NumUser
        rateE(iUser)=real(double(logdet(eye(L)+User(iUser).H_h*Sc{iUser}*User(iUser).H_h')));
    end
    Pow = double(trace(PowTot));
    Rate = sum(rateE);
    tau_next = real(Rate/(Pow+Pc));
    tau_seq(iIteration) = tau_next;
    if(tau_next-tau < epsilon)
        break;
    end
    tau = tau_next;
end
tau_seq(iIteration+1:end) =[];
plot(tau_seq)
xlabel('Iteration count, $n$','Interpreter','latex')
ylabel('$\tau_n$','Interpreter','latex')
saveas(gcf,'../results/convergence.png')
clear
close all
% Enter model parameters

VERSION='8';

m=2; % 2 actions (work and shirk)
w=1; 
% s=1.2;


% Construct state space
nd=15;     % grid size on each dimension
S = gridmake((1:nd)',(1:nd)',(1:nd)');
n    = size(S,1);                    % number of states

% Find in each hierarchy who is J, S and F. Ties are broken randomly
% [ iif,iij,iis,ii2f,ii2j,ii2s,ii3f,ii3j,ii3s, S ] = construct_indices( nd );
iif = (S(:,1)<=S(:,2) & S(:,1)<=S(:,3)); 
iij = ((S(:,1)>=S(:,2) & S(:,1)<=S(:,3)) | (S(:,1)<=S(:,2) & S(:,1)>=S(:,3))); 
iis = (S(:,1)>=S(:,2) & S(:,1)>=S(:,3)); 
ii2f = (S(:,2)<=S(:,3) & S(:,2)<=S(:,1));   % indices where 2 is F
ii2j = ((S(:,2)>=S(:,3) & S(:,2)<=S(:,1)) | (S(:,2)<=S(:,3) & S(:,2)>=S(:,1))); 
ii2s = (S(:,2)>=S(:,3) & S(:,2)>=S(:,1));
ii3f = (S(:,3)<=S(:,2) & S(:,3)<=S(:,1));   % indices where 3 is F
ii3j = ((S(:,3)>=S(:,2) & S(:,3)<=S(:,1)) | (S(:,3)<=S(:,2) & S(:,3)>=S(:,1))); 
ii3s = (S(:,3)>=S(:,2) & S(:,3)>=S(:,1)); 


        
pi_default = zeros(n,1);
pi_default(iif)=0;
pi_default(iij)=1;
pi_default(iis)=2;
pi_default(iif&iij)=.5;
pi_default(iij&iis)=1.5;
pi_default(iij&iif&iis)=1;


% bb=.1:.2:.9;
% dd=.8:.02:.98;

bb=.1:.025:.95;
dd=.8:.005:.98;


d_list=[]; b_list=[];
% Add basic uniform points to the list of nodes
GG=gridmake(dd',bb');
d_list=[d_list; GG(:,1)];
b_list=[b_list; GG(:,2)];
% 
bb2=.1:.2:.9;
dd2=.7:.02:.79;
GG=[];
GG=gridmake(dd2',bb2');
d_list=[d_list; GG(:,1)];
b_list=[b_list; GG(:,2)];
% 
% % Now manually refine the grid to calculate the border more precisely
% 
[d_list_new, b_list_new] = makegrid_in_rec([.92 .2], [.065,.11], 40 );
d_list=[d_list; d_list_new];
b_list=[b_list; b_list_new];
[d_list_new, b_list_new] = makegrid_in_rec([.975 .87], [.01,.05], 10 );
d_list=[d_list; d_list_new];
b_list=[b_list; b_list_new];

%%%%%

eqm_list=ones(size(b_list))*-1;
converged=ones(size(b_list));
SW=ones(size(b_list))*-1;
RECOMP=0;  % if true, forces to recompute v even if there is a saved data
for num_eq=1:length(d_list)
    num_eq_=num_eq;
    b     = b_list(num_eq_);
    delta = d_list(num_eq_);
        
    b_str = num2str(b);
    delta_str = num2str(delta);
    filename = ['saves/' VERSION '/vx_' b_str '_' delta_str '.mat']
    if exist(filename,'file') && ~RECOMP
        LS=load(filename);
        V=LS.V; g=LS.g; gp=LS.gp; gp_stat=LS.gp_stat; pi=LS.pi; 
        conv_flag=LS.conv_flag; X=LS.X; pstat=LS.pstat; pstat2=LS.pstat2;
        VE_cr=LS.VE_cr; 
    else
        % Construct state transition probability matrix


        % cdf evaluated at xt conditional on the current seniority level x
        G = @(xt,x)  max((1-x).^(1/b-1) - (1-xt)^(1/b-1),0) ./ (1-x).^(1/b-1); 

        st = 1/nd;
        % construct individual advancement probability function g
        g = zeros(nd);              % transition probability matrix: advancement only
        ll1 = (((1:nd)-1)./nd)';
        for j=1:nd
            lt1 = (j-1)/nd;
            g(:,j) = G(lt1+st,ll1)-G(lt1,ll1);
        end

        % construct P 

        g_die = zeros(nd);
        g_die(:,1)=1;                   % transition probability matrix: die or leave
        gp = delta*g + (1-delta)*g_die;  % transition probability matrix ...
        % that takes into account death with probability delta

        % next three are auxilliary matrices to be used later in the loop
        % taken out of the loop
        Pii0 = kron(kron(gp,gp),g);    % no one leaves
%               kron(kron(a,b),c): c is own transition probability (l1), 
%                                   b is for l2, a is for l3
        Pii1 = kron(kron(gp,gp),g_die); % 1 leaves
        Pii2 = kron(kron(gp,g_die),g);  % 2 leaves
        Pii3 = kron(kron(g_die,gp),g);  % 3 leaves
        P_out= kron(kron(gp,gp),gp);    % everyone can die (outside observer's persepective)

        P = zeros(m,n,n); %P(1) is for work, P(2) is for shirk
        P(1,1:n,1:n)=Pii0;%kron(kron(gp,gp),g); % opponents either die or advance 
        % P(1) is going to be changing in cycle below: F might be leaving

        % if shirk, start in a new community and advance:         
        gp_stat=stat_distr(gp);
        P(2,1:n,1:n)=repmat(kron(kron(gp_stat',gp_stat'),g(1,:)),n,1);

        % pi's
        pi = pi_default;  % will be modified in the cycle below

        V=zeros(n,1);
        
        f_shirk_ii=zeros(size(S,1),1);
        itern=1;
        gamma=1;      % updating rate parameter
        gamma3=0;
        conv_flag=1;
        while 1        % iterate value and best response to find an equilibrium 
            Reward = zeros(n,m);
            Reward(:,1) = pi.*w;    % work
            Reward(:,2) = pi.*w;    % shirk. Reward is the same because we are looking
                                        % for an eqm where S and J work on path. 
            % Pack model data
%             clear model
            model = struct;
            model.reward = Reward;
            model.discount = delta;
            if norm(sum(squeeze(P(1,:,:)),2)-ones(n,1)) > 1e-4, error 'P is not a prob matrix'; end
            model.transprob = P;

            optset('ddpsolve','algorithm','funcit'); % Value function interation
            optset('ddpsolve','maxit',500);
            optset('ddpsolve','prtiters',0);  % don't print iterations

            % Solve for V
            Vinit = ones(n,1)*w/(1-delta);
            Vold = V;
            [V,X,pstar] = my_ddpsolve(model,Vinit); 
            
            
            % % compute indices about leaving opponents
            f_shirk_ii0 = f_shirk_ii;
            f_shirk_ii = iif&~iij&~iis&(X==2);    % indices in state space where (a) fresman (b) shirks. 
            xx=xor(f_shirk_ii,f_shirk_ii0);
            iis_f_shirks = f2s(S(f_shirk_ii,:),nd); % indices where I am S and F shirks
            iif_=iif;                   % indices where I get zero stage payoff
            iij_=iij|iis_f_shirks;      % indices where I get stage payoff w
            iis_=iis&~iis_f_shirks;     % indices where I get stage payoff 2w
            
            ii2 = opp_shirks(S(f_shirk_ii,:),nd,2);  % inidces where opponent 2 shirks and leaves 
            ii3 = opp_shirks(S(f_shirk_ii,:),nd,3);  % inidces where opponent 3 shirks and leaves 

            
            % % recompute the outside option
            % find stationary distribution 
            P_s=P_out;%kron(kron(gp,gp),gp); 
            P_s(find(f_shirk_ii),:) = Pii1(find(f_shirk_ii),:);
            P_s(find(ii2),:) = Pii2(find(ii2),:);
            P_s(find(ii3),:) = Pii3(find(ii3),:);
            [pstat,ec]=stat_distr(P_s); 
                % pstat is stationary distribution from the outside observer 
                % perspective: if someone leaves, a new person starts from the scratch
                % It is used to calculate pstat2 and later SW
            if ec, error('pstat is not found (while loop)'), end
            % Find stationary distribution in the new community
            ii_1 = (S(:,1)==1);  % indices where player 1 is on the initial position
            pstat2 = pstat(ii_1)/sum(pstat(ii_1));

            % % recompute pi's:
            % seniors who face F who is leaving, get pi=1 instead of
            % pi=2
            pi0=pi;
            
            pi = zeros(n,1);
            pi(iif_)=0;
            pi(iij_)=1;
            pi(iis_)=2;
            pi(iif_&iij_)=.5;
            pi(iij_&iis_)=1.5;
            pi(iij_&iif_&iis_)=1;

            yy=find(pi-pi0);
            yy_fire = yy(find(binornd(1,gamma3,size(yy,1),1)));
            pi(yy_fire) = pi0(yy_fire);  % introduce some randomness to avoid cycles

            
            % % recompute P:
            % opponent F can be leaving on path
            P0=P;
            P(1,1:n,1:n)=Pii0;%kron(kron(gp,gp),g); 
            P(1,find(ii2),:) = Pii2(find(ii2),:);
            P(1,find(ii3),:) = Pii3(find(ii3),:);
            P(1,:,:) = squeeze(P0(1,:,:))*(1-gamma) + squeeze(P(1,:,:))*gamma;
            
            % update the outside option with the new stationary
            % distribution pstat2
            P(2,1:n,1:n)=repmat(kron(pstat2',g(1,:)),n,1);


            eps=1e-4;
            pi_diff = norm(pi-pi0,1);
            P_diff  = sum(sum(sum(abs(P(1,:,:)-P0(1,:,:)))));
            if pi_diff<=0 && (P_diff<1 || norm(V-Vold)<eps)
                break
            else
                itern=itern+1
                pi_diff
                if pi_diff<=5
                    yy=find(pi-pi0);
                    S(yy,:)
                    [pi(yy) pi0(yy)]
                end
                P_diff
                norm(V-Vold)
%                     sum(f_shirk_ii)
                if itern>10, gamma=.8; gamma3=.5; end;
                if itern>20, gamma=.5; end;
                if itern>50, warning('No convergence!'), conv_flag=0; break; end
            end
        end
        
        % Find stationary distribution (will be used to compute v0)
%         [pstat,ec]=stat_distr(pstar);  % wrong formula!
%         if ec, error('pstat is not found'), end
        
        % Find stationary distribution in the new community
%         ii_1 = (S(:,1)==1);
%         sum(ii_1)
%         pstat2 = pstat(ii_1)/sum(pstat(ii_1));
        
         % compute cheat-and-return value function (will be used for no
         % cheat-and-return IC)
         % first compute the state transition probability matrix P_cr
        P_cr = zeros(n,4*n);
        P_cr(1:n,1:n) = delta^2*kron(kron(g,g),g);  % 2&3 survive
        P_cr(1:n,n+1:2*n) = delta*(1-delta)*kron(kron(g,g_die),g);  % 2 dies
        P_cr(1:n,2*n+1:3*n) = delta*(1-delta)*kron(kron(g_die,g),g);  % 3 dies
        P_cr(1:n,3*n+1:4*n) = (1-delta)^2*kron(kron(g_die,g_die),g);  % 2&3 die
        P_cr(find(ii2),n+1:2*n) = P_cr(find(ii2),n+1:2*n) + P_cr(find(ii2),1:n);  % indices where the opponent 2 shirks and leaves
        P_cr(find(ii2),1:n) = zeros(length(find(ii2)),n);
        P_cr(find(ii3),2*n+1:3*n) = P_cr(find(ii3),2*n+1:3*n) + P_cr(find(ii3),1:n);  % indices where the opponent 3 shirks and leaves
        P_cr(find(ii3),1:n) = zeros(length(find(ii3)),n);
        if norm(sum(P_cr,2)-ones(n,1)) > 1e-4, error 'P_cr is not a prob matrix'; end
        
        % second, compute scale of trade vector pi_cr
        pi_cr = zeros(4*n,1);
        pi_cr(1:n) = zeros(n,1); % if both opponents survive, no one works with me
        pi_cr2 = pi;   % this pi is updated in the cycle above to take into account the leaving F
            % but opponent 3 does not work with me:
            pi_cr2(iij&ii3s)=0;
            pi_cr2(iis&ii3j)=max(pi_cr2(iis&ii3j)-1,0);
            pi_cr2(iis&ii3f&~ii3)=max(pi_cr2(iis&ii3f&~ii3)-1,0);
        pi_cr(n+1:2*n) = pi_cr2;
        pi_cr3 = pi;   % this pi is updated in the cycle above to take into account the leaving F
            % but opponent 2 does not work with me:
            pi_cr3(iij&ii2s)=0;
            pi_cr3(iis&ii2j)=max(pi_cr3(iis&ii2j)-1,0);
            pi_cr3(iis&ii2f&~ii2)=max(pi_cr3(iis&ii2f&~ii2)-1,0);        
        pi_cr(2*n+1:3*n) = pi_cr3;
        pi_cr(3*n+1:4*n) = pi;   % this pi is updated in the cycle above to take into account the leaving F
            
        % compute VE_cr 
        VE_cr = P_cr*pi_cr*w + delta*squeeze(P(1,:,:))*V; % P(1,:,:) in the second term because I *return*
        
        % Save data

%         save(filename, 'V','X','delta', 'b','pi','nd','n','w','g','gp','gp_stat','pstat','pstat2','conv_flag','VE_cr');
        parsave(filename,  V , X , delta ,  b , pi , nd , n , w , g , gp , gp_stat , pstat , pstat2 , conv_flag , VE_cr  );
    end
%     continue


    % % Check the eqm conditions

    % Compute Emilio's value function
    % construct P with advancement only

    P_adv=zeros(n,n);
    P_adv(1:n,1:n)=kron(kron(g,g),g);

    VE=P_adv*V;

    v_IF = w/(1-delta);
    v0=kron(pstat2',g(1,:))*V;

    violated=0;

     % Find the range of s s.t. the IC conditions hold and F is not
     % trusted as J
     eps=1e-4;

     % Compare with the utility of IF people
     if v0<v_IF-eps
         disp(['Ins.-free is better off than starting over']);
         violated=bitor(violated,1);
     end

    % Check no cheat-and-leave IC for S and J
     XX = delta*(VE-v0)./pi;
     [minv,~] = min(XX(iis|iij));  % IC for S and J
     s_u=minv+w;  % largest s so that the constraint holds
     if s_u<w-eps
            disp('s_u is less than w!');
            violated=bitor(violated,4);
     end
     
     % Check no cheat-and-retun IC
     XX = delta*(VE-VE_cr)./pi;
      [minv2,~] = min(XX(iis));  % IC for S 
     s_u2=minv2+w;
     if s_u2<w-eps
            disp('s_u2 is less than w!');
            violated=bitor(violated,2);
     end



    eqm_list(num_eq)=violated;
    converged(num_eq)=conv_flag;
    fprintf ('b=%g delta=%g violated=%d\n',b,delta, violated)  % print progress
    SW(num_eq)=pstat'*VE;

end

    

save(['saves/' VERSION '/eqm'],'SW','eqm_list','b_list','d_list','converged');


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

load(['saves/' VERSION '/eqm.mat']);

%% Create homogenous grid 

bb_grid=.1:.01:.9;
dd_grid=.7:.002:.98;
Int = scatteredInterpolant(d_list,b_list,eqm_list,'nearest');
G=gridmake(dd_grid',bb_grid'); 
% next three lines overwrite the loaded data
d_list=G(:,1); 
b_list=G(:,2);
eqm_list = Int(G);

%% Plot raw scattered points

au = unique(eqm_list);
leg_str={};
% e_all=[0 1 2 3 4 6];
for i=1:length(au)
    e=au(i);
    iiu = (eqm_list==e)&(converged==1);

    marker='.';
    switch e
        case 0
            leg_str{i}='Eqm exists [0]';
            col='blue';
            marker='o';
        case 1
            leg_str{i}='(IF) violated [1]';
            col='red';
        case 2
            leg_str{i}='(IC\_AR) violated [2]';
            col='yellow';
        case 3
            leg_str{i}='(IC\_AR)&(IF) violated [3]';
        case 4
            leg_str{i}='(IC\_AL) violated [4]';
            col='magenta';
        case 6
            leg_str{i}='(IC\_AL)&(IC\_AR) violated [6]';
            col='green';
        otherwise
            leg_str{i}=sprintf('[%d]',e);
    end
    if e==0
        scatter(d_list(iiu),b_list(iiu),250,'filled','d');
    else
        scatter(d_list(iiu),b_list(iiu),50,'filled');
    end
    hold on
%         i=i+1;    
end
i=i+1;

% if sum(~converged)>0  % mark points where algorithm didn't converge
%     scatter(d_list(~converged),b_list(~converged),200,'*');
%     leg_str{i}='Didnt converge'; i=i+1;
% end

% d_list_new=[]; b_list_new=[];
% [d_list_new_, b_list_new_] = makegrid_in_rec([.92 .2], [.065,.11], 40 );
% d_list_new=[d_list_new; d_list_new_]; b_list_new=[b_list_new; b_list_new_]; 
% [d_list_new_, b_list_new_] = makegrid_in_rec([.975 .88], [.01,.07], 10 );
% d_list_new=[d_list_new; d_list_new_]; b_list_new=[b_list_new; b_list_new_]; 


% scatter(d_list_new,b_list_new,'filled');
% leg_str{i}='New points'; i=i+1;

hold off
legend(leg_str,'Location','northwestoutside');
xlabel('\delta')
ylabel b
xlim([.8 .98]);
ylim([.1 .9]);
% print2eps(['saves/' VERSION '/pics/grid.eps' ]);


%% interp2 --> scatter
figure

bb_fine=.1:.01:.9;
dd_fine=.7:.005:.98;
eqm_list_exist = double(eqm_list>0);

Int = scatteredInterpolant(d_list,b_list,eqm_list,'nearest');
G=gridmake(dd_fine',bb_fine'); dd_fine_list=G(:,1); bb_fine_list=G(:,2);
eqm_fine_list = Int(G); 


iiu = (abs(eqm_fine_list)<.01)&(dd_fine_list>=.8);
scatter(dd_fine_list(iiu),bb_fine_list(iiu),'filled')
xlim([.8 .98]);
ylim([.1 .9]);

%%

eqm_gridfit=interp2(d_list,b_list,eqm_list_exist,dd_fine,bb_fine,'nearest');
eqm_gridfit=round(eqm_gridfit);
G=gridmake(dd_fine,bb_fine); dd_fine_list=G(:,1); bb_fine_list=G(:,2);
eqm_fine_list = reshape(eqm_gridfit,[length(dd_fine_list) 1]);
iiu = (abs(eqm_fine_list)<.01)&(dd_fine_list>=.8);
scatter(dd_fine_list(iiu),bb_fine_list(iiu))
xlabel('\delta')
ylabel b
print2eps(['saves/' VERSION '/pics/eqm.eps' ]);


%% gridfit --> contourf
figure
bb_fine=.1:.01:.9;
dd_fine=.7:.005:.98;
eqm_list_exist = double(eqm_list==0);
eqm_gridfit=gridfit(d_list,b_list,eqm_list_exist,dd_fine,bb_fine,'smoothness',1);
% surf(dd,bb,eqm_gridfit);
contourf(dd_fine,bb_fine,eqm_gridfit,1)
xlabel('\delta')
ylabel b
print2eps(['saves/' VERSION '/pics/eqm.eps' ]);



%% Plot social welfare 

figure
SW_gridfit=gridfit(d_list,b_list,SW,dd_fine,bb_fine,'smoothness',20);
[c,h]=contourf(dd_fine,bb_fine,SW_gridfit,400);
set(h,'EdgeColor','none');
colorbar
xlabel('\delta')
ylabel b
print2eps(['saves/' VERSION '/pics/SW' ]);
% export_fig(['saves/' VERSION '/pics/SW' ],'-pdf');

%% Combine the existence region and the welfare 

figure
[c,h]=contourf(dd_fine,bb_fine,SW,200);
set(h,'EdgeColor','none');
colorbar
hold on
contour(dd_fine,bb_fine,eqm_gridfit,1,'LineWidth',5)
hold off
print2eps(['saves/' VERSION '/pics/existence-and-SW' ]);


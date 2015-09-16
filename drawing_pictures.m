clear
close all

%save(['saves/' VERSION '/eqm'],'SW','eqm_list','b_list','d_list','converged');

% VERSION='7'
load(['saves/7/eqm.mat']);
L8=load(['saves/8/eqm.mat']);

% substitute data from 8 into 7
[ind7, loc8]=ismember([b_list d_list],[L8.b_list L8.d_list],'rows');
loc8(loc8==0)=[];
SW(ind7)=L8.SW(loc8);
eqm_list(ind7)=L8.eqm_list(loc8);
converged(ind7)=L8.converged(loc8);




%% Create homogenous grid 
% 
bb_grid=.1:.025:.9;
dd_grid=.7:.005:.98;
Int = scatteredInterpolant(d_list,b_list,eqm_list,'nearest');
G=gridmake(dd_grid',bb_grid'); 
% next three lines overwrite the loaded data
d_list=G(:,1); 
b_list=G(:,2);
eqm_list = round(Int(G));


%% Smooth out eqm_list
N=length(d_list);
% convInf = converged;
% convInf(converged==0)=Inf;
% convInf(converged==1)=1;
maxd=max(d_list); mind=min(d_list);
maxb=max(b_list); minb=min(b_list);
rd = (maxd-mind); rb=(maxb-minb);
% dist=(repmat(d_list,1,N)-repmat(d_list'.*convInf',N,1)).^2/rd^2+(repmat(b_list,1,N)-repmat(b_list',N,1)).^2/rb^2;
dist=(repmat(d_list,1,N)-repmat(d_list',N,1)).^2/rd^2+(repmat(b_list,1,N)-repmat(b_list',N,1)).^2/rb^2;

%% 
for rep=1:2
    rep
    eqm_list0=eqm_list;
    for i=1:N
%         num_nb=0;
%         r=0;
%         while num_nb < 5
%             r=r+.01;
%             nb=(dist(i,:)<r);
%             num_nb=sum(nb);
%         end

        [ds,I]=sort(dist(i,:),'ascend');
        nb=I(2:9);

        tbl=tabulate(eqm_list0(nb));
        [perc,maxj]=max(tbl(:,3));
        myc = (tbl(:,1)==eqm_list0(i));

        if ~myc | tbl(myc,2)<=1
            pd = makedist('Multinomial','Probabilities',tbl(:,3)'/100);
            y = random(pd);
            eqm_list(i) = tbl(y,1);
        end
        if perc>= 60
            eqm_list(i)=tbl(maxj,1);
        end
    end
end
eqm_list_back=eqm_list;
%% Manual correction
eqm_list=eqm_list_back;
jj=(b_list<.3 & d_list<.9 & eqm_list==2);
eqm_list(jj)=0;
jj=(b_list<.3 & d_list<=.915 & eqm_list==6);
eqm_list(jj)=4;
jj=(b_list<.3 & d_list>.915 & eqm_list==6);
eqm_list(jj)=2;


%% Plot raw scattered points

au = unique(eqm_list);
% au = flipud(au);
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
        scatter(d_list(iiu),b_list(iiu),50,'filled','o');
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

% j=find(abs(b_list-.225)<.02 & abs(d_list-.894)<.002)
% [ds,I]=sort(dist(j,:),'ascend');
% nb=I(2:9);
% %         tbl=tabulate(eqm_list0(nb));
% %         [perc,maxj]=max(tbl(:,3));
% %         myc = (tbl(:,1)==eqm_list0(i));
% % %         if perc>= 60
% % %             eqm_list(i)=tbl(maxj,1);
% %         if ~myc || tbl(myc,2)<=1
% %             pd = makedist('Multinomial','Probabilities',tbl(:,3)'/100);
% %             y = random(pd);
% %             eqm_list(i) = tbl(y,1);
% %         end
% scatter(d_list(nb),b_list(nb),200,'*');
% leg_str{i}='---'; i=i+1;

hold off
legend(leg_str,'Location','northwestoutside');
xlabel('\delta')
ylabel b
xlim([.8 .98]);
ylim([0 1]);
% print2eps(['saves/' VERSION '/pics/grid.eps' ]);
return

print2eps(['saves/8/pics/grid.eps' ]);

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
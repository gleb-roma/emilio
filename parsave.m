function parsave(varargin)
savefile = varargin{1}; % first input argument
for i=2:nargin
    savevar.(inputname(i)) = varargin{i}; % other input arguments
end
save(savefile,'-struct','savevar')

% function [ ] = parsave( fname,  V , X , delta ,  b , pi , nd , n , w , g , gp , gp_stat , pstat , pstat2 , conv_flag , VE_cr  )
% 
%     save(filename, 'V','X','delta', 'b','pi','nd','n','w','g','gp','gp_stat','pstat','pstat2','conv_flag','VE_cr');
% 
% end


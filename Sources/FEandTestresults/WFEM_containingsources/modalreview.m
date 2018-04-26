function modalreview
% MODALREVIEW provides a really crude review of previous modal
% analysis results.
  
global restart
global K
global Ks
global M
global nodes
global beamprops
global matprops
global points
global elprops
global element
global Fepsn % Forces representing initial strain
global F     % Prescribed forces
global Fnln  % Forces due to changing coordinate systems. 
global lines
global surfs
global bcs   % Boundary Conditions
global slavedofs
global fs
global fms 
global filename


disp('The first 16 natural frequencies are (in Hz):')
format long
fs(1:min(16,length(fs)))
%length(fs)
%xfpause
%sort(abs(eig(K)))
q='f';
answer=input(['Do you want to survey all of the modes (s), pick ' ...
              'modes (p) or quit (q)? '],'s');
while ~strcmp('q',answer)
  if strcmp(answer,'s')
    for i=1:size(fms,2)
      plotgeom('deformed',nodes,lines,surfs,fms(:,i))
      title(['Mode shape ' num2str(i) ', ' num2str(fs(i)) ' Hz.'],'buttondownf','edtext')
      %  disp([num2str(i) 'th natural frequency is ' num2str(fs(i)) ' Hz.'])
      figure(gcf)
      pause
    end
  elseif strcmp(answer,'p')
    question=['What mode would you like to see (0 to goto top menu, ' ...
              num2str(size(fms,2)) ' max)? ']; 
    modenum=input(question);
    while length(modenum)~=1
     modenum=input(['Don''t be greedy. Pick just one mode. '])
    end
       
    while modenum~=0
      if modenum>size(fms,2)
        disp(['What part of ' num2str(size(fms,2)) ' max don''t you ' ...
              'understand?'])
        else
          plotgeom('deformed',nodes,lines,surfs,fms(:,modenum))
          title(['Mode shape ' num2str(modenum) ', ' num2str(fs(modenum)) ' Hz.'],'buttondownf','edtext')
          %  disp([num2str(i) 'th natural frequency is ' num2str(fs(i)) ' Hz.'])
          figure(gcf)
          answer=lower(input(['Would you like to animate this mode? '], ...
                       's'));
          if length(findstr(answer,'y'))>0
            [AZ,EL]=view;T=[AZ,EL];
            animatemode(nodes,lines,surfs,fms(:,modenum),T)
          end
      end
      modenum=input(question);
    end
  end
  answer=input(['Do you want to survey all of the modes (s), pick ' ...
                'modes (p) or quit (q)? '],'s');
end
classdef NR < marmodata.mdbase % vgsaccade.vgsaccade 
  % base class for analysis of the ocular following paradigm with motion
  % cloud stim
methods

  function v = getComplete(d) 
          % true for trials ending in SUCCESS in tarBr and mOnset
          [~,trial,~,state] = d.meta.fixbhv.state; 
 
          v = false([d.numTrials,1]); 
          ix = strcmpi(state,'SUCCESS'); 
          v(trial(ix)) = true;           
          iy = strcmpi(state,'FAIL'); 
          v(trial(iy)) = false;
  end 
end
end
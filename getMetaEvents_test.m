%m = d.meta;
srcNames = m.name2src.keys;
srcNames(cell2mat(m.name2src.values)) = srcNames;

prmNames = m.name2prm.keys;
prmNames(cell2mat(m.name2prm.values)) = prmNames;


id_all = [m.src(:),m.prm(:)];
id = unique(id_all,'rows');

prmList = unique(m.prm(:))';
srcList = unique(m.src(:))';

% for ip = prmList
%     numEventsPrm(ip) = sum(id_all(:,2) == ip);
% end
% [~, prm_sorted] = sort(numEventsPrm, 'descend');
% 
% prmNames(prm_sorted(1:20))
% 
% 
% for is = srcList
%     numEventsSrc(is) = sum(id_all(:,1) == is);
% end
% [~, src_sorted] = sort(numEventsSrc, 'descend');
% srcNames(src_sorted)

%# events per src
for is = srcList
    disp(srcNames{is})
    theseEvents = (id_all(:,1) == is);

    numEventsPrm = [];
    for ip = prmList
        numEventsPrm(ip) = sum(id_all(theseEvents,2) == ip);
    end
    [nEvents_sorted, prm_sorted] = sort(numEventsPrm, 'descend');

    for ipt = 1:10
        disp([prmNames{prm_sorted(ipt)}  ': ' num2str(nEvents_sorted(ipt))]);
    end
    disp('--');
end


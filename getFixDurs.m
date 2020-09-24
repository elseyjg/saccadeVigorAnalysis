function [fixDurByMag]=getFixDurs(TOD,whichAttribute, FixDuration, whichCol,mags, whichAttID)
for ii =1:length(mags)
    temp =[];
    trls = find(TOD(:,whichCol)==mags(ii));
    for jj =1:length(trls)
        idx=whichAttribute(trls(jj),:)==whichAttID;
        if any(idx)
            temp=[temp,FixDuration(trls(jj),idx)];
        end
    end
    fixDurByMag(ii,1:length(temp))=temp;
    nfix(ii,1)= length(temp);
    
end
for ii =1:length(nfix)
    if nfix(ii)<size(fixDurByMag,2)
        fixDurByMag(ii,nfix(ii)+1:end)=nan;
    end
end
end
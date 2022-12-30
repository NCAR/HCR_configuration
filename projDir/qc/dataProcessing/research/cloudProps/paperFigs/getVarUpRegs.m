function [allVars,groups,nums]=getVarUpRegs(plotV,varIn)

minArea=0.1;

allVars=[];
groups={};
nums=[];

varName='upRegsAll';

thisTable=plotV.(varName).socrates.CloudLow;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudLowS'},length(thisVar),1));

thisTable=plotV.(varName).otrec.CloudLow;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudLowO'},length(thisVar),1));

thisTable=plotV.(varName).cset.ConvYoungShallow;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvShallowC'},length(thisVar),1));

thisTable=plotV.(varName).socrates.ConvYoungShallow;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvShallowS'},length(thisVar),1));

thisTable=plotV.(varName).otrec.ConvYoungShallow;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvShallowO'},length(thisVar),1));

thisTable=plotV.(varName).cset.ConvMatureShallow;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvStratShallowC'},length(thisVar),1));

thisTable=plotV.(varName).socrates.ConvMatureShallow;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvStratShallowS'},length(thisVar),1));

thisTable=plotV.(varName).cset.StratShallow;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratShallowC'},length(thisVar),1));

thisTable=plotV.(varName).socrates.StratShallow;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratShallowS'},length(thisVar),1));

thisTable=plotV.(varName).socrates.CloudMid;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudMidS'},length(thisVar),1));

thisTable=plotV.(varName).otrec.CloudMid;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudMidO'},length(thisVar),1));

thisTable=plotV.(varName).socrates.ConvYoungMid;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvMidS'},length(thisVar),1));

thisTable=plotV.(varName).otrec.ConvYoungMid;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvMidO'},length(thisVar),1));

thisTable=plotV.(varName).socrates.StratMid;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratMidS'},length(thisVar),1));

thisTable=plotV.(varName).otrec.StratMid;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratMidO'},length(thisVar),1));

thisTable=plotV.(varName).cset.CloudHigh;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudHighC'},length(thisVar),1));

thisTable=plotV.(varName).socrates.CloudHigh;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudHighS'},length(thisVar),1));

thisTable=plotV.(varName).otrec.CloudHigh;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudHighO'},length(thisVar),1));

thisTable=plotV.(varName).otrec.ConvYoungDeep;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvDeepO'},length(thisVar),1));

thisTable=plotV.(varName).otrec.ConvMatureDeep;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvStratDeepO'},length(thisVar),1));

thisTable=plotV.(varName).otrec.StratDeep;
thisTable(thisTable.area<minArea,:)=[];
thisVar=thisTable.(varIn);
if strcmp(varIn,'cloudAltPerc')
    thisVar=abs(thisVar-100);
end
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratDeepO'},length(thisVar),1));

end
function [allVars,groups,nums]=getVarAll(plotV,varName)
allVars=[];
groups={};
nums=[];

thisVar=plotV.(varName).cset.CloudLow;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudLowC'},length(thisVar),1));

thisVar=plotV.(varName).socrates.CloudLow;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudLowS'},length(thisVar),1));

thisVar=plotV.(varName).otrec.CloudLow;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudLowO'},length(thisVar),1));

thisVar=plotV.(varName).cset.ConvYoungShallow;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvShallowC'},length(thisVar),1));

thisVar=plotV.(varName).socrates.ConvYoungShallow;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvShallowS'},length(thisVar),1));

thisVar=plotV.(varName).otrec.ConvYoungShallow;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvShallowO'},length(thisVar),1));

thisVar=plotV.(varName).cset.ConvMatureShallow;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvStratShallowC'},length(thisVar),1));

thisVar=plotV.(varName).socrates.ConvMatureShallow;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvStratShallowS'},length(thisVar),1));

thisVar=plotV.(varName).otrec.ConvMatureShallow;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvStratShallowO'},length(thisVar),1));

thisVar=plotV.(varName).cset.StratShallow;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratShallowC'},length(thisVar),1));

thisVar=plotV.(varName).socrates.StratShallow;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratShallowS'},length(thisVar),1));

thisVar=plotV.(varName).otrec.StratShallow;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratShallowO'},length(thisVar),1));

thisVar=plotV.(varName).socrates.CloudMid;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudMidS'},length(thisVar),1));

thisVar=plotV.(varName).otrec.CloudMid;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudMidO'},length(thisVar),1));

thisVar=plotV.(varName).socrates.ConvYoungMid;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvMidS'},length(thisVar),1));

thisVar=plotV.(varName).otrec.ConvYoungMid;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvMidO'},length(thisVar),1));

thisVar=plotV.(varName).socrates.StratMid;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratMidS'},length(thisVar),1));

thisVar=plotV.(varName).otrec.StratMid;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratMidO'},length(thisVar),1));

thisVar=plotV.(varName).cset.CloudHigh;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudHighC'},length(thisVar),1));

thisVar=plotV.(varName).socrates.CloudHigh;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudHighS'},length(thisVar),1));

thisVar=plotV.(varName).otrec.CloudHigh;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudHighO'},length(thisVar),1));

thisVar=plotV.(varName).otrec.ConvYoungDeep;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvDeepO'},length(thisVar),1));

thisVar=plotV.(varName).otrec.ConvMatureDeep;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvStratDeepO'},length(thisVar),1));

thisVar=plotV.(varName).otrec.StratDeep;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratDeepO'},length(thisVar),1));



end
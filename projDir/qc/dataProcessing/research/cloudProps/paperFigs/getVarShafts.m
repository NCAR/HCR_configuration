function [allVars,groups,nums]=getVarShafts(plotV,varNameIn)
allVars=[];
groups={};
nums=[];

varName='precShaftsAll';

thisVar=plotV.(varName).cset.ConvYoungShallow.(varNameIn);
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvShallowC'},length(thisVar),1));

thisVar=plotV.(varName).socrates.ConvYoungShallow.(varNameIn);
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvShallowS'},length(thisVar),1));

thisVar=plotV.(varName).otrec.ConvYoungShallow.(varNameIn);
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvShallowO'},length(thisVar),1));

thisVar=plotV.(varName).cset.ConvMatureShallow.(varNameIn);
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvStratShallowC'},length(thisVar),1));

thisVar=plotV.(varName).socrates.ConvMatureShallow.(varNameIn);
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvStratShallowS'},length(thisVar),1));

thisVar=plotV.(varName).otrec.ConvMatureShallow.(varNameIn);
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvStratShallowO'},length(thisVar),1));

thisVar=plotV.(varName).cset.StratShallow.(varNameIn);
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratShallowC'},length(thisVar),1));

thisVar=plotV.(varName).socrates.StratShallow.(varNameIn);
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratShallowS'},length(thisVar),1));

thisVar=plotV.(varName).otrec.StratShallow.(varNameIn);
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratShallowO'},length(thisVar),1));

thisVar=plotV.(varName).socrates.ConvYoungMid.(varNameIn);
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvMidS'},length(thisVar),1));

thisVar=plotV.(varName).otrec.ConvYoungMid.(varNameIn);
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvMidO'},length(thisVar),1));

thisVar=plotV.(varName).socrates.StratMid.(varNameIn);
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratMidS'},length(thisVar),1));

thisVar=plotV.(varName).otrec.StratMid.(varNameIn);
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratMidO'},length(thisVar),1));

thisVar=plotV.(varName).otrec.ConvYoungDeep.(varNameIn);
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvDeepO'},length(thisVar),1));

thisVar=plotV.(varName).otrec.ConvMatureDeep.(varNameIn);
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvStratDeepO'},length(thisVar),1));

thisVar=plotV.(varName).otrec.StratDeep.(varNameIn);
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratDeepO'},length(thisVar),1));



end
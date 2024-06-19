%%  MultiPartialMap    Apply Superoperator Phi on multiple subsystem
%   This function has 4 required arguments:
%     X: a matrix
%     PHI: a superoperator
%     SYS: a subset of 1:numel(DIM) that PHI is acting on.
%     DIM: the array of subsystems [d1,d2,...]
%     
%   PHIX = MultiPartialMap(X,Phi,SYS,DIM) is the operator PHI(X), 
%   in which Phi partially acting on the subsystem SYS.
% 
%   requires: PartialMap.m, PermuteSystems.m
%   author: Bikun Li (bikunli@uchicago.edu)
%   package: QETLAB
%   last updated: June 16, 2024

function [PhiX] = MultiPartialMap(X,Phi,SYS,DIM)
% Apply Superoperator Phi on multiple subsystem
% SYS is a subset of 1:numel(DIM). The original ParitialMap only allows a
% scalar SYS, now this function can adapt to an SYS array

if numel(SYS) == 1
    PhiX = PartialMap(X,Phi,SYS,DIM);
else
    SYS = SYS(:);
    SYS_c = setdiff(1:numel(DIM), SYS); % the complementary system
    PERM = [SYS(:); SYS_c(:)].'; 
    DIM_temp = [prod(DIM(SYS(:).')), DIM(SYS_c(:).')]; 
    X_temp = PermuteSystems(X,PERM,DIM,0,0); % permute all subsystems to the front
    X_temp = PartialMap(X_temp, Phi, 1, DIM_temp); % apply partial map
    PhiX = PermuteSystems(X_temp,PERM,DIM,0,1); % permute all subsystems back  
end
end
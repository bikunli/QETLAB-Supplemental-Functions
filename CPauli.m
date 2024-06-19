%% CPauli    Generates a unitary operator for controlled-Pauli gate 
%   This function has 4 required argument:
%     i: the controlling qubit
%     j: the target qubit
%     p: the type of Pauli operator, e.g. 1 or 'X'
%     DIM: the array of subsystems [2,2,....,2]
% 
%   U = CPauli(i,j,p,DIM) is a unitary operator for controlled-Pauli gate.
%   For example, CPauli(1,2,'Z',[2,2]) is the controlled-Z gate
% 
%   requires: Pauli.m
%   author: Bikun Li (bikunli@uchicago.edu)
%   package: QETLAB
%   last updated: June 16, 2024

function U = CPauli(i,j,p,DIM)
% control site: i
% target site: j
% p: an element in {0,1,2,3} or {'I','X','Y','Z'}
% DIM: the array for the subsystems' dimensionality, e.g. [2,2]

if ~ismember(p,[0,1,2,3]) && ~ismember(p, {'I','X','Y','Z'}) 
    error('p has to be a member of {1,2,3}')
end
if ~ismember(p,[0,1,2,3]) 
    p = find(p=='IXYZ')-1;
end

v_i = zeros(size(DIM));
v_j = zeros(size(DIM));
v_i(i) = 1;
v_j(j) = 1;
U = (Pauli(0*v_i) + Pauli(3*v_i))/2 + (Pauli(0*v_i) - Pauli(3*v_i))/2 * Pauli(p*v_j);
end
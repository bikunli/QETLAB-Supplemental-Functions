%%  CSwap    Generates a unitary operator for controlled-SWAP
%   This function has 4 required arguments:
%     i: the controlling qubit
%     j and k: the qubit pair to be swapped
%     DIM: the array of subsystems [2,2,....,2]
% 
%   CSwap(i,j,k,DIM) is a unitary operator acting on an array of qubits. 
%   i is the controlling qubit, j and k are two qubits to be swapped if
%   qubit i is in |1> state.
% 
%   requires: Pauli.m, Swap.m
%   author: Bikun Li (bikunli@uchicago.edu)
%   package: QETLAB
%   last updated: June 16, 2024

function [U] = CSwap(i,j,k,DIM)

v_i = zeros(size(DIM));
v_j = zeros(size(DIM));
v_k = zeros(size(DIM));
v_i(i) = 1;
v_j(j) = 1;
v_k(k) = 1;
S = sparse(Swap(Pauli(0*v_i),[j,k],DIM,1));
U = (Pauli(0*v_i,1) + Pauli(3*v_i,1))/2 + (Pauli(0*v_i,1) - Pauli(3*v_i,1))/2 *  S;
end
function [A,M]=assembly_surface(Nodes,Elements)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix an r.h.s assembly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% degree of basis functions
%% element-by-element
dof=length(Nodes);

%mátrixok
A=spalloc(dof, dof, dof*6);
M=spalloc(dof, dof, dof*6);
% w=zeros(dof,1);
%% matrices on reference element
forma=[1/3 1/6
       1/6 1/3];
forma_dx=[1  -1
         -1   1];
     
% degree of polinomials on each element
poli_deg=ones(size(Elements,1),1);

%% assembly using a reference element
for i_element=1:size(Elements,1)    
    a=Nodes(Elements(i_element,1),:);
    b=Nodes(Elements(i_element,poli_deg(i_element)+1),:);
    nrm=norm(b-a);
    
    SLoc_mtx=(1/nrm)*forma_dx;
    MLoc_mtx=nrm*forma;
    
    
    A(Elements(i_element,1:2),Elements(i_element,1:2))=A(Elements(i_element,1:2),Elements(i_element,1:2))+SLoc_mtx;
    M(Elements(i_element,1:2),Elements(i_element,1:2))=M(Elements(i_element,1:2),Elements(i_element,1:2))+MLoc_mtx;
    
end
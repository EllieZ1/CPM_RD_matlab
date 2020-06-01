adj_empty = ~cell_mask(up) | ~cell_mask(down) | ~cell_mask(left) | ~cell_mask(right); %cells at the boundry 
adj_full = cell_mask(up) | cell_mask(down) | cell_mask(left) | cell_mask(right); %cells at the boundry 

%finding boundry points 
bndry_cell = cell_mask & adj_empty;
bndry_empty = ~cell_mask & adj_full;
bndry = find( bndry_cell | bndry_empty );

if any(cell_maskp~=cell_mask)
    error('not reseting')
end
ij = bndry(randi(length(bndry)));

r=randi(4);
cell_maskp(ij) = cell_mask(jump(sub2ind([sz,4],ij,r)));% make a new trial configuration

Per=perim(cell_maskp); % perimter
A=nnz(cell_maskp); % area
HA=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; % the hamiltonian after the possible change
dH=HA-H0;
Ncell_mask=squeeze(sum(sum(x))); %for a sanity check 
i=1;diffused=1;rx=0;xi0=0;neg=0;%to call update.m
if getfield(bwconncomp(cell_maskp,4),'NumObjects')==1 %makes sure the cell stays connected 
    grow= cell_maskp(ij) & ~cell_mask(ij);
    shrink= ~cell_maskp(ij) & cell_mask(ij);
    if grow
        dH=dH+B_rho*(RhoRatio(jump(sub2ind([sz,4],ij,r)))-rho_eq)-B_R*(RacRatio(jump(sub2ind([sz,4],ij,r)))-R_eq);
        if rand<exp(-(dH+Hb)/T) %step raises energy boltzman prob forward
            cell_mask=cell_maskp; %changing cell shape
            
            for j=0:(N_species-1) %splitting the molecules with the new lattice
                x(ij+j*sz)=floor(x(jump(sub2ind([sz,4],ij,r))+j*sz)/2);
                x(jump(sub2ind([sz,4],ij,r))+j*sz)=ceil(x(jump(sub2ind([sz,4],ij,r))+j*sz)/2);
            end
            
            
            
            I=[ij jump(sub2ind([sz,4],ij,r))]; %places where molecule number has changed
            
            H0=HA; %changing the hamiltonn to the new one
            
            %similar idea now to the CPM_chem_func dependence tree
            A=nnz(cell_mask);
            cell_inds(1:A)=find(cell_mask);
            [x,A,sz,diffusing_species_sum,alpha_rx,...
            alpha_chem,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
            RacRatio,RbarRatio,I_Ks,N_instantaneous,reaction,ir0,id0,cell_inds,...
            k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
            alpha,PAKtot,rx,diffused,i,I,xi0,neg]=update(x,A,sz,diffusing_species_sum,alpha_rx,...
            alpha_chem,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
            RacRatio,RbarRatio,I_Ks,N_instantaneous,reaction,ir0,id0,cell_inds,...
            k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
            alpha,PAKtot,rx,diffused,i,I,xi0,neg);
            
            alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
            
            
        else
            cell_maskp=cell_mask;
        end
    elseif shrink
        dH=dH-B_rho*(RhoRatio(ij)-rho_eq)+B_R*(RacRatio(ij)-R_eq);
        
        if rand<exp(-(dH+Hb)/T)%step raises energy
            cell_mask=cell_maskp;  %changing cell shape 
            
            neighbors=[jump(ij,1) jump(ij,2) jump(ij,3) jump(ij,4)]; %finding places the molecules will go to
            neighbors=neighbors(find(cell_maskp(neighbors)));
            
            for j=0:(N_species-1) %dumping out molecules from the retracting site
                x(neighbors+j*sz)=x(neighbors +j*sz)+diff(round(linspace(0,x(ij+j*sz),length(neighbors)+1)));
                x(ij+j*sz)=0;
            end
            
            I=[ij neighbors]; %indices where molecule number changed 
            
            
            %Same as grow
            H0=HA;

            A=nnz(cell_mask);
            cell_inds(1:A)=find(cell_mask);
            [x,A,sz,diffusing_species_sum,alpha_rx,...
            alpha_chem,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
            RacRatio,RbarRatio,I_Ks,N_instantaneous,reaction,ir0,id0,cell_inds,...
            k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
            alpha,PAKtot,rx,diffused,i,I,xi0,neg]=update(x,A,sz,diffusing_species_sum,alpha_rx,...
            alpha_chem,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
            RacRatio,RbarRatio,I_Ks,N_instantaneous,reaction,ir0,id0,cell_inds,...
            k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
            alpha,PAKtot,rx,diffused,i,I,xi0,neg);
            
            alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
            
        else
            %reaction doesn't happen
            cell_maskp=cell_mask;
        end
    else
        cell_maskp=cell_mask;
    end
else
    cell_maskp=cell_mask;
end

Ncell_maskp=squeeze(sum(sum(x)));
%sanity checks 
if any(Ncell_mask~=Ncell_maskp)
    errror('molecule loss')
end

if min(cell_mask(:))<0
    error('Oh no! D: (negtive numbers)')
end


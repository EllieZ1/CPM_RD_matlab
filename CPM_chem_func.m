function [x,diffusing_species_sum,D,h,alpha_rx,num_diffuse,...
    alpha_chem,time,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
    RacRatio,RbarRatio,I_Ks,reaction,ij_diffuse,jump,ir0,id0,cell_inds,...
    k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
    alpha,PAKtot,numDiff,numReac] =  CPM_chem_func(x,diffusing_species_sum,D,h,alpha_rx,num_diffuse,...
    alpha_chem,time,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
    RacRatio,RbarRatio,I_Ks,reaction,ij_diffuse,jump,ir0,id0,cell_inds,...
    k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
    alpha,PAKtot,nrx,A,numDiff,numReac)


N_instantaneous=50; % the number of steady reaction itterated at a moment in time 
sz=size(x,1)*size(x,2);

%variables declared so that the c encoder is happy 
i=1;
p=1;
i2=1;
rx=1;
xi0=x(id0+i(1));
I=[i i2];
for kk=1:nrx
    [t1,t2]=meshgrid(ir0,cell_inds(1:A)');
    t3=t1+t2;
    alpha_rx=sum(alpha_chem(t3));
    
    %total propensity to diffuse for all species
    alpha_diff=sum(diffusing_species_sum).*D/(h*h);
    
    
    %total propensity for rxn+diff
    a_total=sum(alpha_diff)+sum(alpha_rx(:));
    
    tau  = (1/a_total)*log(1/rand()); % time increment
    time=time+tau;
    reacted=false;
    diffused=false;
    neg=false;
    RN=rand();

    if RN*a_total <= sum(alpha_diff(:))
        
        diffused=true; 

        %-------------diffusion-----------
        p=find((RN*a_total<=cumsum(alpha_diff)),1); %which protein diffuses
        p=p(1); %for encoder 
        i0=(p-1)*sz; % the vectorized location of the appropriate proteins 

        %adding all previous reactions to the propensity 
        if p==1
            temp=0;
        else
            temp=sum(alpha_diff(1:p-1));
        end
        
        for drx=1:length(num_diffuse) %iterate over possible diffusion directions
            %check if it will diffuse in this direction
            if temp+(D(p)/h^2)*diffusing_species_sum(drx,p)>RN*a_total&&~reacted
                
                %find the point that diffuses
                ii=find((D(p)/h^2)*cumsum(x(i0+ij_diffuse(drx,1:num_diffuse(drx)')))>RN*a_total-temp,1);
                i(1)=ij_diffuse(drx,ii);
                %find the point they diffuses to
                i2=jump(i,drx);
                I=[i i2];
                
                %carry out the diffusion reaction
                x(i0+i) = x(i0+i)-1;
                x(i0+i2) = x(i0+i2)+1;
                %update the sum of diffusing species
                %in case you have diffused to an edge
                diffusing_species_sum(:,p)=diffusing_species_sum(:,p)+(diffuse_mask(:,i2)-diffuse_mask(:,i));
                
                reacted=true;
                neg=any(x(i+i0)<0);
                
            elseif ~reacted
                %keep adding
                temp=temp+(D(p)/h^2)*diffusing_species_sum(drx,p);
            else
                break;
            end
        end
        
        if reacted==false %making sure the code works
            error('Oh no! D: diffusion propensites did not sum correctly')
        end
        numDiff=numDiff+1;%ellie
    else% ---------- Reaction Time -------------
        reacted=false;      
        temp=sum(alpha_diff(:));
        %determing where the reaction happens 
        for rx=1:size(alpha_chem,3)
            i0=(rx-1)*sz;
            if temp+alpha_rx(rx)>=RN*a_total&&~reacted
                ii=find(cumsum(alpha_chem(i0+cell_inds(1:A)))>=RN*a_total-temp,1);
                i(1)=cell_inds(ii);
                I=i;
                xi0=x(id0+i(1));
                reacted = true;
            elseif ~reacted
                temp=temp+alpha_rx(rx);
            end
            if reacted
                break;
            end
        end
        
        if reacted==false
            error('Oh no! D: chemical reaction propensites did not sum correctly')
        end
        numReac=numReac+1;%ellie
        %Inactive rho to active rho
        if rx==1
            x(i+(3-1)*sz) = x(i+(3-1)*sz)+1;
            x(i+(1-1)*sz) = x(i+(1-1)*sz)-1;
        end
        
        %Inactive Rac to active Rac
        if rx==2
            x(i+(2-1)*sz)=x(i+(2-1)*sz)-1;
            x(i+(4-1)*sz)=x(i+(4-1)*sz)+1;
        end
        
        %Active rho to inactive rho
        if rx==3
            x(i+(3-1)*sz) = x(i+(3-1)*sz)-1;
            x(i+(1-1)*sz) = x(i+(1-1)*sz)+1;
        end
        
        %Active Rac to inactive Rac
        if rx==4
            x(i+(2-1)*sz)=x(i+(2-1)*sz)+1;
            x(i+(4-1)*sz)=x(i+(4-1)*sz)-1;
        end
        
        %Unphosphorylated Pax to phosphorylated Pax
        if rx==5
            x(i+(5-1)*sz)=x(i+(5-1)*sz)-1;
            x(i+(6-1)*sz)=x(i+(6-1)*sz)+1;
        end
        
        %Phosphorylated Pax to unphosphorylated
        if rx==6
            x(i+(5-1)*sz)=x(i+(5-1)*sz)+1;
            x(i+(6-1)*sz)=x(i+(6-1)*sz)-1;
        end
        % I moved the following 2 lines into my dependence tree due to
        % complexing (should uncomment if replicating this algorthim)
        %         dxi=xi0-x(id0+i(1));
        %         diffusing_species_sum = diffusing_species_sum - (diffuse_mask(:,i)*dxi);
        %         neg=x(i+(rx-1)*sz)<0;
    end
    
    %-----------recalculating value that would have changed-------------- 
    %the reaction trees are seprate as 2 cells change when diffusion
    %happens
    
    [x,A,sz,diffusing_species_sum,alpha_rx,...
    alpha_chem,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
    RacRatio,RbarRatio,I_Ks,N_instantaneous,reaction,ir0,id0,cell_inds,...
    k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
    alpha,PAKtot,rx,diffused,i,I,xi0,neg]=update(x,A,sz,diffusing_species_sum,alpha_rx,...
    alpha_chem,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
    RacRatio,RbarRatio,I_Ks,N_instantaneous,reaction,ir0,id0,cell_inds,...
    k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
    alpha,PAKtot,rx,diffused,i,I,xi0,neg);
    
    if neg
        error('Oh no! D: (negtive numbers)')
    end
    
end

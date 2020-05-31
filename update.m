function [x,A,sz,diffusing_species_sum,alpha_rx,...
    alpha_chem,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
    RacRatio,RbarRatio,I_Ks,N_instantaneous,reaction,ir0,id0,cell_inds,...
    k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
    alpha,PAKtot,rx,diffused,i,I,xi0,neg]=update(x,A,sz,diffusing_species_sum,alpha_rx,...
    alpha_chem,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
    RacRatio,RbarRatio,I_Ks,N_instantaneous,reaction,ir0,id0,cell_inds,...
    k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
    alpha,PAKtot,rx,diffused,i,I,xi0,neg)
if ~diffused 
    %deal with complexing and decomplexing

    if rx==2||rx==4%Rac
        K_R=(1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(i))*alpha*PAKtot*K_is(i);
        for j=1:N_instantaneous
            if x(i+(7-1)*sz)/(x(i+(7-1)*sz)+x(i+(4-1)*sz)*K_R)>rand()%decomplex
                x(i+(4-1)*sz)=x(i+(4-1)*sz)+1;
                x(i+(7-1)*sz)=x(i+(7-1)*sz)-1;
            elseif x(i+(4-1)*sz)>0
                x(i+(4-1)*sz)=x(i+(4-1)*sz)-1;
                x(i+(7-1)*sz)=x(i+(7-1)*sz)+1;
            end
        end
    elseif rx==5||rx==6%Pax
        K_P=k_G*k_X*k_C*GIT*PIX*K_is(i)*PAKtot*(1+alpha_R*RacRatio(i));
        for j=1:N_instantaneous
            if x(i+(8-1)*sz)/(x(i+(8-1)*sz)+x(i+(6-1)*sz)*K_P)>rand()%decomplex
                x(i+(6-1)*sz)=x(i+(6-1)*sz)+1;
                x(i+(8-1)*sz)=x(i+(8-1)*sz)-1;
            elseif x(i+(6-1)*sz)>0
                x(i+(8-1)*sz)=x(i+(8-1)*sz)+1;
                x(i+(6-1)*sz)=x(i+(6-1)*sz)-1;
            end
        end
    end
    dxi=xi0-x(id0+i(1));
    diffusing_species_sum = diffusing_species_sum - (diffuse_mask(:,i)*dxi);
    neg=x(i+(rx-1)*sz)<0;
end
    
%update Ratios
RacRatio(I)=x(I+(4-1)*sz)./(x(I+(4-1)*sz)+x(I+(2-1)*sz)+x(I+(7-1)*sz));
RbarRatio(I)=x(I+(7-1)*sz)./(x(I+(4-1)*sz)+x(I+(2-1)*sz)+x(I+(7-1)*sz));
RhoRatio(I)=x(I+(3-1)*sz)./(x(I+(3-1)*sz)+x(I+(1-1)*sz));
PaxRatio(I)=x(I+(6-1)*sz)./(x(I+(6-1)*sz)+x(I+(5-1)*sz)+x(I+(8-1)*sz));

RacRatio(isnan(RacRatio))=0;
RbarRatio(isnan(RbarRatio))=0;
RhoRatio(isnan(RhoRatio))=0;
PaxRatio(isnan(PaxRatio))=0;
%update other parameters    
K_is(I)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(I)).*(1+alpha_R*RacRatio(I))+k_G*k_X*GIT*PIX);
K(I)=alpha_R*RacRatio(I).*K_is(I).*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(I));%RbarRatio(I)/gamma;         %changed from paper
I_Ks(I)=I_K*(1-K_is(I).*(1+alpha_R*RacRatio(I)));
reaction(I+(1-1)*sz) = I_rho*(L_R^m./(L_R^m +(RacRatio(I)+RbarRatio(I)).^m));            %From inactive rho to active rho changed from model
reaction(I+(2-1)*sz) = (I_R+I_Ks(I)).*(L_rho^m./(L_rho^m+RhoRatio(I).^m));                %From inactive Rac to active Rac
reaction(I+(5-1)*sz) = B_1*(K(I).^m./(L_K^m+K(I).^m));
[tmp,tmp2]=meshgrid(ir0,I);
tmp3=tmp+tmp2;
alpha_chem(tmp3) = reaction(tmp3).*x(tmp3);


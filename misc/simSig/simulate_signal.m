function dwi=simulate_signal(dwi,wm,gm,csf,prm)

b=prm.b; %b-value


%computing masks for parital volume effect
wm_gm=gm.*wm;
wm_csf=csf.*wm;
wm_gm_csf=wm.*csf.*gm;
gm_csf=csf.*gm;

%WM voxels info for simulating signal
check=find(wm);

%tissue-specific signal-WM
sigWM=exp(-prm.TE/prm.wm_T2)*(1-exp(-prm.TR/prm.wm_T1));
%tissue-specific signal-GM
sigGM=exp(-prm.TE/prm.gm_T2)*(1-exp(-prm.TR/prm.gm_T1)).*exp(diag(-prm.b.*(prm.g'*prm.diff_gm*prm.g)));
%tissue-specific signal-csf
sigCSF=exp(-prm.TE/prm.csf_T2)*(1-exp(-prm.TR/prm.csf_T1)).*exp(diag(-prm.b.*(prm.g'*prm.diff_csf*prm.g)));

for i=1:length(check)
    
[x,y,z] = ind2sub(size(wm),check(i));

%if all the three tissues are present
if wm_gm_csf(x,y,z)~=0
    

dwi(x,y,z,:)=(1-prm.fu_gm)*sigWM*dwi(x,y,z,:) + ...
    reshape( (prm.fu_gm*sigGM)  + (prm.fu_csf*sigCSF),[1,1,1,length(b)]);


%for WM-GM 
elseif wm_gm(x,y,z)~=0


      dwi(x,y,z,:)=(1-prm.fu_gm)*sigWM*dwi(x,y,z,:) + reshape( prm.fu_gm*sigGM,[1,1,1,length(b)]);

%for WM-CSF     
elseif wm_csf(x,y,z)~=0



  dwi(x,y,z,:)=(1-prm.fu_csf)*sigWM*dwi(x,y,z,:) + reshape(prm.fu_csf*sigCSF,[1,1,1,length(b)]);

%only WM    
else
    dwi(x,y,z,:)=sigWM*dwi(x,y,z,:);

    

end

end

%simulating signal for gray-matter
check=find(gm | gm_csf);

for i=1:length(check)

[x,y,z] = ind2sub(size(wm),check(i));
%only GM
if wm_gm_csf(x,y,z)==0 && wm_gm(x,y,z)==0 && wm_csf(x,y,z)==0 && gm_csf(x,y,z)==0 && gm(x,y,z)==1
    
    dwi(x,y,z,:)=sigGM;

%for GM-CSF    
elseif wm_gm_csf(x,y,z)==0 && wm_gm(x,y,z)==0 && wm_csf(x,y,z)==0 && gm_csf(x,y,z)==1
    
    dwi(x,y,z,:)=0.5*sigGM + 0.5*sigCSF;

end
    
end

%simulating signal for CSF
check=find(csf);

for i=1:length(check)

[x,y,z] = ind2sub(size(csf),check(i));
    
if wm_gm_csf(x,y,z)==0 && wm_gm(x,y,z)==0 && wm_csf(x,y,z)==0 && csf(x,y,z)==1
    
    dwi(x,y,z,:)=sigCSF;
    
end
    
end

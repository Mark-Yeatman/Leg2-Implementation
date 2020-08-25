function makeMatlabFunctions()    
    %Which physical parameters to use
    lt = 0;
    ls = 0;
    syms Ms Mf Isz Ifz lt ls lc la rs rfy rfx COPfx g 
    params = [Ms,Mf,Isz,Ifz,lt,ls,lc,la,rs,rfy,rfx,COPfx,g];
    
    x = sym('x',[4,1]);
    
    %convert cos and sin to cosd and sind, so that you're using degrees
    filenames = {'M_matrix.m',...
    'C_matrix.m',...
    'G_matrix.m',...
    'Length.m',...
    'LengthJacobian.m',...
    'LengthJacobianDt.m',...
    'LoadCellJacobian.m',...
    'Theta.m',...
    'ThetaJacobian.m',...
    'ThetaJacobianDt.m'};
    for i=1:length(filenames)        
        fid  = fopen(filenames{i},'r');
        f=fread(fid,'*char')';
        fclose(fid);
        f = strrep(f,'cos(','cosd(');
        f = strrep(f,'sin(','sind(');
        fid  = fopen(filenames{i},'w');
        fprintf(fid,'%s',f);
        fclose(fid);
    end
    
    M_matrix
    C_matrix
    G_matrix
    Length
    LengthJacobian
    LengthJacobianDt
    LoadCellJacobian
    Theta
    ThetaJacobian
    ThetaJacobianDt
    
    folder = '';
    matlabFunction(M,   'File',strcat(folder,'M_func'),'Vars',{x,params})
    matlabFunction(Cmat,'File',strcat(folder,'C_func'),'Vars',{x,params})
    matlabFunction(G,   'File',strcat(folder,'G_func'),'Vars',{x,params})
    
    matlabFunction(L,          'File',strcat(folder,'L_func'),'Vars',{x,params})
    matlabFunction(JacobL,     'File',strcat(folder,'J_L_func'),'Vars',{x,params})
    matlabFunction(JacobLDot, 'File',strcat(folder,'J_dot_L_func'),'Vars',{x,params})
    
    matlabFunction(THETA,      'File',strcat(folder,'Theta_func'),'Vars',{x,params})
    matlabFunction(Jacobtheta, 'File',strcat(folder,'J_Theta_func'),'Vars',{x,params})
    matlabFunction(JacobthetaDot, 'File',strcat(folder,'J_dot_Theta_func'),'Vars',{x,params})
    
    matlabFunction(JC, 'File',strcat(folder,'J_C_func'),'Vars',{x,params})  
end



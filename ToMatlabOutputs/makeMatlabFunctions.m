function makeMatlabFunctions()    
    %Which physical parameters to use
    lt = 0;
    ls = 0;
    syms Ms Mf Isz Ifz lt ls lc la rs rfy rfx COPfx g 
    params = [Ms,Mf,Isz,Ifz,lt,ls,lc,la,rs,rfy,rfx,COPfx,g];
    
    x = sym('x',[4,1]);
    
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



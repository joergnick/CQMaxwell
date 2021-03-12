function [u_magn ] = draw_scatterer( u_ges, scatterer)



n_grid=length(scatterer(:,1));
N=length(u_ges(1,:));

u_magnet=zeros(n_grid^2,N);

%original_ngrid=n_grid/2+0.5;
original_ngrid=n_grid;
u_sq=zeros(original_ngrid,original_ngrid);

u_long=zeros(n_grid^2,1);
u_long_changed=zeros(n_grid^2,1);
for j=1:N
  u_long=u_ges(:,j);
  
    for i=1:original_ngrid
        u_sq(:,i)=u_long((i-1)*original_ngrid+1:i*original_ngrid);
        
    end
%     
%     u_sq=double_resolution(u_sq);
    
    for ind_x=1:n_grid
        for ind_y=1:n_grid
            if scatterer(ind_x,ind_y)==1
            u_sq(ind_x,ind_y)=inf;
            end
        end
    end
    
    
    for i=1:n_grid
       u_long_changed((i-1)*n_grid+1:i*n_grid)=u_sq(:,i);
    end
    
    
    u_magn(:,j)=u_long_changed;
end

    

end


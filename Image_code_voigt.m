f = @(x, y, p) voigt2d_aniso_rot(x, y, p);

%electron repulsion term
vx_rand=rep_std*randn(MAX,1); % Draw MAX samples from a normal distribution with standard deviation rep_std.
vy_rand=rep_std*randn(MAX,1);
vz_rand=rep_std*randn(MAX,1);
Cell_velocity_el=cell(numel(I2_list),1);
vx_pix=cell(numel(I2_list),1);
vy_pix=cell(numel(I2_list),1);

v_rand=sqrt(vx_rand.^2 + vy_rand.^2 + vz_rand.^2 ); % Compute the magnitude of each sampled velocity vector.
for ind = 1:numel(I2_list)
Cell_velocity_el{ind}(:,1)=Cell_velocity{ind}(:,1)+rep_vel*vx_rand./v_rand; % Scale the x component of the unit direction by rep_vel to add the recoil velocity.
Cell_velocity_el{ind}(:,2)=Cell_velocity{ind}(:,2)+rep_vel*vy_rand./v_rand;
vx_pix{ind} = Cell_velocity_el{ind}(:,1)./vel_pix;
vy_pix{ind} = Cell_velocity_el{ind}(:,2)./vel_pix;
end

%image blur term
xgrid=-125:125;
ygrid=-125:125;

[xgird_2, ygrid_2]=meshgrid(xgrid,ygrid);
pFit = [-3.698853876412670	2.554394107105213e+04	0.441357801223283	6.354598584616713	6.584074190809296	8.250163296670557	9.417167979585598	0	0	-0.477546278665701];
% atten_factor = 1;

% pFit([4,6]) = pFit([4,6]) * atten_factorx;
% pFit([5,7]) = pFit([5,7]) * atten_factory;
M=f(xgird_2*atten_factorx, ygrid_2*atten_factory, pFit);
h2d=cell(numel(I2_list),1);
img_matrix=cell(numel(I2_list),1);
for ind = 1:numel(I2_list)
    h2d{ind}=hist3([vx_pix{ind}(:,1),vy_pix{ind}(:,1)],'Ctrs',{xgrid ygrid});
    img_matrix{ind}=conv2(h2d{ind},M,'same');
end

sumx=zeros(numel(xgrid),numel(I2_list));
sumy=zeros(numel(ygrid),numel(I2_list));

for ind=1:numel(I2_list)
sumx(:,ind)=sum(img_matrix{ind},2);
sumy(:,ind)=sum(img_matrix{ind},1);
end
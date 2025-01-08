%% Code for generating cylindrical and concentric cylnder patterns.
 
 % Creat square space (5000 x 5000) unit^2

 clear all
 figure


L = 2500; %units

%xs = [-L L L -L];
%ys = [L L -L -L];
%patch(xs,ys,'w')
axis equal
hold on
axis off


%% Coencentric cylinders n = no. of rings

n =1; %order
d1 = 200; %inner most diameter 
r1 = d1/2;
g_cc = 0; %intra cylinder gap
g =200; inter cylinder gap
t = 100; % Cylinder thickness
d1_e = d1 +2*t; r1_e = d1_e/2; % 
radius_vec = r1*ones(1,2*n);
colors_cell = cell(size(radius_vec));


for ind = 2:2:2*n
     radius_vec(ind) = radius_vec(ind)+t*(ind/2)+g_cc*((ind/2)-1);
end

for ind = 1:2:2*n
    radius_vec(ind) = radius_vec(ind)+t*((ind-1)/2)+g_cc*((ind-1)/2);
end

for ind = 1:length(colors_cell)
    if rem(ind,2)>0
        colors_cell(ind) = {'k'};
    else
        colors_cell(ind) = {'w'};
    end
end

d2 = d1+2*t+2*g_cc; r2 = d2/2;
d2_e = 2*radius_vec(2*n);  r2_e = radius_vec(2*n);





m = floor((L-2*r2_e)/(r2_e+g/2)+1)+4;
Cyl_count = m*m/2+(m-2)*m/2;

x_c_vec_pos = zeros(1,m);
x_c_vec_neg = zeros(1,m);
y_c_vec_pos = zeros(1,m);
y_c_vec_neg = zeros(1,m);


L_eff = (m*d2_e + (m-1)*g)/2;
offset = L-L_eff;

for p = 1:m-1
    x_c_vec_pos(1) = 0*offset+0*r2_e-0*L;
    x_c_vec_neg(1) = 0*offset+0*r2_e-0*L;
    y_c_vec_pos(1) = 0*offset+0*r2_e*sin(pi/3)+0*L*sin(pi/3)-0*L;
    y_c_vec_neg(1) = 0*offset+0*r2_e*sin(pi/3)+0*L*sin(pi/3)-0*L;


    x_c_vec_pos(p+1) = x_c_vec_pos(p)+d2_e+g;
    x_c_vec_neg(p+1) = x_c_vec_neg(p)-d2_e-g;
    y_c_vec_pos(p+1) = y_c_vec_pos(p)+(d2_e+g)*sin(pi/3);
    y_c_vec_neg(p+1) = y_c_vec_neg(p)-(d2_e+g)*sin(pi/3);
end



radius = flip(radius_vec);

%radius = [r1,r1+t1,r2,r2+t2];
u = 0:pi/50: 2*pi;

%colors = {'k','w','k','w'};
colors = colors_cell;
%colors = {'w','k','w','k'};


for y_cent = 1:2:m
    for x_cent = 1:m

     for i = 1:2*n

    x_pos = radius(i)*cos(u)+x_c_vec_pos(x_cent);
    x_neg = radius(i)*cos(u)+x_c_vec_neg(x_cent);
    y_pos = radius (i)*sin(u) +y_c_vec_pos(y_cent);
    y_neg = radius (i)*sin(u) +y_c_vec_neg(y_cent);

    patch(x_pos,y_pos,colors{i});
    patch(x_pos,y_neg,colors{i});
    patch(x_neg,y_pos,colors{i});
    patch(x_neg,y_neg,colors{i});
    hold on
    end
    end
end

%
for y_cent = 2:2:m
    for x_cent = 1:m

     for i = 1:2*n

    x_pos = radius(i)*cos(u)+x_c_vec_pos(x_cent)+r2_e+g/2;
    x_neg = radius(i)*cos(u)+x_c_vec_neg(x_cent)-r2_e-g/2;
    y_pos = radius (i)*sin(u) +y_c_vec_pos(y_cent);
    y_neg = radius (i)*sin(u) +y_c_vec_neg(y_cent);
    patch(x_pos,y_pos,colors{i});
    patch(x_pos,y_neg,colors{i});
    patch(x_neg,y_pos,colors{i});
    patch(x_neg,y_neg,colors{i});
    hold on
    end
    end
end

xlim([-L L]);
ylim([-L L]);
%}

%% Single cylinders plotting

x_c_or = 0;
y_c_or = 0;

radivec = zeros(1,n);
radivec(1) = r1;
for rad = 1:n-1
    radivec(rad+1) = radivec(rad)+t+g_cc;
end


for v = 1:n
    figure
        x1 = (radivec(v)+t)*cos(u)+ x_c_or ;
        y1 = (radivec(v)+t)*sin(u) + y_c_or;
        patch (x1,y1,'k');
        hold on;
        x2 = radivec(v)*cos(u)+ x_c_or ;
        y2 = radivec(v)*sin(u) + y_c_or;
        patch (x2,y2,'w');
        axis equal   %Figure 2-n
xlim([-L L]);
ylim([-L L]);
axis off
end

saveas(figure(1),'figure_con_4_1.jpg');
saveas(figure(2),'figure_con_4_2.jpg');
saveas(figure(3),'figure_con_4_3.jpg');
saveas(figure(4),'figure_con_4_4.jpg');
saveas(figure(5),'figure_con_4_5.jpg');

figure
imshow(imread("figure_con_4_1.jpg"))

%% Code for generating 2nd order fractal arrays.
 
% To Create square space (2L x 2L) unit^2
figure(1)

L = 2500;

xs = [-L L L -L];
ys = [L L -L -L];
patch(xs,ys,'w')
axis equal
hold on
axis off



%% Define input parameters
d1 = 25; % Inner most diameter
t = 5;   %Thikness of the cylinders 

%% Plotter code

r1 = d1/2; 

d2 = 3*(d1+2*t);

r1_e = r1+t;
d1_e = 2*r1_e;
radius1 = [r1+t, r1];
radius2 = [3*r1+4*t, 3*r1+3*t];
u = 0: pi/50:2*pi;
g = 0;
r2_e = 9*r1+13*t;
d2_e = 2*r2_e;


m = 2*ceil(L/d2_e);
m = m+rem(m,2)+4;
Cyl_count = m*m/2+(m-2)*m/2;
x_c_or_vec = zeros(1,m);
y_c_or_vec = zeros(1,m);
L_eff = (m*d2_e + (m-1)*g)/2;
offset = L-L_eff;

% Define the centers of the 3rd order nested arrays
x_c_or_vec(1) = 0;
x_c_or_vec(1) = 0;
for p = 1:m/2
    
    x_c_or_vec(p+1) = x_c_or_vec(p)+d2_e+g;
    y_c_or_vec(p+1) = y_c_or_vec(p)+(d2_e+g)*sin(pi/3);
end

x_c_or_vec(m/2+1:end) = -(x_c_or_vec(1:m/2));
y_c_or_vec(m/2+1:end) = -(y_c_or_vec(1:m/2));




colors = {'k','w','k','w'};

for y_cent = 1:2:m
 for x_cent = 1:m
x_c_or = x_c_or_vec(x_cent);
y_c_or = y_c_or_vec(y_cent);




x_c_vec3 = [0 0 0 0 0 0 0];
y_c_vec3 = [0 0 0 0 0 0 0];

for k = 2:7
x_c_vec3(1) =  x_c_or;
y_c_vec3(1) =  y_c_or;
x_c_vec3(k) = 2*(3*r1+4*t)*cos((k-1)*pi/3) + x_c_or;
y_c_vec3(k) = 2*(3*r1+4*t)*sin((k-1)*pi/3) + y_c_or;
end

  
x_c_vec = [0 0 0 0 0 0 0];
y_c_vec = [0 0 0 0 0 0 0];

for p = 2:7
    x_c_vec(p) = 2*(3*r1+4*t)*cos((p-1)*pi/3);
    y_c_vec(p) = 2*(3*r1+4*t)*sin((p-1)*pi/3);
end

x_nescent_vec1 = [0 0 0 0 0 0];
y_nescent_vec1 = [0 0 0 0 0 0];

% Plot 2nd order outer cylinder

radius3 =[9*r1+13*t, 9*r1+12*t];
for j2 = 1:2
        x = radius3(j2)*cos(u) + x_c_or;
        y = radius3(j2)*sin(u) + y_c_or;
        patch (x,y,colors{j2});
        hold on;
end
 for l = 1:7
    x_c = x_c_or + x_c_vec(l);
    y_c = y_c_or + y_c_vec(l);
    
    
for j1 = 1:2
        x = radius2(j1)*cos(u)+ x_c;
        y = radius2(j1)*sin(u) + y_c;
        patch (x,y,colors{j1});
        hold on;
end

  for ind1 = 1:7
       x_nescent_vec1(ind1+1) =d1_e*cos((ind1-1)*pi/3); 
       x_nescent_vec1(1) =0; 
       y_nescent_vec1(ind1+1) =d1_e*sin((ind1-1)*pi/3); 
       y_nescent_vec1(1) =0; 
   for i = 1:2
        x = radius1(i)*cos(u)+ x_c + x_nescent_vec1(ind1);
        y = radius1(i)*sin(u) + y_c + y_nescent_vec1(ind1);
        patch (x,y,colors{i});
        hold on;
    end
   end
  end
end
    
end

for y_cent = 2:2:m
 for x_cent = 1:m
x_c_or = x_c_or_vec(x_cent)+r2_e+g/2;
y_c_or = y_c_or_vec(y_cent);


x_c_vec3 = [0 0 0 0 0 0 0];
y_c_vec3 = [0 0 0 0 0 0 0];

for k = 2:7
x_c_vec3(1) =  x_c_or;
y_c_vec3(1) =  y_c_or;
x_c_vec3(k) = 2*(3*r1+4*t)*cos((k-1)*pi/3) + x_c_or;
y_c_vec3(k) = 2*(3*r1+4*t)*sin((k-1)*pi/3) + y_c_or;
end

   
x_c_vec = [0 0 0 0 0 0 0];
y_c_vec = [0 0 0 0 0 0 0];

for p = 2:7
    x_c_vec(p) = 2*(3*r1+4*t)*cos((p-1)*pi/3);
    y_c_vec(p) = 2*(3*r1+4*t)*sin((p-1)*pi/3);
end

x_nescent_vec1 = [0 0 0 0 0 0];
y_nescent_vec1 = [0 0 0 0 0 0];

% Plot 2nd order outer cylinder
radius3 =[9*r1+13*t, 9*r1+12*t];

for j2 = 1:2
        x = radius3(j2)*cos(u) + x_c_or;
        y = radius3(j2)*sin(u) + y_c_or;
        patch (x,y,colors{j2});
        hold on;
end
   for l = 1:7
    x_c = x_c_or + x_c_vec(l);
    y_c = y_c_or + y_c_vec(l);
    
    
    for j1 = 1:2
        x = radius2(j1)*cos(u)+ x_c;
        y = radius2(j1)*sin(u) + y_c;
        patch (x,y,colors{j1});
        hold on;
    end

    for ind1 = 1:7
       x_nescent_vec1(ind1+1) =d1_e*cos((ind1-1)*pi/3); 
       x_nescent_vec1(1) =0; 
       y_nescent_vec1(ind1+1) =d1_e*sin((ind1-1)*pi/3); 
       y_nescent_vec1(1) =0; 
        for i = 1:2
        x = radius1(i)*cos(u)+ x_c + x_nescent_vec1(ind1) ;
        y = radius1(i)*sin(u) + y_c + y_nescent_vec1(ind1);
        patch (x,y,colors{i});
        hold on;
        end
    end
  end
 end
end


axis equal
xlim([-L L]);
ylim([-L L]);

%% Single cylinders plotting

x_c_or = 0;
y_c_or = 0;

figure(2)

radius3 =[9*r1+13*t, 9*r1+12*t];
for j2 = 1:2
        x = radius3(j2)*cos(u) + x_c_or;
        y = radius3(j2)*sin(u) + y_c_or;
        patch (x,y,colors{j2});
        hold on;
end
axis equal   %Figure 2
xlim([-L L]);
ylim([-L L]);
axis off

figure(3)

for j1 = 1:2
        x = radius2(j1)*cos(u)+ x_c_or;
        y = radius2(j1)*sin(u) + y_c_or;
        patch (x,y,colors{j1});
        hold on;
end
axis equal   %Figure 3
xlim([-L L]);
ylim([-L L]);
axis off

figure(4)
for i = 1:2
        x = radius1(i)*cos(u)+ x_c_or ;
        y = radius1(i)*sin(u) + y_c_or;
        patch (x,y,colors{i});
        hold on;
end

axis equal   %Figure 4
xlim([-L L]);
ylim([-L L]);
axis off


saveas(figure(1),'figure_nes2_1.jpg');
saveas(figure(2),'figure_nes2_2.jpg');
saveas(figure(3),'figure_nes2_3.jpg');
saveas(figure(4),'figure_nes2_4.jpg');


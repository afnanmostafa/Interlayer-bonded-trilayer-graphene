
clear
fid = fopen('sampletri_diamane.data');
s = textscan(fid,'%f %f %f %f %f','headerlines',14);
fclose(fid);

index=s{1};
atom_type=s{2};
x = s{3};
y = s{4};
z = s{5};

write = true;
%% %%% Section 3: selecting the aligned atoms %%% %%

x_bottom = x(1:length(x)/2);
x_top = x(length(x)/2+1:end);

y_bottom = y(1:length(y)/2);
y_top = y(length(y)/2+1:end);

z_bottom = z(1:length(z)/2);
z_top = z(length(z)/2+1:end);

len = round((max(x)+min(x))/2);
leny = round((max(y)+min(y))/2);

id = index((x>=len & x<(len+50.5)) & (y>=leny & y<(leny+48)));

newX = x(id);
newY = y(id);
newZ = z(id);

plot(newX,newY,'.');

cc  = 1.4;                          %% C-C sp2 bond distance
per_arm = 4.26;                     %% C-C distance in armchair direction
per_zig = 2.46;                     %% C-C distance in zigzag direction
ly = 3*cc;
lx = sqrt(3)*cc;

len = max(newX) - min(newX);
wid = max(newY) - min(newY);

nx = round((len)/per_zig);    %% repetitions in the x direction
ny = ceil((wid)/per_arm);  

% nx = 19;    %% repetitions in the x direction
% ny = 7;    %% repetitions in the y direction
             
cc_height = 3.35;                

total_atoms = length(newX);
if write
    fid = fopen('sam.data','w');
    fprintf(fid,'#graphene\n');
    fprintf(fid,'%g atoms\n\n',total_atoms);
    fprintf(fid,'1 atom types\n\n');
    fprintf(fid,'%g %g xlo xhi\n',min(newX), max(newX)+ 1.25);
    fprintf(fid,'%g %g ylo yhi\n',min(newY),max(newY) + 1.45);
    fprintf(fid,'%g %g zlo zhi\n\n',-2*floor(cc_height*6*3),2*floor(cc_height*6*3));
    fprintf(fid,'Masses\n\n');
    fprintf(fid,'1 12.0107\n');
    
    fprintf(fid, '\n');
    fprintf(fid,'Atoms\n\n');
    
    newID = (1:length(newX));
    for j=1:length(newX)
        fprintf(fid,'%g %g %g %g %g\n',newID(j),1,newX(j),newY(j),newZ(j));
    end
    fclose(fid);
end
hold on
plot((min(newX):max(newX)),min(newY),'r--')


fid = fopen('sam.data');
s = textscan(fid,'%f %f %f %f %f','headerlines',15);
fclose(fid);

index=s{1};
atom_type=s{2};
x = s{3};
y = s{4};
z = s{5};


botC = index((z>-10.56 & z<-10.5));
topC = index((z<-5.6 & z>-6));

%% data file including the H atoms (no C atoms in this data file --> for record)
atomtypes = 1;
fid = fopen('h.data','w');

l=1;
cout = length(botC);
for m=1:cout
    fprintf(fid,'%g %d %g %g %g\n',length(newX)+l,2,x(botC(m)), y(botC(m)),z(botC(m))-0.9);
    l=l+1;
end

clear m l
cout2 = length(topC);
l = 1;
for m=1:cout2
    fprintf(fid,'%g %d %g %g %g\n',length(atom_type)+cout2+l,2,x(topC(m)), y(topC(m)),z(topC(m))+0.9);
    l=l+1;
end
fclose(fid);


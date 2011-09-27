function projectors(filename)

%to change conventions, edit from here %%%%%%%%%
g0_ch = [ 0  0 -1  0;
          0  0  0 -1;
         -1  0  0  0;
          0 -1  0  0];


g1_ch = [ 0  0  0 -i;
          0  0 -i  0;
          0  i  0  0;
          i  0  0  0];
      
g2_ch = [ 0  0  0 -1;
          0  0  1  0;
          0  1  0  0;
         -1  0  0  0];
       
g3_ch = [ 0  0 -i  0;
          0  0  0  i;
          i  0  0  0;
          0 -i  0  0];
       
g4_ch = -g0_ch; % redundand, but practical
       
%u = [1 0 1 0; 0 1 0 1; -1 0 1 0; 0 -1 0 1];
%n=2;

u= eye(4);
n=1;

g0 = (u'*g0_ch*u)/n
g1 = (u'*g1_ch*u)/n       
g2 = (u'*g2_ch*u)/n
g3 = (u'*g3_ch*u)/n
g4 = (u'*g4_ch*u)/n
%to here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

one = eye(4);       
       
g5 = g1*g2*g3*g4;


Proj1  = 0.5*(one+i*g5) * 0.25*((one+g0)*(one-i*g5*g3)).' *(one+i*g5);
Proj2  = 0.5*(one+i*g5) * 0.25*((one+g0)*(one-i*g5*(g1+g2+g3))).' *(one+i*g5);
Proj3  = 0.5*(one+i*g5) * 0.25*((one+g0)*(i*g5*(g1+g2+g3))).' *(one+i*g5);
Proj4  = 0.5*(one+i*g5) * 0.25*((one+g0)*(i*g5*g3)).' *(one+i*g5);

Proj5  = one;
Proj6  = one;
Proj7  = one;
Proj8  = one;
Proj9  = one;
Proj10 = one;
Proj11 = one;
Proj12 = one;

Proj13 = 0.5*(one+i*g5) * 0.5*(one+g0).' *(one+i*g5);
Proj14 = one;
Proj15 = 0.5*(one+i*g5) * 0.25*((one+g0)*(i*g5*g2)).' *(one+i*g5);
Proj16 = 0.5*(one+i*g5) * 0.25*((one+g0)*(i*g5*g1)).' *(one+i*g5);

Proj17 = one;
Proj18 = one;
Proj19 = one;

Proj20 = 0.5*(one+i*g5) * 0.5*(one+(g1+g2+g3)).' *(one+i*g5);
Proj21 = 0.5*(one+i*g5) * 0.25*(one).' *(one+i*g5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% print gamma conventions in gamma_xx.h format into file 'filename'

if (fid = fopen(filename,'w')) == 0,
   printf('Error opening file "%s" for writing\n',filename);
   return;
end;   

fprintf(fid,"/* projectors.h\n * collection of projector matrices\n * ETMC conventions\n *\n * Tomasz Korzec 2009\n **************************************/\n\n");

fprintf("\n// transposed projectors, rotated to TM basis\n\n");

fprintf(fid,'qcd_complex_16 PROJECTOR[22][4][4]=\n  {\n');
printgamma(fid,one);fprintf(fid,',\n\n');
printgamma(fid,Proj1);fprintf(fid,',\n\n');
printgamma(fid,Proj2);fprintf(fid,',\n\n');
printgamma(fid,Proj3);fprintf(fid,',\n\n');
printgamma(fid,Proj4);fprintf(fid,',\n\n');
printgamma(fid,Proj5);fprintf(fid,',\n\n');
printgamma(fid,Proj6);fprintf(fid,',\n\n');
printgamma(fid,Proj7);fprintf(fid,',\n\n');
printgamma(fid,Proj8);fprintf(fid,',\n\n');
printgamma(fid,Proj9);fprintf(fid,',\n\n');
printgamma(fid,Proj10);fprintf(fid,',\n\n');
printgamma(fid,Proj11);fprintf(fid,',\n\n');
printgamma(fid,Proj12);fprintf(fid,',\n\n');
printgamma(fid,Proj13);fprintf(fid,',\n\n');
printgamma(fid,Proj14);fprintf(fid,',\n\n');
printgamma(fid,Proj15);fprintf(fid,',\n\n');
printgamma(fid,Proj16);fprintf(fid,',\n\n');
printgamma(fid,Proj17);fprintf(fid,',\n\n');
printgamma(fid,Proj18);fprintf(fid,',\n\n');
printgamma(fid,Proj19);fprintf(fid,',\n\n');
printgamma(fid,Proj20);fprintf(fid,',\n\n');
printgamma(fid,Proj21);fprintf(fid,'\n');
fprintf(fid,'  };\n\n');

fclose(fid);
return;

function printgamma(fid,g)
fprintf(fid,'    {\n');
fprintf(fid,'      { {.re=%f,.im=%f},{.re=%f,.im=%f},{.re=%f,.im=%f},{.re=%f,.im=%f} },\n',real(g(1,1)),imag(g(1,1)),real(g(1,2)),imag(g(1,2)),real(g(1,3)),imag(g(1,3)),real(g(1,4)),imag(g(1,4)));
fprintf(fid,'      { {.re=%f,.im=%f},{.re=%f,.im=%f},{.re=%f,.im=%f},{.re=%f,.im=%f} },\n',real(g(2,1)),imag(g(2,1)),real(g(2,2)),imag(g(2,2)),real(g(2,3)),imag(g(2,3)),real(g(2,4)),imag(g(2,4)));
fprintf(fid,'      { {.re=%f,.im=%f},{.re=%f,.im=%f},{.re=%f,.im=%f},{.re=%f,.im=%f} },\n',real(g(3,1)),imag(g(3,1)),real(g(3,2)),imag(g(3,2)),real(g(3,3)),imag(g(3,3)),real(g(3,4)),imag(g(3,4)));
fprintf(fid,'      { {.re=%f,.im=%f},{.re=%f,.im=%f},{.re=%f,.im=%f},{.re=%f,.im=%f} }\n',real(g(4,1)),imag(g(4,1)),real(g(4,2)),imag(g(4,2)),real(g(4,3)),imag(g(4,3)),real(g(4,4)),imag(g(4,4)));
fprintf(fid,'    }');
return;

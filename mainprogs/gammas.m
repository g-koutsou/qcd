function gammas(filename)

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
g6 = (g1+i*g2)/2;
g7 = (g1-i*g2)/2;

C = g4*g2;

g5g0 = g5*g0;
g5g1 = g5*g1;
g5g2 = g5*g2;
g5g3 = g5*g3;
g5g4 = g5*g4;
g5g6 = g5*g6;
g5g7 = g5*g7;

Cg0 = C*g0;
Cg1 = C*g1;
Cg2 = C*g2;
Cg3 = C*g3;
Cg4 = C*g4;
Cg5 = C*g5;   

barCg0 = g4*Cg0'*g4;
barCg1 = g4*Cg1'*g4;
barCg2 = g4*Cg2'*g4;
barCg3 = g4*Cg3'*g4;
barCg4 = g4*Cg4'*g4;
barCg5 = g4*Cg5'*g4;    

barg5g0= g4*g0'*g5'*g4;
barg5g1= g4*g1'*g5'*g4;
barg5g2= g4*g2'*g5'*g4;
barg5g3= g4*g3'*g5'*g4;
barg5g4= g4*g4'*g5'*g4;

barg5g5= g4*g5'*g5'*g4;
barg5g6= g4*g6'*g5'*g4;
barg5g7= g4*g7'*g5'*g4;

epg0 = one+g0;
epg1 = one+g1;
epg2 = one+g2;
epg3 = one+g3;
epg4 = one+g4;
epg5 = one+g5;

emg0 = one-g0;
emg1 = one-g1;
emg2 = one-g2;
emg3 = one-g3;
emg4 = one-g4;
emg5 = one-g5;
       
% uncomment the following to check hermiticity and commutation relations
%printf('the following are zero if gammas are hermitian\n');
%g4-g4'
%g1-g1'
%g2-g2'
%g3-g3'
%g5-g5'
%printf('diagonal anticommutation relations - the following should be unit matrices\n');
%g4*g4
%g1*g1
%g2*g2
%g3*g3
%g5*g5
%printf('offdiagonal anticommutation relations - the following should be zero\n');
%g4*g1+g1*g4
%g4*g2+g2*g4
%g4*g3+g3*g4
%g4*g5+g5*g4
%g1*g2+g2*g1
%g1*g3+g3*g1
%g1*g5+g5*g1
%g2*g3+g3*g2
%g2*g5+g5*g2
%g3*g5+g5*g3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% print gamma conventions in gamma_xx.h format into file 'filename'

if (fid = fopen(filename,'w')) == 0,
   printf('Error opening file "%s" for writing\n',filename);
   return;
end;   

fprintf(fid,"/* qcd_gamma.h\n * collection of gamma matrices\n * ETMC conventions\n *\n * Tomasz Korzec 2009\n **************************************/\n\n");


fprintf(fid,'\n int qcd_EPS[6][3]=\n  {\n');
fprintf(fid,'    {0,1,2},\n');
fprintf(fid,'    {2,0,1},\n');
fprintf(fid,'    {1,2,0},\n');
fprintf(fid,'    {2,1,0},\n');
fprintf(fid,'    {0,2,1},\n');
fprintf(fid,'    {1,0,2}\n  };\n\n');

fprintf(fid,'int qcd_SGN_EPS[6]=\n  {\n    +1,+1,+1,-1,-1,-1\n  };\n\n');

fprintf(fid,'qcd_complex_16 qcd_ONE[4][4]=\n ');
printgamma(fid,one);fprintf(fid,';\n\n');


fprintf(fid,'/*\n *\n *    qcd_GAMMA[0,1,2,3,4,5,6,7] -> \\gamma_{0,1,2,3,4,5,+,-}\n *\n */\n');
fprintf(fid,'/*\n *\n *    gamma_4=-gamma_0\n *\n */\n');

fprintf(fid,'qcd_complex_16 qcd_GAMMA[8][4][4]=\n  {\n');
printgamma(fid,g0);fprintf(fid,',\n\n');
printgamma(fid,g1);fprintf(fid,',\n\n');
printgamma(fid,g2);fprintf(fid,',\n\n');
printgamma(fid,g3);fprintf(fid,',\n\n');
printgamma(fid,g4);fprintf(fid,',\n\n');
printgamma(fid,g5);fprintf(fid,',\n\n');
printgamma(fid,g6);fprintf(fid,',\n\n');
printgamma(fid,g7);fprintf(fid,'\n');
fprintf(fid,'  };\n\n');

fprintf(fid,'qcd_complex_16 qcd_CGAMMA[6][4][4]=\n  {\n');
printgamma(fid,Cg0);fprintf(fid,',\n\n');
printgamma(fid,Cg1);fprintf(fid,',\n\n');
printgamma(fid,Cg2);fprintf(fid,',\n\n');
printgamma(fid,Cg3);fprintf(fid,',\n\n');
printgamma(fid,Cg4);fprintf(fid,',\n\n');
printgamma(fid,Cg5);fprintf(fid,'\n');
fprintf(fid,'  };\n\n');

fprintf(fid,'qcd_complex_16 qcd_BAR_CGAMMA[6][4][4]=\n  {\n');
printgamma(fid,barCg0);fprintf(fid,',\n\n');
printgamma(fid,barCg1);fprintf(fid,',\n\n');
printgamma(fid,barCg2);fprintf(fid,',\n\n');
printgamma(fid,barCg3);fprintf(fid,',\n\n');
printgamma(fid,barCg4);fprintf(fid,',\n\n');
printgamma(fid,barCg5);fprintf(fid,'\n');
fprintf(fid,'  };\n\n');

fprintf(fid,'qcd_complex_16 qcd_G5GAMMA[8][4][4]=\n  {\n');
printgamma(fid,g5g0);fprintf(fid,',\n\n');
printgamma(fid,g5g1);fprintf(fid,',\n\n');
printgamma(fid,g5g2);fprintf(fid,',\n\n');
printgamma(fid,g5g3);fprintf(fid,',\n\n');
printgamma(fid,g5g4);fprintf(fid,',\n\n');
printgamma(fid,one);fprintf(fid,',\n\n');
printgamma(fid,g5g6);fprintf(fid,',\n\n');
printgamma(fid,g5g7);fprintf(fid,'\n');
fprintf(fid,'  };\n\n');

fprintf(fid,'qcd_complex_16 qcd_BAR_G5GAMMA[8][4][4]=\n  {\n');
printgamma(fid,barg5g0);fprintf(fid,',\n\n');
printgamma(fid,barg5g1);fprintf(fid,',\n\n');
printgamma(fid,barg5g2);fprintf(fid,',\n\n');
printgamma(fid,barg5g3);fprintf(fid,',\n\n');
printgamma(fid,barg5g4);fprintf(fid,',\n\n');
printgamma(fid,barg5g5);fprintf(fid,',\n\n');
printgamma(fid,barg5g6);fprintf(fid,',\n\n');
printgamma(fid,barg5g7);fprintf(fid,'\n');
fprintf(fid,'  };\n\n');


fprintf(fid,'qcd_complex_16 qcd_ONE_PLUS_GAMMA[6][4][4]=\n  {\n');
printgamma(fid,epg0);fprintf(fid,',\n\n');
printgamma(fid,epg1);fprintf(fid,',\n\n');
printgamma(fid,epg2);fprintf(fid,',\n\n');
printgamma(fid,epg3);fprintf(fid,',\n\n');
printgamma(fid,epg4);fprintf(fid,',\n\n');
printgamma(fid,epg5);fprintf(fid,'\n');
fprintf(fid,'  };\n\n');

fprintf(fid,'qcd_complex_16 qcd_ONE_MINUS_GAMMA[6][4][4]=\n  {\n');
printgamma(fid,emg0);fprintf(fid,',\n\n');
printgamma(fid,emg1);fprintf(fid,',\n\n');
printgamma(fid,emg2);fprintf(fid,',\n\n');
printgamma(fid,emg3);fprintf(fid,',\n\n');
printgamma(fid,emg4);fprintf(fid,',\n\n');
printgamma(fid,emg5);fprintf(fid,'\n');
fprintf(fid,'  };\n\n');


fclose(fid);
return;

function printgamma(fid,g)
fprintf(fid,'    {\n');
fprintf(fid,'      { {.re=%0.1f,.im=%0.1f},{.re=%0.1f,.im=%0.1f},{.re=%0.1f,.im=%0.1f},{.re=%0.1f,.im=%0.1f} },\n',real(g(1,1)),imag(g(1,1)),real(g(1,2)),imag(g(1,2)),real(g(1,3)),imag(g(1,3)),real(g(1,4)),imag(g(1,4)));
fprintf(fid,'      { {.re=%0.1f,.im=%0.1f},{.re=%0.1f,.im=%0.1f},{.re=%0.1f,.im=%0.1f},{.re=%0.1f,.im=%0.1f} },\n',real(g(2,1)),imag(g(2,1)),real(g(2,2)),imag(g(2,2)),real(g(2,3)),imag(g(2,3)),real(g(2,4)),imag(g(2,4)));
fprintf(fid,'      { {.re=%0.1f,.im=%0.1f},{.re=%0.1f,.im=%0.1f},{.re=%0.1f,.im=%0.1f},{.re=%0.1f,.im=%0.1f} },\n',real(g(3,1)),imag(g(3,1)),real(g(3,2)),imag(g(3,2)),real(g(3,3)),imag(g(3,3)),real(g(3,4)),imag(g(3,4)));
fprintf(fid,'      { {.re=%0.1f,.im=%0.1f},{.re=%0.1f,.im=%0.1f},{.re=%0.1f,.im=%0.1f},{.re=%0.1f,.im=%0.1f} }\n',real(g(4,1)),imag(g(4,1)),real(g(4,2)),imag(g(4,2)),real(g(4,3)),imag(g(4,3)),real(g(4,4)),imag(g(4,4)));
fprintf(fid,'    }');
return;

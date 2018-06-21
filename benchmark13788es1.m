function benchmark13788
  
   c=0; #condensation/evaporation
 r=1;#relative humidity
 p=1;#vpres satpres
 p2=0; #vpres satpres over vapour resistance
 t=1;
  plot_options=[c,r,p,p2,t];

  
  %%example C1
Wall=[ 0.01   0.01/0.05     500000;
       0.1    0.1/3         150;
       0.001  1             1000/0.001
       0.012  0.012/0.075   10];


%      %example C2
%      Wall=[  0.01   0.01/0.05     5000/0.01;
%              0.1    .1/3         150;
%              0.001  1            1000/0.001;
%              0.1    .1/3         150;
%              0.012  0.012/0.075   10];
      
      
%DATI CLIMATICI LUOGHI BASSA OCCUPAZIONE
BC=[  -1  .85 20  .39;
      0   .84 20  .40;
      4   .78 20  .44;
      9   .72 20  .49;
      14  .68 22  .54;
      18  .69 24  .58;
      19  .73 24.5  .59;
      19  .75 24.5  .55;
      15  .79 22.5  .55;
      10  .83 20  .50;
      5   .88 20  .45;
      1   .88 20  .41];  


%DATI CLIMATICI ALTA OCCUPAZIONE
%BC=[  -1  .85 20  .49;
%      0   .84 20  .50;
%      4   .78 20  .54;
%      9   .72 20  .59;
%      14  .68 22  .64;
%      18  .69 24  .68;
%      19  .73 24.5  .69;
%      19  .75 24.5  .69;
%      15  .79 22.5  .65;
%      10  .83 20  .60;
%      5   .88 20  .55;
%      1   .88 20  .51];

dl=0.02;
Rs=[0.04, 0.10];



glaser_method(Wall,BC,Rs,dl,plot_options)
  
endfunction
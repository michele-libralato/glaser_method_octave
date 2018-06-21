function benchmark13788
  
 c=0; #condensation/evaporation
 r=0;#relative humidity
 p=1;#vpres satpres
 p2=1; #vpres satpres over vapour resistance
 t=0;
  plot_options=[c,r,p,p2,t];

      %example C4
      Wall=[  0.105   .77206  8;
              0.05    .278    0.05/0.01;
              0.001   1       200;
              0.012   .1304   90;
              0.14    .04     1.4;
              0.001   1       50/.001;
              0.0125  .208     12];
      
      

%DATI CLIMATICI esempio 4
BC=[  8.0  .630 20.0 .4
8.0  .634 20.0 .4
11.0 .644 20.0 .4
15.5 .692 20.0 .4
19.5 .769 20.0 .4
21.5 .841 20.0 .4
26.0 .854 20.0 .4
27.0 .828 20.0 .4
25.0 .802 20.0 .4
20.0 .705 20.0 .4
15.5 .688 20.0 .4
10.5 .664 20.0 .4];

dl=0.14;
Rs=[0.04, 0.13];



glaser_method(Wall,BC,Rs,dl,plot_options)
  
endfunction
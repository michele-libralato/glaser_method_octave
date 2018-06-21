function benchmark13788
  
 c=0; #condensation/evaporation
 r=0;#relative humidity
 p=1;#vpres satpres
 p2=1; #vpres satpres over vapour resistance
 t=0;
  plot_options=[c,r,p,p2,t];

      %example C4
      Wall=[  0.012 0.012/0.015 8  ;
              0.275 0.275/2.500 10 ;
              0.020 0.020/0.500 1  ;
              0.013 0.013/0.023 8];
      
      

%DATI CLIMATICI esempio 4
BC=[  -10 0.95 20 0.49    ;
 -8 0.94 20 0.50    ;
 -5 0.80 20 0.54    ;
 0 0.82 20 0.59     ;
 13 0.68 22 0.64    ;
 18 0.69 24 0.68    ;
 19 0.73 24.5 0.69  ;
 19 0.75 24.5 0.69  ;
 14 0.79 22.5 0.65  ;
 7 0.83 20 0.60     ;
 1 0.88 20 0.55     ;
 -4 0.95 20 0.51];

dl=0.02;
Rs=[0.04, 0.13];



glaser_method(Wall,BC,Rs,dl,plot_options)
  
endfunction
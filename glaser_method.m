function glaser_method(Wall,BC,Rs,dl,plot_options)
%function glaser_method
## input description
## glaser_method(wall,BC,Rs,dl)
## wall=[layer1 thickness, conductivity, permeability;
##       layer2 thickness, conductivity, permeability;
##       ...]; 
## BC=[january Te,RHe,Ti,RHi;
##    february Te,RHe,Ti,RHi;
##      ... ];
## dl=calculation step default dl=0.001
## Rs=[Rse,Rsi]; default Rs=[0.04, 0.13]


%%%example
%Wall=[  1.2000e-001  1.3000e-001  2.0000e+002;
%       1.0000e-003  2.3000e+000  1.5000e+006;
%        4.0000e-002  3.5700e-002  2.0000e+000;
%      2.2000e-001  3.5700e-002  2.0000e+000;
%      4.0000e-002  3.5700e-002  2.0000e+000;
%       3.0000e-002  1.3000e-001  2.0000e+002;
%       5.0000e-002  2.8000e-001  3.2000e-001;
%       3.0000e-002  1.3000e-001  2.0000e+002];
%      
%%BC=[20,.50,0,.20;
%%    10,.80,-10,.30;
%%    10,.80,-10,.30;
%%    10,.80,-10,.30;
%%    10,.80,-10,.30;
%%    10,.80,-10,.30;
%%    10,.80,-10,.30;
%%    10,.80,-10,.30;
%%    10,.80,-10,.30;
%%    10,.80,-10,.30;
%%    10,.80,-10,.30;
%%    5, 1, -10, .5];
%BC=BoundaryConditions('weather_files/AUT_Vienna.Schwechat.110360_IWEC.epw');
%    dl=0.001;
%Rs=[0.04, 0.13];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIVIDE LAYERS IN dl ELEMENTS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
more off;

#PLOT OPTIONS
if !plot_options
 c=0; #condensation/evaporation
 r=0;#relative humidity
 p=0;#vpres satpres
 p2=0; #vpres satpres over vapour resistance
 t=0; #temperatures
else
  c=plot_options(1);
  r=plot_options(2);
  p=plot_options(3);
  p2=plot_options(4);
  t=plot_options(5);
endif
%% elements=[,progressive distance from outside to inside, 
%%            progressive thermal resistance, progressive vapour permeability,
%%          layer number, ];
Rse=Rs(1);
Rsi=Rs(2);

layersNumber=size(Wall,1);
WallThickness=sum(Wall(:,1));
ElementsNumber=sum(ceil(Wall(:,1)/dl))+1
nBC=size(BC,1);
elements(1,1)=0; %first element thickness
elements(1,2)=0; %first element thickness
elements(1,3)=0; %first element thermal resistance
elements(1,4)=0; %first element vapour resistance
elements(1,5)=0; %layer number

beginLayer=2; %first element of the consiered layer
endLayer=1;  %end element of the considered layer (defined later)

wstorage=zeros(ElementsNumber,nBC);
condensate=wstorage;

for layer=1:layersNumber
  LayerThickness=Wall(layer,1);

  endLayer+=ceil(LayerThickness/dl); %end element of the considered layer
  for e=beginLayer:endLayer
      if e==endLayer
        elements(e,1)=LayerThickness-dl*(endLayer-beginLayer);
      else
        elements(e,1)=dl;
      endif
      dlelement=elements(e,1); %element thinckness
      dRt=dlelement/Wall(layer,2);
      dRv=dlelement*Wall(layer,3);
      
      elements(e,2)=elements(e-1,2)+dlelement; %progressive thickness
      elements(e,3)=elements(e-1,3)+dRt; %progressive thermal resistance
      elements(e,4)=elements(e-1,4)+dRv; %progressive vapour resistance
      elements(e,5)=layer; %layer number
    endfor
  beginLayer=endLayer+1;

endfor
if ElementsNumber != size(elements,1)
  disp('WARNING: wrong elements count');
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEMPERATURES, SAT AND VAP PRESSURE %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rttot=Rse+elements(ElementsNumber,3)+Rsi;
Rvtot=elements(ElementsNumber,4);
difference=zeros(ElementsNumber,nBC);

for m=1:nBC
  Te=BC(m,1);
  Ti=BC(m,3);

  RHe=BC(m,2);
  RHi=BC(m,4);

  Pve=RHe*SatPressure(Te);  
  Pvi=RHi*SatPressure(Ti);

  vapPres(1,m)=Pve;
  %vapPres(ElementsNumber,m)=Pvi;
  
  Temperature(1,m)=Te+((Rse)*(Ti-Te))/Rttot;
  %Temperature(ElementsNumber,m)=Te+((Rttot-Rsi)*(Ti-Te))/Rttot; %temperatura superficiale
  
  satPres(1,m)=SatPressure(Temperature(1,m));
  %satPres(ElementsNumber,m)=SatPressure(Temperature(ElementsNumber,m));

  for el=2:ElementsNumber
    Rt=elements(el,3);
    Rv=elements(el,4);
    Temperature(el,m)=Te+((Rse+Rt)*(Ti-Te))/Rttot;
    satPres(el,m)=SatPressure(Temperature(el,m));
    vapPres(el,m)=Pve+(Rv*(Pvi-Pve))/Rvtot;
    if vapPres(el,m)-satPres(el,m)>0.1
      difference(el,m)=vapPres(el,m)-satPres(el,m);
    endif
  endfor
  
%figure
%    plot_title=[ 'vap and sat pression, month ' num2str(m)];
%    plot(elements(:,2),satPres(:,m),elements(:,2),vapPres(:,m))
%    title(plot_title)
%    axis ([0, WallThickness ]);  
 
 endfor
%figure  
%plot(elements(:,2),elements(:,3))
%figure  
%plot(elements(:,2),elements(:,4))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIRST MONTH of CONDENSATION %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
storageswitch=0;
m1=1;
for m=1:nBC
  %m actual month
  %M previous month
  if m>1
    M=m-1;
  else
    M=12;
  endif 
  
  if sum(difference(:,M))==0 && sum(difference(:,m))>0
    m1=m;
    break;
  endif

  endfor
%m1=12

   disp(['initial month ' num2str(m1)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACTUAL VAPPRES %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fixedPoints is a vector with the indexes of 
%of the points in which vapPres is imposed
%vapPres=satPres because vapPres>satPres

months_days=[31 29 31 30 31 30 31 31 30 31 30 31];
condensate=zeros(ElementsNumber,nBC);
wstorage=zeros(ElementsNumber,nBC);
delta0=2e-10;

fixedPoints(1)=1;
fixedPoints(2)=ElementsNumber;

for M=1:nBC
  m=m1+M-1;
  if m>12
    m=m-12;
  endif
    Mm1=m-1;
  if m==1
    Mm1=12;
  endif

  while sum(difference(:,m))>10 || storageswitch==1
    %disp(['inizio ciclo while per ' num2str(m)])
    
    [value,index]=max(difference(:,m));
    if all(index!=fixedPoints)
      fixedPonitsNumber=size(fixedPoints,2);
      fixedPoints(fixedPonitsNumber+1)=index;
     % 'aggiungo fixed points da vapPres'
      index
      m;
      fixedPoints=sort(fixedPoints);
    endif
    fixedPonitsNumber=size(fixedPoints,2)
    
    
    for j=2:fixedPonitsNumber-1
      j
      fixedPonitsNumber
      el=fixedPoints(j)
      elsat1=fixedPoints(j-1)
      elsat2=fixedPoints(j+1)
      vapPres(el,m)=satPres(el,m);

%     figure;plot(elements(:,4),vapPres(:,11),elements(:,4),satPres(:,11))
		  %cambio le pressioni in avanti
      for el1=el+1:elsat2-1
        vapPres(el1,m)=vapPres(el,m)+((elements(el1,4)-elements(el,4))*(vapPres(elsat2,m)-vapPres(el,m)))/(elements(elsat2,4)-elements(el,4));
      endfor
      %cambio le pressioni indietro
      for el0=elsat1+1:el-1
        vapPres(el0,m)=vapPres(elsat1,m)+(elements(el0,4)-elements(elsat1,4))*(vapPres(el,m)-vapPres(elsat1,m))/(elements(el,4)-elements(elsat1,4));
      endfor
    endfor
    %ricalcolo le differenze
    for elm=1:ElementsNumber
      difference(elm,m)=0;
      if vapPres(elm,m)-satPres(elm,m)>0
        difference(elm,m)=vapPres(elm,m)-satPres(elm,m);
      endif
    endfor
%    figure
%    plot(elements(:,2),satPres(:,m),"-",elements(:,2),satPres(:,m),"o",elements(:,2),vapPres(:,m),"-",elements(:,2),vapPres(:,m),"s")

      storageswitch=0;
  endwhile
%ricalcolo RH
  for elm=1:ElementsNumber
    RH(elm,m)=vapPres(elm,m)/satPres(elm,m);
  endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% condensate calculation %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%condensate(el,m) contains the water condensed (+) 
%or the water evaporated (-) for the month m
%wstorage(el,m)


  disp(['inizio calcolo condensa per ' num2str(m)])
 
  for el=1:ElementsNumber-1
%    disp(elements(el+1,4))
%    disp(elements(el,4))
    g0(el,m)=delta0*months_days(m)*24*60*60*(vapPres(el+1,m)-vapPres(el,m))/(elements(el+1,4)-elements(el,4));
  endfor
  
  Ti=BC(m,3);
  Pvi=RHi*SatPressure(Ti);
  g0(ElementsNumber,m)=g0(ElementsNumber-1,m);
%  g0(ElementsNumber,m)=delta0*months_days(m)*24*60*60*(vapPres(ElementsNumber,m)-Pvi)/(elements(ElementsNumber,4));;
  
   
    for el=2:ElementsNumber
    %condensate
      if RH(el,m)==1
        if abs(g0(el,m)-g0(el-1,m))>10e-7
          condensate(el,m)=g0(el,m)-g0(el-1,m);
        endif
      %water storage
        if condensate(el,m)+wstorage(el,Mm1)<0
          disp('moisture content evaporated')
          wstorage(el,m)=0;
        else
          wstorage(el,m)=condensate(el,m)+wstorage(el,Mm1);
          disp('total moisture content')
        endif
      endif
  endfor
%
%  figure
% plot_title=[ 'g0, month ' num2str(m)];
% plot(elements(:,2),g0(:,m))
% title(plot_title)
% axis([0, WallThickness])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set vapour pressure to sat pressure %%%
% in the next month %%%%%%%%%%%%%%%%%%%%%
%if condensate is not dried %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %disp(['inizio calcolo fixedpoints da condensa per ' num2str(m)])
%  c=2;
  clear fixedPoints
  fixedPoints(1)=1;
  fixedPoints(2)=ElementsNumber;
  for el=1:ElementsNumber
    if wstorage(el,m)>0 && all(el!=fixedPoints)
      %'aggiungo un fixed point da condensa'
      fixedPonitsNumber=size(fixedPoints,2)
      el
      m
      fixedPoints(fixedPonitsNumber+1)=el;
      fixedPoints=sort(fixedPoints)
      storageswitch=1;
    endif
  endfor

  
endfor



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT RESULTS  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:nBC
  

 if c==1
 figure
 plot_title=[ 'water stored, month' num2str(m)];
 plot(elements(:,2),wstorage(:,m))
 title(plot_title);
 axis([0,WallThickness]);
 endif
 if r==1
 figure
 plot_title=[ 'Relaive Humidity, month ' num2str(m)];
 plot(elements(:,2),RH(:,m),"o",elements(:,2),RH(:,m),"-")
 title(plot_title)
 axis([0, WallThickness,0,1 ])
 endif
 if p==1
 figure
 plot_title=[ 'Vapour Pressure, month ' num2str(m)];
 plot(elements(:,2),satPres(:,m),"-",elements(:,2),satPres(:,m),"o",elements(:,2),vapPres(:,m),"-",elements(:,2),vapPres(:,m),"s")
 title(plot_title)
 axis ([0, WallThickness ]);
 endif
 if p2==1
 figure
 plot_title=[ 'Vapour Pressure, month ' num2str(m)];
 plot(elements(:,4),satPres(:,m),"-",elements(:,4),satPres(:,m),"o",elements(:,4),vapPres(:,m),"-",elements(:,4),vapPres(:,m),"s")
 title(plot_title)
% axis ([0, WallThickness ]);
 endif
 if t==1
 figure
 plot_title=['Temperature, month ' num2str(m)];
 plot(elements(:,2),Temperature(:,m),'o',elements(:,2),Temperature(:,m),'-')
 title(plot_title)
 axis ([0, WallThickness ]);
 endif

endfor  


for mnt=1:12
results(mnt,:)=[sum(condensate(:,mnt)), sum(abs(wstorage(:,mnt)))];
endfor
'condensate'
results

endfunction






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAT PRESSURE FROM TEMPERATURE %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psat=SatPressure(ta)
 U=2; %U=1 metodo 13788 U=2 metodo ashrae
switch (U)
case 1
 if ta>0
    psat=610.5*exp(17.269*ta/(237.3+ta));
  else
    psat=610.5*exp(21.875*ta/(265.5+ta));
  endif
case 2
T = ta + 273.15;
if ta < 0 
    C1 = -5674.5359;
    c2 = 6.3925247;
    c3 = -0.009677843;
    c4 = 0.00000062215701;
    c5 = 2.0747825E-09;
    c6 = -9.484024E-13;
    c7 = 4.1635019;
else
    C1 = -5800.2206;
    c2 = 1.3914993;
    c3 = -0.048640239;
    c4 = 0.000041764768;
    c5 = -0.000000014452093;
    c6 = 0;
    c7 = 6.5459673;
endif

lnpws = C1 / T + c2 + c3 * T + c4 * T ^ 2 + c5 * T ^ 3 + c6 * T ^ 4 + c7 * log(T);
psat = exp(lnpws);
endswitch
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOUNDARY CONDITIONS           %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [BC]=BoundaryConditions(filename)
  
global wuficases;
global elements;
global resistance;
  datiMediati=mediaDati(filename);
  BC=zeros(12,4);
  %Pvap esterne
  for i=1:12
    te=datiMediati(i,1);
    RHe=datiMediati(i,3)*.01;
    ti=ContIntTemp(te);
    RHi=IntRelHumidity(te,ti,RHe);
    estPsat=SatPressure(te)*RHe;
    intPsat=SatPressure(ti)*RHi;
    
    BC(i,:)=[te RHe ti RHi];
%    BC(i,:)=[te RHe];
  endfor
  
endfunction




function RHi=IntRelHumidity(te,ti,RHe)

   c=1;
   class_p=1360;
   class_pm=200;
   switch c
  
  %Internal Vapour Pressure - Continental and Tropical Climates
  
  %Normal Occupancy
   case 1
   if te<-10
    RHi=.35;
  elseif te>-10&&te<20
    RHi=.35+(te+10)*.01;
  else
    RHi=.65;
  endif
  
  %High Occupancy
  
   case 2
   if te<-10
    RHi=.4;
  elseif te>-10&&te<20
    RHi=.4+(te+10)*.01;
  else
    RHi=.70;
  endif
  
  %maritime climates
  
  case 3
  if te<0
    delta_p=class_p;
  elseif te<20 && te>0
	delta_p=class_p-te*(class_p-class_pm)/20;
	else
	delta_p=class_pm;
	endif
  
  delta_p;
  pe=SatPressure(te)*RHe;
  pi=pe+delta_p;
  RHi=pi/SatPressure(ti);
  
  
  endswitch
endfunction




function [datiMediati]=mediaDati(filename)

global wuficases;
global elements;
global resistance;

  datiMediati=zeros(12,3);
  [date, hour, dry_bulb_temperature,dew_point_temperature,relative_humidity]=EPWopen(filename);
  Length=size(dry_bulb_temperature);
  datiMediatiTemp=zeros(12,4);
  for i=1:Length
    mese=date(i,2);
    datiMediatiTemp(mese,1) += dry_bulb_temperature(i);
    datiMediatiTemp(mese,2) += dew_point_temperature(i);
    datiMediatiTemp(mese,3) += relative_humidity(i);
    datiMediatiTemp(mese,4)++;
  endfor
  for i=1:12
    datiMediati(i,1)=datiMediatiTemp(i,1)/datiMediatiTemp(i,4);
    datiMediati(i,2)=datiMediatiTemp(i,2)/datiMediatiTemp(i,4);
    datiMediati(i,3)=datiMediatiTemp(i,3)/datiMediatiTemp(i,4);
  endfor
  %disp(datiMediati);
endfunction





function ti = ContIntTemp(te)
  %Internal Temperature- Continental and Tropical Climates
  if te<10
    ti=20;
  elseif te>10&&te<20
    ti=20+.5*(te-10);
  else
    ti=25;
  endif 
endfunction


function [date, hour, dry_bulb_temperature,dew_point_temperature,relative_humidity]=EPWopen(filename)
  
global wuficases;
global elements;
global resistance;
  
  disp('------------- EPW open -----------');
  disp(filename);
  fid=fopen(filename,'r'); % Open TMY file
  if fid == -1
    disp(['file not found: ',filename]);
  else
    tmy = csvread (filename);      %read the file 
  endif
  % file header information
  
  date = zeros(8760,3); 
  hour = zeros(8760,1); 
  dry_bulb_temperature = ones(8760,1); 
  dew_point_temperature = ones(8760,1); 
  relative_humidity = ones(8760,1);
  
  for i=9:8768 
    date(i-8,1:3) = tmy(i,1:3);
    hour(i-8) = tmy(i,4);
    dry_bulb_temperature(i-8) = (tmy(i,7));
    dew_point_temperature(i-8) = (tmy(i,8));
    relative_humidity(i-8) = (tmy(i,9));
   endfor  
 %disp([ date(1) hour(1) dry_bulb_temperature(1) dew_point_temperature(1) relative_humidity(1)]);
endfunction
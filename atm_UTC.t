// include atm_UTC.t;
// transformation des données en données UTC



macro asos4{

dir = "Donnees/";

mimin=0,24*60-1;

mP1=mimin*0;
mP2=mimin*0;
mP3=mimin*0;
nP=mP1;

// on cumule 
plotreset;
RAD="64060";
//STA="K0J4";
STA="PHOG";
//STA="PHLI";
AN="2002";

//for mo=1 to 12 {
//file =file_fullname(dir+RAD+STA+AN+str("%02d",mo)+".dat");
file1="PHOG2002.dat";
_read = ( "skipbrokenlines", "skipemptylines" );

_info off;
%readasos4[file1];

max(P1);
max(P2);
max(P3);
min(P1);
min(P2);
min(P3);
t=1,size(P1);
//write("PHOG2002bis.dat",t,33.8637526*P1,33.8637526*P2,33.8637526*P3);
//ce bout de code permet de detecter une erreur faite en nettoyant les data
/*
min_idx=1;
for j=1 to size(P3){
if (P3[j]<P3[min_idx])  then
{

min_idx=j;
};
};
*/
%acump3;
//};

mP1=mP1/nP;
mP2=mP2/nP;
mP3=mP3/nP;
ihg2mbar=33.8637526;

/*
// PNG
gnuplot;
#set terminal png;
set output "PatmUTC.png";
set xlabel "temps(heures UTC)";
set ylabel "Pression(mbar)";
set xrange[0:24];
end;
//plot(mimin/60, mP1*ihg2mbar, "w l title'Patm1'");
//replot(mimin/60, mP2*ihg2mbar, "w l title'Patm2'");
//replot(mimin/60, mP3*ihg2mbar, "w l title'Patm3'");
*/

mbP1=mP1*ihg2mbar;
mbP2=mP2*ihg2mbar;
mbP3=mP3*ihg2mbar;

//plot(mimin/60, mbP1-mbP3, "notitle w l");

// ajustement 
//mbP0=(mbP1+mbP2+mbP3)/3; 

mbP0=(mbP1+mbP2)/2;

mbP0=(mbP1+mbP2)/2;

P1off=sum(mbP1-mbP0)/size(mbP1);
P2off=sum(mbP2-mbP0)/size(mbP1);
P3off=sum(mbP3-mbP0)/size(mbP1);

xmbP1=mbP1-P1off;
xmbP2=mbP2-P2off;
xmbP3=mbP3-P3off;

//plot(mimin/60, xmbP1, "notitle w l");
//replot(mimin/60, xmbP2, "notitle w l");
//replot(mimin/60, xmbP3, "notitle w l");

//plot(mimin/60, xmbP2-xmbP1, "notitle w l");
//replot(mimin/60, xmbP3-xmbP1, "notitle w l");

// sauvegarde 
fout=STA+AN+"UTC.dat";
fmt="%6.0f %10.4f  %10.4f  %10.4f\n"; 
//write(fout,fmt,mimin,xmbP1,xmbP2,xmbP3);


// analyse de Fourier

mday=mimin/60/24;

// on ramene les temps des données a 1 
tM=mday;
//PM=mbP0;
vnumR DP;
read("DP1UTC.dat",DP);
PM=DP;
%four[tM,PM];

n1=size(xmbP1);
n2=size(xmbP2);
n3=size(xmbP3);

DP1 =xmbP1[2:n1]-xmbP1[1:n1-1];
DP2 =xmbP2[2:n2]-xmbP2[1:n2-1];
DP3 =xmbP3[2:n3]-xmbP3[1:n3-1];

};

macro readasos4[file]{
//PHLI

vnumR WBAN,Yr,Mo,Dy,hlt,mlt,hutc,mutc,P1,P2,P3,T1,T2;
vnumR pl,ps,sfr;
// ch14,ch24,ch3
form="%5g%4s%4s%4g%2g%2g%2g%2g%2g%2g %37s  %g  %g %g";

read(file,form, WBAN,ch14,ch24,Yr,Mo,Dy,hlt,mlt,hutc,mutc, sfr,P1,P2,P3);

nn=size(P1);
nm=1,nn;
//plot (nm,P1, "notitle w dots lw 2");
//replot (nm,P2, "notitle w  dots lw 2");
//replot (nm,P3, "notitle w  dots lw 2");

};

macro acump3{
nutc=hutc*60+mutc;


for k=1 to size(P1) {
    kk=k$
    nP[nutc[k]+1]=nP[nutc[k]+1]+1$
    mP1[nutc[k]+1]= mP1[nutc[k]+1]+ P1[k]$
    mP2[nutc[k]+1]= mP2[nutc[k]+1]+ P2[k]$
    mP3[nutc[k]+1]= mP3[nutc[k]+1]+ P3[k]$
   };
};

macro four[tM,PM]{

// analyse de Fourier
// on reechantillonne 
tf=0,2^12-1;
// on ramene les temps Fourier a  1
tF=tf/2^12;
// on reechantillonne 
PMF=interpol(SPLINE, tM,PM,tF);
// Elle retourne un vecteur de rÈels des valeurs interpolÈes aux temps tF ‡ partir de (tM,PM)
//plot (tM,PM, "notitle w l");
//replot (tF,PMF, "notitle w l");



// analyse de Fourier 

fmin=0;
fmax=20;
_naf_dtour=1;
%plotfftR[tF,PMF,_naf_dtour,fmin,fmax,0];
 gnuplot;
 #set terminal png;
 #set output "SPECTRE_FFT_PHOG2002.png";
 set xlabel "frequence(1/jour)";
 set ylabel "Amplitude(mbar)";
 set xrange[0:24];
 set xtics 1
 end;
plot(_ftf, 2*abs(_fza),"  w impulses lw 2 title'FFT'");

// reconstruction 
NT=10;
N1=int(size(_ftf)/2)-NT+1;
N2=int(size(_ftf)/2)+1+NT;

sP=sertrig(_fza[N1:N2],_ftf[N1:N2],T);

tt=0,1439;
ttt=tt/1440*2*pi;
vP=evalnum(sP,REAL,(T,ttt));
// vP est  un vecteur numÈrique de rÈels 

mPM=sum(PM)/size(PM);


//plot("set term X11 1",tt,vP-mPM,"notitle w l");
//replot ("set term X11 1",tt,PM-mPM, "notitle w l");

//plot("set term X11 2",tt,PM-vP,"notitle w l");

// recuperation des termes 
zap=_fza[N1:N2];
frp=_ftf[N1:N2];
//writes(zap[1:24]);
argp=arg(zap);
modp=abs(zap);

mop0=abs(zap[NT+1:NT+1]);
mop=modp[NT+2:size(zap)];
arp=argp[NT+2:size(zap)];
fr=frp[NT+2:size(zap)];

darp=arp/pi*180;

writes("  %12.6f\n",mop0);
writes("+ %12.6f *cos(%2.0f * t + %10.3f)\n",2*mop,abs(fr),darp);

};


macro plotfftR[_tt,_xx,_dtour,_fmin,_fmax,_IW] {
vnumR _ftf;
vnumC _fza;

_yy=_xx*0;

_step = (_tt[2]-_tt[1]);
fft(_xx,_yy,_fza,_ftf,_step,_tt[1],_dtour,_IW);

// divise par 1440 
ifft(_fza,Tx,Ty,_step,_tt[1]);
write("IFFT.dat",Tx,Ty);
strange=str(_fmin)+":"+str(_fmax);
gnuplot;
set xrange [@strange@]
end;
// Plot de la frequence par rapport au log de amplitude trouvee
//plot(_ftf,log10(abs(_fza))," notitle w l");
};



%asos4;



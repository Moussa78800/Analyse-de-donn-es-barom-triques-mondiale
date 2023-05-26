

macro conversion{

vnumR temps,P1,P2,P3;
file="K6R62018.dat";
//read(file, temps,P1,P2,P3);
%readasos4[file];
mimin=0,24*60-1;

mP1=mimin*0;
mP2=mimin*0;

nP=mP1;


%acump3;



mP1=mP1/nP;
mP2=mP2/nP;
ihg2mbar=33.8637526;

mbP1=mP1*ihg2mbar;
mbP2=mP2*ihg2mbar;

mbP0=(mbP1+mbP2)/2;

P1off=sum(mbP1-mbP0)/size(mbP1);
P2off=sum(mbP2-mbP0)/size(mbP1);

xmbP1=mbP1-P1off;
xmbP2=mbP2-P2off;
n=size(xmbP2);
DxmbP1=xmbP2-sum(xmbP2)/n;
vnumR IFFTR,IFFTI;
gnuplot;
#set terminal png;
#set output "ST+TFINV.png";
set xlabel "Temps(Heures locales)";
set ylabel "Pression(mbar)";
set xrange[0:24];
set xtics 1;
end;
write("DP1UTC.dat",DxmbP1);
plot(mimin/60,DxmbP1,"w l lw 2/1 lt 3 title'ST'");
/*
min_idx=1;
for j=1 to size(DxmbP1){
if (DxmbP1[j]>DxmbP1[min_idx])  then
{

min_idx=j;
};
};
*/
vnumR IFFTR,IFFTI;
//read("IFFT.dat",IFFTR,IFFTI);
Pas=1,4096;
//replot(Pas/170,IFFTR," w l lw 1 lt 2 title'TFINV'");

};
macro acump3{
nutc1=hutc*60+mutc;
nutc=hlt*60+mlt;
//nutc1+int(19.18387179970457/3600);
for k=1 to size(P1) {
    kk=k$
    nP[nutc[k]+1]=nP[nutc[k]+1]+1$
    mP1[nutc[k]+1]= mP1[nutc[k]+1]+ P1[k]$
    mP2[nutc[k]+1]= mP2[nutc[k]+1]+ P2[k]$
    
   };
};

macro readasos4[file]{

vnumR WBAN,Yr,Mo,Dy,hlt,mlt,hutc,mutc,P1,P2,P3,T1,T2;
vnumR pl,ps,sfr;
// ch14,ch24,ch3
form="%5g%4s%4s%4g%4g%2g%2g%2g%2g %41s  %g  %g ";
_read = ( "skipbrokenlines", "skipemptylines" );
_info off;
read(file,form, WBAN,ch14,Yr,Mo,Dy,hlt,mlt,hutc,mutc, sfr,P1,P2);
max(P1);
max(P2);
min(P1);
min(P2);
/*
min_idx=1;
for j=1 to size(P2){
if (P2[j]>P2[min_idx])  then
{

min_idx=j;
};
};
*/
};

%conversion;

// Code pour la superposition d1 signal barométrique replié sur 1440 minutes, i.e 24 heures avec une approximation d1e somme de fourier limitée à N Harmoniques

//Implementation macro pour le calcul des coefficients de Fourier
//Implementation de la macro pour faire approximation et superposition avec le signal barométrique original( de départ)
macro ASF[N]{
// Periode 
T=1440;
t=0,T-1,1;
vnumC c,Coeffs;
resize(c,N);
resize(Coeffs,N);
vnumR DP;
resize(DP,size(t));
read("DP1UTC.dat",DP);

//Boucle for pour calculer les coefficients de Fourier

for j=1 to N
{
c=sum(DP*exp(-2*i*pi*j*t/size(t)))/T$
Coeffs[j]=c$
};
write("Coeffs.dat",Coeffs);
write("2ABSCoeffs.dat",2*abs(Coeffs));
// Les amplitudes des oscillations 
//writes(2*abs(Coeffs));
Freq=1,24$

// On peut traçer le spectre egalement
gnuplot;
#set terminal png;
#set output "SpectreAFourier.png"
set xlabel "Temps solaire(h)"
set ylabel "Patm(mbar)";
set xrange[0:24];
set xtics 1;
end;

// Traçage du Spectre de Frequence
plot(Freq, 2*abs(Coeffs)," notitle w impulses lw 3");

vnumR CReal,CImg,Y_r;
resize(Y_r,size(t));
resize(CReal,N);
resize(CImg,N);
read("Coeffs.dat",CReal,CImg);
//writes(CReal,CImg);
//writes(2*abs(CReal+i*CImg));
Y_r=Y_r+CReal[1]/size(t);

//Boucle for pour calculer les valeurs pour la reconstitution

for j=1 to N
{
Y_r=Y_r+2*(CReal[j]*cos(2*pi*j*t/T)-CImg[j]*sin(2*pi*j*t/T))$
};
//writes(Y_r);
//Affichage du Signal original superposer au signal reconstruit par Somme de Fourier
gnuplot;
#set terminal png;
#set output "Superposition.png"
set xlabel "Temps solaire(h)"
set ylabel "Patm(mbar)";
set xrange[0:24];
set xtics 1;
end;
//plot(t/60,Y_r,"w l lw 3 title'5 Harmoniques'");
//replot(t/60,DP,"w l lw 2 title'Data'");


};

// reconstitution du signal barometrique de départ avec 10 harmoniques obtenues par analyse de Fourier.
%ASF[5];


//reconstitution par Analyse en Frequence

macro SAF{
/*
Cette macro trace le spectre de frequence obtenu par une analyse en frequence
*/

vnumR A, F;
resize(A,3);
resize(F,3);
A[1]=6.4415899827929457E-01;
A[2]=6.3670963094595445E-01;
A[3]=1.5614485696110281E-01;
F[1]=1.8179214096567367E+00;
F[2]=-1.8052391485785577E+00;
F[3]=3.5808427818290484E-01;
gnuplot;
#set terminal png;
#set output "SpectreAFreq.png"
set xlabel "Temps solaire(h)"
set ylabel "Patm(mbar)";
set xrange[0:24];
set xtics 1;
end;
plot(F,2*A,"w impulses lw 3 title'Spectre'");
};

//%SAF;

macro RAF[N]{
T=1440;
t=0,T-1,1;
vnumC c,Coeffs;
resize(c,N);
resize(Coeffs,N);
vnumR DP;
resize(DP,size(t));
read("DP1UTC.dat",DP);

//Boucle for pour calculer les coefficients 

for j=1 to N
{
c=sum(DP*exp(-i*2*pi*j*t/size(t)))/T$
Coeffs[j]=c$
};
write("CoeffsAF.dat", Coeffs);



vnumR CReal,CImg,Y_r;
resize(Y_r,size(t));
resize(CReal,N);
resize(CImg,N);
read("CoeffsAF.dat",CReal,CImg);
//writes(CReal,CImg);
//writes(2*abs(CReal+i*CImg));
Y_r=Y_r+CReal[1]/size(t);

//Boucle for pour calculer les valeurs pour la reconstitution

for j=1 to N
{
Y_r=Y_r+2*(CReal[j]*cos(2*pi*j*t/T)-CImg[j]*sin(2*pi*j*t/T))$
};
//writes(Y_r);
//Affichage du Signal original superposer au signal reconstruit par Somme de Fourier
gnuplot;
#set terminal png;
#set output "Superposition.png"
set xlabel "Temps solaire(h)"
set ylabel "Patm(mbar)";
set xrange[0:24];
set xtics 1;
end;
//plot(t/60,Y_r,"w l lw 3 title'Signal reconstruit avec 5 harmoniques'");
replot(t/60,DP,"w l lw 2 title'Signal original'");

};

//%RAF[5];

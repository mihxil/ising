/* ------------------------------------------------------
   Ising Model 
   test-version september 1998

   based on:
   opgave bij MOSI
   Michiel Meeuwissen
   november 1996
   
   ------------------------------------------------------
*/    


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <iostream.h>

const hrandmax = RAND_MAX /2;
enum  boolean {False,True};

// Constanten van het Model:
const     unsigned maxN=511;  
   // maximaal te kiezen N, want er is geen moeite gedaan om geheugen 
   // dynamisch te reserveren.

// Globale instellingen met default waarden:
float    J=1;        // koppeling
float    B=0;        // magneetveld
unsigned Nx=10, Ny = 10;
unsigned long nit;     // number of iterations for one calculation
unsigned long Ncalc;   // number of calculations

boolean  terugleggen=True; // zonder nog niet geimplementeerd
char     upc='|';
char     downc='-';

// Voor en random getal van 0 tot en met n:
unsigned rndom(unsigned);

// Voor het rooster:

// Een paar types:
const int sint = sizeof(int) * 8,  //bits per int (that many spins we can store in a int).
  Nint = maxN/sint + 1;   // this many integers we need to store our spins.
typedef unsigned rij[Nint];
typedef rij  rooster[maxN];

enum spin {down, up};

// Een paar functies
void setrooster(rooster,int,int, spin);
spin getrooster(rooster,int,int);
void initrooster(rooster);
void initrooster(rooster,spin);
void printrooster(rooster);
unsigned aantalup(rooster);


// It would be nice to make our 'rooster' a class.
// TODO I'll do that later.


/* ---------------------------
   Metropolis algoritme.
   'Ising Model'
   
   zie ook: Reif, F. - Fundamentals of Statistical and Thermal Physics 
                     - McGraw-Hill 1985
*/

inline float fspin(spin m) { return m==up ? 1 : -1; }
// translates a spin to a float. TODO Clearly it would be nice to make a spin 
// a class also such that we can define some operators...


void itereer(rooster r)
{ 
  unsigned i,j;
  // kies willekeurige lokatie (metropolis):
  
  if(terugleggen) {  i=rndom(Nx-1);  j=rndom(Ny-1); }
  else {   }  

  // uitrekenen energie:

  float E;
  float sij; 
  sij = fspin(getrooster(r,i,j));


  E= ( J * ( - fspin( getrooster(r,(i-1+Nx)%Nx,j)        )
           - fspin( getrooster(r,(i+1)%Nx,  j)        )
           - fspin( getrooster(r,i,        (j+1)%Ny)  )
           - fspin( getrooster(r,i,        (j-1+Ny)%Ny))
         )
        -B)*sij;
    
  // toekennen nieuwe waarde:
  float p;
  p=(float)rand()/RAND_MAX;
  // printf("E:%f  %f  %f\n",E,1.0/(1+exp(-2*E)),p);
  if(p>=1.0/(1+exp(+2*E))){ setrooster(r,i,j,(sij==1 ? down : up)); }
                     
  nit++;
}

/* ------------------------------------------------------------------------
   Swendsen / Wang 
   bond algoritme
   
  zie ook Phys. Rev. Letters {\bf 58} (1987) 86

*/

void verander(rooster r, int i, int j, rooster hulp, rooster b[2], spin R)
// Deze functie is in 4 richtingen recursief.
{
  setrooster(r,i,j,R);
  setrooster(hulp,i,j,up);
  if((getrooster(b[0],i,j)==up) && (getrooster(hulp,(i+1)%Nx,j)==down)) verander(r,(i+1)%Nx,j,hulp,b,R);
  if((getrooster(b[0],(i-1+Nx)%Nx,j)==up) && (getrooster(hulp,(i-1+Nx)%Nx,j)==down)) verander(r,(i-1+Nx)%Nx,j,hulp,b,R);
  if((getrooster(b[1],i,j)==up) && (getrooster(hulp,i,(j+1)%Ny)==down)) verander(r,i,(j+1)%Ny,hulp,b,R);
  if((getrooster(b[1],i,(j-1+Ny)%Ny)==up) && (getrooster(hulp,i,(j-1+Ny)%Ny)==down)) verander(r,i,(j-1+Ny)%Ny,hulp,b,R);
}


void itereerbonds(rooster r)
{ 
  rooster b[2]; // horizontale en verticale bonds; 
  rooster hulp; // hierin worden de al veranderde roosterpunten bij gehouden.

  spin p;
  int i,j;
  // het maken van de bonds:
  for(i=0; i<Nx; i++)
    for(j=0; j<Ny; j++)
      { 
	p=getrooster(r,i,j);
	if(p==getrooster(r,(i+1)%Nx,j))
	  if(((float)rand()/RAND_MAX)<=(1-exp(-2*J))) setrooster(b[0],i,j,up);
	  else setrooster(b[0],i,j,down);
	else setrooster(b[0],i,j,down);
	if(p==getrooster(r,i,(j+1)%Ny))
	  if(((float)rand()/RAND_MAX)<=(1-exp(-2*J))) setrooster(b[1],i,j,up);
	  else               setrooster(b[1],i,j,down);
	else setrooster(b[1],i,j,down);
      }
  // clusters omswitchen (of niet):
  initrooster (hulp, down);

  unsigned nc=0; // aantal cluster wordt bijgehouden (maar dat is nergens goed voor)
  float R;
  for(i=0; i<Nx; i++)
    for(j=0; j<Ny; j++)
      if(getrooster(hulp,i,j)==down)
	{ 
	  nc++;
	  R=(float)rand()/RAND_MAX;
	  verander(r,i,j,hulp,b, R < 1/(1+exp(-4*B)) ? up : down);
	}     // tis me nog niet goed duidelijk waarom hier 4 moet staan
  // printf("nc: %i  ",nc);
  nit++;
}

/* -------------------------------------------------------------------
   MAIN

*/


int main(int argc, char *argv[])
{
  float deJ=J,deB=B;
  float minJ=0, minB=0;
  boolean bond=False;   // 'Bond' of 'Metropolis' algoritme 
  boolean tekenr=False; // Of het rooster getekend moet worden.
  boolean schrijfnit=True; // Of tijdens de berekening het iteratienummer zichtbaar moet zijn. 
  boolean Sweep=False; // Of B of J moeten worden gevarieerd
  boolean SweepoverJ;  // Of dat B of J is.
  boolean calcQ=False;

  unsigned aantal=100; //  aantal waarden welke B en J kunnen aannemen.  
  boolean eenheidmaxnitNN=False; // de eenheid waarin maxnit werd opgegeven.
                                 // dit kan N^2 worden.
  long unsigned maxnit=100; // hoeveel iteraties er worden gedaan.
  int seed=1;
  char *commentaar=NULL;

  unsigned Ncalc = 20; // hoevaak de berekening herhaald wordt om Q te kunnen berekenen.  
  
  // Via de commandline kan veel worden ingesteld:
  if (argc>1)
     { 
       if(argv[1][0]=='h')
       {
	 printf(" Ising Model -- September 1998 \n Michiel Meeuwissen \n"
		" gebruik %s <opties> \n\n"
		" met opties: \n"
		" J??? : specificeer  waarde van J (%f)\n"          
		" B??? : specificeer een B-veld (%f)\n"
		" S??  : sweep: SJ?? over J; SB?? over B; ??:aantal (%u).\n"
		" b??  : bij sweep de startwaarde: bJ?? van J; bB?? van B (%f,%f)\n"
		" N??  : specificeer grootte van rooster (Nx.Ny) (%u, %u) \n"
		" n??  : specificeer aantal iteraties (1 sweep n=NxN) (%u)\n"
		" T    : teken het eerste en laatste rooster \n"
		" c??  : te gebruiken up en down karakter (%c,%c)\n"
		" C??? : extra te schrijven commentaar ()\n"
		" s??  : specificeer random seed (%i) \n"
		" t    : schrijf tussendoor niet nit en aantalup \n"
		" m?   : methode: mm: Metropolis - volledig random \n" 
		"                 mt: Metropolis - zonder terugleggen (niet geimplementeerd)\n"
		"                 mb: Swendsen/Wang - bond-algoritme\n"
                " Q??  : Berekenen amplitude ratio. Daarvoor wordt de berekening ?? keer herhaald (%u).\n"

		,argv[0],J,B,aantal,minJ,minB,Nx, Ny, maxnit,upc,downc,seed, Ncalc);
	 exit(1);
       }
       for (int arg=1; arg<argc; arg++)
         { 
	   if (argv[arg][0] == 'J') sscanf(argv[arg]+1,"%f",&deJ);
           if (argv[arg][0] == 'B') sscanf(argv[arg]+1,"%f",&deB);
           if (argv[arg][0] == 'N')
	     {
	       switch(argv[arg][1])
		 {
		 case 'x': sscanf(argv[arg]+2,"%u",&Nx); break;
		 case 'y': sscanf(argv[arg]+2,"%u",&Ny); break;
                 default:  sscanf(argv[arg]+1,"%u",&Nx); Ny = Nx; break;
		 }
	     }
           if (argv[arg][0] == 'r') sscanf(argv[arg]+1,"%i",&seed);
           if (argv[arg][0] == 'T') tekenr=True;
           if (argv[arg][0] == 'c') 
               if(argv[arg][1]!='\0') 
                 { 
		   upc=argv[arg][1];
                   if(argv[arg][2]!='\0') downc=argv[arg][2];
                 }  
           if (argv[arg][0]=='C') commentaar=argv[arg]+1;      
           if (argv[arg][0] == 'S')
              { 
		if(argv[arg][1]=='J') 
		  {
		    Sweep=True; SweepoverJ=True;
		    if(argv[arg][2]!='\0')
		      sscanf(argv[arg]+2,"%u",&aantal);
		  }
                if(argv[arg][1]=='B') 
		  {
		    Sweep=True; SweepoverJ=False;
		    if(argv[arg][2]!='\0')
		      sscanf(argv[arg]+2,"%u",&aantal);
		  }
              }  
           if (argv[arg][0] == 'b')
              { 
		if((argv[arg][1] == 'J') && (argv[arg][2] != '\0') ) sscanf(argv[arg]+2,"%f",&minJ);
                if((argv[arg][1] == 'B') && (argv[arg][2] != '\0') ) sscanf(argv[arg]+2,"%f",&minB);
              }                       
           if (argv[arg][0] == 'm') 
              { 
		if(argv[arg][1] == 'm') {bond=False; terugleggen=True;}
                if(argv[arg][1] == 't') {bond=False; terugleggen=False;} // dit werkt nog niet
                if(argv[arg][1] == 'b') {bond=True;}              
              }  
           if (argv[arg][0] == 't') schrijfnit=False;
           if (argv[arg][0] == 'n') 
                if (argv[arg][1] == 'N') { eenheidmaxnitNN=True;
                                           sscanf(argv[arg]+2,"%lu",&maxnit);
                                         }  
                else sscanf(argv[arg]+1,"%lu",&maxnit);
           if (argv[arg][0] == 's') sscanf(argv[arg]+1,"%i",&seed);
	   if (argv[arg][0] == 'Q')  
	     {
	       calcQ = True;
               if(argv[arg][1] != '\0') sscanf(argv[arg]+1, "%u", &Ncalc);
	     }

         }
     }
     
  if (Nx>maxN) { printf("Nx kan maximaal %i zijn. (of compileer opnieuw)\n",maxN); return 1;}
  if (Ny>maxN) { printf("Ny kan maximaal %i zijn. (of compileer opnieuw)\n",maxN); return 1;}
  if (!bond && !terugleggen) {printf("sorry, 'met terugleggen' werkt nog niet...."); return 1;}
  if (eenheidmaxnitNN) maxnit = maxnit * Nx * Ny;
  printf("# meth: %s  J:%f  B:%f  Nx:%u Ny:%u maxnit:%lu  seed:%i\n",
	 bond ? "Bond" : "Metropolis", deJ,deB,Nx,Ny,maxnit,seed);
  if (commentaar!=NULL) printf("# %s\n",commentaar);              
  J=deJ; B=deB;

  rooster hetrooster;
  // rooster hulprooster; // voor in 'zonder terugleggen' 
  // initrooster(hulprooster,down);
  // unsigned long ne=0;

  srand(seed);
  
  nit=0;
  long unsigned  i;
  unsigned nup;





  if(calcQ)
    {
      float sumM2;
      float sumM4;     
      float M, M2;
       
      sumM2 = 0;
      sumM4 = 0; 

      for(unsigned ncalc = 0; ncalc < Ncalc; ncalc++)
        {
	  initrooster(hetrooster);
	  for(i=1; i<=maxnit; i++)
	    { 
	    if(bond) itereerbonds(hetrooster);
	    else     itereer(hetrooster);
	    }
	  nup = aantalup(hetrooster);
	  
	  M = (float)2 * nup - Nx*Ny;
          M2 = M*M;
          sumM2 += M2;
          sumM4 += M2*M2;          
	//  cout << nup << sumM2 << "\n";
	}

      printf(" <M^2> = %e \n <M^4> = %e \n Q = %e \n", (float)sumM2/Ncalc, (float)sumM4/Ncalc, (float)sumM2*sumM2 / (Ncalc * sumM4));
      return 0;
    }


 
  for (int k=0; k<=aantal; k++)
    {
      if(Sweep)
	{ 
	  if(SweepoverJ) J=minJ+(float)k*(deJ-minJ)/aantal;
	  else   B=minB+(float)k*(deB-minB)/aantal;
	}
      else k=aantal; 
      initrooster(hetrooster);
      if (tekenr) printrooster(hetrooster);
      for(i=1; i<=maxnit; i++)
	{ 
	  if(schrijfnit)
	    { 
	      nup=aantalup(hetrooster);
	      printf(" nit: %lu up: %u down: %u  \r",nit, nup, Nx*Ny -nup);
	    }  
	  if(bond) itereerbonds(hetrooster);
	  else     itereer(hetrooster);
	}
      if (tekenr){printf("\n"); printrooster(hetrooster);}
      nup=aantalup(hetrooster);
      if (Sweep) printf("%e   %e   %u   %u  \n", J, B, nup, Nx*Ny -nup);
      else  printf("up: %u  down: %u                   \n",nup,Nx*Ny-nup);
    }



  return 0;
}

/* --------------------------------------------------------------------
  Implementatie van een paar 'rooster' functies  en andere

  TODO This of course should be methods of some class.
       (at least a namespace)
*/

// geeft een random getal van 0 tot en met n;
unsigned rndom(unsigned n)
{ 
  unsigned h1, h2,n1;
  n1=n+1;
  h1=(RAND_MAX / n1) * n1;
  h2=RAND_MAX;
  while(h2>=h1) h2=rand();
  return (h2 % n1);
}


void setrooster(rooster r, int i, int j, spin v)
{ 
  if(v==up) { r[i][j/sint] |= (1<<(j%sint) );}
  else { r[i][j/sint] = ~( ~(r[i][j/sint]) | (1<<(j%sint)) ); }
}

spin getrooster(rooster r, int i, int j)
{ 
  return (spin)((r[i][j/sint] & (1 << (j%sint))) >> (j%sint));
}

void initrooster(rooster r)
// initialiseert het rooster r met willekeurige waarden.
{
  for(int i=0; i<Nx; i++)     
    for(int j=0; j<Ny; j++)
      setrooster(r,i,j,(spin)rndom(1));
  return;  
}
void initrooster(rooster r, spin v)
// initialiseert het rooster r met v
{
  for(int i=0; i<Nx; i++)
        for(int j=0; j<Ny; j++)
            setrooster(r,i,j,v);
  return;
}

void printrooster(rooster r)
{
  for(int i=0; i<Nx; i++)
      {
        for(int j=0; j<Ny; j++) printf("%c",getrooster(r,i,j)==up ? upc : downc); 
        printf("\n");
      }
  return;
}

unsigned  aantalup(rooster r)
{ 
  unsigned hulp=0;
  for(int i=0; i<Nx; i++)
    for (int j=0; j<Ny; j++)
      hulp+=(int)getrooster(r,i,j);
  return hulp;
}

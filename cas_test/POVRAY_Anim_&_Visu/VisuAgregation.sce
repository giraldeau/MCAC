clear
cd 'C:\Users\dlca\Desktop\DLCA_sous_Qt\POVRAY_Anim_&_Visu'
exec('DLCA2POV.sci'); //procedure de d'export pour POV
exec('lecture.sce'); //Charge la conversion C++ -> sci

//Export  l'aggregat "AggSelect" vers POV ou DDSCAT
//Si AggSelect=0 alors on exporte l'ensemble des agregats
//si AggSelect!=0 et AxeProj=0 alors on conserve l'orientation initiale de l'agrégat

//si AggSelect!=0 et AxeProj=i alors on aligne l'agrégat suivant son ième axe principale d'inertie
//si AggSelect!=0 alors tetha représente l'angle de rotation autour de l'axe vertical du plan d'observation au microscope

nomfichier='params.txt';
NumFichierAgg=994;

AggSelect=3;
AxeProj=0;
tetha=0*%pi/180;

cd 'C:\Users\dlca\Desktop\DLCA_sous_Qt'

[temps,NbAgg,N,T,P,Fv,X,L,rho,ModeDistrib,Dpmodal,sigmaDp,Chemin]=LectureInfoGenerales(nomfichier,NumFichierAgg);
if AggSelect>NbAgg then 
  AggSelect=NbAgg;
end

[temps,NbAgg,RgTab,NpTab,DmTab,lpmTab,deltatTab,RmaxTab,XTab,YTab,ZTab,N,T,P,Fv,X,L,rho,ModeDistrib,Dpmodal,sigmaDp,Chemin]=LectureAgg(nomfichier,NumFichierAgg);
[temps,NbAgg,Label,RpTab,Posi,N,T,P,Fv,X,L,rho,ModeDistrib,Dpmodal,sigmaDp,Chemin]=LectureSphere(nomfichier,NumFichierAgg);
RpTab=RpTab*1E9;
Posi=Posi*1E9;
Dpmodal=Dpmodal;
RgTab=RgTab*1E9;


//Commande : function imprima(Nom_suie,Posi,RpTab,Label,L,AggSelect,AxeProj,tetha,comment)

//comment=sprintf('Fichier N°%1.5d LabelAgregat N°%1.5d  Rg=%1.5f nm',NumFichierAgg,AggSelect,RgTab(AggSelect));
comment='';
imprima('test-pov',Posi,RpTab,Label,L,AggSelect,AxeProj,tetha,comment)

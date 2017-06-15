//Script permettant de générer les fichiers POV nécéssaires à la creation d'une animation


clear

nomfichier='params.txt';  //Nom du fichier contenant les paramètres d'execution du code DLCA
NbdeltatT=2000;             //Nombre d'images souhaitées dans l'animation


exec('DLCA2POV.sci'); //procedure de d'export pour POV
exec('lecture.sce'); //Charge la conversion C++ -> sci


//Récupère le fichier de donnees i et imprime avec l'index indi
function ImprimePOV(i,NomSortie,L)  //Récupère le fichier de donnees i et imprime avec l'index indi
    //Lecture du fichier modele pour calcul POV
    Nom_Model_POV="agregat-model.pov";
    f1=mopen(Nom_Model_POV,'r');
    Nblignes=0;
    while meof(f1)==0 then 
      Nblignes=Nblignes+1;
      Fichier(Nblignes)=mgetl(f1,1);
    end
    mclose(f1);
    //On adapte les lignes à modifier  
    Tmp=msprintf('#include ""agregat-DLCA-%1.6d.inc""',i); 
    Fichier(2)=Tmp;
      
    Tmp=msprintf('#declare L =%f;',L);
    Fichier(7)=Tmp;
      
    //On ecrit le fichier
    fout=mopen(NomSortie,'w');
    mfprintf(fout,'%s\n',Fichier);
    mclose(fout); 
endfunction
  
  
//Lecture du fichier contenant les paramètres d'execution 
[temps,NbAgg,N,T,P,Fv,X,L,rho,ModeDistrib,Dpmodal,sigmaDp,Chemin]=LectureInfoGenerales(nomfichier,0);
infodir=ls(Chemin+'\Sphere*.txt');
NbFich=size(infodir,'*')-1;
mkdir(Chemin+'\AnimPOV\');

//Lecture des fichiers de donnees pour trier les fichiers dans 
//l'ordre croissant du temps

for i=1:NbFich   
  //Tmp=msprintf('%1.5d.txt',i);
  [temps,NbAgg,N,T,P,Fv,X,L,rho,ModeDistrib,Dpmodal,sigmaDp,Chemin]=LectureInfoGenerales(nomfichier,i);
  Classement(i)=temps;
  NAgg(i)=NbAgg;
end

[ordonne,indiceclass]=gsort(Classement,'g','i');
plot2d(Classement(indiceclass),NAgg(indiceclass))
tmax=ordonne($);

//On balaye le temps à intervalle de temps régulier pour observer l'agglomeration
deltaT=tmax/NbdeltatT;
temps=0;
indice=0;
for i=1:NbFich
 // if Classement(indiceclass(i))>temps then
    indice=indice+1;
    printf('image N° : %d\n',indice)
    NomSortie1=Chemin+'\AnimPOV\'+msprintf('animation-%0.4d.pov',indice);
    ImprimePOV(indiceclass(i),NomSortie1,L);
    [temps,NbAgg,Label,RpTab,Posi,N,T,P,Fv,X,L,rho,ModeDistrib,Dpmodal2,sigmaDp,Chemin]=LectureSphere(nomfichier,indiceclass(i));
    RpTab=RpTab*1E9;
    Posi=Posi*1E9;    
    NomSortie2=Chemin+'\AnimPOV\'+msprintf('agregat-DLCA-%1.6d',indiceclass(i));
    imprima(NomSortie2,Posi,RpTab,Label,L,1,0,0*%pi/180)
    temps=temps+deltaT;
//end
end





  
  
  


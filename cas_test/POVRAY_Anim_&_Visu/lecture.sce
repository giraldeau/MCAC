////////////////////////////////////////////////////////////
//Fonction qui retourne toutes les informations nécéssaires
////////////////////////////////////////////////////////////
 function [temps,NbAgg,RgTab,NpTab,DmTab,lpmTab,deltatTab,RmaxTab,XTab,YTab,ZTab,Label,RpTab,Posi,N,T,P,Fv,X,L,rho,ModeDistrib,Dpmodal,sigmaDp,Chemin]=Lecture(nomfichier,NumFichier)
   //Lecture des conditons génériques du calcul
   fichierparam = mopen(nomfichier,"rb");
   [nb,N,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,T,Comment]=mfscanf(fichierparam,"%lf\t%s\n");
   [nb,Dpmodal,Comment]=mfscanf(fichierparam,"%e\t%s\n");
   [nb,sigmaDp,Comment]=mfscanf(fichierparam,"%e\t%s\n");
   [nb,Fv,Comment]=mfscanf(fichierparam,"%e\t%s\n");Fv=Fv*1E-6;
   [nb,P,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,rho,Comment]=mfscanf(fichierparam,"%e\t%s\n");
   [nb,ModeDistrib,Comment]=mfscanf(fichierparam,"%d\t%s\n"); //Mode:1=Normal,0=LogNormal
   [nb,deltasauve,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,Chemin,Comment]=mfscanf(fichierparam,"%s\t%s\n");
   X=(N*%pi/6/Fv*exp(9/2*(log(sigmaDp))^2))^(1/3);
   L=X*Dpmodal;                //Largeur de la boite en nm
   mclose(fichierparam);
   
   
   //Lecture du fichier de données N° NumFichier associé aux agrégats
    fichierData = mopen(Chemin+'\Agg'+sprintf('%0.5d\n',NumFichier)+'.txt',"rb");
    for i=1:5, Comment=mgetl(fichierData,1); end  //on saute 5 lignes
   [nb,NbAgg,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,temps,Comment]=mfscanf(fichierparam,"%e\t%s\n");temps=temps*1E-6;
    Comment=mgetl(fichierData,1);
    for i=1:NbAgg
      [nb,Val1,Val2,Val3,Val4,Val5,Val6,Val7,Val8,Val9]=mfscanf(fichierparam,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"); 
      RgTab(i)=Val1*1E-9;
      NpTab(i)=Val2;
      DmTab(i)=Val3*1E-9;
      lpmTab(i)=Val4*1E-9;
      deltatTab(i)=Val5*1E-6;
      RmaxTab(i)=Val6*1E-9;
      XTab(i)=Val7*1E-9;
      YTab(i)=Val8*1E-9;
      ZTab(i)=Val9*1E-9;
    end
    mclose(fichierData);
    
    
   //Lecture du fichier de données N° NumFichier associé aux à la description des sphères
    fichierData = mopen(Chemin+'\Sphere'+sprintf('%0.5d\n',NumFichier)+'.txt',"rb");
    for i=1:5, Comment=mgetl(fichierData,1); end  //on saute 5 lignes
   [nb,NbAgg,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,temps,Comment]=mfscanf(fichierparam,"%e\t%s\n");temps=temps*1E-6;
    Comment=mgetl(fichierData,1);
    for i=1:N
      [nb,Val1,Val2,Val3,Val4,Val5,Val6,Val7,Val8,Val9]=mfscanf(fichierparam,"%d\t%lf\t%lf\t%lf\t%lf\n"); 
      Label(i)=Val1;            //Label de l'agrégat d'appartenance
      RpTab(i)=Val2*1E-9;            //Rayon de chaque sphérules en m
      Posi(i,1)=Val3*1E-9;           //Posi : positions x,y et z de chaque sphérules en m
      Posi(i,2)=Val4*1E-9;
      Posi(i,3)=Val5*1E-9;
    end
    mclose(fichierData);
endfunction
  
  
  
  
  
  
  ////////////////////////////////////////////////////////////
//Fonction qui retourne toutes les informations liées aux agrégats
////////////////////////////////////////////////////////////
 function [temps,NbAgg,RgTab,NpTab,DmTab,lpmTab,deltatTab,RmaxTab,XTab,YTab,ZTab,N,T,P,Fv,X,L,rho,ModeDistrib,Dpmodal,sigmaDp,Chemin]=LectureAgg(nomfichier,NumFichier)
   //Lecture des conditons génériques du calcul
   fichierparam = mopen(nomfichier,"rb");
   [nb,N,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,T,Comment]=mfscanf(fichierparam,"%lf\t%s\n");
   [nb,Dpmodal,Comment]=mfscanf(fichierparam,"%e\t%s\n");
   [nb,sigmaDp,Comment]=mfscanf(fichierparam,"%e\t%s\n");
   [nb,Fv,Comment]=mfscanf(fichierparam,"%e\t%s\n");Fv=Fv*1E-6;
   [nb,P,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,rho,Comment]=mfscanf(fichierparam,"%e\t%s\n");
   [nb,ModeDistrib,Comment]=mfscanf(fichierparam,"%d\t%s\n"); //Mode:1=Normal,0=LogNormal
   [nb,deltasauve,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,Chemin,Comment]=mfscanf(fichierparam,"%s\t%s\n");
   X=(N*%pi/6/Fv*exp(9/2*(log(sigmaDp))^2))^(1/3);
   L=X*Dpmodal;                //Largeur de la boite en nm
   mclose(fichierparam);
   
   
   
   //Lecture du fichier de données N° NumFichier associé aux agrégats
    fichierData = mopen(Chemin+'\Agg'+sprintf('%0.5d\n',NumFichier)+'.txt',"rb");
    for i=1:5, Comment=mgetl(fichierData,1); end  //on saute 5 lignes
   [nb,NbAgg,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,temps,Comment]=mfscanf(fichierparam,"%e\t%s\n");temps=temps*1E-6;
    Comment=mgetl(fichierData,1);
    for i=1:NbAgg
      [nb,Val1,Val2,Val3,Val4,Val5,Val6,Val7,Val8,Val9]=mfscanf(fichierparam,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"); 
      RgTab(i)=Val1*1E-9;
      NpTab(i)=Val2;
      DmTab(i)=Val3*1E-9;
      lpmTab(i)=Val4*1E-9;
      deltatTab(i)=Val5*1E-6;
      RmaxTab(i)=Val6*1E-9;
      XTab(i)=Val7*1E-9;
      YTab(i)=Val8*1E-9;
      ZTab(i)=Val9*1E-9;
    end
    mclose(fichierData);
endfunction

  
  
  
  
  
////////////////////////////////////////////////////////////
//Fonction qui retourne toutes les informations liées aux Sphérules
////////////////////////////////////////////////////////////
 function [temps,NbAgg,Label,RpTab,Posi,N,T,P,Fv,X,L,rho,ModeDistrib,Dpmodal,sigmaDp,Chemin]=LectureSphere(nomfichier,NumFichier)
   //Lecture des conditons génériques du calcul
   fichierparam = mopen(nomfichier,"rb");
   [nb,N,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,T,Comment]=mfscanf(fichierparam,"%lf\t%s\n");
   [nb,Dpmodal,Comment]=mfscanf(fichierparam,"%e\t%s\n");
   [nb,sigmaDp,Comment]=mfscanf(fichierparam,"%e\t%s\n");
   [nb,Fv,Comment]=mfscanf(fichierparam,"%e\t%s\n");Fv=Fv*1E-6;
   [nb,P,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,rho,Comment]=mfscanf(fichierparam,"%e\t%s\n");
   [nb,ModeDistrib,Comment]=mfscanf(fichierparam,"%d\t%s\n"); //Mode:1=Normal,0=LogNormal
   [nb,deltasauve,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,Chemin,Comment]=mfscanf(fichierparam,"%s\t%s\n");
   X=(N*%pi/6/Fv*exp(9/2*(log(sigmaDp))^2))^(1/3);
   L=X*Dpmodal;                //Largeur de la boite en nm
   mclose(fichierparam);

   //Lecture du fichier de données N° NumFichier associé aux à la description des sphères
    fichierData = mopen(Chemin+'\Sphere'+sprintf('%0.5d\n',NumFichier)+'.txt',"rb");
    for i=1:5, Comment=mgetl(fichierData,1); end  //on saute 5 lignes
   [nb,NbAgg,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,temps,Comment]=mfscanf(fichierparam,"%e\t%s\n");temps=temps*1E-6;
    Comment=mgetl(fichierData,1);
    for i=1:N
      [nb,Val1,Val2,Val3,Val4,Val5,Val6,Val7,Val8,Val9]=mfscanf(fichierparam,"%d\t%lf\t%lf\t%lf\t%lf\n"); 
      Label(i)=Val1;            //Label de l'agrégat d'appartenance
      RpTab(i)=Val2*1E-9;            //Rayon de chaque sphérules en m
      Posi(i,1)=Val3*1E-9;           //Posi : positions x,y et z de chaque sphérules en m
      Posi(i,2)=Val4*1E-9;
      Posi(i,3)=Val5*1E-9;
    end
    mclose(fichierData);    
endfunction
  
  
  
  ////////////////////////////////////////////////////////////
//Fonction qui retourne toutes les informations générales
////////////////////////////////////////////////////////////
 function [temps,NbAgg,N,T,P,Fv,X,L,rho,ModeDistrib,Dpmodal,sigmaDp,Chemin]=LectureInfoGenerales(nomfichier,NumFichier)
   //Lecture des conditons génériques du calcul
   fichierparam = mopen(nomfichier,"rb");
   [nb,N,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,T,Comment]=mfscanf(fichierparam,"%lf\t%s\n");
   [nb,Dpmodal,Comment]=mfscanf(fichierparam,"%e\t%s\n");
   [nb,sigmaDp,Comment]=mfscanf(fichierparam,"%e\t%s\n");
   [nb,Fv,Comment]=mfscanf(fichierparam,"%e\t%s\n");Fv=Fv*1E-6;
   [nb,P,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,rho,Comment]=mfscanf(fichierparam,"%e\t%s\n");
   [nb,ModeDistrib,Comment]=mfscanf(fichierparam,"%d\t%s\n"); //Mode:1=Normal,0=LogNormal
   [nb,deltasauve,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,Chemin,Comment]=mfscanf(fichierparam,"%s\t%s\n");
   X=(N*%pi/6/Fv*exp(9/2*(log(sigmaDp))^2))^(1/3);
   L=X*Dpmodal;                //Largeur de la boite en nm
   mclose(fichierparam);

   //Lecture du fichier de données N° NumFichier associé à la description des sphères
    fichierData = mopen(Chemin+'\Sphere'+sprintf('%0.5d\n',NumFichier)+'.txt',"rb");
    for i=1:5, Comment=mgetl(fichierData,1); end  //on saute 5 lignes
   [nb,NbAgg,Comment]=mfscanf(fichierparam,"%d\t%s\n");
   [nb,temps,Comment]=mfscanf(fichierparam,"%e\t%s\n");temps=temps*1E-6;
    mclose(fichierData);    
endfunction
  
  

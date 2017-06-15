
//Export  l'aggregat "AggSelect" vers POV ou DDSCAT
//Si AggSelect=0 alors on exporte l'ensemble des agregats
//si AggSelect!=0 et AxeProj=0 alors on conserve l'orientation initiale de l'agrégat

//si AggSelect!=0 et AxeProj=i alors on aligne l'agrégat suivant son ième axe principale d'inertie
//si AggSelect!=0 alors tetha représente l'angle de rotation autour de l'axe vertical du plan dobservation au microscope

function imprima(Nom_suie,Posi,RpTab,Label,L,AggSelect,AxeProj,tetha,comment)
  Select=(Label==AggSelect);
  if AggSelect==0 then 
    [Np,C]=size(Posi);
    PosiSelect=Posi-ones(Np,1)*[L/2 L/2 L/2]; //On recadre l'objet autour de son centre de masse
    RpSelect=RpTab;
  else
    PosiSelect=Posi(Select,:);
    [Np,C]=size(PosiSelect);
    RpSelect=RpTab(Select);
    Vol=0;xm=0;ym=0;zm=0;
    for j=1:Np
      xm=xm+PosiSelect(j,1)*(RpSelect(j))^3;
      ym=ym+PosiSelect(j,2)*(RpSelect(j))^3;
      zm=zm+PosiSelect(j,3)*(RpSelect(j))^3;
      Vol=Vol+(RpSelect(j))^3;
    end
    Posim=[xm ym zm]/Vol; //Coordonnée du centre de masse
    PosiSelect=PosiSelect-ones(Np,1)*Posim; //On recadre l'objet autour de son centre de masse
    
    if AxeProj>0 then  //On calcule la matrice d'inertie
      A=0;B=0;C=0;D=0;E=0;F=0;
      for j=1:Np
        x=PosiSelect(j,1);
        y=PosiSelect(j,2);
        z=PosiSelect(j,3);
        A=A+(y^2+z^2)*(RpSelect(j))^3;
        B=B+(x^2+z^2)*(RpSelect(j))^3;
        C=C+(x^2+y^2)*(RpSelect(j))^3;
        D=D+y*z*(RpSelect(j))^3;
        E=E+x*z*(RpSelect(j))^3;
        F=F+x*y*(RpSelect(j))^3;
      end
      MatIn=[A -F -E;-F B -D;-E -D C]/Vol;
      //Diagonalisation et calcul des vecteurs d'inertie
      [Propres,MatriceDiag]=spec(MatIn);
      [val,indice]=sort(diag(MatriceDiag));
      Vect1=Propres(:,indice(3)); //Vecteur principal d'inertie
      Vect2=Propres(:,indice(2)); //Vecteur secondaire d'inertie
      Vect3=Propres(:,indice(1)); //Vecteur terciaire
      
      if AxeProj==1 then
        Rotation=[Vect2 Vect3 Vect1];
      elseif AxeProj==2 then
        Rotation=[Vect1 Vect3 Vect2];
      else
        Rotation=[Vect1 Vect2 Vect3];
      end
      Rotation=1/Rotation;  

      //On tourne maintenant l'agrégat pour qu'il se presente suivant l'ace principal selectionne
      for j=1:Np
        PosiSelect(j,:)=(Rotation*PosiSelect(j,:)')';
      end
    end
    //On va maintenant faire pivoter l'agrégat de l'angle AngleRot autour de l'axe y
    phi=%pi/2;
    Rotation=[cos(tetha) 0 sin(tetha);0 1 0;-sin(tetha) 0 cos(tetha)];
    for j=1:Np
      PosiSelect(j,:)=(Rotation*PosiSelect(j,:)')';
    end
  end


  //generation d'un script pour POV
  Nom = [Nom_suie+'.inc'];
  final = mopen(Nom,'w'); //ouverture du fichier de donnees à creer pour POVRAY
  mfprintf(final,'#declare agreg=union{\n');
  //on copie les coordonnees de points
  for n=1:Np //NbMonoMax
    mfprintf(final,'sphere { <%-12.5E, %-12.5E, %-12.5E>, %-12.5E} \n',PosiSelect(n,1),PosiSelect(n,2),PosiSelect(n,3),RpSelect(n));
  end
  mfprintf(final,'}\n');
  mclose(final);

 // //generation d'un script pour DDSCAT
//  Nom = [Nom_suie+'.tab'];
//  final = mopen(Nom,'w'); //ouverture du fichier de donnees à creer pour DDSCAT
//  mfprintf(final,'%d\n',Np);  //on indique le nombre de monomères
//  mfprintf(final,'%f\n',(mean(RpSelect.^3)^(1/3)));  //on indique le rayon moyen d'ordre 3
//  mfprintf(final,'%d\t%f\t5.0\t0.0\n',Np,(mean(RpSelect.^3)^(1/3))*1E-3);  //on indique Np, le rayon moyen d'ordre 3 en µm
//  mfprintf(final,'\n'); //quelques infos
//  mfprintf(final,'\n'); 
//  //on copie les coordonnées de points
//  for n=1:Np
//    mfprintf(final,'%-12.5E %-12.5E %-12.5E %-12.5E \n',PosiSelect(n,1),PosiSelect(n,2),PosiSelect(n,3),RpSelect(n));
//  end
//  mclose(final);
     
endfunction





  
    
  

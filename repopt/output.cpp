#include <iostream>
#include <string>
#include <iomanip>
#include "auxiliary.hpp"
#include "allequations.hpp"

using namespace std;

extern Allequations allequations;
extern double gtol;
extern double expandR,deltar;
extern double ga_pmut,ga_pcross,gauss_dev;
extern int    ilmsfit,idecompose,nreplicate,dftbout_new;
extern int    ga_popsize,ga_ngen,ga_scoref,ga_flushf;
extern int    score_type,preserved_num,destroy_num,popsizemin,seed;
extern bool   ga,read_spline,fsmooth_spline,aa_spline;
extern bool   fderivative1st,fderivative2nd,fderivative3rd;
extern bool   runtest,grid_update;
extern sddh   ddh;

int size0,size1,size2,size3;

void Allequations::writeout() const {
  cout<<std::fixed;

  cout << "_____________________________________________________________________________________" << endl << endl;
  cout<<"\n### interpreted  repopt input ###\n";

  cout<<"\n$system:\n";
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  dftb_version    "<<dftbversion<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  ilmsfit         "<<ilmsfit<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  idecompose      "<<idecompose<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  nreplicate      "<<nreplicate<<endl;
  cout<<"$end\n";

  cout<<"\n$genetic_algorithm:\n";
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  popsizemax      "<<ga_popsize<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  ngen            "<<ga_ngen<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  preserved_num   "<<preserved_num<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  destroy_num     "<<destroy_num<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  popsizemin      "<<popsizemin<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  scoref          "<<ga_scoref<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  flushf          "<<ga_flushf<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  score_type      "<<score_type<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  seed            "<<seed<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  pmut            "<<ga_pmut<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  pcross          "<<ga_pcross<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  dev             "<<gauss_dev<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  tol             "<<gtol<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  deltar          "<<deltar<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  ga              "<<ga<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  runtest         "<<runtest<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  read_spline     "<<read_spline<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "  grid_update     "<<grid_update<<endl;
  cout<<"$end\n";

  cout<<"\n$element_types:# element name, eatom" << endl;
  for (int ieat=0;ieat<veatom.size();ieat++)
    cout<<"  "<<fixed
        <<left<<setw(10)<< velem[ieat]<< " " 
        <<right<<setw(10)<<setprecision(6)<<veatom[ieat] << endl;
  cout<<"$end\n";
 
  cout<<"\n$potentials:# potentials to optimize, max_expand, file of knot-vector, max step size, order of spline, constrain level, aa" << endl;
  size0=size1=size2=size3=0;
  for (int ipot=0; ipot<vpot.size(); ipot++ ){
    if(size0<vpot[ipot].potname.size()) size0=vpot[ipot].potname.size()+5;
    if(size1<vpot[ipot].gridname.size()) size1=vpot[ipot].gridname.size()+5;
  }
  for (int ipot=0; ipot<vpot.size(); ipot++ ) 
    cout<<"  "<<fixed
        <<left<<setw(size0)<<vpot[ipot].potname<< " " 
        <<right<<setw(10)<<setprecision(3)<<vpot[ipot].max_expand*AA_Bohr<<" " 
        <<right<<setw(size1)<<setprecision(3)<<vpot[ipot].gridname<<" " 
        <<right<<setw(10)<<setprecision(3)<<vpot[ipot].max_step*AA_Bohr<<" " 
        <<right<<setw(5)<<setprecision(3)<<vpot[ipot].ordspl<<" "
        <<right<<setw(5)<<setprecision(3)<<vpot[ipot].smooth_order<<" "
        <<right<<setw(5)<<setprecision(3)<<vpot[ipot].aa << endl;
  cout<<"$end\n";

  cout<<"\n$compounds:# geometries, Ediss(kcal/mol), eweight, fweight, dftbinp, forceinput" << endl;
  size0=size1=size2=size3=0;
  for (int imol=0; imol<vmol.size(); imol++ ){
    if(size0<vmol[imol].name.size()) size0=vmol[imol].name.size()+5;
    if(size1<vmol[imol].dftbin.size()) size1=vmol[imol].dftbin.size()+5;
  }
  for (int imol=0; imol<vmol.size(); imol++ )
    cout<<"  "
        <<fixed<<left<<setw(size0) << vmol[imol].name<<" "
        <<setprecision(3)<<right<<setw(10)<<-vmol[imol].ebind*kcal_H<<" "
        <<setprecision(2)<<right<<setw(5)<< vmol[imol].eweight<<" " 
        <<setprecision(2)<<right<<setw(5)<< vmol[imol].fweight<<" "
        <<setprecision(2)<<left<<setw(size1)<< vmol[imol].dftbin<<" "
        << vmol[imol].finp << endl;
  cout<<"$end\n";

 cout<<"\n$definition_reactions:# abbreviation, filename, dftbinp" << endl;
  size0=size1=size2=size3=0;
  for (int imol=0; imol<vreamol.size(); imol++){
    if(size0<vreamol[imol].abbr.size()) size0=vreamol[imol].abbr.size()+5;
    if(size1<vreamol[imol].name.size()) size1=vreamol[imol].name.size()+5;
    if(size2<vreamol[imol].dftbin.size()) size2=vreamol[imol].dftbin.size()+5;
  }
  for (int imol=0; imol<vreamol.size(); imol++)
    cout<<"  "
        <<fixed<<left<<setw(size0)<<vreamol[imol].abbr<<" "
        <<left<<setw(size1)<<vreamol[imol].name<<" "
        <<left<<setw(size2)<<vreamol[imol].dftbin << endl;
  cout<<"$end\n";

  cout<<"\n$reactions:#  reaction, \" -> \", Reactionenergy [kcal/mol], reaweight" << endl;
  for (int irea=0; irea<vrea.size(); irea++)
    cout<<""<<vrea[irea].reastr << endl;
  cout<<"$end\n";

  cout << "_____________________________________________________________________________________" << endl << endl;
  if ( vexclpot.size() > 0 ) {
    cout << "WARNING: The following potentials were not asked to optimize in the input! " << endl;
    for (int i=0; i<vexclpot.size(); i++){
      cout<< "  " <<setw(6)<< vexclpot[i].potname 
          << "  " <<setw(6)<< vexclpot[i].ndist 
          << "  " <<setw(6)<< vexclpot[i].mindist*AA_Bohr 
          << "  " << vexclpot[i].molname << endl;
    }
    cout << endl << endl;
  }
     
  if ( nzl > 0 ) {
    cout << "WARNING: There are not enough conditions to solve the equation system." << endl;
    cout << "  "<< nzl << " zero-lines added" << endl;
  }
  cout << endl << endl;

  cout << "knot-vector for each potential (Angstrom)" << endl << endl;
  for (int ipot=0; ipot<vpot.size(); ipot++ ) {
    cout << vpot[ipot].potname << "     ";
    for (int ir=0; ir<vpot[ipot].nknots; ir++)  
      cout << vpot[ipot].vr[ir]*AA_Bohr << " ";
    cout << endl; 
  }
  cout << endl << endl; 

  cout << "_____________________________________________________________________________________" << endl << endl;
  cout << "Distances in Angstrom between all atoms in all molecules appearing in the input (duplicates are removed)" 
       << endl << endl << fixed;
  for (int ipot=0;ipot<vpot.size();ipot++){
    for (int idist=0;idist<vpot[ipot].vmoldist.size();idist++){
      cout<<fixed<< vpot[ipot].potname<<"  "<<setprecision(3)
          <<setw(10)<<vpot[ipot].vmoldist[idist].dist*AA_Bohr<<" "
          <<vpot[ipot].vmoldist[idist].molname << endl;
    }
    cout << endl << endl; 
  }

  cout << "_____________________________________________________________________________________" << endl << endl;
  cout << "equation-matrix:"  << endl;
  cout << "number of rows: " << nrows << endl;
  cout << "                                      number of spline-equations:     " << setw(10) << nspleq << endl;
  cout << "                                      number of energy-equations:     " << setw(10) << neeq   << endl;
  cout << "                                      number of force-equations:      " << setw(10) << nfeq   << endl;
  cout << "                                      number of reaction-equations:   " << setw(10) << nreaeq << endl;
  cout << "                                      number of zero-lines:           " << setw(10) << nzl    << endl;
  cout << "number of cols: " << ncols << endl;
  cout << "                                      number of spline-coefficients:  " << setw(10) << nallequations << endl;
  cout << "                                      number of elements:             " << setw(10) << nelem     << endl;
  cout << "                                      number of eatom to be fitted:   " << setw(10) << natfit    << endl;
  cout << endl << endl;

  cout << "_____________________________________________________________________________________" << endl << endl;
  if (!runtest){
    int ires=0;
    cout << "residual:" << endl << endl;

    cout << "    component          MSE       MUE       RMS        in kcal/mol:  MSE    MUE    RMS" << endl << endl << fixed;
    if (neeq>0)      cout << setprecision(6)<<"    energies        " 
                          << setw(10) << reseS/neeq
                          << setw(10) << reseU/neeq     
                          << setw(10) << sqrt(rese2/neeq) 
                          << "              "<< setprecision(1) 
                          << setw(7) << reseS/neeq*kcal_H     
                          << setw(7) << reseU/neeq*kcal_H     
                          << setw(7) << sqrt(rese2/neeq)*kcal_H << endl;
    if (nfeq>0)      cout << setprecision(6)<< "    forces          " 
                          << setw(10) << resfS/nfeq            
                          << setw(10) << resfU/nfeq     
                          << setw(10) << sqrt(resf2/nfeq) << endl;
    if (nreaeq>0)    cout << setprecision(6)<< "    reactions       " 
                          << setw(10) << resreaS/nreaeq        
                          << setw(10) << resreaU/nreaeq 
                          << setw(10) << sqrt(resrea2/nreaeq) 
                          << "              "<< setprecision(1) 
                          << setw(7) << resreaS/nreaeq*kcal_H 
                          << setw(7) << resreaU/nreaeq*kcal_H 
                          << setw(7) << sqrt(resrea2/nreaeq)*kcal_H << endl;
    cout << "    ----------------" << endl;
    cout << "    total           " << setprecision(6) 
         << setw(10) << restotS/(nrows-nspleq) 
         << setw(10) << restotU/(nrows-nspleq) 
         << setw(10) << sqrt(restot2/(nrows-nspleq)) << endl << endl;

    cout << "_____________________________________________________________________________________" << endl << endl;
    ires=nspleq;
    if (neeq>0) cout << "energy equations (residual in kcal/mol)" << endl << fixed << setprecision(1);
    for (int i=0; i<vmol.size(); i++) {
      if (vmol[i].eweight!=0){
        cout<<fixed
            <<right<<setw(10)<<setprecision(3)<<-vres[ires]*kcal_H<<" " 
            <<left<<vmol[i].name << endl; 
        ires++;
      }
    } 
   
    int iconstraineda;
    bool constrained;

    cout << "_____________________________________________________________________________________" << endl << endl;
    if (nfeq>0) cout << endl << "force equations (residual in atomic units)" << endl << fixed << setprecision(6);
    for (int i=0; i<vmol.size(); i++) {
      if (vmol[i].fweight!=0){
        cout << vmol[i].name << endl; 
        for (int j=0; j<vmol[i].natom; j++){
          if(vmol[i].ncontraineda>0){
            constrained=false; 
            for(iconstraineda=0;iconstraineda<vmol[i].ncontraineda;iconstraineda++){
              if(j==vmol[i].contraineda[iconstraineda]) constrained=true;
            }
            if(!constrained) continue;
          }
          cout<<fixed<<left<<setw(3)<< vmol[i].atomname[j] <<" "
              <<right<<setprecision(6)<<setw(12)<<vres[ires] 
              <<right<<setprecision(6)<<setw(12)<<vres[ires+1] 
              <<right<<setprecision(6)<<setw(12)<<vres[ires+2] << endl; 
          ires+=3;
        }
      }
    }
    if (nreaeq>0) cout << endl << "reaction equations (residual in kcal/mol)" << endl << fixed << setprecision(1);
    for (int i=0; i<nreaeq; i++){
      cout<<fixed<<right<<setw(10)<<setprecision(3)<< vres[ires]*kcal_H << " " 
          << vrea[i].reastr << endl; 
      ires++;
    }

    cout << "_____________________________________________________________________________________" << endl << endl;

    if (natfit>0) cout << "Atom  Energy(Hartree) " << endl;
    else          cout << "Atomic energies were not fitted." << endl;
    for (int i=0; i<natfit; i++){
      cout << velem[nelem-natfit+i] << "    " << setprecision(12) << vunknown[vunknown.size()-natfit+i] << endl;
    }

    cout << "_____________________________________________________________________________________" << endl << endl;
    
    if(ddh.td3){
      cout<<"$d3:\n"; 
      for(int i=0;i<ddh.d3.size();i++){ 
        cout.precision(ddh.d3[i].precision) ;
        cout<<ddh.d3[i].name<<" "<<ddh.d3[i].min<<" "<<ddh.d3[i].value<<" "<<ddh.d3[i].max<<" "<<ddh.d3[i].precision<<endl;
      }
      cout<<"$end\n"; 
    }
    if(ddh.tdamph){
      cout<<"$damph:\n"; 
      cout.precision(ddh.damph.precision);
        cout<<ddh.damph.min<<" "<<ddh.damph.value<<" "<<ddh.damph.max<<" "<<ddh.damph.precision<<endl;
      cout<<"$end\n"; 
    }
    if(ddh.thubbardderivs){
      cout<<"$hubbardderivs:\n"; 
      for(int i=0;i<ddh.hubbardderivs.size();i++){ 
        cout.precision(ddh.hubbardderivs[i].precision);
        cout<<ddh.hubbardderivs[i].name<<" "<<ddh.hubbardderivs[i].min<<" "<<ddh.hubbardderivs[i].value<<" "<<ddh.hubbardderivs[i].max<<" "<<ddh.hubbardderivs[i].precision<<endl;
      }
      cout<<"$end\n"; 
    }
    if(ddh.tdamphver2){
      cout<<"$damphver2:\n"; 
      for(int i=0;i<ddh.damphver2.size();i++){ 
        cout.precision(ddh.damphver2[i].precision);
        cout<<ddh.damphver2[i].name<<" "<<ddh.damphver2[i].min<<" "<<ddh.damphver2[i].value<<" "<<ddh.damphver2[i].max<<" "<<ddh.damphver2[i].precision<<endl;
      }
      cout<<"$end\n"; 
    }

    cout << "_____________________________________________________________________________________" << endl << endl;

    int icoeff=0;
    cout << "calculated potentials: " << endl << endl;
    //for each potential
    for (int i=0; i < vpot.size(); i++){
      // header and first three lines
      cout << vpot[i].potname  << "   r1,r2,Allequations ( " << vpot[i].ordspl << " th order splines)"  << endl;
      cout << "Spline" << endl;
      cout << setiosflags(ios::showpoint) << resetiosflags(ios::scientific);
      cout << vpot[i].nknots-1 << "  " << setprecision(6) << vpot[i].vr[vpot[i].nknots-1] << endl;
      cout << setprecision(15) << scientific << vpot[i].expA << setw(23) << vpot[i].expB << setw(24) << vpot[i].expC << endl;
      // r1, r2, a, b, c, d, e
      for (int j=0; j<vpot[i].nknots-1; j++){
        cout << setiosflags(ios::showpoint) << resetiosflags(ios::scientific);
        cout << setprecision(6) << setw(8) << vpot[i].vr[j] << "  " << setw(8) << vpot[i].vr[j+1] << "   ";
        for (int k=0; k<=vpot[i].ordspl; k++){
          cout << setprecision(15) << scientific << setw(23) << vunknown[icoeff] << " ";
          icoeff++;
        }
        cout << endl;
      }
      cout << endl;
    }

    cout << endl << endl << "** repopt normal termination **" << endl << endl;
    //ofstream fdist("dist.repopt");
    //for (int i=0; i < vmol.size(); i++){
    //  for (int iat1=0; iat1 < vmol[i].natom -1; iat1++){
    //    for (int iat2=iat1+1; iat2 < vmol[i].natom; iat2++){
    //      fdist << fixed << setprecision(3) << vmol[i].dist(iat1,iat2) << "  " << vmol[i].dist(iat1,iat2)*AA_Bohr
    //            << "  " << vmol[i].atomname[iat1] << vmol[i].atomname[iat2] << "  "
    //            << vmol[i].name << endl;
    //    }
    //  }
    //}
    //fdist.close();
  }else{
        //cout<<fixed<<left<<setw(3)<< vmol[i].atomname[j] <<" "
        //    <<right<<setprecision(6)<<setw(12)<<vres[ires] 
        //    <<right<<setprecision(6)<<setw(12)<<vres[ires+1] 
        //    <<right<<setprecision(6)<<setw(12)<<vres[ires+2] << endl; 

    double mse,mad,rmse,max,size,tmp; 
    double msre,mard,rmsre,maxr,tmpr; 
    if (neeq>0){ 
      cout << fixed << "\nEnergy equations (residual in kcal/mol)" << endl;
      cout << fixed <<right<< setw(10) << "Ref" <<right<< setw(10) << "Cal" <<right<< setw(10) << "Error" << endl;
      cout << fixed << setprecision(2)<<right;
      mse=mad=rmse=max=size=msre=mard=rmsre=maxr=tmpr=0.0;
      for (int imol=0; imol < allequations.vmol.size(); imol++){
        if (allequations.vmol[imol].eweight!=0){
          tmp=(allequations.vmol[imol].ebind-allequations.vmol[imol].teel)*kcal_H;
          tmpr=-(allequations.vmol[imol].ebind-allequations.vmol[imol].teel)/allequations.vmol[imol].ebind;
          cout<<fixed<<right<<setprecision(2)<<setw(10)<<-allequations.vmol[imol].ebind*kcal_H
                     <<right<<setprecision(2)<<setw(10)<<-allequations.vmol[imol].teel*kcal_H
                     <<right<<setprecision(2)<<setw(10)<<tmp<<"  "
                     <<allequations.vmol[imol].name <<endl; 
          mse+=tmp;
          rmse+=tmp*tmp;
          mad+=abs(tmp);
          if(abs(tmp)>max) max=abs(tmp);

          msre+=tmpr;
          rmsre+=tmpr*tmpr;
          mard+=abs(tmpr);
          if(abs(tmpr)>maxr) maxr=abs(tmpr);

          size+=1.0;
        }
      }
      mse/=size;
      rmse/=size;
      mad/=size;  
      rmse=sqrt(rmse);

      msre/=size;
      rmsre/=size;
      mard/=size;  
      rmsre=sqrt(rmsre);

      msre*=100.0;
      rmsre*=100.0;
      mard*=100.0;  
      maxr*=100.0;  

      cout<<"#\n";
      cout<<fixed<<left<<setw(10)<<"#"<<right<<setw(10)<<"MSE"<<setw(10)<<"MAD"<<setw(10)<<"RMSE"<<setw(10)<<"MAX"<<"\n";
      cout<<fixed<<left<<setw(10)<<"#"<<right<<setprecision(2)<<setw(10)<<mse<<setprecision(2)<<setw(10)<<mad<<setprecision(2)<<setw(10)<<rmse<<setprecision(2)<<setw(10)<<max<<"\n";
      cout<<"#\n";

      cout<<"#\n";
      cout<<fixed<<left<<setw(10)<<"#"<<right<<setw(10)<<"MSRE"<<setw(10)<<"MARD"<<setw(10)<<"RMSRE"<<setw(10)<<"MAXR"<<"\n";
      cout<<fixed<<left<<setw(10)<<"#"<<right<<setprecision(2)<<setw(10)<<msre<<setprecision(2)<<setw(10)<<mard<<setprecision(2)<<setw(10)<<rmsre<<setprecision(2)<<setw(10)<<maxr<<"\n";
      cout<<"#\n";
    }
    if(nreaeq>0){
      cout << fixed << "\nReaction equations (residual in kcal/mol)" << endl;
      cout << fixed <<right<< setw(10) << "Ref" <<right<< setw(10) << "Cal" <<right<< setw(10) << "Error" << endl;
      cout << fixed << setprecision(2)<<right;
      mse=mad=rmse=max=size=msre=mard=rmsre=maxr=tmpr=0.0;
      for (int irea=0; irea < allequations.vrea.size(); irea++ ){
        if (allequations.vrea[irea].reaweight!=0){
          tmp=(allequations.vrea[irea].treaE-allequations.vrea[irea].reaE)*kcal_H;
          tmpr=(allequations.vrea[irea].treaE-allequations.vrea[irea].reaE)/allequations.vrea[irea].reaE;
          cout<<fixed<<right<<setprecision(2)<<setw(10)<<allequations.vrea[irea].reaE*kcal_H
                     <<right<<setprecision(2)<<setw(10)<<allequations.vrea[irea].treaE*kcal_H
                     <<right<<setprecision(2)<<setw(10)<<tmp<<"  "
                     <<allequations.vrea[irea].reastr<<endl; 
          mse+=tmp;
          rmse+=tmp*tmp;
          mad+=abs(tmp);
          if(abs(tmp)>max) max=abs(tmp);

          msre+=tmpr;
          rmsre+=tmpr*tmpr;
          mard+=abs(tmpr);
          if(abs(tmpr)>maxr) maxr=abs(tmpr);

          size+=1.0;
        }
      }
      mse/=size;
      rmse/=size;
      mad/=size;  
      rmse=sqrt(rmse);

      msre/=size;
      rmsre/=size;
      mard/=size;  
      rmsre=sqrt(rmsre);

      msre*=100.0;
      rmsre*=100.0;
      mard*=100.0;  
      maxr*=100.0;  

      cout<<"#\n";
      cout<<fixed<<setprecision(2);
      cout<<fixed<<left<<setw(10)<<"#"<<right<<setw(10)<<"MSE"<<setw(10)<<"MAD"<<setw(10)<<"RMSE"<<setw(10)<<"MAX"<<"\n";
      cout<<fixed<<left<<setw(10)<<"#"<<right<<setprecision(2)<<setw(10)<<mse<<setprecision(2)<<setw(10)<<mad<<setprecision(2)<<setw(10)<<rmse<<setprecision(2)<<setw(10)<<max<<"\n";
      cout<<"#\n";

      cout<<"#\n";
      cout<<fixed<<left<<setw(10)<<"#"<<right<<setw(10)<<"MSRE"<<setw(10)<<"MARD"<<setw(10)<<"RMSRE"<<setw(10)<<"MAXR"<<"\n";
      cout<<fixed<<left<<setw(10)<<"#"<<right<<setprecision(2)<<setw(10)<<msre<<setprecision(2)<<setw(10)<<mard<<setprecision(2)<<setw(10)<<rmsre<<setprecision(2)<<setw(10)<<maxr<<"\n";
      cout<<"#\n";
    }
  }
}


#include <iostream>
#include <string>
#include <iomanip>
#include "auxiliary.hpp"
#include "allequations.hpp"

using namespace std;

extern Allequations allequations;
extern int ga_popsize, ga_ngen, ga_scoref, ga_flushf;
extern int score_type,preserved_num, seed;
extern double gtol;
extern double ga_pmut, ga_pcross, gauss_dev;
extern double deltar;
extern bool ga,read_spline,fsmooth_spline,aa_spline,fderivative1st,fderivative2nd,fderivative3rd;
extern bool runtest;
extern sddh ddh;

void Allequations::writeout() const {
  cout<<std::fixed;

  cout << "** repopt **" << endl << endl;

  if ( vexclpot.size() > 0 ) {
    cout << "WARNING: The following potentials were not asked to optimize in the input! " << endl;
    cout << "         Repulsive energy and forces resulting from these potentials are calculated from potentials" << endl;
    cout << "         which are defined in the dftbinp for each certain molecule and added to the corresponding electronic energy/force." << endl;
    cout << "         potential, number of tried access (energy and force eq.), minimum distance (A, a.u.) in molecule" << endl << endl;
    for (int i=0; i<vexclpot.size(); i++){
      cout << "         " << vexclpot[i].potname << "  " << vexclpot[i].ndist << "  " << vexclpot[i].mindist*AA_Bohr 
           << "  " << vexclpot[i].mindist << "  " << vexclpot[i].molname << endl;
    }
    cout << endl << endl;
  }
     
  if ( nzl > 0 ) {
    cout << "WARNING: There are not enough conditions to solve the equation system." << endl;
    cout << "         "<< nzl << " zero-lines added"                                << endl;
  }
  cout << endl << endl;

  cout<<fixed<<setprecision(6)<<left<<setw(30) << "###interpreted input###"      << endl << endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "dftbversion "  << dftbversion << endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "popsizemax "   <<ga_popsize   <<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "ngen "         <<ga_ngen      <<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "preserved_num "<<preserved_num<<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "scoref "       <<ga_scoref    <<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "flushf "       <<ga_flushf    <<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "score_type "   <<score_type   <<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "seed "         <<seed         <<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "pmut "         <<ga_pmut      <<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "pcross "       <<ga_pcross    <<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "gauss_dev "    <<gauss_dev    <<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "tol "          <<gtol         <<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "deltar "       <<deltar       <<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "ga "           <<ga           <<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "runtest "      <<runtest      <<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << "read_spline "  <<read_spline  <<endl;
  cout<<fixed<<setprecision(6)<<left<<setw(30) << endl << endl;

  cout << "atom_types:# atom type, eatom" << endl;
  for (int ieat=0;ieat<veatom.size();ieat++) 
    cout<<fixed<<left<<setw(10)<< velem[ieat]<< " "
        <<right<<setw(10)<<setprecision(6)<<veatom[ieat] << endl;
  cout << endl << endl;
 
  cout << "potentials:# potentials to optimize, max_expand, file of knot-vector, order of spline, constrain, aa" << endl;
  for (int ipot=0; ipot<vpot.size(); ipot++ ) 
    cout<<fixed<<left<<setw(10)<<vpot[ipot].potname<< " " 
        <<right<<setw(10)<<setprecision(3)<<vpot[ipot].max_expand*AA_Bohr<<" " 
        <<right<<setw(20)<<setprecision(3)<<vpot[ipot].gridname<<" " 
        <<right<<setw(5)<<setprecision(3)<<vpot[ipot].ordspl<<" "
        <<right<<setw(5)<<setprecision(3)<<vpot[ipot].smooth_order<<" "
        <<right<<setw(5)<<setprecision(3)<<vpot[ipot].aa << endl;
  cout << endl << endl;

  cout << "compounds:# geometries, Ediss(kcal/mol), eweight, fweight, dftbinp, forceinput" << endl;
  for (int imol=0; imol<vmol.size(); imol++ )
    cout<<fixed<<left<<setw(60) << vmol[imol].name<<" "
        <<setprecision(3)<<right<<setw(10)<<-vmol[imol].ebind*kcal_H<<" "
        <<setprecision(2)<<right<<setw(5)<< vmol[imol].eweight<<" " 
        <<setprecision(2)<<right<<setw(5)<< vmol[imol].fweight<<" "
        <<setprecision(2)<<left<<setw(30)<< vmol[imol].dftbin<<" "
        << vmol[imol].finp << endl;
  cout << endl << endl;

  cout << "definition_reactions:# abbreviation, filename, dftbinp" << endl;
  for (int imol=0; imol<vreamol.size(); imol++)
    cout<<fixed<<left<<setw(30)<<vreamol[imol].abbr<<" "
        <<left<<setw(60)<<vreamol[imol].name<<" "
        <<left<<vreamol[imol].dftbin << endl;
  cout << endl << endl;

  cout << "reactions:#  reaction, \" -> \", Reactionenergy [kcal/mol], reaweight" << endl;
  for (int irea=0; irea<vrea.size(); irea++)
    cout << vrea[irea].reastr << endl;
  cout << endl << endl;

  cout << "_____________________________________________________________________________________" << endl << endl;

  cout << "knot-vector for each potential (au)" << endl << endl;
  for (int ipot=0; ipot<vpot.size(); ipot++ ) {
    cout << vpot[ipot].potname << "     ";
    for (int ir=0; ir<vpot[ipot].nknots; ir++)  
      cout << vpot[ipot].vr[ir] << " ";
    cout << endl; 
  }
  cout << endl << endl; 

  cout << "_____________________________________________________________________________________" << endl << endl;

  cout << "Distances in a.u. and AA between all atoms in all molecules appearing in the input (duplicates are removed)" 
       << endl << endl << fixed;
  for (int ipot=0;ipot<vpot.size();ipot++){
    for (int idist=0;idist<vpot[ipot].vmoldist.size();idist++){
      cout<<fixed<< vpot[ipot].potname<<" "
          <<setprecision(3)<<setw(10)<<vpot[ipot].vmoldist[idist].dist<<" " 
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
    if (neeq>0)      cout << setprecision(6) << "    energies        " << setw(10) 
                          << reseS/neeq     << setw(10) << reseU/neeq     << setw(10) << sqrt(rese2/neeq) 
                          << "              " << setprecision(1) << setw(7) 
                          << reseS/neeq*kcal_H     << setw(7) << reseU/neeq*kcal_H     << setw(7) << sqrt(rese2/neeq)*kcal_H << endl;
    if (nfeq>0)      cout << setprecision(6) << "    forces          " << setw(10)
                          << resfS/nfeq     << setw(10) << resfU/nfeq     << setw(10) << sqrt(resf2/nfeq) << endl;
    if (nreaeq>0)    cout << setprecision(6) << "    reactions       " << setw(10)
                          << resreaS/nreaeq << setw(10) << resreaU/nreaeq << setw(10) << sqrt(resrea2/nreaeq) 
                          << "              " << setprecision(1) << setw(7)
                          << resreaS/nreaeq*kcal_H << setw(7) << resreaU/nreaeq*kcal_H << setw(7) << sqrt(resrea2/nreaeq)*kcal_H << endl;
    cout << "    ----------------" << endl;
    cout << "    total           " << setprecision(6) << setw(10) 
         << restotS/(nrows-nspleq) << setw(10) << restotU/(nrows-nspleq) << setw(10) << sqrt(restot2/(nrows-nspleq)) << endl << endl << endl;

    ires=nspleq;

    if (neeq>0) cout << "energy equations (residual in kcal/mol)" << endl << fixed << setprecision(1);
    for (int i=0; i<vmol.size(); i++) {
      if (vmol[i].eweight!=0){
        cout<<fixed<<right<<setw(10)<<setprecision(3)<<vres[ires]*kcal_H<<" " 
            <<left<<vmol[i].name << endl; 
        ires++;
      }
    } 
   
    int iconstraineda;
    bool constrained;

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

    if ( vexclpot.size() > 0 ) {
      cout << "WARNING: The following potentials were not asked to optimize in the input! " << endl;
      cout << "         Repulsive energy and forces resulting from these potentials are calculated from potentials" << endl;
      cout << "         which are defined in the dftbinp for each certain molecule and added to the corresponding electronic energy/force." << endl;
      cout << endl << endl;
      cout << "         potential, number of tried access (energy and force eq.), minimum distance (A, a.u.) in molecule" << endl << endl;
      for (int i=0; i<vexclpot.size(); i++){
      cout << "         " << vexclpot[i].potname << "  " << vexclpot[i].ndist << "  " << vexclpot[i].mindist*AA_Bohr 
           << "  " << vexclpot[i].mindist << "  " << vexclpot[i].molname << endl;
      }
      cout << endl << endl;
    }
       
    if ( nzl > 0 ) {
      cout << "WARNING: There are not enough conditions to solve the equation system." << endl;
      cout << "         " << nzl << " zero-lines added"                        << endl << endl;
    }

    //cout << endl << endl << "** repopt normal termination **" << endl << endl;

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

    double mse,mad,max,size,tmp; 
    if (neeq>0) cout << "energy equations (residual in kcal/mol)" << endl << fixed << setprecision(2)<<right;
    mse=mad=max=size=0.0;
    for (int imol=0; imol < allequations.vmol.size(); imol++){
      if (allequations.vmol[imol].eweight!=0){
        tmp=(allequations.vmol[imol].ebind-allequations.vmol[imol].teel)*kcal_H;
        //cout<<allequations.vmol[imol].teel<<endl;
        cout<< setw(10)<<-allequations.vmol[imol].ebind*kcal_H<<setw(10)<<-allequations.vmol[imol].teel*kcal_H<<setw(10)<<tmp<<"  "<< allequations.vmol[imol].name <<endl; 
        mse+=tmp;
        mad+=abs(tmp);
        if(abs(tmp)>max) max=abs(tmp);
        size+=1.0;
      }
    }
    mse/=size;
    mad/=size;  
    cout<<"#\n";
    cout<<left<<setw(10)<<"#"<<right<<setw(10)<<"MSE"<<setw(10)<<"MAD"<<setw(10)<<"MAX"<<"\n";
    cout<<left<<setw(10)<<"#"<<right<<setw(10)<<mse<<setw(10)<<mad<<setw(10)<<max<<"\n";
    cout<<"#\n";
    if (nreaeq>0) cout << endl << "reaction equations (residual in kcal/mol)" << endl << fixed << setprecision(1);
    mse=mad=max=size=0.0;
    for (int irea=0; irea < allequations.vrea.size(); irea++ ){
      if (allequations.vrea[irea].reaweight!=0){
        tmp=(allequations.vrea[irea].treaE-allequations.vrea[irea].reaE)*kcal_H;
        cout<< setw(10)<<allequations.vrea[irea].reaE*kcal_H<<setw(10)<<allequations.vrea[irea].treaE*kcal_H<<setw(10)<<tmp<<"  "<<allequations.vrea[irea].reastr<<endl; 
        mse+=tmp;
        mad+=abs(tmp);
        if(abs(tmp)>max) max=abs(tmp);
        size+=1.0;
      }
    }
    mse/=size;
    mad/=size;  
    cout<<"#\n";
    cout<<left<<setw(10)<<"#"<<right<<setw(10)<<"MSE"<<setw(10)<<"MAD"<<setw(10)<<"MAX"<<"\n";
    cout<<left<<setw(10)<<"#"<<right<<setw(10)<<mse<<setw(10)<<mad<<setw(10)<<max<<"\n";
    cout<<"#\n";
  }

}


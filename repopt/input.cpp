#include <iostream>
#include <string>
#include <iomanip>
#include "tools.hpp"
#include "auxiliary.hpp"
#include "allequations.hpp"

using namespace std;
extern double gtol;
extern double min_step01;
extern double expandR,deltar;
extern double ga_pmut,ga_pcross,gauss_dev;
extern double s1,s2,s3,s4,s5,s6,s7,s8,s9; 
extern int    ilmsfit,idecompose,nreplicate,dftbout_new;
extern int    ga_popsize,ga_ngen,ga_scoref,ga_flushf;
extern int    score_type,preserved_num,destroy_num,popsizemin,seed;
extern bool   ga,read_spline,fsmooth_spline,aa_spline;
extern bool   fderivative1st,fderivative2nd,fderivative3rd;
extern bool   runtest,grid_update;
extern bool   runxtb,runmopac;
extern string xtbarg,xtbarg2,mopacarg;
extern sddh   ddh;

void Allequations::readinp(const string inputfile){
  ifstream infile;
  string str,dftbin,potname,filename,finp,stmp;
  int tsmooth,npot,nadd,nmol,nabbr,nrea,derivsmooth=0,neighsmooth=1;
  int ordspl=0;
  double eweight,fweight;
  double ediss=0,eatom=0,psmooth=0;
  string stemp,stemp1;
  char cline[512],ctemp0[512],ctemp[512],ctemp1[512],ctemp2[512],ctemp3[512],ctemp4[512];
  int i,j,k,itmp,ntmp;
  double fmin,ftmp,fmax;
  ifnotfile(inputfile.c_str());
  infile.open(inputfile.c_str(),ios::in);
 
  if(!infile.is_open()){
    cout<<"unable to open "<<inputfile<<" file.\n";
    exit(1);
  }else{
    ctemp[0]=ctemp0[0]=ctemp1[0]=ctemp2[0]=ctemp3[0]=ctemp4[0]='\0';
    while(infile.getline(cline,512)){
      remove_comment(cline);
      stemp=cline;
      if(stemp.find("$system:")!=string::npos){
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          if(stemp.find("dftb_version ")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); dftbversion=ctemp2;}
          if(stemp.find("dftbout_new ")!=string::npos) sscanf(cline,"%s %d",ctemp1,&dftbout_new);
          if(stemp.find("ilmsfit ")!=string::npos) sscanf(cline,"%s %d",ctemp1,&ilmsfit);
          if(stemp.find("idecompose ")!=string::npos) sscanf(cline,"%s %d",ctemp1,&idecompose);
          if(stemp.find("nreplicate ")!=string::npos) sscanf(cline,"%s %d",ctemp1,&nreplicate);
          if(stemp.find("xtbarg ")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); xtbarg=ctemp2; runxtb=true;}
          if(stemp.find("xtbarg2 ")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); xtbarg2=ctemp2; runxtb=true;}
          if(stemp.find("mopacarg ")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); mopacarg=ctemp2; runmopac=true;}
        }
      }else if(stemp.find("$variables:")!=string::npos){
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          if(stemp.find("s1 ")!=string::npos) sscanf(cline,"%s %le",ctemp1,&s1);
          if(stemp.find("s2 ")!=string::npos) sscanf(cline,"%s %le",ctemp1,&s2);
          if(stemp.find("s3 ")!=string::npos) sscanf(cline,"%s %le",ctemp1,&s3);
          if(stemp.find("s4 ")!=string::npos) sscanf(cline,"%s %le",ctemp1,&s4);
          if(stemp.find("s5 ")!=string::npos) sscanf(cline,"%s %le",ctemp1,&s5);
          if(stemp.find("s6 ")!=string::npos) sscanf(cline,"%s %le",ctemp1,&s6);
          if(stemp.find("s7 ")!=string::npos) sscanf(cline,"%s %le",ctemp1,&s7);
          if(stemp.find("s8 ")!=string::npos) sscanf(cline,"%s %le",ctemp1,&s8);
          if(stemp.find("s9 ")!=string::npos) sscanf(cline,"%s %le",ctemp1,&s9);
          if(stemp.find("min_step01 ")!=string::npos) sscanf(cline,"%s %le",ctemp1,&min_step01);
        }
      }else if(stemp.find("$genetic_algorithm:")!=string::npos){
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          if(stemp.find("popsizemax ")!=string::npos) sscanf(cline,"%s %d",ctemp1,&ga_popsize);
          if(stemp.find("ngen ")!=string::npos) sscanf(cline,"%s %d",ctemp1,&ga_ngen);
          if(stemp.find("preserved_num ")!=string::npos) sscanf(cline,"%s %d",ctemp1,&preserved_num);
          if(stemp.find("destroy_num ")!=string::npos) sscanf(cline,"%s %d",ctemp1,&destroy_num);
          if(stemp.find("popsizemin ")!=string::npos) sscanf(cline,"%s %d",ctemp1,&popsizemin);
          if(stemp.find("scoref ")!=string::npos) sscanf(cline,"%s %d",ctemp1,&ga_scoref);
          if(stemp.find("flushf ")!=string::npos) sscanf(cline,"%s %d",ctemp1,&ga_flushf);
          if(stemp.find("score_type ")!=string::npos) sscanf(cline,"%s %d",ctemp1,&score_type);
          if(stemp.find("seed ")!=string::npos) sscanf(cline,"%s %d",ctemp1,&seed);
          if(stemp.find("pmut ")!=string::npos) sscanf(cline,"%s %le",ctemp1,&ga_pmut);
          if(stemp.find("pcross ")!=string::npos) sscanf(cline,"%s %le",ctemp1,&ga_pcross);
          if(stemp.find("dev ")!=string::npos) sscanf(cline,"%s %le",ctemp1,&gauss_dev);
          if(stemp.find("tol ")!=string::npos) sscanf(cline,"%s %le",ctemp1,&gtol);
          if(stemp.find("deltar ")!=string::npos) sscanf(cline,"%s %le",ctemp1,&deltar);
          if(stemp.find("ga ")!=string::npos){ sscanf(cline,"%s %d",ctemp1,&itmp); ga = itmp;}
          if(stemp.find("runtest ")!=string::npos){ sscanf(cline,"%s %d",ctemp1,&itmp); runtest = itmp;}
          if(stemp.find("read_spline ")!=string::npos){ sscanf(cline,"%s %d",ctemp1,&itmp); read_spline = itmp;}
          if(stemp.find("grid_update ")!=string::npos){ sscanf(cline,"%s %d",ctemp1,&itmp); grid_update = itmp;}
        }
      }else if(stemp.find("$element_types:")!=string::npos){
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          sscanf(cline,"%s %le",ctemp1,&eatom);
          str=ctemp1;
          for (int ielem=0; ielem < velem.size(); ielem++){
            if ( str == velem[ielem] ) {
              cerr << endl << "ERROR: " << endl << endl;
              cerr << "in " << inputfile << ": there are two energies for atom: " << str << endl
                   << endl << "\" " << str << "\" was read in twice. Change this in the input and rerun." 
                   << "exit repopt" << endl << endl;
              exit(1);
            }
          }
          velem.push_back(str);
          veatom.push_back(eatom);
          eatom=0;
        }
      }else if(stemp.find("$d3:")!=string::npos){
        ddh.td3=true;
        ntmp=0;
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#') continue;
          sscanf(cline,"%s %le %le %le %d",ctemp1,&fmin,&ftmp,&fmax,&itmp);
          ddh.d3.resize(ntmp+1);
          ddh.d3[ntmp].name=ctemp1;
          ddh.d3[ntmp].min=fmin;
          ddh.d3[ntmp].value=ftmp;
          ddh.d3[ntmp].max=fmax;
          ddh.d3[ntmp].precision=itmp;
          ntmp++;
        }
      }else if(stemp.find("$damph:")!=string::npos){
        ddh.tdamph=true;
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          sscanf(cline,"%le %le %le %d",&fmin,&ftmp,&fmax,&itmp);
          ddh.damph.min=fmin;
          ddh.damph.value=ftmp;
          ddh.damph.max=fmax;
          ddh.damph.precision=itmp;
        }
      }else if(stemp.find("$hubbardderivs:")!=string::npos){
        ddh.thubbardderivs=true;
        ddh.thirdorderfull=false;
        ntmp=0;
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue ;
          sscanf(cline,"%s %le %le %le %d",ctemp1,&fmin,&ftmp,&fmax,&itmp);
          ddh.hubbardderivs.resize(ntmp+1);
          ddh.hubbardderivs[ntmp].name=ctemp1;
          ddh.hubbardderivs[ntmp].min=fmin;
          ddh.hubbardderivs[ntmp].value=ftmp;
          ddh.hubbardderivs[ntmp].max=fmax;
          ddh.hubbardderivs[ntmp].precision=itmp;
          ntmp++;
        }
      }else if(stemp.find("$hubbardderivsfull:")!=string::npos){
        ddh.thubbardderivs=true;
        ddh.thirdorderfull=true;
        ntmp=0;
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue ;
          sscanf(cline,"%s %le %le %le %d",ctemp1,&fmin,&ftmp,&fmax,&itmp);
          ddh.hubbardderivs.resize(ntmp+1);
          ddh.hubbardderivs[ntmp].name=ctemp1;
          ddh.hubbardderivs[ntmp].min=fmin;
          ddh.hubbardderivs[ntmp].value=ftmp;
          ddh.hubbardderivs[ntmp].max=fmax;
          ddh.hubbardderivs[ntmp].precision=itmp;
          ntmp++;
        }
      }else if(stemp.find("$damphver2:")!=string::npos){
        ddh.tdamphver2=true;
        ntmp=0;
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue ;
          sscanf(cline,"%s %le %le %le %d",ctemp1,&fmin,&ftmp,&fmax,&itmp);
          ddh.damphver2.resize(ntmp+1);
          ddh.damphver2[ntmp].name=ctemp1;
          ddh.damphver2[ntmp].min=fmin;
          ddh.damphver2[ntmp].value=ftmp;
          ddh.damphver2[ntmp].max=fmax;
          ddh.damphver2[ntmp].precision=itmp;
          ntmp++;
        }
      }else if(stemp.find("$potentials:")!=string::npos){
        npot=0;
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          vpot.resize(npot+1);
          sscanf(cline,"%s%le%s%le%d%d%d",ctemp1,&vpot[npot].max_expand,ctemp2,&vpot[npot].max_step,&ordspl,&vpot[npot].smooth_order,&itmp);
          vpot[npot].aa = itmp;
          potname=ctemp1;
          if (potname.length() != 4) {
            cerr << endl << "ERROR: " << endl << endl;
            cerr << "in " << inputfile << ": name of potential has to be four characters long, e.g. \'c_c_\', \'clcl\'" 
                 << endl << "\" " << potname << "\" was read instead. Change this in the input and rerun." 
                 << endl << endl  << "exit repopt" << endl << endl;
            exit(1);
          }
          vpot[npot].max_expand=vpot[npot].max_expand/AA_Bohr;
          vpot[npot].max_step=vpot[npot].max_step/AA_Bohr;
          filename=ctemp2;
          if (filename.substr(0,8)!="autogrid")  ifnotfile(filename);
          if (ordspl == 0) {
            cerr << endl << "ERROR: " << endl << endl;
            cerr << "in " << inputfile << ": order of spline of potential \"" << potname << "\" with corresponding grid-filename \"" 
                 << filename << "\" should be greater than zero!" 
                 << endl << endl << "exit repopt" << endl << endl; 
            exit(1);
          }
          psmooth=0; derivsmooth=0; neighsmooth=1;
          vpot[npot].init(potname,filename,ordspl,"no",psmooth,derivsmooth,neighsmooth);   
          ordspl=0; 
          npot++;
        }
      }else if(stemp.find("$additional_potentials:")!=string::npos){
        nadd=0;
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          vadd.resize(nadd+1);
          sscanf(cline,"%s%le%d%le%le",ctemp1,&vadd[nadd].dist,&vadd[nadd].deriv,&vadd[nadd].val,&vadd[nadd].weight);
          potname=ctemp1;
          if (potname.length() != 4) {
            cerr << endl << "ERROR: " << endl << endl;
            cerr << "in " << inputfile << ": name of potential has to be four characters long, e.g. \'c_c_\', \'ClCl\'" 
                 << endl << "\" " << potname << "\" was read instead. Change this in the input and rerun." 
                 << endl << endl << "exit repopt" << endl << endl;
            exit(1);
          }
          vadd[nadd].potname = potname;
          //vadd[nadd].dist = vadd[nadd].dist / AA_Bohr;
          vadd[nadd].dist = vadd[nadd].dist;
          if (vadd[nadd].deriv > 3 ) {
            cerr << "ERROR: " << endl << endl;
            cerr << "additional equations can only be added for 0st,1st,2nd or 3rd dervative (0,1,2,3)!" 
                 << endl << endl << "exit repopt" << endl << endl;
            exit(1);
          }
          nadd++;
        }
      }else if(stemp.find("$compounds:")!=string::npos){
        nmol=0;
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          vmol.resize(nmol+1);
          sscanf(cline,"%s%le%le%le%s%s",ctemp1,&ediss,&eweight,&fweight,ctemp2,ctemp3);
          str=ctemp1;
          ifnotfile(str.c_str());
          dftbin=ctemp2;
          if(!runxtb) ifnotfile(dftbin.c_str());
          finp=ctemp3;
          if (finp != "0") ifnotfile(finp.c_str());
          vmol[nmol].init(str,ediss,eweight,fweight,dftbin," ",dftbversion,finp);
          vmol[nmol].ncontraineda = 0;
          ediss=0;
          nmol++;
        }
      }else if(stemp.find("$definition_reactions:")!=string::npos){
        nabbr=0;
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          vreamol.resize(nabbr+1);
          sscanf(cline,"%s%s%s",ctemp1,ctemp2,ctemp3);
          str=ctemp1;
          filename=ctemp2;
          ifnotfile(filename.c_str());
          dftbin=ctemp3;
          if(!runxtb) ifnotfile(dftbin.c_str());
          for (i=0; i<vmol.size(); i++){
            if (filename == vmol[i].name){
              vreamol[nabbr] = vmol[i];
              vreamol[nabbr].abbr = str;
              break;
            }
          }
          if (i==vmol.size()) vreamol[nabbr].init(filename,0,0,0,dftbin,str,dftbversion,"0");
          nabbr++;
        }
      }else if(stemp.find("$reactions:")!=string::npos){
        nrea=0;
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          vrea.resize(nrea+1);
          str=cline;
          vrea[nrea].init(str,vreamol);
          nrea++;
        }
      }else continue;
    }
    infile.close();
  }
  // check if a_b_ potential and b_a_ potential was given to optimize, if yes then only one will be considered
  for (int ipot=0; ipot < vpot.size()-1; ipot++){
    for (int jpot=ipot+1; jpot < vpot.size(); jpot++){  
      if ( vpot[ipot].potname == vpot[jpot].potname ){
        cout << "CAUTION: The potentials \"" << vpot[ipot].potname << "\" and \"" << vpot[jpot].potname
             << "\" will not be treated seperately. " << endl 
             << "Only the input data from potential \"" << vpot[ipot].potname 
             << "\" with the knot-vetor file \"" << vpot[ipot].gridname << "\" will be regarded in this repopt." << endl;
        vpot.erase(vpot.begin()+jpot);
      }
    }   
  }
}


void Allequations::rereadmol(const string inputfile){

  ifstream infile;
  string str, dftbin, potname, filename,finp,stmp;
  int tsmooth,npot,nadd,nmol,nabbr,nrea,derivsmooth=0,neighsmooth=1;
  int ordspl=0;
  double ediss=0,eatom=0,psmooth=0;
  double eweight,fweight;
  
  string stemp,stemp1;
  char cline[512],ctemp0[512],ctemp[512],ctemp1[512],ctemp2[512],ctemp3[512],ctemp4[512];
  int i,j,k,itmp,ntmp;
  double ftmp;
  ifnotfile(inputfile.c_str());
  infile.open(inputfile.c_str(),ios::in);
 
  if(!infile.is_open()){
    cout<<"unable to open "<<inputfile<<" file.\n";
    exit(1);
  }else{
    ctemp[0]=ctemp0[0]=ctemp1[0]=ctemp2[0]=ctemp3[0]=ctemp4[0]='\0';
    while(infile.getline(cline,512)){
      remove_comment(cline);
      stemp=cline;
      if(stemp.find("$compounds:")!=string::npos){
        nmol=0;
        vmol.resize(nmol);
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          vmol.resize(nmol+1);
          sscanf(cline,"%s%le%le%le%s%s",ctemp1,&ediss,&eweight,&fweight,ctemp2,ctemp3);
          str=ctemp1;
          ifnotfile(str.c_str());
          dftbin=ctemp2;
          if(!runxtb) ifnotfile(dftbin.c_str());
          finp=ctemp3;
          if (finp != "0") ifnotfile(finp.c_str());
          vmol[nmol].init(str,ediss,eweight,fweight,dftbin," ",dftbversion,finp);
          vmol[nmol].ncontraineda = 0;
          ediss=0;
          nmol++;
        }
      }else if(stemp.find("$definition_reactions:")!=string::npos){
        nabbr=0;
        vreamol.resize(nabbr);
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          vreamol.resize(nabbr+1);
          sscanf(cline,"%s%s%s",ctemp1,ctemp2,ctemp3);
          str=ctemp1;
          filename=ctemp2;
          ifnotfile(filename.c_str());
          dftbin=ctemp3;
          if(!runxtb) ifnotfile(dftbin.c_str());
          for (i=0; i<vmol.size(); i++){
            if (filename == vmol[i].name){
              vreamol[nabbr] = vmol[i];
              vreamol[nabbr].abbr = str;
              break;
            }
          }
          if (i==vmol.size()) vreamol[nabbr].init(filename,0,0,0,dftbin,str,dftbversion,"0");
          nabbr++;
        }
      }else if(stemp.find("$reactions:")!=string::npos){
        nrea=0;
        vrea.resize(nrea);
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          vrea.resize(nrea+1);
          str=cline;
          vrea[nrea].init(str,vreamol);
          nrea++;
        }
      }else continue;
    }
    infile.close();
  }
}

void Allequations::sort_inputdist(){
//set dummy vpot[ipot].vmoldist[0].dist =-1
  for (int ipot=0;ipot<vpot.size();ipot++){
    vpot[ipot].init_vmoldist();
  }
// use vmol to fill vpot[ipot].vmoldist
  for (int imol=0;imol<vmol.size();imol++){
    for (int iat1=0;iat1<vmol[imol].natom;iat1++){
      for (int iat2=iat1+1;iat2<vmol[imol].natom;iat2++){
        for (int ipot=0;ipot<vpot.size();ipot++){
          if(vmol[imol].atomname[iat1]+vmol[imol].atomname[iat2]==vpot[ipot].potname || 
          vmol[imol].atomname[iat2]+vmol[imol].atomname[iat1]==vpot[ipot].potname){
            vpot[ipot].insert_moldist(vmol[imol].name,vmol[imol].dist(iat1,iat2));     
            break;
          }
        }
      }
    }
  }

// use vreamol to fill vpot[ipot].vmoldist
  for (int imol=0;imol<vreamol.size();imol++){
    for (int iat1=0;iat1<vreamol[imol].natom;iat1++){
      for (int iat2=iat1+1;iat2<vreamol[imol].natom;iat2++){
        for (int ipot=0;ipot<vpot.size();ipot++){
          if(vreamol[imol].atomname[iat1]+vreamol[imol].atomname[iat2]==vpot[ipot].potname || 
           vreamol[imol].atomname[iat2]+vreamol[imol].atomname[iat1]==vpot[ipot].potname){
            vpot[ipot].insert_moldist(vreamol[imol].name,vreamol[imol].dist(iat1,iat2));     
            break;
          }
        }
      }
    }
  }

// delete dummy vpot[ipot].vmoldist[firstelement]
  for (int ipot=0;ipot<vpot.size();ipot++){
    vpot[ipot].remove_dummy_vmoldist();
  }

}

void Allequations::get_autogrids(){
  for (int ipot=0;ipot<vpot.size();ipot++){
    if (vpot[ipot].gridname.substr(0,8)=="autogrid")  vpot[ipot].get_autogrid();
  }
}


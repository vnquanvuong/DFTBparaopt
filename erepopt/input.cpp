#include<iostream>
#include<string>
#include<iomanip>
#include"erepobj.hpp"
#include"ddh.hpp"

using namespace std;
extern vector<string> addHamiltonian;
extern int ga_popsize, ga_ngen, ga_scoref, ga_flushf;
extern int preserved_num, seed;
extern double gtol;
extern double ga_pmut, ga_pcross;
extern bool ga,readr,restart;
extern bool runtest, skfclean;
extern int cpu_number;
extern int power;
extern sddh ddh;
extern string popfinalfile, popinitialfile;
extern string scratch;

void remove_comment(char *cstr);

void Erepobj::readinp(const string inputfile){

  addskf=addrep=false;
  ngrid=100;
  dgrid=0.1;
  omega=0.3;
  onsiten=onsitep=false;
  score_type=2;
  fit_type=0;
  outfilename="erepopt.log";
  ifstream infile;
  string stemp,stemp1,xyzfile,dosfile;
  char cline[512],ctemp[512],ctemp0[512],ctemp1[512],ctemp2[512],ctemp3[512],ctemp4[512];
  int i,j,k,i1,i2,i3,i4,i5,itmp,ntmp;
  double d1,d2,d3,d4,fmin,ftmp,fmax;
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
          sscanf(cline,"%s",ctemp1);
          if(stemp.find("dftbversion")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); dftbversion=ctemp2;}
          if(stemp.find("scratchfolder")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2);  scratch=ctemp2;}
          if(stemp.find("popfinalfile")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2);  popfinalfile=ctemp2;}
          if(stemp.find("popinitialfile")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); popinitialfile=ctemp2;}
          if(stemp.find("skgen")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); skgen=ctemp2;}
          if(stemp.find("onecent")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); onecent=ctemp2;}
          if(stemp.find("twocent")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); twocent=ctemp2;}
          if(stemp.find("outfile")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); outfilename=ctemp2;}
          if(stemp.find("libdir")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); libdir=ctemp2;}
          if(stemp.find("libadddir")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); libadddir=ctemp2;addskf=true;}
          if(stemp.find("repadddir")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); repadddir=ctemp2;addrep=true;}
          if(stemp.find("omega")!=string::npos){ sscanf(cline,"%s %le",ctemp1,&omega); lc=true;}
          if(stemp.find("dgrid")!=string::npos){ sscanf(cline,"%s %le",ctemp1,&dgrid);}
          if(stemp.find("ngrid")!=string::npos){ sscanf(cline,"%s %d",ctemp1,&ngrid);}
          if(stemp.find("power")!=string::npos){ sscanf(cline,"%s %d",ctemp1,&power);}
          if(stemp.find("nthreads")!=string::npos){ sscanf(cline,"%s %d",ctemp1,&cpu_number);}
          if(stemp.find("grids")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); fgrid=ctemp2;}
          if(stemp.find("dftb_inp")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); fdftb_inp=ctemp2;}
          if(stemp.find("rep.in")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); frep_in=ctemp2;}
          if(stemp.find("gasrepfit")!=string::npos){ sscanf(cline,"%s %s",ctemp1,ctemp2); gasrepfit=ctemp2;}
          if(stemp.find("onsiten")!=string::npos){ sscanf(cline,"%s %le",ctemp1,&esn); onsiten=true;}
          if(stemp.find("onsitep")!=string::npos){ sscanf(cline,"%s %le",ctemp1,&esp); onsitep=true;}
          if(stemp.find("skfclean")!=string::npos){ sscanf(cline,"%s %d",ctemp1,&itmp); skfclean= itmp;}
        }
      }else if(stemp.find("addHamiltonian:")!=string::npos){
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          addHamiltonian.push_back(stemp);
        }
      }else if(stemp.find("$genetic_algorithm:")!=string::npos){
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          if(stemp.find("popsize")!=string::npos) sscanf(cline,"%s %d",ctemp1,&ga_popsize);
          if(stemp.find("ngen")!=string::npos) sscanf(cline,"%s %d",ctemp1,&ga_ngen);
          if(stemp.find("preserved_num")!=string::npos) sscanf(cline,"%s %d",ctemp1,&preserved_num);
          if(stemp.find("scoref")!=string::npos) sscanf(cline,"%s %d",ctemp1,&ga_scoref);
          if(stemp.find("flushf")!=string::npos) sscanf(cline,"%s %d",ctemp1,&ga_flushf);
          if(stemp.find("score_type")!=string::npos) sscanf(cline,"%s %d",ctemp1,&score_type);
          if(stemp.find("fit_type")!=string::npos) sscanf(cline,"%s %d",ctemp1,&fit_type);
          if(stemp.find("seed")!=string::npos) sscanf(cline,"%s %d",ctemp1,&seed);
          if(stemp.find("pmut")!=string::npos) sscanf(cline,"%s %le",ctemp1,&ga_pmut);
          if(stemp.find("pcross")!=string::npos) sscanf(cline,"%s %le",ctemp1,&ga_pcross);
          if(stemp.find("tol")!=string::npos) sscanf(cline,"%s %le",ctemp1,&gtol);
          if(stemp.find("ga")!=string::npos){ sscanf(cline,"%s %d",ctemp1,&itmp); ga = itmp;}
          if(stemp.find("runtest")!=string::npos){ sscanf(cline,"%s %d",ctemp1,&itmp); runtest = itmp;}
          if(stemp.find("readr")!=string::npos){ sscanf(cline,"%s %d",ctemp1,&itmp); readr= itmp;}
          if(stemp.find("restart")!=string::npos){ sscanf(cline,"%s %d",ctemp1,&itmp); restart= itmp;}
        }
      }else if(stemp.find("$atomic_energy:")!=string::npos){
        ntmp=0;
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          atomicenergy.resize(ntmp+1);
          atomicname.resize(ntmp+1);
          sscanf(cline,"%s%le",ctemp1,&d1);
          atomicname[ntmp]=ctemp1;
          atomicenergy[ntmp]=d1;
          ntmp++;
        }
      }else if(stemp.find("$element_types:")!=string::npos){
        nelem=0;
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          velem.resize(nelem+1);
          sscanf(cline,"%s%s%d",ctemp1,ctemp2,&ntmp);
          if(ntmp==0){
            velem[nelem].radius.resize(2);
            sscanf(cline,"%s%d%d%le%le%le%d%le%le%le%d",ctemp1,&velem[nelem].optype,&velem[nelem].lmax,
                   &velem[nelem].radius[0].minr,&velem[nelem].radius[0].r,&velem[nelem].radius[0].maxr,&velem[nelem].radius[0].precision,
                   &velem[nelem].radius[1].minr,&velem[nelem].radius[1].r,&velem[nelem].radius[1].maxr,&velem[nelem].radius[1].precision);
            velem[nelem].name=ctemp1;
            
          }else if(ntmp==1){
            velem[nelem].radius.resize(3);
            sscanf(cline,"%s%d%d%le%le%le%d%le%le%le%d%le%le%le%d",ctemp1,&velem[nelem].optype,&velem[nelem].lmax,
                   &velem[nelem].radius[0].minr,&velem[nelem].radius[0].r,&velem[nelem].radius[0].maxr,&velem[nelem].radius[0].precision,
                   &velem[nelem].radius[1].minr,&velem[nelem].radius[1].r,&velem[nelem].radius[1].maxr,&velem[nelem].radius[1].precision,
                   &velem[nelem].radius[2].minr,&velem[nelem].radius[2].r,&velem[nelem].radius[2].maxr,&velem[nelem].radius[2].precision);
            velem[nelem].name=ctemp1;
            
          }else if(ntmp==2){
            velem[nelem].radius.resize(4);
            sscanf(cline,"%s%d%d%le%le%le%d%le%le%le%d%le%le%le%d%le%le%le%d",ctemp1,&velem[nelem].optype,&velem[nelem].lmax,
                   &velem[nelem].radius[0].minr,&velem[nelem].radius[0].r,&velem[nelem].radius[0].maxr,&velem[nelem].radius[0].precision,
                   &velem[nelem].radius[1].minr,&velem[nelem].radius[1].r,&velem[nelem].radius[1].maxr,&velem[nelem].radius[1].precision,
                   &velem[nelem].radius[2].minr,&velem[nelem].radius[2].r,&velem[nelem].radius[2].maxr,&velem[nelem].radius[2].precision,
                   &velem[nelem].radius[3].minr,&velem[nelem].radius[3].r,&velem[nelem].radius[3].maxr,&velem[nelem].radius[3].precision);
            velem[nelem].name=ctemp1;
          }else if(ntmp==3){
            velem[nelem].radius.resize(5);
            sscanf(cline,"%s%d%d%le%le%le%d%le%le%le%d%le%le%le%d%le%le%le%d%le%le%le%d",ctemp1,&velem[nelem].optype,&velem[nelem].lmax,
                   &velem[nelem].radius[0].minr,&velem[nelem].radius[0].r,&velem[nelem].radius[0].maxr,&velem[nelem].radius[0].precision,
                   &velem[nelem].radius[1].minr,&velem[nelem].radius[1].r,&velem[nelem].radius[1].maxr,&velem[nelem].radius[1].precision,
                   &velem[nelem].radius[2].minr,&velem[nelem].radius[2].r,&velem[nelem].radius[2].maxr,&velem[nelem].radius[2].precision,
                   &velem[nelem].radius[3].minr,&velem[nelem].radius[3].r,&velem[nelem].radius[3].maxr,&velem[nelem].radius[3].precision,
                   &velem[nelem].radius[4].minr,&velem[nelem].radius[4].r,&velem[nelem].radius[4].maxr,&velem[nelem].radius[4].precision);
            velem[nelem].name=ctemp1;
          }
          nelem++;
        }
        for(i=0;i<velem.size();i++){ 
          for(j=0; j<velem[i].lmax+2; j++){
            if(velem[i].radius[j].precision==0) velem[i].radius[j].delta=1.0;
            else if(velem[i].radius[j].precision==1) velem[i].radius[j].delta=0.1;
            else if(velem[i].radius[j].precision==2) velem[i].radius[j].delta=0.01;
            else if(velem[i].radius[j].precision==3) velem[i].radius[j].delta=0.001;
            else if(velem[i].radius[j].precision==4) velem[i].radius[j].delta=0.0001;
            else if(velem[i].radius[j].precision==5) velem[i].radius[j].delta=0.00001;
          }
        }
      }else if(stemp.find("$d3:")!=string::npos){
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
        if(ntmp>0) ddh.td3=true;
        for(i=0;i<ddh.d3.size();i++){ 
          if     (ddh.d3[i].precision==0) ddh.d3[i].delta=1.0;
          else if(ddh.d3[i].precision==1) ddh.d3[i].delta=0.1;
          else if(ddh.d3[i].precision==2) ddh.d3[i].delta=0.01;
          else if(ddh.d3[i].precision==3) ddh.d3[i].delta=0.001;
          else if(ddh.d3[i].precision==4) ddh.d3[i].delta=0.0001;
          else if(ddh.d3[i].precision==5) ddh.d3[i].delta=0.00001;
        }
      }else if(stemp.find("$damph:")!=string::npos){
        ntmp=0;
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
          ntmp++;
        }
        if(ntmp>0) ddh.tdamph=true;
        if     (ddh.damph.precision==0) ddh.damph.delta=1.0;
        else if(ddh.damph.precision==1) ddh.damph.delta=0.1;
        else if(ddh.damph.precision==2) ddh.damph.delta=0.01;
        else if(ddh.damph.precision==3) ddh.damph.delta=0.001;
        else if(ddh.damph.precision==4) ddh.damph.delta=0.0001;
        else if(ddh.damph.precision==5) ddh.damph.delta=0.00001;
      }else if(stemp.find("$hubbardderivs:")!=string::npos){
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
        if(ntmp>0) ddh.thubbardderivs=true;
        for(i=0;i<ddh.hubbardderivs.size();i++){ 
          if     (ddh.hubbardderivs[i].precision==0) ddh.hubbardderivs[i].delta=1.0;
          else if(ddh.hubbardderivs[i].precision==1) ddh.hubbardderivs[i].delta=0.1;
          else if(ddh.hubbardderivs[i].precision==2) ddh.hubbardderivs[i].delta=0.01;
          else if(ddh.hubbardderivs[i].precision==3) ddh.hubbardderivs[i].delta=0.001;
          else if(ddh.hubbardderivs[i].precision==4) ddh.hubbardderivs[i].delta=0.0001;
          else if(ddh.hubbardderivs[i].precision==5) ddh.hubbardderivs[i].delta=0.00001;
        }
      }else if(stemp.find("$hubbardderivsfull:")!=string::npos){
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
        if(ntmp>0) ddh.thubbardderivs=true;
        for(i=0;i<ddh.hubbardderivs.size();i++){ 
          if     (ddh.hubbardderivs[i].precision==0) ddh.hubbardderivs[i].delta=1.0;
          else if(ddh.hubbardderivs[i].precision==1) ddh.hubbardderivs[i].delta=0.1;
          else if(ddh.hubbardderivs[i].precision==2) ddh.hubbardderivs[i].delta=0.01;
          else if(ddh.hubbardderivs[i].precision==3) ddh.hubbardderivs[i].delta=0.001;
          else if(ddh.hubbardderivs[i].precision==4) ddh.hubbardderivs[i].delta=0.0001;
          else if(ddh.hubbardderivs[i].precision==5) ddh.hubbardderivs[i].delta=0.00001;
        }
      }else if(stemp.find("$damphver2:")!=string::npos){
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
        if(ntmp>0) ddh.tdamphver2=true;
        for(i=0;i<ddh.damphver2.size();i++){ 
          if     (ddh.damphver2[i].precision==0) ddh.damphver2[i].delta=1.0;
          else if(ddh.damphver2[i].precision==1) ddh.damphver2[i].delta=0.1;
          else if(ddh.damphver2[i].precision==2) ddh.damphver2[i].delta=0.01;
          else if(ddh.damphver2[i].precision==3) ddh.damphver2[i].delta=0.001;
          else if(ddh.damphver2[i].precision==4) ddh.damphver2[i].delta=0.0001;
          else if(ddh.damphver2[i].precision==5) ddh.damphver2[i].delta=0.00001;
        }
      }else if(stemp.find("$hubbards:")!=string::npos){
        ntmp=0;
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue ;
          sscanf(cline,"%s %s %le %le %le %d",ctemp1,ctemp2,&fmin,&ftmp,&fmax,&itmp);
          ddh.hubbards.resize(ntmp+1);
          ddh.hubbards[ntmp].name=ctemp1;
          ddh.hubbards[ntmp].type=ctemp2;
          ddh.hubbards[ntmp].min=fmin;
          ddh.hubbards[ntmp].value=ftmp;
          ddh.hubbards[ntmp].max=fmax;
          ddh.hubbards[ntmp].precision=itmp;
          ntmp++;
        }
        if(ntmp>0) ddh.thubbards=true;
        for(i=0;i<ddh.hubbards.size();i++){ 
          if     (ddh.hubbards[i].precision==0) ddh.hubbards[i].delta=1.0;
          else if(ddh.hubbards[i].precision==1) ddh.hubbards[i].delta=0.1;
          else if(ddh.hubbards[i].precision==2) ddh.hubbards[i].delta=0.01;
          else if(ddh.hubbards[i].precision==3) ddh.hubbards[i].delta=0.001;
          else if(ddh.hubbards[i].precision==4) ddh.hubbards[i].delta=0.0001;
          else if(ddh.hubbards[i].precision==5) ddh.hubbards[i].delta=0.00001;
        }
      }else if(stemp.find("$vorbes:")!=string::npos){
        ntmp=0;
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue ;
          sscanf(cline,"%s %s %le %le %le %d",ctemp1,ctemp2,&fmin,&ftmp,&fmax,&itmp);
          ddh.vorbes.resize(ntmp+1);
          ddh.vorbes[ntmp].name=ctemp1;
          ddh.vorbes[ntmp].type=ctemp2;
          ddh.vorbes[ntmp].min=fmin;
          ddh.vorbes[ntmp].value=ftmp;
          ddh.vorbes[ntmp].max=fmax;
          ddh.vorbes[ntmp].precision=itmp;
          ntmp++;
        }
        if(ntmp>0) ddh.tvorbes=true;
        for(i=0;i<ddh.vorbes.size();i++){ 
          if     (ddh.vorbes[i].precision==0) ddh.vorbes[i].delta=1.0;
          else if(ddh.vorbes[i].precision==1) ddh.vorbes[i].delta=0.1;
          else if(ddh.vorbes[i].precision==2) ddh.vorbes[i].delta=0.01;
          else if(ddh.vorbes[i].precision==3) ddh.vorbes[i].delta=0.001;
          else if(ddh.vorbes[i].precision==4) ddh.vorbes[i].delta=0.0001;
          else if(ddh.vorbes[i].precision==5) ddh.vorbes[i].delta=0.00001;
        }
      }else if(stemp.find("$compounds:")!=string::npos){
        ncompound=0;
        while(infile.getline(cline,512)){
          remove_comment(cline);
          stemp=cline;
          if(stemp.find("$end")!=string::npos) break;
          if(sscanf(cline,"%s",ctemp)<0) continue;
          if(ctemp[0]=='#' ) continue;
          vmol.resize(ncompound+1);
          sscanf(cline,"%s%d%le%le%le%s%le",ctemp1,&i1,&d1,&d2,&d3,ctemp2,&d4);
          xyzfile=ctemp1;
          ifnotfile(xyzfile.c_str());
          dosfile=ctemp2;
          ifnotfile(dosfile.c_str());
          vmol[ncompound].init (xyzfile, i1, d1, d2, d3, dosfile, d4);
          ncompound++;
        }
      }
    }
    infile.close();
  }
//fout.open(outfilename.c_str());
//fout<<std::fixed;
}

void Erepobj::remove_comment(char *cstr){
  int i,ipos;
  string sstr;
  sstr=cstr;
  std::string::size_type pos = sstr.find('#');
  if (pos!= std::string::npos){
    ipos=pos;
    for(i=ipos;i<=sstr.length();i++){
      cstr[i]='\0';
    }
  }
}

void Erepobj::ifnotfile(const string filename){
  ifstream file;

  file.open(filename.c_str());
  if ( !file ) {
    cerr << endl << "ERROR: " << endl << endl;
    cerr << filename << ": file not found" 
         << endl << "exit erepobj" << endl << endl;
    exit(1);
  }
  file.close();
}


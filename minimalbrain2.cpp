/***********************************************************************************
   Created: 2016-09-28
   Modified: 2016-12-12
   Authors: Anders Lansner, Örjan Ekeberg
************************************************************************************/

#include "math.h"
#include <vector>
#include <string>
#include <mpi.h>

#include "base/source/Globals.h"
#include "base/source/Pop.h"
#include "base/source/SensorPop.h" 
#include "base/source/MotorPop.h" 
#include "base/source/SensoryMotorPop.h" 
#include "base/source/PrjD.h"
#include "base/source/Prj.h"
#include "base/source/Prj11.h"
#include "base/source/Logger.h"
#include "utils/Utils.h"
#include "utils/Patio.h"
#include "utils/Parseparam.h"

using namespace std;

int nrank = 8,H = 8,U = 88/H,dn = 40,nldot = 0,ctlmode = 0,trnpat = -1,stoprun = 0;
float taum = 0.002,wtagain = 8.,musicdelay = 0.010,momobgain = 1.,momowegain = 1.,momowigain = 14.,
    semobgain = 1.,semowgain = 1.,tauzi = 0.004,tauzj = 0.004,taue = 0.004,taup = 10.,noise = 0.,
    lgbias = 0.,dmax = 16,da = tauzi,dq = 2.,motaua = 1.,moadampl = 0.,learntime = 4.;
bool dolog = true,doprn = true,externctl;
int selgilogstep = 0,seactlogstep = 1,mowsulogstep = 0,modsulogstep = 0,moactlogstep = 1,moadalogstep = 0,
    ctllgilogstep = 0;

string paramfile = "params_3.par",traininfile = "bwv772_100Hz.dat";

SensorPop *sepop,*ctlpop;
MotorPop *mopop;
PrjD *momoprj1,*momoprj2;
Prj11 *semoprj;

vector<int> tapsvec;
vector<float> tauzivec;

void parseparams(string paramfile) {
    Parseparam *pp = new Parseparam(paramfile);

    pp->postparam("nrank",&nrank,Parseparam::Int);
    pp->postparam("H",&H,Parseparam::Int);
    pp->postparam("U",&U,Parseparam::Int);
    pp->postparam("taum",&taum,Parseparam::Float);
    pp->postparam("noise",&noise,Parseparam::Float);
    pp->postparam("wtagain",&wtagain,Parseparam::Float);
    pp->postparam("lgbias",&lgbias,Parseparam::Float);

    pp->postparam("dolog",&dolog,Parseparam::Bool);
    pp->postparam("doprn",&doprn,Parseparam::Bool);
    pp->postparam("selgilogstep",&selgilogstep,Parseparam::Int);
    pp->postparam("seactlogstep",&seactlogstep,Parseparam::Int);
    pp->postparam("mowsulogstep",&mowsulogstep,Parseparam::Int);
    pp->postparam("modsulogstep",&modsulogstep,Parseparam::Int);
    pp->postparam("moactlogstep",&moactlogstep,Parseparam::Int);
    pp->postparam("moadalogstep",&moadalogstep,Parseparam::Int);
    pp->postparam("ctllgilogstep",&ctllgilogstep,Parseparam::Int);

    pp->postparam("tauzi",&tauzi,Parseparam::Float);
    pp->postparam("tauzj",&tauzj,Parseparam::Float);
    pp->postparam("taue",&taue,Parseparam::Float);
    pp->postparam("taup",&taup,Parseparam::Float);
    pp->postparam("momobgain",&momobgain,Parseparam::Float);
    pp->postparam("momowegain",&momowegain,Parseparam::Float);
    pp->postparam("momowigain",&momowigain,Parseparam::Float);
    pp->postparam("semobgain",&semobgain,Parseparam::Float);
    pp->postparam("semowgain",&semowgain,Parseparam::Float);
    pp->postparam("dn",&dn,Parseparam::Int);
    pp->postparam("da",&da,Parseparam::Float);
    pp->postparam("dmax",&dmax,Parseparam::Float);
    pp->postparam("dq",&dq,Parseparam::Float);
    pp->postparam("nldot",&nldot,Parseparam::Int);
    pp->postparam("ctlmode",&ctlmode,Parseparam::Int);
    pp->postparam("externctl",&externctl,Parseparam::Int);
    pp->postparam("traininfile",&traininfile,Parseparam::String);
    pp->postparam("trnpat",&trnpat,Parseparam::Int);
    pp->postparam("learntime",&learntime,Parseparam::Float);
    pp->postparam("motaua",&motaua,Parseparam::Float);
    pp->postparam("moadampl",&moadampl,Parseparam::Float);

    pp->doparse();
}

void mkdelay(float dmin,float dmax,int n,float tauzi,bool doprn = false) {
    if (dmax<=dmin)
	Utils::mpierror("main::mkdelay","Illegal dmax<=dmin");
    if (n<=0)
	Utils::mpierror("main::mkdelay","Illegal n<=0");
    tapsvec = vector<int>(n);
    tauzivec = vector<float>(n);
    float ds = (dmax - dmin)/n;
    for (int i=0; i<n; i++) {
	tapsvec[i] = (dmin + (i+1)*ds + Globals::_dt/2)/Globals::_dt;
	tauzivec[i] = tauzi;
    }
    if (doprn and ISROOT) {
	for (int i=0; i<n; i++)
	    fprintf(stderr,"%.4f %.4f\n",tapsvec[i]*Globals::_dt,tauzivec[i]);
    }
}

void setuplogging() {
    new PopLogger(sepop,"selgi.log",Pop::LGI,selgilogstep);
    new PopLogger(sepop,"seact.log",Pop::ACT,seactlogstep);

    new PopLogger(mopop,"mowsu.log",Pop::WSU,mowsulogstep);
    new PopLogger(mopop,"modsu.log",Pop::DSU,modsulogstep);
    new PopLogger(mopop,"moact.log",Pop::ACT,moactlogstep);
    new PopLogger(mopop,"moada.log",Pop::ADA,moadalogstep);

    new PopLogger(ctlpop,"ctlgi.log",Pop::LGI,ctllgilogstep);
    new PopLogger(ctlpop,"ctact.log",Pop::ACT,ctllgilogstep);

}

void doprinting() {
    new PrjPrinter(semoprj,"semobj.bin",Prj::BJ);
    new PrjPrinter(semoprj,"semowij.bin",Prj::WIJ);
    new PrjPrinter(momoprj1,"momobj1.bin",Prj::BJ);
    new PrjPrinter(momoprj1,"momowij1.bin",Prj::WIJ);
    new PrjPrinter(momoprj2,"momobj2.bin",Prj::BJ);
    new PrjPrinter(momoprj2,"momowij2.bin",Prj::WIJ);
}

#define RESETDELAYBUF

enum BrainMode { RESET = 0, LEARN, RECALL } ;
BrainMode _brainmode = LEARN;

void prnbrainmode(BrainMode brainmode) {
    if (ISROOT) {
	string modestr;
	switch (brainmode) {
	case RESET: modestr = "RESET"; break;
	case LEARN: modestr = "LEARN"; break;
	case RECALL: modestr = "RECALL"; break;
	}
	fprintf(stderr," [%.2f %s %.2f] ",Globals::_simtime,modestr.c_str(),lgbias);
    }
}

void setmode(BrainMode brainmode) {
    if (Globals::_simtime<0.01) brainmode = LEARN;
    if (brainmode!=_brainmode) prnbrainmode(brainmode);
    switch (brainmode) {
    case RESET:
	semoprj->setprn(0.); semoprj->setbgain(0.); semoprj->setwgain(0.);
	momoprj1->setprn(0.); momoprj1->setbgain(0.); momoprj1->setwegain(0.); momoprj1->setwigain(0.);
	momoprj2->setprn(0.); momoprj2->setbgain(0.); momoprj2->setwegain(0.); momoprj2->setwigain(0.);
#ifdef RESETDELAYBUF
	if (brainmode!=_brainmode) {
	    momoprj1->resetdelaybuf();
	    momoprj2->resetdelaybuf();
	}
#endif
	mopop->setlgbias(0.);
	mopop->setadapt(motaua,0.);
	break;
    case LEARN:
	semoprj->setprn(1.); semoprj->setbgain(1.); semoprj->setwgain(1.);
	momoprj1->setprn(1.); momoprj1->setbgain(0.); momoprj1->setwegain(0.); momoprj1->setwigain(0.);
	momoprj2->setprn(1.); momoprj2->setbgain(0.); momoprj2->setwegain(0.); momoprj2->setwigain(0.);
	mopop->setlgbias(0.);
	mopop->setadapt(motaua,0.);
	break;
    case RECALL:
	semoprj->setprn(0.); semoprj->setbgain(1.); semoprj->setwgain(1.);
	momoprj1->setprn(0.); momoprj1->setbgain(momobgain/2);
	momoprj1->setwegain(momowegain); momoprj1->setwigain(momowigain);
	momoprj2->setprn(0.); momoprj2->setbgain(momobgain/2);
	momoprj2->setwegain(momowegain); momoprj2->setwigain(momowigain);
	mopop->setlgbias(lgbias);
	mopop->setadapt(motaua,moadampl);
	break;
    }
    _brainmode = brainmode;
}

/* _brainmode is the enabled mode */
void proctl() {
    int brainmode;
    float lgb;
    if (Globals::getmpirank()==ctlpop->getrank0()) {
	int dolea = (ctlpop->getact()[0]>0.5);
 	int dores = (ctlpop->getact()[1]>0.5);
 	int dostop = (ctlpop->getact()[3]>0.5);
	lgb = 6 * ctlpop->getact()[4] - 3.;
	lgbias = lgb;
	if (dostop>0) {
	    stoprun = 1;
	    if (ISROOT) fprintf(stderr,"Stopping!\n");
	} else if (dores>0)
	    brainmode = (int)RESET;
	else if (dolea>0)
	    brainmode = (int)LEARN;
	else
	    brainmode = (int)RECALL;
    }
    MPI_Bcast(&stoprun,1,MPI_INT,ctlpop->getrank0(),Globals::_comm_world);
    MPI_Bcast(&brainmode,1,MPI_INT,ctlpop->getrank0(),Globals::_comm_world);
    MPI_Bcast(&lgbias,1,MPI_FLOAT,ctlpop->getrank0(),Globals::_comm_world);
    setmode((BrainMode)brainmode);

}

int main(int argc, char** argv) {
   
  Globals::init(argc,argv);
  Globals::setdt(0.001);

  parseparams(paramfile);

  sepop = new SensorPop("sensoryinput",musicdelay,nrank,H,U,taum); 
  sepop->setnormtype(Pop::NONE);
  sepop->setwtagain(1.);

  ctlpop = new SensorPop("controlinput",musicdelay,1,1,5,Globals::_dt); 

  mopop = new MotorPop("motoroutput",musicdelay,nrank,H,U,taum); 
  mopop->setnormtype(Pop::NONE);
  mopop->setwtagain(wtagain);
  mopop->setiwgain(0.);
    
  mkdelay(0,dmax/2,dn/2,tauzi);
  momoprj1 = new PrjD(mopop,mopop,tapsvec,tauzivec,tauzj,taue,taup);
  momoprj1->setselfconn(true,true,true);
  momoprj1->reinit();

  mkdelay(dmax/2,dmax,dn/2,tauzi);
  momoprj2 = new PrjD(mopop,mopop,tapsvec,tauzivec,tauzj,taue,taup);
  momoprj2->setselfconn(true,true,true);
  momoprj2->reinit();

  semoprj = new Prj11(sepop,mopop,Globals::_dt,Globals::_dt,taue,10.);
  vector<float> bj(H*U,-3.5); // Something strange with this -3.5 -- 5.
  vector<vector<float> > wij(H*U,vector<float>(H*U,0.));
  for (int u=0; u<H*U; u++) wij[u][u] = 5;
  semoprj->initbw(bj,wij);

  if (ISROOT) std::cerr << "music simtime = " << Globals::_musicstoptime << std::endl;

  int n = tapsvec.size();

  float k = (float)n/16.;
  momowegain /= k; momowigain /= k;

  MPI_Barrier(Globals::_comm_world);

  if (ISROOT) {
      fprintf(stderr,"H = %d U = %d taum = %.3f ",H,U,taum);
      fprintf(stderr,"taup = %.3f sec semobgain = %.4f semowgain = %.4f ",
	      taup,semobgain,semowgain,n);
      fprintf(stderr,"momobgain = %.4f momowegain = %.4f momowigain = %.4f n = %d\n",
	      momobgain,momowegain,momowigain,n);
  }
  MPI_Barrier(Globals::_comm_world);
  
  if (dolog) setuplogging();

  if (ISROOT) std::cerr << "Setup done!" << std::endl;

  sepop->setiwgain(1.);

  momoprj1->setbgain(0.);
  momoprj1->setwegain(0.);
  momoprj1->setwigain(0.);
  momoprj1->setprn(1.);

  momoprj2->setbgain(0.);
  momoprj2->setwegain(0.);
  momoprj2->setwigain(0.);
  momoprj2->setprn(1.);

  semoprj->setbwgain(semobgain); semoprj->setwgain(semowgain);
  semoprj->setprn(1.);


  if (ISROOT) std::cerr << "music simtime = " << Globals::_musicstoptime << " sec" << std::endl;

  Globals::setnldot(nldot);

  float simtime; int inon;

  Globals::start();

  /* Read from file here */
  
  Patio *patio = new Patio(sepop);
  int idur;
  if (patio->checkfile(traininfile)) {
      if (ISROOT) fprintf(stderr,"Training with infile %s\n",traininfile.c_str());
      patio->readfile(traininfile,-1,false);
      idur = (Globals::_musictimestep + Globals::_dt/2)/Globals::_dt;
      for (int p=0; p<min(trnpat,patio->getnpat()); p++) {
	  sepop->setlgi(patio->gethypcolpat(p),true);
	  for (int s=0; s<idur; s++) Globals::updstateall(false);
      }
      if (ISROOT) fprintf(stderr," Done!\n");
  }

  while ((stoprun==0) and (Globals::_musicruntime->time() < Globals::_musicstoptime - 3)) {

      simtime = Globals::_simtime;

      if (not externctl) {

	  switch (ctlmode) {
	  case 0:
	      setmode(RESET);
	      break;
	  case 1:
	      setmode(LEARN); 
	      break;
	  case 2:
	      setmode(RECALL);
	      break;
	  case 10:
	      inon = sepop->inon();
	      if (inon==0 or simtime>learntime) {
		  setmode(RECALL);
		  //	      if (ISROOT) fprintf(stderr,"R:%5.3f:%d\n",simtime,inon);
	      } else {
		  setmode(LEARN);
		  //	      if (ISROOT) fprintf(stderr,"L:%5.3f:%d\n",simtime,inon);
	      }
	      break;
	  case 101:
	      if (simtime>0. and simtime<learntime) setmode(LEARN);
	      if (simtime>learntime and simtime<7.) setmode(RECALL);
	      if (simtime>7. and simtime<10.) setmode(RESET);
	      if (simtime>10. and simtime<16.) setmode(RECALL);
	      if (simtime>16. and simtime<17.) setmode(RESET);
	      if (simtime>17. and simtime<21.) setmode(LEARN);
	      if (simtime>21. and simtime<24.) setmode(RECALL);
	      if (simtime>24. and simtime<26.) setmode(RESET);
	      if (simtime>26. and simtime<100.) setmode(RECALL);
	      break;
	  case 102:
	      if (simtime<learntime) setmode(LEARN);
	      if (simtime>3. and simtime<5.) setmode(RESET);
	      if (simtime>5.) setmode(RECALL);
	      break;
	  case 200:
	      if (simtime>0. and simtime<10.) setmode(LEARN);
	      if (simtime>10. and simtime<100.) setmode(RECALL);
	      break;
	  case 201:
	      if (simtime>0. and simtime<25. and simtime<learntime) setmode(LEARN);
	      if ((simtime>25. or simtime>learntime) and simtime<35.5) setmode(RECALL);
	      if (simtime>35.5 and simtime<35.6) setmode(RESET);
	      if (simtime>36. and simtime<100.) setmode(RECALL);
	      break;
	  case 202:
	      if (simtime>0. and simtime<10. and simtime<learntime) setmode(LEARN);
	      if ((simtime>10. or simtime>learntime) and simtime<13.5) setmode(RECALL);
	      if (simtime>13.5 and simtime<13.6) setmode(RESET);
	      if (simtime>13.6 and simtime<100.) setmode(RECALL);
	      break;
	  case 203:
	      if (simtime>0. and simtime<2.1 and simtime<learntime) setmode(LEARN);
	      if (simtime>2.1 and simtime<4.2) setmode(RESET);
	      if (simtime>4.2) setmode(RECALL);
	      break;
	  default: Utils::mpierror("main::mkdelay","Illegal ctlmode");
	  }
      } else proctl();

      Globals::updstateall();
  }

  Globals::stop();

  MPI_Barrier(Globals::_comm_world);
    
  if (doprn) doprinting();

  Globals::report();
    
  Globals::finish();

  return 0;
}

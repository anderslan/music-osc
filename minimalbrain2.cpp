/***********************************************************************************
   Created: 2016-09-28
   Modified: 2016-11-24
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
#include "utils/Parseparam.h"

using namespace std;

int nrank = 8,H = 8,U = 88/H,dn = 40,nstep = 100,idur = 50,nldot = 0;
float taum = 0.002,wtagain = 8.,musicdelay = 0.010,momobgain = 1.,momowegain = 1.,momowigain = -14.,
    semobgain = 1.,semowgain = 1.,tauzi = 0.004,tauzj = 0.004,taue = 0.004,taup = 10.,noise = 0.,
    lgbias = 0.,dmax = 16,da = tauzi,dq = 2.,time1 = 1.e6,time2 = 1.e6,time3 = 1.e6,time4 = 1.e6,
    time5 = 1.e6;
bool dolog = true;

string paramfile = "params_3.par";

SensorPop *sepop;
MotorPop *mopop;
PrjD *momoprj;
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
    pp->postparam("nstep",&nstep,Parseparam::Int);
    pp->postparam("time1",&time1,Parseparam::Float);
    pp->postparam("time2",&time2,Parseparam::Float);
    pp->postparam("time3",&time3,Parseparam::Float);
    pp->postparam("time4",&time4,Parseparam::Float);
    pp->postparam("time5",&time5,Parseparam::Float);

    pp->doparse();
}

void mkdelay(float d0,float dmax,float a,int n,float q) {
    float dt = Globals::_dt,c;
    tapsvec = vector<int>(n);
    tauzivec = vector<float>(n);
    if (n==1) {
	c = 0;
	a = d0;
    } else if (a==0) {
	a = (dmax - d0)/(n-1);
	c = 0;
    } else if (dmax<d0 + n*a) {
	if (ISROOT) fprintf(stderr,"dmax = %f d0 = %f n = %d a = %f\n",dmax,d0,n,a);
	Utils::mpierror("main","Illegal dmax<d0+n*a");
    } else
	//	c = 1./n * log((dmax - d0)/d0);
	c = 0.;
    for (int i=0; i<n; i++) {
	tapsvec[i] = ((a*i + d0*exp(c*i)) + dt/2)/dt;
	tauzivec[i] = q * (a + d0*c*exp(c*i));
    }
//     if (ISROOT) {
// 	fprintf(stderr,"a = %.4f c = %.4f q = %.4f\n",a,c,q);
// 	for (int i=0; i<n; i++)
// 	    fprintf(stderr,"%.4f %.4f\n",tapsvec[i]*Globals::_dt,tauzivec[i]);
//     }
}

void setuplogging() {
    new PopLogger(sepop,"selgi.log",Pop::LGI);
    new PopLogger(sepop,"seact.log",Pop::ACT);

    new PopLogger(mopop,"mowsu.log",Pop::WSU);
    new PopLogger(mopop,"moact.log",Pop::ACT);
}

void doprinting() {
  new PrjPrinter(semoprj,"semobj.bin",Prj::BJ);
  new PrjPrinter(semoprj,"semowij.bin",Prj::WIJ);
  new PrjPrinter(momoprj,"momobj.bin",Prj::BJ);
  new PrjPrinter(momoprj,"momowij.bin",Prj::WIJ);
}

enum BrainMode { RELAX, LEARN, RECALL, COMPLETE } ;

void setmode(BrainMode brainmode) {
    switch (brainmode) {
	break;
    case RELAX:
	semoprj->setprn(0.);
	momoprj->setprn(0.);
	momoprj->setbgain(0.); momoprj->setwegain(0.); momoprj->setwigain(0.);
	break;
    case LEARN:
	semoprj->setprn(1.);
	momoprj->setprn(1.);
	momoprj->setbgain(0.); momoprj->setwegain(0.); momoprj->setwigain(0.);
	break;
    case RECALL:
	semoprj->setprn(0.);
	momoprj->setprn(0.);
	momoprj->setbgain(momobgain); momoprj->setwegain(momowegain); momoprj->setwigain(momowigain);
    }
}

int main(int argc, char** argv) {
   
  Globals::init(argc,argv);
  Globals::setdt(0.001);

  parseparams(paramfile);

  sepop = new SensorPop("sensoryinput",musicdelay,nrank,H,U,taum); 
  sepop->setnormtype(Pop::NONE);
  sepop->setwtagain(1.);

  mopop = new MotorPop("motoroutput",musicdelay,nrank,H,U,taum); 
  mopop->setnormtype(Pop::NONE);
  mopop->setwtagain(wtagain);
  mopop->setiwgain(0.);
    
  mkdelay(tauzi,dmax,da,dn,dq);

  int n = tapsvec.size();

  //  momoprj = new Prj(mopop,mopop,tauzi,tauzj,taue,taup);
  momoprj = new PrjD(mopop,mopop,tapsvec,tauzivec,tauzj,taue,taup);
  momoprj->setselfconn(true,true,true);
  momoprj->reinit();

  semoprj = new Prj11(sepop,mopop,tauzi,tauzj,taue,taup);
  semoprj->reinit();

  if (ISROOT) std::cerr << "music simtime = " << Globals::_musicstoptime << std::endl;
  if (ISROOT) std::cerr << "music timestep = " << Globals::_musictimestep << std::endl;

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

  momoprj->setbgain(0.);
  momoprj->setwegain(0.);
  momoprj->setwigain(0.);
  momoprj->setprn(1.);

  semoprj->setbwgain(semobgain); semoprj->setwgain(semowgain);
  semoprj->setprn(1.);


  if (ISROOT) std::cerr << "music simtime = " << Globals::_musicstoptime << " sec" << std::endl;
  if (ISROOT) std::cerr << "music timestep = " << Globals::_musictimestep << " sec" << std::endl;
  idur = (Globals::_musictimestep+Globals::_dt)/Globals::_dt;

  Globals::setnldot(nldot);

  Globals::start();

  for (double simtime=0.; Globals::_musicruntime->time() < Globals::_musicstoptime;
       simtime+=Globals::_musictimestep) {

      if (simtime>0. and simtime<4.) setmode(LEARN);
      if (simtime>4. and simtime<7.) setmode(RECALL);
      if (simtime>7. and simtime<10.) setmode(RELAX);
      if (simtime>10. and simtime<17.) setmode(RECALL);
      if (simtime>17. and simtime<21.) setmode(LEARN);
      if (simtime>21. and simtime<24.) setmode(RECALL);
      if (simtime>24. and simtime<25.) setmode(RELAX);
      if (simtime>25. and simtime<100.) setmode(RECALL);

      for (int step=0; step<idur; step++) Globals::updstateall();
  }

  Globals::stop();

  MPI_Barrier(Globals::_comm_world);
    
  doprinting();

  Globals::report();
    
  Globals::finish();

  return 0;
}

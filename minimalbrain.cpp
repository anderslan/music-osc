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
#include "base/source/Prj.h"
#include "base/source/PrjD.h"
#include "base/source/Logger.h"
#include "utils/Utils.h"
#include "utils/RndGen.h"
#include "utils/Parseparam.h"
#include "utils/Patio.h"

using namespace std;

int nrank = 8,H = 8,U = 88/H,dn = 40,nstep = 100,idur = 50,nldot = 0;
float taum = 0.002,wtagain = 8.,musicdelay = 0.010,bgain = 1.,wegain = 1.,wigain = 14.,
    tauzi = 0.004,tauzj = 0.004,taue = 0.004,taup = 10.,noise = 0.,lgbias = 0.,
    dmax = 16,da = tauzi,dq = 2.;

string paramfile = "params_3.par";

SensoryMotorPop *semopop;
PrjD *semoprj;

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

    pp->postparam("tauzi",&tauzi,Parseparam::Float);
    pp->postparam("tauzj",&tauzj,Parseparam::Float);
    pp->postparam("taue",&taue,Parseparam::Float);
    pp->postparam("taup",&taup,Parseparam::Float);
    pp->postparam("bgain",&bgain,Parseparam::Float);
    pp->postparam("wegain",&wegain,Parseparam::Float);
    pp->postparam("wigain",&wigain,Parseparam::Float);
    pp->postparam("dn",&dn,Parseparam::Int);
    pp->postparam("da",&da,Parseparam::Float);
    pp->postparam("dmax",&dmax,Parseparam::Float);
    pp->postparam("dq",&dq,Parseparam::Float);
    pp->postparam("nldot",&nldot,Parseparam::Int);
    pp->postparam("nstep",&nstep,Parseparam::Int);

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
// 	printf("a = %.4f c = %.4f q = %.4f\n",a,c,q);
// 	for (int i=0; i<n; i++)
// 	    printf("%.4f %.4f\n",tapsvec[i]*Globals::_dt,tauzivec[i]);
//     }
}

void setuplogging() {
    new PopLogger(semopop,"semolgi.log",Pop::LGI);
    new PopLogger(semopop,"semowsu.log",Pop::WSU);
    new PopLogger(semopop,"semodsu.log",Pop::DSU);
    new PopLogger(semopop,"semoact.log",Pop::ACT);
    new PrjLogger(semoprj,"semoprn.log",Prj::PRN);
}

void doprinting() {
}

float calciwgain(float inent,float inlike) {
    float iwgain;
    if (ISROOT)
	fprintf(stderr,"%8d %15.6f %15.6f\n",Globals::_simstep,inent,inlike);
    //    if (inent>0) iwgain = 0.; else iwgain =  1. - exp(inlike);
    return iwgain;
}

int main(int argc, char** argv) {
   
  Globals::init(argc,argv);
  Globals::setdt(0.001);

<<<<<<< HEAD
  SensorPop* sensorpop = new SensorPop("sensorpool", 0.010, 1, 1, 88, 0.010, 0.0); 
  MotorPop* motorpop = new MotorPop("motorpool", 0.010, 1, 1, 88, 0.010, 0.0); 
    
  Globals::start();

=======
  parseparams(paramfile);

  semopop = new SensoryMotorPop("sensoryinput","motoroutput",musicdelay,nrank,H,U,taum); 
  semopop->setnormtype(Pop::NONE);
  semopop->setwtagain(wtagain);
    
  mkdelay(tauzi,dmax,da,dn,dq);

  int n = tapsvec.size();

  semoprj = new PrjD(semopop,semopop,tapsvec,tauzivec,tauzj,taue,taup);
  semoprj->setselfconn(true,true,true);
  semoprj->reinit();

  MPI_Barrier(Globals::_comm_world);
  if (ISROOT) {
      fprintf(stderr,"H = %d U = %d taum = %.3f ",H,U,taum);
      fprintf(stderr,"taup = %.3f sec bgain = %.4f wegain = %.4f wigain = %.4f n = %d\n",
	      taup,bgain,wegain,wigain,n);
  }
  MPI_Barrier(Globals::_comm_world);
  
  setuplogging();

  if (ISROOT) std::cerr << "Setup done!" << std::endl;

  semoprj->setwegain(wegain);
  semoprj->setwigain(wigain);

  float inent,inlike,inon,igain = 1.,wgain = 0.,prn;
  semopop->setigain(1.);
  semopop->setwgain(0.);
  semoprj->setprn(1.);

  if (ISROOT) std::cerr << "Running with music input " << std::endl;

  // Main simulation loop
  while (Globals::_musicruntime->time() < Globals::_musicstoptime)
    Globals::updstateall();

  if (ISROOT)
    std::cerr << "Done!" << std::endl;

  if (ISROOT) std::cerr << "music simtime = " << Globals::_musicstoptime << std::endl;
  if (ISROOT) std::cerr << "music timestep = " << Globals::_musictimestep << std::endl;
  idur = (Globals::_musictimestep+Globals::_dt)/Globals::_dt;
  if (ISROOT) std::cerr << "idur = " << idur << std::endl;

  Globals::setnldot(nldot);

  Globals::start();

  for (int simt=0; Globals::_musicruntime->time() < Globals::_musicstoptime; simt++) {
      inent = semopop->inentropy();
      inlike = semopop->inlikelihood();
      inon = semopop->inon();
      if (ISROOT)
	  fprintf(outf,"%8d %15.6f %15.6f %4d\n",Globals::_simstep,inent,inlike,inon);
      if (inon>0) {
	  wgain = 0.;
	  prn = 1.;
      } else {
	  wgain = 1.;
	  prn = 0.;
      }
      if (Globals::_simstep>nstep) { prn = 0.; igain = 0.; wgain = 1.; }
      MPI_Bcast(&prn,1,MPI_FLOAT,0,Globals::_comm_world);
      MPI_Bcast(&wgain,1,MPI_FLOAT,0,Globals::_comm_world);
      semopop->setigain(igain);
      semopop->setwgain(wgain);
      semoprj->setprn(prn);
      for (int step=0; step<idur; step++) Globals::updstateall();
  }

  Globals::stop();

  fclose(outf);
    
  MPI_Barrier(Globals::_comm_world);
    
  doprinting();

  Globals::report();
    
  Globals::finish();

  return 0;
}

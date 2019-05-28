#include<sys/types.h>
#include<dirent.h>
#include<unistd.h>

string findFilename(string folder, string time, string ch);
int findAllFields(string folder, string fields[100]);
void getFieldValues(string name, int &E1, int &E2, int &E3);

double getMaximum(TGraph *g1);
double getMinimum(TGraph *g1);
void getMeanAndERR(TFile *fin, double &mean, double &err, bool isCathode);

void checkGridTransparency(string infolder="/unix/dune/purity/CERN/2019/Liquid/PrM2/Day4/AllFibres_wFilters_newResistors/"){

  string outname="GridTransparency";

  string fields[100];
  int nfields = findAllFields(infolder, fields);

  int E1[100], E2[100], E3[100];

  double QK[100], QA[100], R[100], x[100];
  double QKerr[100], QAerr[100], Rerr[100], xerr[100];

  double R40[100], x40[100], R40err[100], x40err[100];
  int count40=0;
  

  for (int ifield=0; ifield<nfields; ifield++){
    //    cout << fields[ifield] << endl;
    getFieldValues(fields[ifield], E1[ifield], E2[ifield], E3[ifield]);

    string fnameA = infolder + "/" + findFilename(infolder, fields[ifield], "ch3");
    string fnameK = infolder + "/" + findFilename(infolder, fields[ifield], "ch4");

    cout << fnameA << endl;
    
    TFile *fA = new TFile(fnameA.c_str(), "read");
    getMeanAndERR(fA, QA[ifield], QAerr[ifield], false);
    // int nA     = gA->GetN();
    // double *xA = gA->GetX();
    // double *yA = gA->GetY();    
    // QA[ifield]=0;
    // for (int i=0; i<nA; i++){
    //   if (yA[i]>QA[ifield] && xA[i]>0.02E-3){
    // 	QA[ifield]=yA[i];
    //   }
    // }
    


    TFile *fK = new TFile(fnameK.c_str(), "read");
    getMeanAndERR(fK, QK[ifield], QKerr[ifield], true);
    // TGraph *gK = (TGraph*)fK->Get("justAvg");    
    // int nK     = gK->GetN();
    // double *xK = gK->GetX();
    // double *yK = gK->GetY();
    // QK[ifield]=0;
    // for (int i=0; i<nK; i++){
    //   if (yK[i]<QK[ifield] && xK[i]>0.02E-3){
    // 	QK[ifield]=yK[i];
    //   }
    // }


    R[ifield]  = -QA[ifield]/QK[ifield];
    x[ifield]  = E2[ifield]*1.0/E1[ifield];
    xerr[ifield] = 0;
    Rerr[ifield] = R[ifield]*TMath::Sqrt(TMath::Power(QAerr[ifield]/QA[ifield], 2)+TMath::Power(QKerr[ifield]/QK[ifield], 2));

    cout << E2[ifield] <<  " " << E1[ifield] << " " << x[ifield] << " " << R[ifield] << endl;

    if (E1[ifield]==40){
      R40[count40]=R[ifield];
      x40[count40]=x[ifield];
      R40err[count40]=Rerr[ifield];
      x40err[count40]=xerr[ifield];
      count40++;
    }

    delete fA;
    delete fK;



  }

  TCanvas *c = new TCanvas("c");
  //  TGraph *g = new TGraph (nfields, x, R);
  //  TGraph *g = new TGraph (count40, x40, R40);
  TGraphErrors *g = new TGraphErrors (count40, x40, R40, x40err, R40err);
  g->GetYaxis()->SetRangeUser(0, 1.5);
  g->SetMarkerStyle(20);
  g->SetTitle("Grid transparency check;E2/E1;-QA/QK");
  g->Draw("Ape");

  c->Print(Form("plots/%s.png", outname.c_str()));


}


void getFieldValues(string name, int &E1, int &E2, int &E3){

  string pattern1="Field_";
  string pattern2="Vcm";

  if (name.find(pattern1) != std::string::npos) {
    std::string::size_type i = name.find(pattern1);
    name.erase(i, pattern1.length());
    i = name.find(pattern2);
    name.erase(i);
    // cout << name << endl;
    sscanf(name.c_str(), "%d.%d.%d", &E1, &E2, &E3);
    // cout << E1 << " " << E2 << " "  << E3 << endl;
  }
}

string findFilename(string folder, string field, string ch){

  DIR *dp;
  dirent *d;
  int count = 0;
  int secs0, secs1;
  int h, m;
  string pattern1="."+ch+".traces_averages.root";
  if((dp = opendir(folder.c_str())) == NULL)
    perror("opendir");
  
  while((d = readdir(dp)) != NULL)
    {
      if(!strcmp(d->d_name,".") || !strcmp(d->d_name,".."))
	continue;

      std::string temps = d->d_name;
      if (temps.find(pattern1) != std::string::npos) {
	if (temps.find(field) != std::string::npos) {
	  //	  cout << temps << endl;
	  return temps;
	}
      } 
      
    }
  return "nothing";
  
}


int findAllFields(string folder, string fields[100]){

  DIR *dp;
  dirent *d;
  int count = 0;
  int secs0, secs1;
  int h, m;
  string pattern1=".ch4.traces_averages.root";
  string pattern2="FibreIn_";
  if((dp = opendir(folder.c_str())) == NULL)
    perror("opendir");
  
  while((d = readdir(dp)) != NULL)
    {
      if(!strcmp(d->d_name,".") || !strcmp(d->d_name,".."))
	continue;

      std::string temps = d->d_name;
      if (temps.find(pattern1) != std::string::npos) {
	std::string::size_type i = temps.find(pattern1);
	temps.erase(i, pattern1.length());
	fields[count] = temps;
	
	count++;
      } 
      
    }
  return count;
  
  
}


void getMeanAndERR(TFile *fin, double &mean, double &err, bool isCathode){

  double somek[10];
  mean=0.;
  err=0.;

  int ngraphs=10;
  for (int i=0; i<ngraphs; i++){
    TGraph *gt = (TGraph*)fin->Get(Form("avg100/gavg100_%i", i));
    if (isCathode) somek[i] = getMinimum(gt);
    else somek[i] = getMaximum(gt);
    delete gt;
    mean += somek[i]/ngraphs;
  }

  for (int i=0; i<ngraphs; i++){
    err += (somek[i]-mean)*(somek[i]-mean);
  }
  err = TMath::Sqrt(err)/ngraphs;  

}

double getMaximum(TGraph *g1){

  int n1     = g1->GetN();
  double *x1 = g1->GetX();
  double *y1 = g1->GetY();
  double k1=0.;
  for (int i=0; i<n1; i++){
    if (y1[i]>k1 && x1[i]>0.02E-3){
      k1=y1[i];
    }
  }

  return k1;
}


double getMinimum(TGraph *g1){

  int n1     = g1->GetN();
  double *x1 = g1->GetX();
  double *y1 = g1->GetY();
  double k1=0.;
  for (int i=0; i<n1; i++){
    if (y1[i]<k1 && x1[i]>0.02E-3){
      k1=y1[i];
    }
  }

  return k1;
}



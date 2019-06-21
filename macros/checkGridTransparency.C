#include<sys/types.h>
#include<dirent.h>
#include<unistd.h>

string findFilename(string folder, string time, string ch);
int findAllFields(string folder, string fields[100]);
void getFieldValues(string name, int &E1, int &E2, int &E3);

void checkGridTransparency(string infolder="/unix/dune/purity/CERN/2019/Liquid/PrM2/Day4/AllFibres_wFilters_newResistors/"){

  string outname="GridTransparency";

  string fields[100];
  int nfields = findAllFields(infolder, fields);

  int E1[100], E2[100], E3[100];

  double QK[100], QA[100], R[100], x[100];

  double R40[100], x40[100];
  int count40=0;
  

  for (int ifield=0; ifield<nfields; ifield++){
    //    cout << fields[ifield] << endl;
    getFieldValues(fields[ifield], E1[ifield], E2[ifield], E3[ifield]);

    string fnameA = infolder + "/" + findFilename(infolder, fields[ifield], "ch3");
    string fnameK = infolder + "/" + findFilename(infolder, fields[ifield], "ch4");

    cout << fnameA << endl;
    
    TFile *fA = new TFile(fnameA.c_str(), "read");
    TGraph *gA = (TGraph*)fA->Get("justAvg");

    int nA     = gA->GetN();
    double *xA = gA->GetX();
    double *yA = gA->GetY();
    
    QA[ifield]=0;
    for (int i=0; i<nA; i++){
      if (yA[i]>QA[ifield] && xA[i]>0.02E-3){
	QA[ifield]=yA[i];
      }
    }

    TFile *fK = new TFile(fnameK.c_str(), "read");
    TGraph *gK = (TGraph*)fK->Get("justAvg");
    
    int nK     = gK->GetN();
    double *xK = gK->GetX();
    double *yK = gK->GetY();
    QK[ifield]=0;
    for (int i=0; i<nK; i++){
      if (yK[i]<QK[ifield] && xK[i]>0.02E-3){
	QK[ifield]=yK[i];
      }
    }

    //    QA[ifield] = TMath::MaxElement(gA->GetN(), gA->GetY());
    //    QK[ifield] = TMath::MinElement(gK->GetN(), gK->GetY());
    R[ifield]  = -QA[ifield]/QK[ifield];
    x[ifield]  = E2[ifield]*1.0/E1[ifield];

    cout << E2[ifield] <<  " " << E1[ifield] << " " << x[ifield] << " " << R[ifield] << endl;

    if (E1[ifield]==40){
      R40[count40]=R[ifield];
      x40[count40]=x[ifield];
      count40++;
    }

    delete fA;
    delete fK;



  }

  TCanvas *c = new TCanvas("c");
  //  TGraph *g = new TGraph (nfields, x, R);
  TGraph *g = new TGraph (count40, x40, R40);
  g->GetYaxis()->SetRangeUser(0, 1.5);
  g->SetMarkerStyle(20);
  g->SetTitle("Grid transparency check;E2/E1;-QA/QK");
  g->Draw("Ap");

  c->Print(Form("%s.png", outname.c_str()));


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

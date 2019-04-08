#include<sys/types.h>
#include<dirent.h>
#include<unistd.h>

string findFilename(string folder, string time, string ch);
int findAllTimes(string folder, string times[100], int minutes[100]);


void getPhotodiodeVsTime(){

  //  string basename = "/unix/dune/purity/CERN/2019/Vacuum/PrM1/FinalConfiguration/Fibre_cleanSilvsTime/";
  string basename = "/unix/dune/purity/CERN/2019/Vacuum/PrM2/FinalConfiguration_newPreampD_vsTime/";

  cout << basename << endl;
  string times[100];
  int minutes[100];
  int ntimes = findAllTimes(basename, times, minutes);
  

  double peak[3][100], integral[3][100], timesSec[3][100], der[3][100];
  
  string channels[3] = {"ch1", "ch3", "ch4"};

  TGraph *gp[3], *gi[3];
  double mixmax;

  for (int ich=0; ich<3; ich++){
    for (int it=0; it<ntimes; it++){
      
      string filename = findFilename(basename, times[it], channels[ich]);
      
      TFile *fin = new TFile((basename+"/"+filename).c_str(), "read");
      TGraph *g = (TGraph*)fin->Get("justAvg");
      mixmax = TMath::Abs(TMath::MaxElement(g->GetN(), g->GetY()));
      if (TMath::Abs(TMath::MinElement(g->GetN(), g->GetY())) > mixmax) mixmax = TMath::Abs(TMath::MinElement(g->GetN(), g->GetY()));
      peak[ich][it] = mixmax;
      integral[ich][it] = TMath::Abs(g->Integral());
      timesSec[ich][it] = it*30;
      if (it>0) der[ich][it-1] = 1-(peak[ich][it]/peak[ich][it-1]);
      //    cout << times[ich][it] << " " << peak[it] << " " << integral[it] << endl;
      fin->Close();
      
    }
    gp[ich] = new TGraph(ntimes, timesSec[ich], peak[ich]);
    gi[ich] = new TGraph(ntimes-1, timesSec[ich], der[ich]);

  }

  gp[1]->SetLineColor(kBlue);
  gp[2]->SetLineColor(kRed);

  gi[1]->SetLineColor(kBlue);
  gi[2]->SetLineColor(kRed);
  
  TLegend *leg = new TLegend(0.6, 0.6, 0.89, 0.89);
  leg->AddEntry(gp[0], "Photodiode peak", "l");
  leg->AddEntry(gp[1], "Anode peak", "l");
  leg->AddEntry(gp[2], "Cathode peak", "l");

  TCanvas *cp = new TCanvas("cp");
  gp[0]->SetTitle("Voltage peak;Time [s];Peak [V]");
  gp[0]->GetYaxis()->SetRangeUser(0, 3.);
  gp[0]->Draw("Al");
  gp[1]->Draw("l");
  gp[2]->Draw("l");
  leg->Draw();
  cp->Print("VoltagePeak.png");


  TCanvas *ci = new TCanvas("ci");
  gi[0]->SetTitle("Voltage peak relative ratio;Time [s];Relative ratio");  gi[0]->Draw("Al");
  gi[1]->Draw("l");
  gi[2]->Draw("l");
  leg->Draw();
  ci->Print("VoltagePeak_FractionalIncrease.png");

}



string findFilename(string folder, string time, string ch){

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
	if (temps.find(time) != std::string::npos) {
	  //	  cout << temps << endl;
	  return temps;
	}
      } 
      
    }
  return "nothing";
  
}



int findAllTimes(string folder, string times[100], int minutes[100]){

  DIR *dp;
  dirent *d;
  int count = 0;
  int secs0, secs1;
  int h, m;
  string pattern1=".ch1.traces_averages.root";
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
	cout << " " << __LINE__ << " " << temps << endl;
	times[count] = temps;
	if (count==0) minutes[count]=0;
	else {
	  sscanf(times[count].c_str(), "%d.%d", &h, &m);
	  secs0 = h*60 + m;
	  sscanf(times[count-1].c_str(), "%d.%d", &h, &m);
	  secs1 = h*60 + m;
	  minutes[count] = secs0-secs1 + minutes[count-1];
	  cout << times[count] << " " << minutes[count] << endl;
	}
	count++;
      } 
      
    }
  return count;
  
  
}


#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

int findFile(string folder, int ch, string &filename);

int getField(string folder, string &field);

void quickPlotAverages(string folder="/data/PurityMonitor/GasTests/Run040/", int whichPrM=1){


  if (whichPrM!=1 && whichPrM!=2) {
    cout << "Which PrM are you using? it's either 1 or 2!" << endl;
    return;
  }

  int channels[2];
  if (whichPrM==1){
    channels[0] = 1;
    channels[1] = 2;
  } else {
    channels[0] = 4;
    channels[1] = 5;
  }

  TGraph *g[2];

  double min=0.;
  double max=0.;
  double tmin[2] = {0., 0.};
  double tmax[2] = {0., 0.};

  // Loop over cathode and anode
  for (int ich=0; ich<2; ich++){

    string filename;
    int ret = findFile(folder, channels[ich], filename);
    cout << filename << endl;
    TFile *fin = new TFile(Form("%s/%s", folder.c_str(), filename.c_str()));

    g[ich] = (TGraph*)fin->Get("justAvg");

    tmax[ich] = TMath::MaxElement(g[ich]->GetN(), g[ich]->GetY());
    tmin[ich] = TMath::MinElement(g[ich]->GetN(), g[ich]->GetY());
    

    cout << tmax[ich] << " " << tmin[ich] << endl;

    delete fin;

  }

  max = tmax[0] ? tmax[1] : (tmax[0] > tmax[1]);

  if (tmin[0]<tmin[1]) min = tmin[0];
  else min = tmin[1];

  //  min = tmin[0] ? tmin[1] : (tmin[0] < tmin[1]);

  string field;

  getField(folder, field);
  TCanvas *c = new TCanvas("c");

  g[0]->SetTitle(Form("%s Vcm;Time [ns];Amplitude [mV]", field.c_str()));

  cout << min << " " << max << endl;

  g[0]->GetYaxis()->SetRangeUser(min*1.1, max*1.1);

  g[0]->Draw("Al");

  g[1]->Draw("l");

  c->Print(Form("%s/JustAverages_PrM%d.png", folder.c_str(), whichPrM));
  c->Print(Form("%s/JustAverages_PrM%d.pdf", folder.c_str(), whichPrM));
  


}




int findFile(string folder, int ch, string &filename){

  DIR *dp;
  dirent *d;
  int count = 0;
  string pattern1=Form("ch%d", ch);
  if((dp = opendir(folder.c_str())) == NULL)
    perror("opendir");
  
  while((d = readdir(dp)) != NULL)
    {
      if(!strcmp(d->d_name,".") || !strcmp(d->d_name,".."))
	continue;
      
      std::string temps = d->d_name;
      //      cout << temps << " " << pattern1 << endl;
      if (temps.find(pattern1) != std::string::npos) {
	filename += temps;
	//	cout << "ONE WAY " << filename << endl;
	return 1;
      } 
      
    }
  return 0;
  
  
}



int getField(string folder, string &field){

  FILE *fp;
  fp = fopen(Form("%s/Elog.txt", folder.c_str()), "r");
  
  char str[1000];
  fgets(str, 10000, fp);

  string s = Form("%s", str);
  cout << s << endl;

  std::string ext("Vcm");

  if (s.find(ext) != std::string::npos) {

    int place = s.find(ext);
      // if so then strip them off                                                                                                        
      s = s.substr(0, place);
      //      cout << s << endl;
    }

  //  cout << s << endl;

  field+=s;

  fclose(fp);

  return 1;
}

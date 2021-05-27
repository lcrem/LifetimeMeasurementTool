#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <ctime>
#include <string>

#include "TSystem.h"

using namespace std;

int ID;

void strreplace(string &str, string oldStr, string newStr);

void getSsoCookie(int ID);

int main(){

  getSsoCookie(ID);    

  string runSummary = " ";
  time_t now = time(0);
  string date = ctime(&now);

  string datadir1="/data/PurityMonitor/Filling/";
  string datadir2="/data/PurityMonitor/Filling/";
  string datadir3="/data/PurityMonitor/Filling/";
   
  int run=0;
  FILE *frun;
  if (( frun = fopen(Form("%s/LastRun", datadir1.c_str()), "r")) ){
    fscanf(frun,"%d", &run);
    printf("The last run number was %d\n", run);
    fclose(frun);
  }


  datadir1+=Form("Run%03d", run-2);
  datadir2+=Form("Run%03d", run-1);
  datadir3+=Form("Run%03d", run);

  // cout << datadir1 << endl;
  // cout << datadir2 << endl;
  cout << datadir3 << endl;

  // FILE *ftime1 = fopen(Form("%s/Timestamp.txt", datadir1.c_str()), "r");
  // int inttime1;
  // fscanf(ftime1, "%d", &inttime1);
  // cout << Form("%s/Timestamp.txt", datadir1.c_str())  << " " << inttime1 << endl;
  // std::time_t timestamp1 = inttime1;
  // string dateTaken1 = ctime(&timestamp1);
  // fclose(ftime1);
  // runSummary += Form("Purity monitor run %d on %s \n", run-2, dateTaken1.c_str());

  // std::ifstream fprm1_1(Form("%s/PrM1_avgLifeInfo.txt", datadir1.c_str()));
  // std::string prm1info_1( (std::istreambuf_iterator<char>(fprm1_1) ),
  //                      (std::istreambuf_iterator<char>()    ) );

  // std::ifstream fprm2_1(Form("%s/PrM2_avgLifeInfo.txt", datadir1.c_str()));
  // std::string prm2info_1( (std::istreambuf_iterator<char>(fprm2_1) ),
  //                      (std::istreambuf_iterator<char>()    ) );

  
  // std::ifstream felog1(Form("%s/Elog.txt", datadir1.c_str()));
  // std::string elog1( (std::istreambuf_iterator<char>(felog1) ),
  // 			(std::istreambuf_iterator<char>()    ) );

  // runSummary +=  prm1info_1 + " \n" + prm2info_1 + " \n";


  FILE *ftime2 = fopen(Form("%s/Timestamp.txt", datadir2.c_str()), "r");
  int inttime2;
  fscanf(ftime2, "%d", &inttime2);
  cout << Form("%s/Timestamp.txt", datadir2.c_str())  << " " << inttime2 << endl;
  std::time_t timestamp2 = inttime2;
  string dateTaken2 = ctime(&timestamp2);
  fclose(ftime2);
  // runSummary += Form("Purity monitor run %d on %s \n", run-1, dateTaken2.c_str());

  std::ifstream fprm1_2(Form("%s/PrM1_avgLifeInfo.txt", datadir2.c_str()));
  std::string prm1info_2( (std::istreambuf_iterator<char>(fprm1_2) ),
                       (std::istreambuf_iterator<char>()    ) );

  std::ifstream fprm2_2(Form("%s/PrM2_avgLifeInfo.txt", datadir2.c_str()));
  std::string prm2info_2( (std::istreambuf_iterator<char>(fprm2_2) ),
                       (std::istreambuf_iterator<char>()    ) );

  
  std::ifstream felog2(Form("%s/Elog.txt", datadir2.c_str()));
  std::string elog2( (std::istreambuf_iterator<char>(felog2) ),
			(std::istreambuf_iterator<char>()    ) );

  // runSummary += prm1info_2 + " \n" + prm2info_2 + "\n" ;

  FILE *ftime3 = fopen(Form("%s/Timestamp.txt", datadir3.c_str()), "r");
  int inttime3;
  fscanf(ftime3, "%d", &inttime3);
  cout << Form("%s/Timestamp.txt", datadir3.c_str())  << " " << inttime3 << endl;
  std::time_t timestamp3 = inttime3;
  string dateTaken3 = ctime(&timestamp3);
  fclose(ftime3);
  runSummary += Form("Purity monitor run %d on %s \n", run, dateTaken3.c_str());

  std::ifstream fprm1_3(Form("%s/PrM1_avgLifeInfo.txt", datadir3.c_str()));
  std::string prm1info_3( (std::istreambuf_iterator<char>(fprm1_3) ),
                       (std::istreambuf_iterator<char>()    ) );

  std::ifstream fprm2_3(Form("%s/PrM2_avgLifeInfo.txt", datadir3.c_str()));
  std::string prm2info_3( (std::istreambuf_iterator<char>(fprm2_3) ),
                       (std::istreambuf_iterator<char>()    ) );

  
  std::ifstream felog3(Form("%s/Elog.txt", datadir3.c_str()));
  std::string elog3( (std::istreambuf_iterator<char>(felog3) ),
			(std::istreambuf_iterator<char>()    ) );

  runSummary +=  prm1info_3 + " \n" + prm2info_3;

  cout << runSummary << endl;

  // strreplace(runSummary, "<", "&#60;");
  // strreplace(runSummary, ">", "&#62;");
  // strreplace(runSummary, "&", "&#38;");
  // strreplace(runSummary, "'", "&#39;");
  // strreplace(runSummary, "\"", "&#34;");
  // strreplace(runSummary, "\n", "&#x0d; &#x0a;");
 strreplace(runSummary, "\n", "\n \n");
  cout << runSummary << endl;


  string stdoutput,stderror;
  int returnCode;

  string IDstring = to_string(ID);
  // DebugTN("------------------- Uploading the summary into the e-log -----------------------");
  //   
  //  returnCode = system( ("curl -s -H \"Content-Type: application/xml\" --cookie /home/lindac/credentials/.secrets/ssocookie"+IDstring+".txt --cookie-jar /home/lindac/credentials/.secrets/ssocookie"+IDstring+".txt -X POST --data-ascii '<input_message><author>Purity Monitors DAQ pc</author><subject>Purity Monitor Measurements at "+date+"</subject><message_type>Automatic</message_type><systems_affected><count>1</count><system_affected>Purity Monitors</system_affected></systems_affected><entry_type>Information</entry_type><body>"+runSummary+"</body></input_message>' https://pddpelog.web.cern.ch/elisa.api/api/ProtoDUNE-DP/messages").c_str()); //,stdoutput,stderror);

  string subject = "Purity Monitor Measurements at "+date;
  // string attachments1 = Form("%s/PurityMonitor1_filAvg_subNoise.png", datadir1.c_str());
  // string attachments2 = Form("%s/PurityMonitor2_filAvg_subNoise.png", datadir1.c_str());
  // string attachments3 = Form("%s/PurityMonitor1_filAvg_subNoise.png", datadir2.c_str());
  // string attachments4 = Form("%s/PurityMonitor2_filAvg_subNoise.png", datadir2.c_str());
  string attachments5 = Form("%s/PurityMonitor1_filAvg_subNoise.png", datadir3.c_str());
  string attachments6 = Form("%s/PurityMonitor2_filAvg_subNoise.png", datadir3.c_str());

  returnCode = system( ("python $HOME/python/elisa_client_api/elisa_insert.py -o /home/lindac/credentials/.secrets/ssocookie"+IDstring+".txt -s https://pddp-elog.cern.ch -j \""+subject+"\" -y Automatic -b \""+runSummary+"\" -a \"Purity Monitors DAQ pc\" -e \"Purity Monitors\" -m \""+attachments5+"\" -m \""+attachments6+"\" -k \"ProtoDUNE-DP\"").c_str());

  //  returnCode = system( ("python $HOME/python/elisa_client_api/elisa_insert.py -o /home/lindac/credentials/.secrets/ssocookie"+IDstring+".txt -s https://pddpelog.web.cern.ch -j \""+subject+"\" -y Automatic -b \""+runSummary+"\" -a \"Purity Monitors DAQ pc\" -e \"Purity Monitors\" -m \"plots/OverlayLatestLifetime.png\"  -m \"plots/OverlayLatestOverLifetime.png\" -m \""+attachments3+"\"  -m \""+attachments4+"\" -m \""+attachments5+"\"  -m \""+attachments6+"\" -k \"ProtoDUNE-DP\"").c_str());

  //curl -s -H \"Content-Type: application/xml\" --cookie /home/lindac/credentials/.secrets/ssocookie"+IDstring+".txt --cookie-jar /home/lindac/credentials/.secrets/ssocookie"+IDstring+".txt -X POST --data-ascii '<input_message><author>Purity Monitors DAQ pc</author><subject>Purity Monitor Measurements at "+date+"</subject><message_type>Automatic</message_type><systems_affected><count>1</count><system_affected>Purity Monitors</system_affected></systems_affected><entry_type>Information</entry_type><body>"+runSummary+"</body></input_message>' https://pddpelog.web.cern.ch/elisa.api/api/ProtoDUNE-DP/messages").c_str()); //,stdoutput,stderror);

  // DebugN("This is the standard output, ", stdoutput);
  // DebugN("This is the standard error, ", stderror);

  if(returnCode)
    {
      // DebugTN("Return code: " + returnCode);

      
    }
}


void getSsoCookie(int ID)
{
  string secretsPath = "/home/lindac/credentials/.secrets/";
  //  string cmd = "kinit -kt /home/lindac/credentials/np02prm.keytab np02prm";
  string cmd = "kinit -kt /home/lindac/credentials/lcremone.keytab lcremone@CERN.CH";
  system(cmd.c_str());

  string IDstring = to_string(ID);
  system(("rm -f "+secretsPath+"ssocookie"+IDstring+".txt").c_str());
  //  system(("cern-get-sso-cookie --krb -r -u   https://pddpelog.web.cern.ch/elisa.api -o "+secretsPath+"ssocookie"+IDstring+".txt").c_str());
  system(("cern-get-sso-cookie --krb -r -u  https://pddp-elog.cern.ch/elisa -o "+secretsPath+"ssocookie"+IDstring+".txt").c_str());
}


void strreplace(string &line, string oldString, string newString){
  
  const size_t oldSize = oldString.length();

  // do nothing if line is shorter than the string to find
  if( oldSize > line.length() ) return;

  const size_t newSize = newString.length();
  for( size_t pos = 0; ; pos += newSize ) {
    // Locate the substring to replace
    pos = line.find( oldString, pos );
    if( pos == string::npos ) return;
    if( oldSize == newSize ) {
      // if they're same size, use std::string::replace
      line.replace( pos, oldSize, newString );
    } else {
      // if not same size, replace by erasing and inserting
      line.erase( pos, oldSize );
      line.insert( pos, newString );
    }
  }

  // std::string::size_type pos = 0u;
  // while((pos = str.find(oldStr, pos)) != std::string::npos){
  //   str.replace(pos, oldStr.length(), newStr);
  //   pos += newStr.length();
  // }
}

#include"libs.hh"
#include <TMathBase.h>
#define LENGTH(x) (sizeof x / sizeof *x)
#define NDETS LENGTH(foot_id)
using namespace std;
namespace{//don't change the order!!!
    int foot_id[] = {15, 16, 1, 2, 13, 4, 11, 6, 7, 12, 9, 10};
}
Double_t pedestal[NDETS][640];
Double_t sigma[NDETS][640];
const int Nbins = 2000; const double ymin = -1000; const double ymax = 2000;
const double THRESHOLD=10;
typedef std::pair<int, double> strip_data;
typedef std::vector<strip_data> cluster; 
typedef std::vector<cluster> foot_data; 

int     event_monitor = 1;//to display events one by one

bool is_good_strip(UInt_t det, UInt_t strip)
{
    switch(det){
        case 1:
            if(
                    strip==449 || strip==456 || strip==507 || strip==509 || 
                    strip==510 || strip==508 || strip==511 || strip==512 ||
                    strip==94  || strip==506 || strip==513 || strip==514 ||
                    strip==94  || strip==64  || strip==65  || strip==129 ||
                    strip==257 || strip==526 || strip==577 || strip==582 ||
                    strip==527 || strip==529 || strip==528 || strip==530 ||
                    strip==531 || strip==532 || strip==533 || strip==534 ||
                    strip==93  || strip==629 || strip==630 || strip==627 ||
                    strip==632 || strip==630 || strip==450 || strip==455 ||
                    strip>450 
              ) return false;
        case 2:
            if(
                    strip==343 || strip==229 || strip ==230
              ) return false;
        case 7:
            if(
                    strip==597 || strip==598 || strip==344 || strip==345 ||
                    strip==304 || strip==305 || strip==596 || strip==599
              ) return false;
        case 11:
            if(
                    strip==30 || strip==31 || strip==32 || strip==33
              ) return false;
        case 12:
            if(
                    strip==232
              ) return false;
        case 15:
            if(
                    strip==204 || strip==203
              ) return false;
    }
    if((strip%64)>62 || (strip%64)<3) return false;
    return true;
}

void Check_Strip(UInt_t number, double energy, foot_data& fdata)
{
    //cout << "\n\n---------------------------------------";
    //cout << "\nChecking strip: " << number << "\t Energy: " << energy << endl;
    strip_data strip = std::make_pair(number,energy);
    cluster clust;
    if(fdata.size()==0)//no cluster yet, create new
    {
        clust.push_back(strip);
        fdata.push_back(clust);
        //cout << "\n\t New cluster is created for this strip";
        return;
    }
    cluster    this_clust = fdata.back();
    strip_data this_strip = this_clust.back();
    if(abs(strip.first-this_strip.first)<2)//neighbour found 
    {
        //cout << "\n\tStrip belong to exisitng cluster! Adding it...";
        fdata.back().push_back(strip);
        return;
    }
    else
    {
        //cout << "\n\tStrip is a new cluster! Making it...";
        clust.clear();
        clust.push_back(strip);
        fdata.push_back(clust);
        return;
    }
}

void FOOT_ana(int firstEvent, int max_events, TString filename = "")
{
    TApplication* theApp = new TApplication("App", 0, 0);
    TFile *datafile = new TFile(filename,"read");
    if(!datafile){ cout << "\n--ERROR in in the input file name!";  return; }

    //---------- Definition of main  histograms ----------------
    TH2F * h2_peds_raw[NDETS];//container of raw pedestal data 
    TH2F * h2_cal_coarse[NDETS];
    TH2F * h2_cal_fine[NDETS];
    TH1F * h1_signal_sum[NDETS];
    TH1F * h1_single_event[NDETS];
    TH1D * h1_peds[NDETS];
    TH1D * h1_sigma[NDETS];
    for(int i=0; i<NDETS; i++){
        h2_cal_coarse[i]  = new TH2F(Form("h2_cal_coarse_FOOT%d",foot_id[i]),Form("h2_cal_coarse_FOOT%d",foot_id[i]),640,1,641,Nbins,ymin,ymax);
        h2_cal_fine[i]    = new TH2F(Form("h2_cal_fine_FOOT%d",foot_id[i]),Form("h2_cal_fine_FOOT%d",foot_id[i]),640,1,641,Nbins,ymin,ymax);
        h2_peds_raw[i]    = new TH2F(Form("h%d",foot_id[i]),Form("h%d",foot_id[i]),640,1,641,1500,0,1500);
        h1_signal_sum[i]  = new TH1F(Form("h_signal_sum_FOOT%d",foot_id[i]),Form("hsignal_sum_FOOT%d",foot_id[i]),80,0,40);

        h1_single_event[i]= new TH1F(Form("h1_event_FOOT%d",foot_id[i]),Form("h1_event_FOOT%d",foot_id[i]),640,1,641);
        h1_single_event[i]->SetBarWidth(1);
        if(i%2) h1_single_event[i]->SetFillColor(kBlue);
        else h1_single_event[i]->SetFillColor(kRed);
    }
    TH2I * h2_mul1_vs_mul15 = new TH2I("h2_mul1_vs_mul15","h2_mul1_vs_mul15",20,0,20,20,0,20);
    TH2I * h2_mul2_vs_mul16 = new TH2I("h2_mul2_vs_mul16","h2_mul2_vs_mul16",20,0,20,20,0,20);
    TCanvas * canvas_single_event = new TCanvas("canvas_single_event","canvas_single_event",1200,900);
    canvas_single_event->Divide(4,3);

    //-------- Define tree data for all foot detectors ----------
    TTree *tree = (TTree*)datafile->Get("h101");
    UInt_t TRIGGER;
    UInt_t FOOT[NDETS];//numbers of foot detectors 
    UInt_t FOOTE[NDETS][640];
    UInt_t FOOTI[NDETS][640];
    UInt_t TPATv[1];
    tree->SetBranchAddress("TRIGGER",&TRIGGER);
    tree->SetBranchAddress("TPATv",TPATv);
    for(int i=0; i<NDETS; i++){
        TString bname   = Form("FOOT%d",  foot_id[i]);
        TString bname_E = Form("FOOT%dE", foot_id[i]);
        TString bname_I = Form("FOOT%dI", foot_id[i]);
        tree->SetBranchAddress(bname.Data(),&FOOT[i]);
        tree->SetBranchAddress(bname_E.Data(),FOOTE[i]);
        tree->SetBranchAddress(bname_I.Data(),FOOTI[i]);
    }
    //-------- Get starting set of pedestal data for all dets -------
    int Nevents = tree->GetEntries();
    if(max_events>0) Nevents = max_events; 
    int     stat=0;
    for(int ev=firstEvent; ev<firstEvent+Nevents; ev++)
    {
        tree->GetEntry(ev);
        //if((TPATv[0] & 4) !=4 ) continue;
        if((TPATv[0] & 8)==8) continue;//p2pv
        //if((TPATv[0] & 16)==16) continue;//OR
        //if(TRIGGER!=3) continue;
        for(int f=0; f<NDETS; f++){
            for(int  j=0 ; j<640 ; j++){
                h2_peds_raw[f]->Fill(FOOTI[f][j],FOOTE[f][j]); 
            }
        }
        stat++;
        //if(stat==1000) break;
    }
    //---------  Slicing, fitting, saving pedestals params ---------
    TF1* foo = new TF1("foo","gaus",0,1000);
    for(int i=0; i<NDETS; i++)
    {
        h2_peds_raw[i]->FitSlicesY(foo,1,640,0,"QNR",0);
        h1_peds[i]   = (TH1D*)gDirectory->Get(Form("h%d_1",foot_id[i]))->Clone(Form("h1_peds_%d",  foot_id[i]));
        h1_sigma[i]  = (TH1D*)gDirectory->Get(Form("h%d_2",foot_id[i]))->Clone(Form("h1_sigma_%d",foot_id[i]));
        for(int j=0; j<640; j++){
            pedestal[i][j] = h1_peds[i]->GetBinContent(j+1);
            sigma[i][j]    = h1_sigma[i]->GetBinContent(i+1);
            cout << "\nFOOT: " << foot_id[i] << "\tStrip: " << j << "\t Ped: " << pedestal[i][j] << "\tSig: " << sigma[i][j];
        }
    }
    //--------- Analyzing  data ------------
    double  mean_ssd = 0;
    double  asic_offset[10];
    double  signal = 0;
    double  signal_sum = 0;
    int     counter_asic =0;

    foot_data data_left_arm[4];//collection of clusters from all foots
    foot_data data_right_arm[4];//collection of clusters from all foots
    foot_data data_central_arm[4];//collection of clusters from all foots

    for(int ev=firstEvent; ev<firstEvent+Nevents; ev++)
    {
        cout << "\r-- Event # : " << ev << flush;
        tree->GetEntry(ev);
        //if((TPATv[0] & 4)==4) continue;//p2pv
        if((TPATv[0] & 8)==8) continue;//p2pv
        //if((TPATv[0] & 16)==16) continue;//OR
        for(int f=0; f<NDETS; f++)//loop over all foots
        {
            if(f<4) data_central_arm[f].clear();
            else if(f>3 && f<8) data_right_arm[f-4].clear();
            else data_left_arm[f-8].clear(); 

            if(event_monitor) h1_single_event[f]->Reset("ICESM");//reset individual events 

            //--------  Global base line correction in every FOOT in this event ---------
            mean_ssd=0; stat=0;
            for(int i=0; i<640; i++){
                if(!is_good_strip(foot_id[f],FOOTI[f][i])) continue; 
                signal = FOOTE[f][i] - pedestal[f][i];
                stat++;    mean_ssd += signal;
            }
            mean_ssd = mean_ssd/stat;
            if(mean_ssd>20){
                cout << "\n--[WARNING]: In FOOT " << foot_id[f] << "Mean ssd = " << mean_ssd << endl;
                //continue;
            }
            //------------ Calculating fine baseline correction for individual asics ---------
            stat=0; counter_asic=0;
            for(int i=0; i<10; i++){  asic_offset[i]=0.; }//reset asic baselines
            for(int i=0; i<640; i++)
            {
                signal = FOOTE[f][i] - pedestal[f][i] - mean_ssd;
                if(fabs(signal) < (3 * sigma[f][i]) && //ignore possible hit candidates
                        is_good_strip(foot_id[f],FOOTI[f][i]))//and bad strips
                {
                    stat++;
                    asic_offset[counter_asic] += signal;
                }
                if((FOOTI[f][i]%64)==0){//switch to next asic
                    asic_offset[counter_asic] /= stat;
                    counter_asic++;  stat=0;
                }
            }
            //-------- Apply baseline correction,fill histograms and cluster data ------------
            counter_asic=0; signal_sum=0;
            for(int i=0; i<640; i++)
            {
                if((FOOTI[f][i]%64) == 1 && FOOTI[f][i]>1) counter_asic++;
                if(!is_good_strip(foot_id[f],FOOTI[f][i])) continue; 
                signal = FOOTE[f][i] - pedestal[f][i] - mean_ssd - asic_offset[counter_asic];
                h2_cal_coarse[f]->Fill(FOOTI[f][i], (FOOTE[f][i] - pedestal[f][i]));
                h2_cal_fine[f]->Fill(FOOTI[f][i], signal);
                if(event_monitor) h1_single_event[f]->Fill(FOOTI[f][i], signal);
                signal_sum+=signal;
                if(signal>THRESHOLD) 
                {
                    if(f<4) Check_Strip(FOOTI[f][i], signal, data_central_arm[f]);
                    else if(f>3 && f<8) Check_Strip(FOOTI[f][i], signal, data_right_arm[f-4]);
                    else Check_Strip(FOOTI[f][i], signal, data_left_arm[f-8]);
                }
            }
            h1_signal_sum[f]->Fill(signal_sum);
        }//end of loop over detectors
        if( 
                data_left_arm[0].size() ==0 || data_left_arm[1].size() ==0 ||
                data_left_arm[2].size() ==0 || data_left_arm[3].size() !=1 ||
                data_right_arm[0].size()==0 || data_right_arm[1].size()==0 ||
                data_right_arm[2].size()==0 || data_right_arm[3].size()!=1)
            continue;
        
        h2_mul1_vs_mul15->Fill(data_central_arm[0].size(), data_central_arm[2].size());
        h2_mul2_vs_mul16->Fill(data_central_arm[1].size(), data_central_arm[3].size());

        if(event_monitor) //plot individual events
        {
            double ymin_new = ymin; double ymax_new=ymax;
            for(int f=0; f<NDETS; f++)
            {
                canvas_single_event->cd(f+1);
                if(foot_id[f]!=1 && foot_id[f]!=2 && foot_id[f]!=15 && foot_id[f]!=16 ){
                    ymin_new = -50.; ymax_new = 100.;
                }
                h1_single_event[f]->GetYaxis()->SetRangeUser(ymin_new, ymax_new);
                h1_single_event[f]->Draw("HISTO B");
                TText *t = new TText(320,10,Form("event: %d",ev));
                t->Draw();
                for(int i_asic=1; i_asic<10; i_asic++)
                {
                    TLine* l = new TLine(64.5*i_asic,ymin_new,64.5*i_asic,ymax_new);
                    l->Draw("same");
                    l->SetLineStyle(7);
                    l->SetLineWidth(1);
                    l->SetLineColor(13);
                }
                gPad->Update();
                canvas_single_event->Update();
            }
            theApp->Run(kTRUE);
        }
    }//end of eventloop
    TCanvas* canvas_coarse = new TCanvas("coarse","coarse",1200,900);
    TCanvas* canvas_fine   = new TCanvas("fine","fine",1200,900);
    canvas_coarse->Divide(4,3);
    canvas_fine->Divide(4,3);
    for(int f=0; f<NDETS; f++)
    {
        canvas_coarse->cd(f+1);
        gPad->SetLogz();
        h2_cal_coarse[f]->Draw("colz");

        canvas_fine->cd(f+1);
        gPad->SetLogz();
        h2_cal_fine[f]->Draw("colz");
    }
    TCanvas* canvas_mul = new TCanvas("canvas_mul","canvas_mul",1000,500);
    canvas_mul->Divide(2,1);
    canvas_mul->cd(1);
    gPad->SetLogz();
    h2_mul1_vs_mul15->Draw("colz");
    canvas_mul->cd(2);
    gPad->SetLogz();
    h2_mul2_vs_mul16->Draw("colz");

    TCanvas* canvas_signal_sum = new TCanvas("canvas_signal_sum","canvas_signal_sum",1200,900);
    canvas_signal_sum->Divide(4,3);
    for(int f=0; f<NDETS; f++)
    {
        gPad->SetLogz();
        canvas_signal_sum->cd(f+1);
        h1_signal_sum[f]->Draw();
    }
    theApp->Run();
    return;
}

int main(Int_t argc, Char_t* argv[])
{
    gRandom = new TRandom3();
    gRandom->SetSeed(0);
    gROOT->Macro("rootlogon.C");
    gStyle->SetPalette(kRainBow);
    //TString filename = Form("/u/land/r3broot/202205_s522/foot/macros/util/rootfiles/physics_run130.root");
    //TString filename = Form("../main0110_0001.root");
    TString filename = Form("../main0131_0041.root");
    FOOT_ana(100,1e5,filename);
    //FOOT_ana(0,-1,filename);
    return 0;
}

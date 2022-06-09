#include"libs.hh"
#include <TMathBase.h>
#define LENGTH(x) (sizeof x / sizeof *x)
#define NDETS LENGTH(foot_id)
using namespace std;

namespace{//don't change the order!!!
    int foot_id[] = {15, 16, 1, 2, 13, 4, 11, 6, 7, 12, 9, 10};
}
//Histogram binning and ranges
int Nbins = 2000; double ymin = -1000; double ymax = 2000;

//Clustering data structures
typedef std::pair<int, double> strip_data;
typedef std::vector<strip_data> cluster; 
typedef std::vector<cluster> foot_data; 

const double    FOOT_LENGTH   = 96.;//mm
const double    THRESHOLD     = 10;
const double    MIN_DISTANCE  = 1;//mm
const double    MAX_THRESHOLD = 80;
const int       EVENT_MONITOR = 0;//run with dsiplaying single events

struct track{ TVector3 point[4]; };

std::vector<track> tracks_left; 
std::vector<track> tracks_right; 


struct vertex
{
    Double_t dca;
    TVector3 vertexPos;
    TVector3 normalVec;
};

std::vector<vertex> vertex_array;

vertex get_vertex(TVector3 pos1, TVector3 pos2, TVector3 slope1, TVector3 slope2)
{
    Double_t t = 100000;
    Double_t v = 100000;
    TVector3 normalVec = slope2.Cross(slope1).Unit();

    TVector3 posDifference = pos2 - pos1;
    t = (-posDifference.Dot(slope2) + posDifference.Dot(slope1) * slope2.Dot(slope1)) / (1 - (slope2.Dot(slope1)) * (slope2.Dot(slope1)));
    v = (posDifference.Dot(slope1) - (posDifference.Dot(slope2) * (slope2.Dot(slope1)))) / (1 - (slope2.Dot(slope1)) * (slope2.Dot(slope1)));

    TVector3 vertex1Global = pos2 + t * slope2;
    TVector3 vertex2Global = pos1 + v * slope1;

    TVector3 Vertex_calc;
    TVector3 Vertex_dir_helper;

    Vertex_dir_helper = vertex2Global - vertex1Global;
    Vertex_dir_helper.SetX(Vertex_dir_helper.X() / 2);
    Vertex_dir_helper.SetY(Vertex_dir_helper.Y() / 2);
    Vertex_dir_helper.SetZ(Vertex_dir_helper.Z() / 2);

    Vertex_calc = vertex1Global + Vertex_dir_helper;
    vertex recVertex;
    recVertex.vertexPos = Vertex_calc;
    recVertex.dca = Vertex_dir_helper.Mag();
    recVertex.normalVec = normalVec;

    return recVertex;
}

//list of bad or dead strips
bool is_good_strip(UInt_t det, UInt_t strip)
{
    switch(det){
        case 1:
            if(
                    strip==449 || strip==456 || strip==507 || strip==509 || 
                    strip==510 || strip==508 || strip==511 || strip==512 ||
                    strip==94  || strip==506 || strip==513 || strip==514 ||
                    strip==64  || strip==65  || strip==129 ||
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

//Functions for clustering
double get_cog(cluster c)//calculate center of gravity of a cluster
{
    double value = 0;  double esum = 0.;
    for(auto & s: c){
        value  += (s.first * s.second);
        esum += s.second;
    }
    value /= esum;//center of gravity
    return value;
}

double get_esum(cluster c)//calculate cluster sum
{
    double esum = 0.;
    for(auto & s: c){
        esum += s.second;
    }
    return esum;
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

//Main analysis function
void FOOT_ana(int firstEvent, int max_events, TString filename = "")
{
    TApplication* theApp = new TApplication("App", 0, 0);

    TFile *datafile = new TFile(filename,"read");
    if(!datafile){ cout << "\n--ERROR in in the input file name!";  return; }

    Double_t pedestal[NDETS][640];
    Double_t sigma[NDETS][640];

    //---------- Definition of main  histograms ----------------
    TH2F * h2_peds_raw[NDETS];//container of raw pedestal data 
    TH2F * h2_cal_coarse[NDETS];
    TH2F * h2_cal_fine[NDETS];
    TH1F * h1_single_event[NDETS];
    TH1D * h1_peds[NDETS];
    TH1D * h1_sigma[NDETS];
    TH1D * h1_cluster_e[NDETS];
    for(int i=0; i<NDETS; i++){
        h2_cal_coarse[i]   = new TH2F(Form("h2_cal_coarse_FOOT%d",foot_id[i]),Form("h2_cal_coarse_FOOT%d",foot_id[i]),640,1,641,Nbins,ymin,ymax);
        h2_cal_fine[i]     = new TH2F(Form("h2_cal_fine_FOOT%d",foot_id[i]),Form("h2_cal_fine_FOOT%d",foot_id[i]),640,1,641,Nbins,ymin,ymax);
        h2_peds_raw[i]     = new TH2F(Form("h%d",foot_id[i]),Form("h%d",foot_id[i]),640,1,641,1500,0,1500);
        h1_single_event[i] = new TH1F(Form("h1_event_FOOT%d",foot_id[i]),Form("h1_event_FOOT%d",foot_id[i]),640,1,641);
        h1_cluster_e[i] = new TH1D(Form("h1_cluster_e_FOOT%d",foot_id[i]),Form("h1_cluster_e_FOOT%d",foot_id[i]),1000,-100,200);
        h1_single_event[i]->SetBarWidth(1);
        if(i%2) h1_single_event[i]->SetFillColor(kBlue);
        else h1_single_event[i]->SetFillColor(kRed);
    }
    TH2I * h2_mul1_vs_mul15 = new TH2I("h2_mul1_vs_mul15","h2_mul1_vs_mul15",20,0,20,20,0,20);
    TH2I * h2_mul2_vs_mul16 = new TH2I("h2_mul2_vs_mul16","h2_mul2_vs_mul16",20,0,20,20,0,20);

    TCanvas * canvas_single_event = new TCanvas("canvas_single_event","canvas_single_event",1200,900);
    canvas_single_event->Divide(4,3);

    TCanvas* canvas_positions = new TCanvas("canvas_positions","canvas_positions",2000,1000);
    canvas_positions->Divide(2,1);
    TH2F * h2_x_vs_z   = new TH2F("h2_x_vs_z","h2_x_vs_z",1000,-50,250,1000,-150,150);
    TH2F * h2_y_vs_x   = new TH2F("h2_y_vs_x","h2_y_vs_x",1000,-150,150,1000,-150,150);

    //Reconstructed vertex
    TCanvas* canvas_vertex = new TCanvas("canvas_vertex","canvas_vertex",1500,1500);
    canvas_vertex->Divide(2,2);
    TH1F * h1_vertex_min   = new TH1F("h1_vertex_min", "h1_vertex_min",1000,0,100);//min distance between tracks
    TH2F * h2_vertex_XY    = new TH2F("h2_vertex_XY",  "h2_vertex_XY", 1000,-300,300,1000,-300,300);//vertex XY profile
    TH2F * h2_vertex_ZX    = new TH2F("h2_vertex_ZX",  "h2_vertex_ZX", 1000,-200,400,1000,-300,300);//vertex ZX profile
    TH2F * h2_vertex_ZY    = new TH2F("h2_vertex_ZY",  "h2_vertex_ZY", 1000,-200,400,1000,-300,300);//vertex ZY profile

    //Angular correlations
    TCanvas* canvas_angcorr = new TCanvas("canvas_angcorr","canvas_angcorr",1000,1000);
    canvas_angcorr->Divide(2,2);
    TH2F * h2_theta1_vs_theta2    = new TH2F("h2_theta1_vs_theta2","h2_theta1_vs_theta2", 100, 0,100,100,0,100);
    TH2F * h2_phi1_vs_phi2  = new TH2F("h2_phi1_vs_phi2","h2_phi1_vs_phi2", 720, -360, 360, 720, -360, 360);
    TH1F * h1_opang         = new TH1F("h1_opang","h1_opang", 180, 0, 180);


    //------ Define tree data for all foot detectors ----------
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

    //----- Get starting set of pedestal data for all dets ------
    int Nevents = tree->GetEntries();
    if(max_events>0) Nevents = max_events; 
    int     stat=0;
    for(int ev=firstEvent; ev<firstEvent+Nevents; ev++)
    {
        tree->GetEntry(ev);
        if((TPATv[0] & 4) !=4 ) continue;
        //if((TPATv[0] & 8)==8) continue;//p2pv
        //if((TPATv[0] & 16)!=16) continue;//OR
        //if(TRIGGER!=3) continue;
        for(int f=0; f<NDETS; f++){
            for(int  j=0 ; j<640 ; j++){
                h2_peds_raw[f]->Fill(FOOTI[f][j],FOOTE[f][j]); 
            }
        }
        stat++;
        if(stat==1e4) break;
    }

    //-------  Slicing, fitting, saving pedestals params -------
    TF1* foo = new TF1("foo","gaus",0,1000);
    for(int i=0; i<NDETS; i++){
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

    foot_data data_left_arm[4];//collection of clusters from left arm: 7, 12, 9, 10 
    foot_data data_right_arm[4];//collection of clusters from right arm: 13, 4, 11, 6
    foot_data data_central_arm[4];//collection of clusters inbeam dets: 15, 16, 1, 2

    for(int ev=firstEvent; ev<firstEvent+Nevents; ev++)
    {
        cout << "\r-- Event # : " << ev << flush;
        tree->GetEntry(ev);
        if((TPATv[0] & 4)!=4) continue;//p2p
        //if((TPATv[0] & 8)!=8) continue;//p2pv
        //if((TPATv[0] & 16)!=16) continue;//OR

        if(EVENT_MONITOR){
            canvas_positions->cd(1); gPad->Clear(); 
            canvas_positions->cd(2); gPad->Clear(); 
        }


        for(int f=0; f<NDETS; f++)//loop over all foots
        {
            //------- Clear data containers ---------
            if(f<4)             data_central_arm[f].clear();
            else if(f>3 && f<8) data_right_arm[f-4].clear();
            else                data_left_arm[f-8].clear(); 
            if(EVENT_MONITOR){
                h1_single_event[f]->Reset("ICESM");//reset individual events 
            }

            //--------  Global base line correction in every FOOT in this event ---------
            mean_ssd=0; stat=0;
            for(int i=0; i<640; i++){
                if(!is_good_strip(foot_id[f],FOOTI[f][i])) continue; 
                signal = FOOTE[f][i] - pedestal[f][i];
                stat++;    mean_ssd += signal;
            }
            mean_ssd /= stat;

            //if(mean_ssd>10){
            //cout << "\n--[WARNING]: In FOOT " << foot_id[f] << "Mean ssd = " << mean_ssd << endl;
            //continue;
            //}

            //------------ Calculating fine baseline correction for individual asics ---------
            stat=0; counter_asic=0;
            for(int i=0; i<10; i++){  asic_offset[i]=0.; }//reset asic baselines
            for(int i=0; i<640; i++)
            {
                signal = FOOTE[f][i] - pedestal[f][i] - mean_ssd;
                if(fabs(signal) < (4 * sigma[f][i]) && //ignore possible hit candidates
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

            //-------- Apply baseline correction, fill cluster data and histos ------------
            counter_asic=0; signal_sum=0;
            for(int i=0; i<640; i++){
                if((FOOTI[f][i]%64) == 1 && FOOTI[f][i]>1) counter_asic++;
                if(!is_good_strip(foot_id[f],FOOTI[f][i])) continue; 
                signal = FOOTE[f][i] - pedestal[f][i] - mean_ssd - asic_offset[counter_asic];
                h2_cal_coarse[f]->Fill(FOOTI[f][i], (FOOTE[f][i] - pedestal[f][i]));
                h2_cal_fine[f]->Fill(FOOTI[f][i], signal);
                if(EVENT_MONITOR) h1_single_event[f]->Fill(FOOTI[f][i], signal);
                if(signal>THRESHOLD){
                    if(f<4) Check_Strip(FOOTI[f][i], signal, data_central_arm[f]);
                    else if(f>3 && f<8) Check_Strip(FOOTI[f][i], signal, data_right_arm[f-4]);
                    else Check_Strip(FOOTI[f][i], signal, data_left_arm[f-8]);
                }
            }
        }//------ end of loop over detectors

        if( 
                data_left_arm[0].size() ==0 || data_left_arm[1].size() ==0 ||
                data_left_arm[2].size() ==0 || data_left_arm[3].size() ==0 ||
                data_right_arm[0].size()==0 || data_right_arm[1].size()==0 ||
                data_right_arm[2].size()==0 || data_right_arm[3].size()==0)
            continue;

        h2_mul1_vs_mul15->Fill(data_central_arm[0].size(), data_central_arm[2].size());
        h2_mul2_vs_mul16->Fill(data_central_arm[1].size(), data_central_arm[3].size());

        //-------- Analyse clusters and get 2p tracks ---------

        //------------- Right arm -----------------------
        double pos[4];
        track t;
        tracks_right.clear();
        TVector3 temp_vec(0., 0., 1.);
        TVector3 beam_vec(0., 0., 1.);// nominal beam axis

        //Define two X and Z points for plane 1 and 4
        double plane1_x[2];
        double plane1_z[2];

        double plane4_x[2];
        double plane4_z[2];

        //central points
        plane1_x[0] = -40.;
        plane1_z[0] =  50.;

        plane4_x[0] = -61.;
        plane4_z[0] =  88.;

        //Second points on 15 deg line:
        plane1_x[1] = plane1_x[0] + 10.* TMath::Sin(15. * TMath::Pi() / 180.);
        plane1_z[1] = plane1_z[0] + 10.* TMath::Cos(15. * TMath::Pi() / 180.);

        plane4_x[1] = plane4_x[0] + 10.* TMath::Sin(15. * TMath::Pi() / 180.);
        plane4_z[1] = plane4_z[0] + 10.* TMath::Cos(15. * TMath::Pi() / 180.);

        //Line parameters for layer 1 and 4
        double plane1_slope  = (plane1_x[1] - plane1_x[0]) / (plane1_z[1] - plane1_z[0]);
        double plane1_offset =  plane1_x[0] - plane1_slope * plane1_z[0];

        double plane4_slope  = (plane4_x[1] - plane4_x[0]) / (plane4_z[1] - plane4_z[0]);
        double plane4_offset =  plane4_x[0] - plane4_slope * plane4_z[0];

        double track_slope, track_offset;


        for(auto & c0: data_right_arm[0]){
            pos[0] = get_cog(c0) * FOOT_LENGTH/640. - 0.5*FOOT_LENGTH;

            for(auto & c1: data_right_arm[1]){
                pos[1] = get_cog(c1) * FOOT_LENGTH/640. - 0.5*FOOT_LENGTH;

                for(auto & c2: data_right_arm[2]){
                    pos[2] = get_cog(c2) * FOOT_LENGTH/640. - 0.5*FOOT_LENGTH;

                    for(auto & c3: data_right_arm[3]){
                        pos[3] = get_cog(c3) * FOOT_LENGTH/640. - 0.5*FOOT_LENGTH;

                        //if( 
                        //        get_esum(c0)<THRESHOLD || get_esum(c0)>MAX_THRESHOLD ||
                        //        get_esum(c1)<THRESHOLD || get_esum(c1)>MAX_THRESHOLD ||
                        //        get_esum(c2)<THRESHOLD || get_esum(c2)>MAX_THRESHOLD ||
                        //        get_esum(c3)<THRESHOLD || get_esum(c3)>MAX_THRESHOLD
                        //  ) continue; 


                        //2nd layer
                        t.point[1].SetX(pos[1] * TMath::Sin(15. * TMath::Pi() / 180.) - 44.);   
                        t.point[1].SetY(0);   
                        t.point[1].SetZ(pos[1] * TMath::Cos(15. * TMath::Pi() / 180.) + 57.);   

                        //3rd layer
                        t.point[2].SetX(pos[2] * TMath::Sin(15. * TMath::Pi() / 180.) - 59.);   
                        t.point[2].SetY(0);   
                        t.point[2].SetZ(pos[2] * TMath::Cos(15. * TMath::Pi() / 180.) + 77.);   

                        track_slope  = (t.point[2].X() - t.point[1].X())/(t.point[2].Z() - t.point[1].Z());  
                        track_offset = (t.point[2].X() - track_slope*t.point[2].Z());

                        //1st layer
                        t.point[0].SetY(pos[0]*(-1));   
                        t.point[0].SetZ((plane1_offset - track_offset) / (track_slope - plane1_slope));//extrapolated 
                        t.point[0].SetX(track_slope * t.point[0].Z() + track_offset);//extrapolated   

                        //4th layer
                        t.point[3].SetY(pos[3]);   
                        t.point[3].SetZ((plane4_offset - track_offset) / (track_slope - plane4_slope));   
                        t.point[3].SetX(track_slope * t.point[3].Z() + track_offset);//extrapolated   

                        //Y position in 2 middle plains
                        track_slope  = (t.point[3].Y() - t.point[0].Y())/(t.point[3].X() - t.point[0].X());  
                        track_offset = (t.point[3].Y() - track_slope*t.point[3].X());

                        //2nd layer
                        t.point[1].SetY(track_offset + track_slope * t.point[1].X());   

                        //3rd layer
                        t.point[2].SetY(track_offset + track_slope * t.point[2].X());   

                        temp_vec.SetXYZ(
                                t.point[3].X() - t.point[0].X(),
                                t.point[3].Y() - t.point[0].Y(),
                                t.point[3].Z() - t.point[0].Z());

                        if(
                                temp_vec.Angle(beam_vec)*180./3.14 > 85. || 
                                temp_vec.Angle(beam_vec)*180./3.14 < 20.
                          ) continue;

                        tracks_right.push_back(t);

                        h2_x_vs_z->Fill(t.point[0].Z(),t.point[0].X());
                        h2_x_vs_z->Fill(t.point[1].Z(),t.point[1].X());
                        h2_x_vs_z->Fill(t.point[2].Z(),t.point[2].X());
                        h2_x_vs_z->Fill(t.point[3].Z(),t.point[3].X());

                        h2_y_vs_x->Fill(t.point[0].X(),t.point[0].Y());
                        h2_y_vs_x->Fill(t.point[1].X(),t.point[1].Y());
                        h2_y_vs_x->Fill(t.point[2].X(),t.point[2].Y());
                        h2_y_vs_x->Fill(t.point[3].X(),t.point[3].Y());

                    }
                }
            }
        }

        //------------- Left arm ----------------------
        tracks_left.clear();
        //Define two X and Z points for plane 1 and 4
        //central points
        plane1_x[0] =  40.;
        plane1_z[0] =  50.;

        plane4_x[0] =  61.;
        plane4_z[0] =  88.;

        //Second points on 15 deg line:
        plane1_x[1] = plane1_x[0] + 10.* TMath::Sin(15. * TMath::Pi() / 180.);
        plane1_z[1] = plane1_z[0] - 10.* TMath::Cos(15. * TMath::Pi() / 180.);

        plane4_x[1] = plane4_x[0] + 10.* TMath::Sin(15. * TMath::Pi() / 180.);
        plane4_z[1] = plane4_z[0] - 10.* TMath::Cos(15. * TMath::Pi() / 180.);

        //Line parameters for layer 1 and 4
        plane1_slope  = (plane1_x[1] - plane1_x[0]) / (plane1_z[1] - plane1_z[0]);
        plane1_offset =  plane1_x[0] - plane1_slope * plane1_z[0];

        plane4_slope  = (plane4_x[1] - plane4_x[0]) / (plane4_z[1] - plane4_z[0]);
        plane4_offset =  plane4_x[0] - plane4_slope * plane4_z[0];

        
        //Analyse

        for(auto & c0: data_left_arm[0]){
            pos[0] = get_cog(c0) * FOOT_LENGTH/640. - 0.5*FOOT_LENGTH;

            for(auto & c1: data_left_arm[1]){
                pos[1] = get_cog(c1) * FOOT_LENGTH/640. - 0.5*FOOT_LENGTH;

                for(auto & c2: data_left_arm[2]){
                    pos[2] = get_cog(c2) * FOOT_LENGTH/640. - 0.5*FOOT_LENGTH;

                    for(auto & c3: data_left_arm[3]){
                        pos[3] = get_cog(c3) * FOOT_LENGTH/640. - 0.5*FOOT_LENGTH;

                        //if( 
                        //        get_esum(c0)<THRESHOLD || get_esum(c0)>MAX_THRESHOLD ||
                        //        get_esum(c1)<THRESHOLD || get_esum(c1)>MAX_THRESHOLD ||
                        //        get_esum(c2)<THRESHOLD || get_esum(c2)>MAX_THRESHOLD ||
                        //        get_esum(c3)<THRESHOLD || get_esum(c3)>MAX_THRESHOLD
                        //  ) continue; 


                        //2nd layer
                        t.point[1].SetX(pos[1] * TMath::Sin(15. * TMath::Pi() / 180.) + 44.);   
                        t.point[1].SetY(0);   
                        t.point[1].SetZ((-1)*pos[1] * TMath::Cos(15. * TMath::Pi() / 180.) + 57.);   

                        //3rd layer
                        t.point[2].SetX(pos[2] * TMath::Sin(15. * TMath::Pi() / 180.) + 59.);   
                        t.point[2].SetY(0);   
                        t.point[2].SetZ((-1)*pos[2] * TMath::Cos(15. * TMath::Pi() / 180.) + 77.);   

                        track_slope  = (t.point[2].X() - t.point[1].X())/(t.point[2].Z() - t.point[1].Z());  
                        track_offset = (t.point[2].X() - track_slope*t.point[2].Z());

                        //1st layer
                        t.point[0].SetY(pos[0]);   
                        t.point[0].SetZ((plane1_offset - track_offset) / (track_slope - plane1_slope));//extrapolated 
                        t.point[0].SetX(track_slope * t.point[0].Z() + track_offset);//extrapolated   

                        //4th layer
                        t.point[3].SetY((-1)*pos[3]);   
                        t.point[3].SetZ((plane4_offset - track_offset) / (track_slope - plane4_slope));   
                        t.point[3].SetX(track_slope * t.point[3].Z() + track_offset);//extrapolated   

                        //Y position in 2 middle plains
                        track_slope  = (t.point[3].Y() - t.point[0].Y())/(t.point[3].X() - t.point[0].X());  
                        track_offset = (t.point[0].Y() - track_slope*t.point[0].X());

                        //2nd layer
                        t.point[1].SetY(track_offset + track_slope * t.point[1].X());   

                        ////3rd layer
                        t.point[2].SetY(track_offset + track_slope * t.point[2].X());   

                        temp_vec.SetXYZ(
                                t.point[3].X()-t.point[0].X(),
                                t.point[3].Y()-t.point[0].Y(),
                                t.point[3].Z()-t.point[0].Z());

                        if(
                                temp_vec.Angle(beam_vec)*180./3.14 > 85. || 
                                temp_vec.Angle(beam_vec)*180./3.14 < 20.
                          ) continue;

                        tracks_left.push_back(t);

                        h2_x_vs_z->Fill(t.point[0].Z(),t.point[0].X());
                        h2_x_vs_z->Fill(t.point[1].Z(),t.point[1].X());
                        h2_x_vs_z->Fill(t.point[2].Z(),t.point[2].X());
                        h2_x_vs_z->Fill(t.point[3].Z(),t.point[3].X());

                        h2_y_vs_x->Fill(t.point[0].X(),t.point[0].Y());
                        h2_y_vs_x->Fill(t.point[1].X(),t.point[1].Y());
                        h2_y_vs_x->Fill(t.point[2].X(),t.point[2].Y());
                        h2_y_vs_x->Fill(t.point[3].X(),t.point[3].Y());

                    }
                }
            }
        }

        //if(tracks_right.size()>10 || tracks_left.size()>10) continue;
        if(tracks_right.size()==0 || tracks_left.size()==0) continue;
        
        for(auto & c0: data_right_arm[0]){
            h1_cluster_e[4]->Fill( get_esum(c0) );
        }
        for(auto & c0: data_right_arm[1]){
            h1_cluster_e[5]->Fill( get_esum(c0) );
        }
        for(auto & c0: data_right_arm[2]){
            h1_cluster_e[6]->Fill( get_esum(c0) );
        }
        for(auto & c0: data_right_arm[3]){
            h1_cluster_e[7]->Fill( get_esum(c0) );
        }
        for(auto & c0: data_left_arm[0]){
            h1_cluster_e[8]->Fill( get_esum(c0) );
        }
        for(auto & c0: data_left_arm[1]){
            h1_cluster_e[9]->Fill( get_esum(c0) );
        }
        for(auto & c0: data_left_arm[2]){
            h1_cluster_e[10]->Fill( get_esum(c0) );
        }
        for(auto & c0: data_left_arm[3]){
            h1_cluster_e[11]->Fill( get_esum(c0) );
        }
      
        //Eventually calculating vertex and min distance with any left-right combinations
        TVector3 pos_left;
        TVector3 vec_left;
        TVector3 vec_left_good;

        TVector3 pos_right;
        TVector3 vec_right;
        TVector3 vec_right_good;

        vertex this_vertex;
        vertex min_vertex;

        
        min_vertex.dca = 10000;//starting value
        int n_good_pairs = 0;
        vertex_array.clear();
        for(auto & tl: tracks_left)
        {
            for(auto & tr: tracks_right)
            {
                pos_left.SetXYZ(tl.point[0].X(), tl.point[0].Y(), tl.point[0].Z()); 
                pos_right.SetXYZ(tr.point[0].X(), tr.point[0].Y(), tr.point[0].Z()); 

                vec_left.SetXYZ( 
                        tl.point[3].X() - tl.point[0].X(),
                        tl.point[3].Y() - tl.point[0].Y(),
                        tl.point[3].Z() - tl.point[0].Z()
                        );

                vec_right.SetXYZ( 
                        tr.point[3].X() - tr.point[0].X(),
                        tr.point[3].Y() - tr.point[0].Y(),
                        tr.point[3].Z() - tr.point[0].Z()
                        );

                vec_left  = vec_left.Unit();
                vec_right = vec_right.Unit();

                this_vertex = get_vertex(pos_left,pos_right,vec_left,vec_right);
                //if(this_vertex.dca<40 && this_vertex.dca>50) continue;

                //if(
                //        this_vertex.vertexPos.X() < -30 || this_vertex.vertexPos.X()>30 || 
                //        this_vertex.vertexPos.Y() < -30 || this_vertex.vertexPos.Y()>30 || 
                //        this_vertex.vertexPos.Z() < -90 || this_vertex.vertexPos.Z()>90
                //  ) continue;



                if(this_vertex.dca<min_vertex.dca)
                {
                    min_vertex.dca = this_vertex.dca;
                    min_vertex.vertexPos = this_vertex.vertexPos;
                    min_vertex.normalVec = this_vertex.normalVec;

                    vec_left_good = vec_left;
                    vec_right_good = vec_right;
                }

                if(this_vertex.dca<MIN_DISTANCE)
                {   
                    n_good_pairs++;
                    vertex_array.push_back(this_vertex);
                }
                //h1_vertex_min->Fill(this_vertex.dca);
                //h2_vertex_XY->Fill(this_vertex.vertexPos.X(), this_vertex.vertexPos.Y());
                //h2_vertex_ZX->Fill(this_vertex.vertexPos.Z(), this_vertex.vertexPos.X());
                //h2_vertex_ZY->Fill(this_vertex.vertexPos.Z(), this_vertex.vertexPos.Y());
            }
        }
        
        if(n_good_pairs!=1) continue;
        h1_vertex_min->Fill(min_vertex.dca);
        //if(min_vertex.dca<MIN_DISTANCE)
        cout << "\nNumber of good pairs in this event: " << n_good_pairs << endl;
        for(auto & v: vertex_array)
        {

            //   h2_vertex_XY->Fill(min_vertex.vertexPos.X(), min_vertex.vertexPos.Y());
            //   h2_vertex_ZX->Fill(min_vertex.vertexPos.Z(), min_vertex.vertexPos.X());
            //   h2_vertex_ZY->Fill(min_vertex.vertexPos.Z(), min_vertex.vertexPos.Y());

            h2_vertex_XY->Fill(v.vertexPos.X(), v.vertexPos.Y());
            h2_vertex_ZX->Fill(v.vertexPos.Z(), v.vertexPos.X());
            h2_vertex_ZY->Fill(v.vertexPos.Z(), v.vertexPos.Y());

            h2_theta1_vs_theta2->Fill(vec_left_good.Theta()*180./3.14,vec_right_good.Theta()*180./3.14);
            h2_phi1_vs_phi2->Fill(vec_left_good.Phi()*180./3.14,vec_right_good.Phi()*180./3.14);
            h1_opang->Fill(vec_left_good.Angle(vec_right_good)*180./3.14);
        }

        if(EVENT_MONITOR) //plot individual events
        {
            TLine* lt;
            TMarker * m1;
            TMarker * m2;
            TMarker * m3;
            TMarker * m4;
            double k, b;

            cout << "\nNumber of track candidates in right (WX) arm: \t" << tracks_right.size();
            cout << "\nNumber of track candidates in left (Messel) arm: \t" << tracks_left.size() << endl;
            canvas_positions->cd(1);
            h2_x_vs_z->Draw("colz");

            TBox * lh2_box = new TBox(-28.3,-25.,21.7,25);
            lh2_box->SetFillColorAlpha(kBlue, 0.25);
            lh2_box->Draw("same");

            canvas_positions->cd(2);
            h2_y_vs_x->Draw("colz");
            TEllipse *lh2_circ = new TEllipse(0.,0.,20,20);
            lh2_circ->SetFillColorAlpha(kBlue, 0.25);
            lh2_circ->Draw("same");

            //----- left arm tracks ---------
            for(auto & t: tracks_left)
            {
                //Paramters of the track line

                canvas_positions->cd(1);
                k = (t.point[2].X() - t.point[1].X())/(t.point[2].Z() - t.point[1].Z());
                b = (t.point[2].X() - k* t.point[2].Z());

                lt = new TLine(t.point[3].Z(), t.point[3].X(), -50,k*(-50.)+b);
                lt->SetLineWidth(1);
                lt->SetLineColor(kBlack);
                lt->Draw();

                m1 = new TMarker(t.point[0].Z(),t.point[0].X(),8);
                m1->SetMarkerColor(kGray);
                m1->Draw();

                m2 = new TMarker(t.point[1].Z(),t.point[1].X(),8);
                m2->SetMarkerColor(kRed);
                m2->Draw();

                m3 = new TMarker(t.point[2].Z(),t.point[2].X(),8);
                m3->SetMarkerColor(kRed);
                m3->Draw();

                m4 = new TMarker(t.point[3].Z(),t.point[3].X(),8);
                m4->SetMarkerColor(kGray);
                m4->Draw();               


                canvas_positions->cd(2);
                k = (t.point[3].Y() - t.point[0].Y())/(t.point[3].X() - t.point[0].X());
                b = (t.point[3].Y() - k* t.point[3].X());

                lt = new TLine(t.point[3].X(), t.point[3].Y(), -100, k*(-100.)+b);
                lt->SetLineWidth(1);
                lt->SetLineColor(kBlack);
                lt->Draw();

                m1 = new TMarker(t.point[0].X(),t.point[0].Y(),8);
                m1->SetMarkerColor(kRed);
                m1->Draw();

                m2 = new TMarker(t.point[1].X(),t.point[1].Y(),8);
                m2->SetMarkerColor(kGray);
                m2->Draw();

                m3 = new TMarker(t.point[2].X(),t.point[2].Y(),8);
                m3->SetMarkerColor(kGray);
                m3->Draw();

                m4 = new TMarker(t.point[3].X(),t.point[3].Y(),8);
                m4->SetMarkerColor(kRed);
                m4->Draw();
            }
            //----- right arm tracks ---------
            for(auto & t: tracks_right)
            {
                //Paramters of the track line
                canvas_positions->cd(1);
                k = (t.point[2].X() - t.point[1].X())/(t.point[2].Z() - t.point[1].Z());
                b = (t.point[2].X() - k* t.point[2].Z());

                lt = new TLine(t.point[3].Z(), t.point[3].X(), -50, k*(-50.)+b);
                lt->SetLineWidth(1);
                lt->SetLineColor(kBlack);
                lt->Draw();

                m1 = new TMarker(t.point[0].Z(),t.point[0].X(),8);
                m1->SetMarkerColor(kGray);
                m1->Draw();

                m2 = new TMarker(t.point[1].Z(),t.point[1].X(),8);
                m2->SetMarkerColor(kRed);
                m2->Draw();

                m3 = new TMarker(t.point[2].Z(),t.point[2].X(),8);
                m3->SetMarkerColor(kRed);
                m3->Draw();

                m4 = new TMarker(t.point[3].Z(),t.point[3].X(),8);
                m4->SetMarkerColor(kGray);
                m4->Draw();

                canvas_positions->cd(2);
                k = (t.point[3].Y() - t.point[0].Y())/(t.point[3].X() - t.point[0].X());
                b = (t.point[3].Y() - k* t.point[3].X());

                lt = new TLine(t.point[3].X(), t.point[3].Y(), 100, k*(100.)+b);
                lt->SetLineWidth(1);
                lt->SetLineColor(kBlack);
                lt->Draw();

                m1 = new TMarker(t.point[0].X(),t.point[0].Y(),8);
                m1->SetMarkerColor(kRed);
                m1->Draw();

                m2 = new TMarker(t.point[1].X(),t.point[1].Y(),8);
                m2->SetMarkerColor(kGray);
                m2->Draw();

                m3 = new TMarker(t.point[2].X(),t.point[2].Y(),8);
                m3->SetMarkerColor(kGray);
                m3->Draw();

                m4 = new TMarker(t.point[3].X(),t.point[3].Y(),8);
                m4->SetMarkerColor(kRed);
                m4->Draw();

            }

            canvas_positions->Modified();
            canvas_positions->Update();

            //------- Plotting actual FOOT data for every detector ---------
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
                    TLine* l = new TLine(64*i_asic,ymin_new,64*i_asic,ymax_new);
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
    TCanvas* canvas_esum   = new TCanvas("esum","esum",1200,900);
    canvas_coarse->Divide(4,3);
    canvas_fine->Divide(4,3);
    canvas_esum->Divide(4,3);
    for(int f=0; f<NDETS; f++)
    {
        canvas_coarse->cd(f+1);
        gPad->SetLogz();
        h2_cal_coarse[f]->Draw("colz");

        canvas_fine->cd(f+1);
        gPad->SetLogz();
        h2_cal_fine[f]->Draw("colz");
    
        canvas_esum->cd(f+1);
        h1_cluster_e[f]->Draw();
    }

    TCanvas* canvas_mul = new TCanvas("canvas_mul","canvas_mul",1000,500);
    canvas_mul->Divide(2,1);
    canvas_mul->cd(1);
    gPad->SetLogz();
    h2_mul1_vs_mul15->Draw("colz");
    canvas_mul->cd(2);
    gPad->SetLogz();
    h2_mul2_vs_mul16->Draw("colz");

    canvas_positions->cd();
    h2_x_vs_z->Draw("colz");

    canvas_vertex->cd(1);
    h1_vertex_min->Draw();
    canvas_vertex->cd(2);
    h2_vertex_XY->Draw("colz");
    canvas_vertex->cd(3);
    h2_vertex_ZX->Draw("colz");
    canvas_vertex->cd(4);
    h2_vertex_ZY->Draw("colz");
    
    canvas_angcorr->cd(1);
    h2_theta1_vs_theta2->Draw("colz");
    canvas_angcorr->cd(2);
    h2_phi1_vs_phi2->Draw("colz");
    canvas_angcorr->cd(3);
    h1_opang->Draw();


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
    //TString filename = Form("../main0131_0041.root");
    TString filename = Form("../main0131_0001.root");
    FOOT_ana(0,-1,filename);
    //FOOT_ana(0,-1,filename);
    return 0;
}

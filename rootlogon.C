void rootlogon()
{
	gROOT->SetStyle("Default");

	//gROOT->SetStyle("Plain");
	//gStyle->SetPalette(1);
	gStyle->SetDrawBorder(0);

	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasBorderSize(0);
	gStyle->SetCanvasColor(kWhite);
	gStyle->SetPadColor(kWhite);
	//gStyle->SetPadColor(kBlack);
	gStyle->SetCanvasDefH(900);
	gStyle->SetCanvasDefW(1300);
	gStyle->SetFrameFillColor(kWhite);
	//gStyle->SetFrameFillColor(kBlack);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameBorderSize(0);

	gStyle->SetOptFit(0);
	//gStyle->SetHatchesLineWidth(1);
	//gStyle->SetHatchesSpacing(1);

	gStyle->SetStatColor(kWhite);
	gStyle->SetStatBorderSize(1);
	//gStyle->SetStatFont(42);
	gStyle->SetStatFont(132);
	gStyle->SetOptStat(0);
	//gStyle->SetOptStat(1110000);
	//gStyle->SetOptStat(1111111);
	gStyle->SetHistFillColor(kWhite);
	//gStyle->SetOptTitle(0);
	//gStyle->SetTextFont(42);//Helvetica
	gStyle->SetTextFont(132);//Times
	//gStyle->SetTitleFont(42,"xyz");//Helvetica
	gStyle->SetTitleFont(132,"xyz");//Helvetica
	//gStyle->SetTitleFont(42,"p");//Helvetica
	gStyle->SetTitleFillColor(kWhite);
	gStyle->SetTitleBorderSize(1);
	//gStyle->SetTitleAlign(5);
	gStyle->SetTitleW(0.7);
	gStyle->SetTitleH(0.075);
	gStyle->SetTitleX(0.15);
	gStyle->SetTitleY(1.00);
	gStyle->SetTitleTextColor(kBlue);

	gStyle->SetFuncColor(kRed);
	gStyle->SetFuncWidth(3);
	//gStyle->SetTitleFont(132,"xyz");//Times
	gStyle->SetTitleXSize(0.06);
	gStyle->SetTitleXOffset(1.);
	gStyle->SetTitleYSize(0.06);
	gStyle->SetTitleYOffset(1.15);
	//gStyle->SetLabelFont(42,"xyz");//Helvetica
	//gStyle->SetLabelFont(42,"p");//Helvetica
	gStyle->SetLabelFont(132,"xyz");//Times
	gStyle->SetLabelSize(0.04,"XYZ");

	//gStyle->SetPadGridX(1);
	//gStyle->SetPadGridY(1);

	gStyle->SetLegendBorderSize(1);
	//gStyle->SetLegendFillColor(kWhite);
	//gStyle->SetLegendFont(42);

	gStyle->SetPadColor(kWhite);
	gStyle->SetPadBorderSize(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadBottomMargin(0.13);
	gStyle->SetPadTopMargin(0.08);
	gStyle->SetPadLeftMargin(0.13);
	gStyle->SetPadRightMargin(0.13);

        //gStyle->SetPadBottomMargin(0.05);
	//gStyle->SetPadTopMargin(0.001);
	//gStyle->SetPadLeftMargin(0.15);
	//gStyle->SetPadRightMargin(0.001);



	gStyle->SetGridColor(kWhite);
	//gStyle->SetGridStyle(3);
	//gStyle->SetGridWidth(3);
	//gStyle->SetPadTickX(1);
	//gStyle->SetPadTickY(1);

	//Color Palette
	const Int_t NRGBs = 5;
	//const Int_t NRGBs = 4;
	//const Int_t NCont = 50;
	const Int_t NCont = 99;
	//const Int_t NCont = 20;

	//black, blue, magenta,red, yellow, white
	Double_t stops[NRGBs] = { 0.00, 0.15,  0.3,  0.7,   1.00};
	Double_t red[NRGBs]   = { 0.00, 0.00, 1.00, 1.00,  1.00}; 
	Double_t green[NRGBs] = { 0.00, 0.24, 0.00, 0.89,  1.00};
	Double_t blue[NRGBs]  = { 0.00, 1.00, 0.00, 0.01,  1.00};
	
	//black, blue, magenta,red, yellow
	//Double_t stops[NRGBs] = { 0.0, 0.3,  0.6, 1.00};
	//Double_t red[NRGBs]   = { 0.0, 0.00, 1.00, 1.00}; 
	//Double_t green[NRGBs] = { 0.0, 0.24, 0.00, 0.89};
	//Double_t blue[NRGBs]  = { 0.0, 1.00, 0.00, 0.01};
	
	//black, blue, magenta,red, yellow
	//Double_t stops[NRGBs] = { 0.0, 0.3,  0.6, 1.00};
	//Double_t red[NRGBs]   = { 1.0, 0.00, 1.00, 1.00}; 
	//Double_t green[NRGBs] = { 1.0, 0.24, 0.00, 0.89};
	//Double_t blue[NRGBs]  = { 1.0, 1.00, 0.00, 0.01};

	////The best one!!! black, blue, magenta,red, yellowi, white
        //Double_t stops[NRGBs] = { 0., 0.25,  0.7, 1.00};//best Nature physics SRC
        //Double_t stops[NRGBs] = { 0., 0.2,  0.5, 1.00};//best Nature physics SRC

        //Double_t stops[NRGBs] = { 0.0, 0.1,  0.4, 1.00};
        //Double_t red[NRGBs]   = { 1.0, 0.00, 1.00, 1.00}; 
        //Double_t green[NRGBs] = { 1.0, 0.24, 0.00, 0.89};
        //Double_t blue[NRGBs]  = { 1.0, 1.00, 0.00, 0.01};


        Int_t NewPalette = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
        gStyle->SetNumberContours(NCont);


        //const Int_t NRGBs = 3;
        //const Int_t NCont = 255;
        //Double_t stops[NRGBs] = { 0.00, 0.50, 1.00};
        //Double_t green[NRGBs]  = { 0.00, 0.00, 1.00};
        //Double_t orange[NRGBs] = { 0.00, 1.00, 0.00};
        //Double_t violet[NRGBs]  = { 1.00, 0.00, 0.00};

        //Int_t NewPalette = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
        //Int_t NewPalette = TColor::CreateGradientColorTable(NRGBs, stops, green, orange, violet, NCont);
        //gStyle->SetNumberContours(NCont);

        //const int NRGBs = 7;
        //const int NCont = 255;
        //Double_t stops[NRGBs] = { 0.00, 0.10, 0.25, 0.45, 0.60, 0.75, 1.00 };
        //Double_t red[NRGBs]   = { 1.00, 0.00, 0.00, 0.00, 0.97, 0.97, 0.10 };
        //Double_t green[NRGBs] = { 1.00, 0.97, 0.30, 0.40, 0.97, 0.00, 0.00 };
        //Double_t blue[NRGBs]  = { 1.00, 0.97, 0.97, 0.00, 0.00, 0.00, 0.00 };
        //Int_t NewPalette = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
        //gStyle->SetNumberContours(NCont);





        //Color Palette
        //const Int_t NRGBs = 5;
        //const Int_t NCont = 255;
        //Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
        //Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
        //Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
        //Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
        //Int_t NewPalette = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
        //gStyle->SetNumberContours(NCont);

        //Int_t palette[13] = {321,320,420,420,520,520,620,620,720,720,820,920,950}; 
        //gStyle->SetPalette(13,palette);

        //Black and white Palette
        //Int_t myGreyPalette[10] = {0, 19, 18, 17, 16, 15, 14, 13, 12, 1};
        //gStyle->SetPalette(10,myGreyPalette);


        gStyle->UseCurrentStyle();
        gROOT->ForceStyle();

        printf("\nCustom styles are loaded! \n\n");
        return;
}


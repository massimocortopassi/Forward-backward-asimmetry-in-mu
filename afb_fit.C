#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <algorithm>
#include <TStyle.h>
#include <sstream>
#include <iomanip>

void afb_fit() {

    // --- CONFIGURAZIONE DELLO STILE DEL BOX STATISTICHE ---
    gStyle->SetOptStat(0);    // Disabilita COMPLETAMENTE il box statistiche di default
    gStyle->SetOptFit(0);     // Disabilita i parametri del fit di default

    const char* files[] = {
        //"zdec_da9300.root"
        //"zdec_da9125_quart.root"
        //"zdec_da8950.root"
        //"zdec_hadr.root"
        //"zdec_taus.root"
        "zdec_muon.root"
    };

    int num_files = sizeof(files) / sizeof(files[0]);

    TH1D* h = new TH1D(
        "h",
        ";cos#theta*; Events",
        45, -0.9, 0.9
    );

    // --- NUOVO: CONTATORE EVENTI INIZIALI ---
    Long64_t total_initial_events = 0;

    for (int k = 0; k < num_files; k++) {

        TFile* f = TFile::Open(files[k]);
        if (!f || f->IsZombie()) continue;

        TTree* tree = (TTree*)f->Get("h100");
        if (!tree) continue;

        Int_t Nch;
        Int_t Nplanes[500];
        Float_t Cthr;
        Int_t Qcha[500];
        Int_t Nhits[500];
        Float_t pcha[500][3];

        // Nuova variabile: parametro d'impatto trasverso
        Float_t D0[500];
        Float_t Z0[500];
        // Variabili energetiche
        Float_t Ehc[500];
        Float_t Eec[500];
        Float_t Elep;

        tree->SetBranchAddress("Nhits", Nhits);
        tree->SetBranchAddress("Nch", &Nch);
        tree->SetBranchAddress("Nplanes", Nplanes);
        tree->SetBranchAddress("Cthr", &Cthr);
        tree->SetBranchAddress("Qcha", Qcha);
        tree->SetBranchAddress("pcha", pcha);

        // Collegamento della nuova branch D0
        tree->SetBranchAddress("D0", D0);
        tree->SetBranchAddress("Z0", Z0);
        // Collegamento delle branch energetiche
        tree->SetBranchAddress("Ehc", Ehc);
        tree->SetBranchAddress("Eec", Eec);
        tree->SetBranchAddress("Elep", &Elep);

        Long64_t nev = tree->GetEntries();

        // --- NUOVO: AGGIORNO IL CONTATORE E STAMPO ---
        total_initial_events += nev;
        std::cout << "Apertura file: " << files[k] << " - Eventi iniziali contenuti: " << nev << std::endl;
        // ---------------------------------------------

        for (Long64_t i = 0; i < nev; i++) {

            tree->GetEntry(i);

            // Taglio sull'energia nel centro di massa (entro 1 GeV dal picco della Z)
            if (std::abs(Elep - 91.19) >= 1.0) continue;

            // 1. Taglio sul numero di tracce cariche (esattamente 2)
            if (Nch != 2) continue;

            double p1 = sqrt(pcha[0][0] * pcha[0][0] + pcha[0][1] * pcha[0][1] + pcha[0][2] * pcha[0][2]);
            double p2 = sqrt(pcha[1][0] * pcha[1][0] + pcha[1][1] * pcha[1][1] + pcha[1][2] * pcha[1][2]);

            //Controllo sulle singole tracce nel laboratorio
            double cos_traccia1 = pcha[0][2] / p1;
            double cos_traccia2 = pcha[1][2] / p2;

            if (std::abs(cos_traccia1) >= 0.9 || std::abs(cos_traccia2) >= 0.9) continue;

            // 0. Taglio sul numero di hits
            if (Nhits[0] < 4 && Nhits[1] < 4) continue;

            // 2. Taglio di qualità sui piani del rivelatore
            // Almeno UNA traccia deve avere > 3 hits negli ultimi piani dell'HCAL
            if (Nplanes[0] <= 3 && Nplanes[1] <= 3) continue;

            // 3. Taglio sul parametro d'impatto (D0 < 0.5 cm per entrambe le tracce)
            if (fabs(D0[0]) >= 0.5 || fabs(D0[1]) >= 0.5) continue;

            if (fabs(Z0[0]) >= 4 || fabs(Z0[1]) >= 4) continue;


            // --------------------------------------------

            // 5. Taglio geometrico sull'asse di thrust????????????????????????????????????????'
            //if (fabs(Cthr) >= 0.9) continue;

            if (Qcha[0] == Qcha[1]) continue;


            // --- NUOVO TAGLIO DI ACOLLINEARITÀ ---
            // 1. Calcolo del prodotto scalare tra i vettori tridimensionali dei momenti
            double dot_product =
                pcha[0][0] * pcha[1][0]
                +
                pcha[0][1] * pcha[1][1]
                +
                pcha[0][2] * pcha[1][2];

            // 2. Coseno dell'angolo tra le due tracce
            double cos_theta12 = dot_product / (p1 * p2);

            // 3. Calcolo di eta in gradi (usando -cos_theta12)
            double acollinearity =
                acos(std::clamp(-cos_theta12, -1.0, 1.0)) * (180.0 / TMath::Pi());

            // 4. Applicazione del taglio
            if (acollinearity >= 23.0) continue;
            // ------------------------------------

            double pmax =
                std::max(p1, p2);

            double pmin =
                std::min(p1, p2);

            if (
                pmax <= 35
                ||
                pmin <= 22
                )
                continue;

            // svoglere megio l'arcocoseno
            double theta1 =
                acos(
                    std::clamp(
                        pcha[0][2] / p1,
                        -1.0,
                        1.0
                    )
                );

            double theta2 =
                acos(
                    std::clamp(
                        pcha[1][2] / p2,
                        -1.0,
                        1.0
                    )
                );

            if (Qcha[0] > 0) {
                std::swap(theta1, theta2);
            }

            double den = cos(0.5 * (theta1 - TMath::Pi() + theta2));

            if (fabs(den) < 1e-8) continue;

            double CthrStar =
                cos(0.5 * (theta1 + TMath::Pi() - theta2)) / den;

            CthrStar =
                std::clamp(
                    CthrStar,
                    -1.0,
                    1.0
                );

            //Controllo nel Centro di Massa per il Fit
            if (CthrStar <= -0.9 || CthrStar >= 0.9) continue;

            h->Fill(CthrStar);
        }

        f->Close();
    }



    TF1* fit = new TF1(
        "fit",
        "[0]*(1+x*x+(8./3.)*[1]*x)",
        -0.9, 0.9
    );

    fit->SetParNames("C", "A_FB");

    fit->SetParameters(
        h->GetMaximum() / 2,
        0.0
    );

    h->Fit(
        fit,
        "RM"
    );


    // --- CREAZIONE DEL CANVAS E DIVISIONE IN PAD ---
    // Aumentiamo un po' l'altezza del canvas per far spazio ai residui
    TCanvas* c = new TCanvas("c", "Asimmetria Forward-Backward", 800, 800);

    // Pad superiore per il fit
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.02); // Rimuove il margine inferiore per incollarlo al pad sotto
    pad1->Draw();

    // Pad inferiore per i residui
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
    pad2->SetTopMargin(0.02);    // Rimuove il margine superiore
    pad2->SetBottomMargin(0.3);  // Lascia spazio per i titoli dell'asse X
    pad2->Draw();

    // ---------------------------------------------------------
    // --- DISEGNO NEL PAD SUPERIORE (GRAFICO ORIGINALE) ---
    pad1->cd();

    h->SetMinimum(0.0);
    h->SetMaximum(h->GetMaximum() * 1.3);

    // Rimuoviamo il titolo dell'asse X qui, lo metteremo nel pad inferiore
    h->SetTitle("#sqrt{s} = 91.19 GeV;; Events");
    h->GetXaxis()->SetLabelSize(0); // Nasconde i numeri sull'asse X del pad superiore

    h->Draw("HIST E");
    fit->Draw("same");

    // Crea la label
    TPaveText* label = new TPaveText(0.15, 0.70, 0.45, 0.89, "NDC");
    label->SetFillColor(0);
    label->SetTextFont(42);
    label->SetTextSize(0.04);
    label->SetBorderSize(0);

    std::stringstream ss;
    ss << "Entries: " << std::fixed << std::setprecision(0) << h->GetEntries();
    label->AddText(ss.str().c_str());

    ss.str("");
    ss << "#chi^{2}/ndf = " << std::fixed << std::setprecision(2)
        << fit->GetChisquare() << "/" << fit->GetNDF()
        << " = " << fit->GetChisquare() / fit->GetNDF();
    label->AddText(ss.str().c_str());

    ss.str("");
    ss << "C = " << std::fixed << std::setprecision(1)
        << fit->GetParameter(0) << " #pm " << fit->GetParError(0);
    label->AddText(ss.str().c_str());

    ss.str("");
    ss << "A_{FB} = " << std::fixed << std::setprecision(4)
        << fit->GetParameter(1) << " #pm " << fit->GetParError(1);
    label->AddText(ss.str().c_str());

    label->Draw("same");

    // ---------------------------------------------------------
    // (RESIDUI NORMALIZZATI) ---
    pad2->cd();

    // Cloniamo l'istogramma per mantenere la stessa binning
    TH1D* h_res = (TH1D*)h->Clone("h_res");
    h_res->Reset(); // Svuotiamo il contenuto
    h_res->SetTitle("");

    // Calcolo dei residui normalizzati
    for (int i = 1; i <= h->GetNbinsX(); i++) {
        double data = h->GetBinContent(i);
        double err = h->GetBinError(i);
        double x = h->GetBinCenter(i);
        double fit_val = fit->Eval(x);

        if (err > 0) {
            h_res->SetBinContent(i, (data - fit_val) / err);
            h_res->SetBinError(i, 0.0); // Nessun errore sulle barre Y per chiarezza visiva
        }
    }

    // Configurazione assi per il pad inferiore (devono essere più grandi perché il pad è piccolo)
    h_res->GetYaxis()->SetTitle("(Data-Fit)/#sigma");
    h_res->GetYaxis()->SetNdivisions(505); // Configura 5 divisioni principali
    h_res->GetYaxis()->SetTitleSize(0.10);
    h_res->GetYaxis()->SetLabelSize(0.08);
    h_res->GetYaxis()->SetTitleOffset(0.4);

    h_res->GetXaxis()->SetTitle("cos#theta*");
    h_res->GetXaxis()->SetTitleSize(0.12);
    h_res->GetXaxis()->SetLabelSize(0.10);

    // Range asse Y per i residui (tipicamente tra -4 e +4 sigma)
    h_res->SetMinimum(-4.5);
    h_res->SetMaximum(4.5);

    // Stile dei punti
    h_res->SetMarkerStyle(20);
    h_res->SetMarkerSize(0.6);
    h_res->SetLineColor(kBlack);
    h_res->Draw("P"); // Disegna solo i punti

    // Linea di riferimento a Y=0 (fit perfetto)
    TLine* line = new TLine(-0.9, 0, 0.9, 0);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");

    // Stampa su terminale (invariata)
    std::cout
        << "chi2 = " << fit->GetChisquare()
        << "\nndf = " << fit->GetNDF()
        << "\nchi2/ndf = " << fit->GetChisquare() / fit->GetNDF()
        << "\nA_FB = " << fit->GetParameter(1)
        << " +- " << fit->GetParError(1)
        << std::endl;
}
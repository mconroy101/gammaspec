#include<iostream>
#include <fstream>

TGraph *xcom_dy = new TGraph("xcom_dy.csv");
const double PI = TMath::Pi();
TRandom3 *rndm = new TRandom3();

void test() {
    cout<<xcom_dy->Eval(0.1)<<endl;
}

double foil_calc_sq(double side_l, double thickness, double det_dist, double det_R, int N, bool attenuate = false, double density = 8.55, double E=100.10595) {

    int counts = 0;
    double mu_p = xcom_dy->Eval(E/1000);
    
    for(int i=0;i<N;i++) {
        double start[3]={side_l*rndm->Rndm() - side_l/2, side_l*rndm->Rndm() - side_l/2, thickness*rndm->Rndm()};
        double theta = 2*PI*rndm->Rndm(); // azimuthal angle
        double phi = acos(-1*rndm->Rndm()); // polar angle for downward hemisphere
        
        double dir[3] = {sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi)};
        double scale_factor = abs((det_dist+start[2])/dir[2]);
        double det_vec[3] = {0,0,0};
        for (int v=0;v<3;v++) {
            det_vec[v] = start[v] + scale_factor*dir[v];
        }
        if ((det_vec[0]*det_vec[0] + det_vec[1]*det_vec[1]) > det_R*det_R) { // miss
            continue;
        }
        else if (attenuate==false) {
            counts += 1;
        }
        else { // Check if attenuated
            double distance_sq = 0;
            double scale_factor = abs(start[2]/dir[2]);
            for (int u=0;u<3;u++) {
                double end = ((start[u] + scale_factor*dir[u]) - start[u]);
                distance_sq += end*end;
            }
            double survival_prop = exp(-mu_p*density*sqrt(distance_sq));
            if (survival_prop > rndm->Rndm()) {
                counts += 1;
            }
            
        }

    }

    return counts;
}

double point_calc(double det_dist=2, double det_R=2.335, int N=1000) {

    int counts = 0;
    
    for(int i=0;i<N;i++) {
        double start[3]={0, 0, 0};
        double theta = 2*PI*rndm->Rndm(); // azimuthal angle
        double phi = acos(-rndm->Rndm()); // polar angle for downward hemisphere
        
        double dir[3] = {sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi)};
        double scale_factor = abs((det_dist+start[2])/dir[2]);
        double det_vec[3] = {0,0,0};
        for (int v=0;v<3;v++) {
            det_vec[v] = start[v] + scale_factor*dir[v];
        }
        if ((det_vec[0]*det_vec[0] + det_vec[1]*det_vec[1]) > det_R*det_R) { // miss
            continue;
        }
        else {
            counts += 1;
        }

    }

    return counts;
}

void calc_correction(double distance) {
    int N = 5E7; // Number of iterations
    double L = 2.5; // Square foil side length
    double T = 0.0125; // Foil thickness
    double R = 2.465; // Detector crystal radius

    double N_foil = foil_calc_sq(L, T, distance, R, N);
    cout << "Foil counts: " << (N_foil) << endl;
    double N_point = point_calc(distance, R, N);
    cout << "Point counts: "<< (N_point) << endl;
    double correction = N_foil/N_point;
    cout << (correction) << endl;

}
